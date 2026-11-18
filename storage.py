from __future__ import annotations

from collections import Counter, OrderedDict
from contextlib import closing
from dataclasses import dataclass
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import re
import sqlite3
import time
from urllib.parse import unquote

import h5py
import numpy as np
from ase import Atoms
from ase.data import chemical_symbols, covalent_radii
from ase.io.extxyz import read_xyz
from ase.neighborlist import neighbor_list
from rdkit import Chem

from hdf_layout import HdfMapping, HdfMappingNeeded, Layout, file_geometries, scale_for


HDF_INDEX_VERSION = 2
ELEMENTS = frozenset(chemical_symbols[1:])
PROPERTIES = re.compile(
    rb"""(?:^|\s)Properties\s*=\s*(?:"([^"]*)"|'([^']*)'|\{([^}]*)\}|\[([^\]]*)\]|([^\s]+))"""
)


@dataclass(frozen=True)
class Source:
    path: Path
    cache_path: Path
    mapping: HdfMapping | None = None

    @classmethod
    def open(cls, path, cache_dir, mapping=None):
        path = Path(path).resolve(strict=True)
        if not path.is_file():
            raise ValueError("Выберите файл.")
        if path.suffix.lower() not in {".hdf5", ".h5", ".extxyz"}:
            raise ValueError("Поддерживаются .hdf5, .h5 и .extxyz.")
        return cls(path, cls.index_path(path, cache_dir, mapping), mapping)

    @staticmethod
    def index_path(path, cache_dir, mapping=None):
        stat = path.stat()
        identity = (str(path), stat.st_size, stat.st_mtime_ns,
                    stat.st_ctime_ns, stat.st_dev, stat.st_ino)
        signature = ((1, *identity) if path.suffix.lower() == ".extxyz"
                     else (HDF_INDEX_VERSION, *identity, mapping))
        digest = sha256(repr(signature).encode("utf-8")).hexdigest()
        return cache_dir / f"{digest}.index"

    def check_unchanged(self):
        if self.cache_path != self.index_path(self.path, self.cache_path.parent, self.mapping):
            raise ValueError("Файл изменился. Откройте его снова.")

    @property
    def is_xyz(self):
        return self.path.suffix.lower() == ".extxyz"


@dataclass
class Frame:
    xyz: np.ndarray
    elements: tuple
    bonds: tuple
    fragments: tuple
    description: str
    inferred_bonds: bool = False

    @property
    def size(self):
        return self.xyz.nbytes + len(self.elements) * 80 + len(self.bonds) * 120


class Index:
    def __init__(self, path):
        self.path = path

    def connect(self):
        return closing(sqlite3.connect(self.path.resolve().as_uri() + "?mode=ro",
                                       uri=True, timeout=0.2))

    def stats(self):
        with self.connect() as db:
            complete, frames = db.execute("SELECT complete, frames FROM state").fetchone()
            groups = db.execute("SELECT COALESCE(MAX(id), 0) FROM formulas").fetchone()[0]
        return complete, frames, groups

    def formulas(self, start, size):
        with self.connect() as db:
            return db.execute(
                "SELECT id, formula, frames FROM formulas WHERE id >= ? ORDER BY id LIMIT ?",
                (start, size),
            ).fetchall()

    def find(self, prefix):
        with self.connect() as db:
            return db.execute(
                "SELECT id FROM formulas WHERE formula GLOB ? ORDER BY formula LIMIT 1",
                (prefix + "*",),
            ).fetchone()

    def location(self, group_id, frame_number):
        with self.connect() as db:
            row = db.execute(
                """SELECT address, last_frame, count, original FROM spans
                   WHERE formula_id = ? AND last_frame >= ? ORDER BY last_frame LIMIT 1""",
                (group_id, frame_number),
            ).fetchone()
        if row is None:
            raise ValueError("Кадр ещё не проиндексирован.")
        address, last, count, original = row
        offset = frame_number - (last - count) - 1
        return address, offset, original + offset


def formula(symbols):
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(counts.keys() - set(order))
    else:
        order = sorted(counts)
    return "".join(symbol + (str(counts[symbol]) if counts[symbol] > 1 else "")
                   for symbol in order)


def atom_count(line):
    try:
        count = int(line.removeprefix(b"\xef\xbb\xbf"))
    except ValueError:
        raise ValueError("Ожидалось число атомов в начале кадра.") from None
    if count < 1:
        raise ValueError("Кадр должен содержать хотя бы один атом.")
    return count


def species_columns(header):
    match = PROPERTIES.search(header)
    if not match:
        return 0, None
    fields = next(value for value in match.groups() if value is not None).split(b":")
    if len(fields) % 3:
        raise ValueError("Некорректное поле Properties в заголовке кадра.")
    offset = 0
    species = atomic_number = None
    for name, _, raw_width in zip(fields[::3], fields[1::3], fields[2::3]):
        width = int(raw_width)
        if width < 1:
            raise ValueError("Некорректная ширина поля Properties.")
        if name == b"species" and width == 1:
            species = offset
        elif name == b"Z" and width == 1:
            atomic_number = offset
        offset += width
    if species is None and atomic_number is None:
        raise ValueError("В Properties отсутствует species или Z.")
    return species, atomic_number


def xyz_records(path, cancelled):
    symbol_bytes = {symbol.encode("ascii"): symbol for symbol in ELEMENTS}
    with path.open("rb") as stream:
        number = 0
        while not cancelled.is_set():
            address = stream.tell()
            line = stream.readline()
            if not line:
                return
            if not line.strip():
                continue
            number += 1
            try:
                count = atom_count(line)
                header = stream.readline()
                if not header:
                    raise ValueError("Отсутствует строка заголовка.")
                species, atomic_number = species_columns(header)
                column = species if species is not None else atomic_number
                elements = []
                for atom in range(count):
                    if atom % 4096 == 0 and cancelled.is_set():
                        return
                    fields = stream.readline().split(None, column + 1)
                    if len(fields) <= column:
                        raise ValueError(f"Неполная строка атома {atom + 1}.")
                    token = fields[column]
                    if species is None:
                        z = int(token)
                        if not 0 < z < len(chemical_symbols):
                            raise ValueError(f"Неизвестный атомный номер {z}.")
                        elements.append(chemical_symbols[z])
                    else:
                        try:
                            elements.append(symbol_bytes[token.capitalize()])
                        except KeyError:
                            raise ValueError(f"Неизвестный элемент {token!r}.") from None
            except ValueError as error:
                raise ValueError(f"Кадр {number}: {error}") from error
            yield formula(elements), address, 1


def hdf_topology(group):
    smiles_data = group["smiles"]
    smiles = smiles_data.asstr()[(0,) * smiles_data.ndim]
    params = Chem.SmilesParserParams()
    params.removeHs = False
    molecule = Chem.MolFromSmiles(smiles, params)
    coordinates = group["conformations"]
    if molecule is None or coordinates.ndim != 3 or coordinates.shape[-1] != 3:
        raise ValueError(f"{group.name}: некорректные SMILES или форма координат.")
    indices = [atom.GetAtomMapNum() - 1 for atom in molecule.GetAtoms()]
    if sorted(indices) != list(range(coordinates.shape[1])):
        raise ValueError(f"{group.name}: нумерация атомов в SMILES и координатах расходится.")
    return molecule, indices, scale_for(coordinates)


def geometry_address(layout, first):
    return json.dumps((layout.coordinates.name, layout.elements.name,
                       layout.offsets.name if layout.offsets is not None else None,
                       layout.scale, layout.transposed, layout.per_frame_elements, first),
                      ensure_ascii=False, separators=(",", ":"))


def hdf_records(source, cancelled):
    with h5py.File(source.path, "r") as file:
        for layout in file_geometries(file, cancelled, source.mapping):
            if isinstance(layout, h5py.Group):
                molecule, _, _ = hdf_topology(layout)
                count = layout["conformations"].shape[0]
                if count:
                    yield formula(atom.GetSymbol() for atom in molecule.GetAtoms()), layout.name, count
                continue
            if layout.offsets is None and not layout.per_frame_elements:
                yield formula(layout.symbols(0)), geometry_address(layout, 0), layout.count
                continue
            for frame in range(layout.count):
                if cancelled.is_set():
                    return
                yield formula(layout.symbols(frame)), geometry_address(layout, frame), 1

def build_index(source, destination, cancelled, progress):
    records = xyz_records(source.path, cancelled) if source.is_xyz else hdf_records(source, cancelled)
    with closing(sqlite3.connect(destination)) as db:
        db.executescript("""
            PRAGMA journal_mode = WAL;
            PRAGMA synchronous = NORMAL;
            CREATE TABLE formulas (id INTEGER PRIMARY KEY, formula TEXT UNIQUE, frames INTEGER);
            CREATE TABLE spans (
                formula_id INTEGER, last_frame INTEGER, count INTEGER, address, original INTEGER,
                PRIMARY KEY (formula_id, last_frame)
            ) WITHOUT ROWID;
            CREATE TABLE state (complete INTEGER, frames INTEGER);
            INSERT INTO state VALUES (0, 0);
        """)
        total = 0
        last_notify = time.monotonic()
        with closing(records) as stream:
            for name, address, count in stream:
                group, last = db.execute(
                    """INSERT INTO formulas (formula, frames) VALUES (?, ?)
                       ON CONFLICT(formula) DO UPDATE SET frames = frames + excluded.frames
                       RETURNING id, frames""",
                    (name, count),
                ).fetchone()
                db.execute("INSERT INTO spans VALUES (?, ?, ?, ?, ?)",
                           (group, last, count, address, total + 1))
                first = total == 0
                total += count
                if first or time.monotonic() - last_notify >= 0.35:
                    db.execute("UPDATE state SET frames = ?", (total,))
                    db.commit()
                    progress()
                    last_notify = time.monotonic()
        if cancelled.is_set():
            return False
        source.check_unchanged()
        if total == 0:
            if source.is_xyz:
                raise ValueError("В файле нет кадров.")
            raise HdfMappingNeeded(
                "Не удалось подобрать координаты и элементы. Укажите поля HDF5 вручную.")
        db.execute("UPDATE state SET complete = 1, frames = ?", (total,))
        db.commit()
        db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    return True


def read_xyz_atoms(path, address):
    with path.open("rb") as stream:
        stream.seek(int(address))
        start = stream.readline()
        count = atom_count(start)
        lines = [start, stream.readline()]
        lines.extend(stream.readline() for _ in range(count))
    return next(read_xyz(StringIO(b"".join(lines).decode("utf-8-sig")), index=0))


def hdf_frame_layout(file, address, offset):
    if not address.startswith("["):
        return None, offset
    coordinates, elements, offsets, scale, transpose, per_frame, first = json.loads(address)
    layout = Layout(file[coordinates], file[elements], file[offsets] if offsets else None,
                    scale, transpose, per_frame)
    return layout, first + offset


def load_frame(source, location):
    address, offset, original = location
    source.check_unchanged()
    if source.is_xyz:
        return inferred_frame(read_xyz_atoms(source.path, address), f"Кадр файла № {original:,}")

    with h5py.File(source.path, "r") as file:
        layout, frame = hdf_frame_layout(file, address, offset)
        if layout is not None:
            group = layout.coordinates.parent.name.strip("/")
            description = f"Кадр файла № {original:,}"
            if group:
                parts = " · ".join(unquote(part) for part in group.split("/"))
                description = f"{parts} · {description}"
            atoms = Atoms(layout.symbols(frame), positions=layout.positions(frame))
            return inferred_frame(atoms, description)
        group = file[address]
        molecule, indices, scale = hdf_topology(group)
        xyz = np.asarray(group["conformations"][offset], dtype=float) * scale
    ensure_coordinates(xyz)
    elements = [""] * len(indices)
    fragments = [0] * len(indices)
    for component, fragment in enumerate(Chem.GetMolFrags(molecule)):
        for atom_index in fragment:
            index = indices[atom_index]
            elements[index] = molecule.GetAtomWithIdx(atom_index).GetSymbol()
            fragments[index] = component
    bonds = tuple((indices[bond.GetBeginAtomIdx()], indices[bond.GetEndAtomIdx()])
                  for bond in molecule.GetBonds())
    return Frame(xyz, tuple(elements), bonds, tuple(fragments),
                 f"{address.lstrip('/')} · кадр группы № {offset + 1:,}")


def inferred_frame(atoms, description):
    xyz = atoms.get_positions()
    ensure_coordinates(xyz)
    # Use the coordinates as listed, including for inputs that specify a periodic cell.
    atoms.set_pbc(False)
    starts, ends = neighbor_list("ij", atoms, 1.12 * covalent_radii[atoms.numbers],
                                 max_nbins=max(1, len(atoms)))
    bonds = tuple((int(a), int(b)) for a, b in zip(starts, ends) if a < b)
    return Frame(xyz, tuple(atoms.get_chemical_symbols()), bonds, (), description, True)


def ensure_coordinates(xyz):
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("В кадре отсутствуют конечные координаты атомов.")


class FrameCache:
    def __init__(self, capacity=32 * 1024 * 1024):
        self.frames = OrderedDict()
        self.capacity = capacity
        self.size = 0

    def get(self, key):
        if key in self.frames:
            self.frames.move_to_end(key)
            return self.frames[key]
        return None

    def add(self, key, frame):
        if key in self.frames:
            self.size -= self.frames.pop(key).size
        if frame.size > self.capacity:
            return
        self.frames[key] = frame
        self.size += frame.size
        while self.size > self.capacity or len(self.frames) > 7:
            self.size -= self.frames.popitem(last=False)[1].size
