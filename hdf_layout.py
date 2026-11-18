from __future__ import annotations

from dataclasses import dataclass
from fnmatch import fnmatchcase
from functools import lru_cache
import posixpath
import re
from threading import Event

import h5py
import numpy as np
from ase.data import chemical_symbols


class HdfMappingNeeded(ValueError):
    pass


@dataclass(frozen=True)
class HdfMapping:
    coordinates: str
    elements: str
    offsets: str = ""
    units: str = ""


COORDINATES = {"coords", "coordinates", "coord", "coordinate", "positions", "position",
               "xyz", "conformations", "conformation", "geometries", "geometry",
               "cartesiancoordinates", "pos", "r"}
ELEMENT_FIELDS = ("atomicnumbers", "z", "numbers", "nuclearcharges", "symbols",
                  "elementsymbols", "chemicalsymbols", "elements", "species", "atomicspecies")
OFFSETS = ("atomoffsets", "frameoffsets", "atomptr", "frameptr", "offsets")
UNITS = {"nanometers": 10.0, "nanometer": 10.0, "nm": 10.0,
         "angstrom": 1.0, "angstroms": 1.0, "ang": 1.0, "å": 1.0, "a": 1.0,
         "bohr": 0.529177210903, "a0": 0.529177210903, "au": 0.529177210903}
SYMBOLS = frozenset(chemical_symbols[1:])


def normalize(key):
    return re.sub(r"[^a-z0-9]+", "", key.casefold())


def text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def scale_for(coordinates, units=""):
    if units:
        return UNITS[units.casefold()]
    node = coordinates
    keys = ("units", "unit", "coordinates_unit", "coordinate_unit",
            "position_unit", "positions_unit", "length_unit")
    while True:
        for key in keys:
            if key not in node.attrs:
                continue
            value = text(node.attrs[key]).strip().casefold()
            if value in UNITS:
                return UNITS[value]
            match = re.search(r"(?:coord(?:inate)?s?|positions?|length)\s*[:=]\s*([^;,\s]+)", value)
            if match:
                value = match[1]
                if value in UNITS:
                    return UNITS[value]
            if key not in ("units", "unit") or node is coordinates:
                raise HdfMappingNeeded(f"{coordinates.name}: неизвестная единица длины {value!r}.")
        if node.name == "/":
            return 1.0
        node = node.parent


def coordinate_shape(dataset):
    if (not isinstance(dataset, h5py.Dataset) or dataset.shape is None
            or dataset.dtype.kind not in "fiu" or dataset.ndim not in (2, 3)
            or not all(dataset.shape)):
        return None
    if dataset.shape[-1] == 3:
        return dataset.shape[-2], False
    if dataset.shape[-2] == 3:
        return dataset.shape[-1], True
    return None


def compatible(dataset, shapes):
    return (isinstance(dataset, h5py.Dataset) and dataset.shape in shapes
            and (dataset.dtype.kind in "iu" or h5py.check_string_dtype(dataset.dtype) is not None))


def hierarchy(root, cancelled, ancestors=frozenset()):
    if cancelled.is_set():
        return
    marker = hash(root.id)
    if marker in ancestors:
        return
    yield root
    ancestry = ancestors | {marker}
    for name in root:
        if cancelled.is_set():
            return
        if not isinstance(root.get(name, getlink=True), h5py.HardLink):
            continue
        if root.get(name, getclass=True) is h5py.Group:
            yield from hierarchy(root[name], cancelled, ancestry)


def named_fields(group, names):
    for key in group:
        alias = normalize(key)
        if alias not in names or not isinstance(group.get(key, getlink=True), h5py.HardLink):
            continue
        obj = group[key]
        if isinstance(obj, h5py.Group) and isinstance(obj.get("value", getlink=True), h5py.HardLink):
            obj = obj["value"]
        if isinstance(obj, h5py.Dataset):
            yield obj, alias


def element_fields(group):
    return tuple(named_fields(group, ELEMENT_FIELDS))


def find_elements(coordinates, shapes, fields):
    group = coordinates.parent
    while True:
        candidates = [(ELEMENT_FIELDS.index(alias), dataset)
                      for dataset, alias in fields(group)
                      if compatible(dataset, shapes)]
        if candidates:
            candidates.sort(key=lambda pair: pair[0])
            if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
                raise HdfMappingNeeded(f"{coordinates.name}: найдено несколько полей элементов.")
            return candidates[0][1]
        if group.name == "/":
            return None
        group = group.parent


def mapped_field(group, name):
    return group.file.get(posixpath.normpath(posixpath.join(group.name, name))) if name else None


def valid_offsets(dataset, atom_count):
    return (isinstance(dataset, h5py.Dataset) and dataset.ndim == 1 and dataset.size >= 2
            and dataset.dtype.kind in "iu" and int(dataset[0]) == 0 and int(dataset[-1]) == atom_count)


@dataclass(frozen=True)
class Layout:
    coordinates: h5py.Dataset
    elements: h5py.Dataset
    offsets: h5py.Dataset | None
    scale: float
    transposed: bool
    per_frame_elements: bool

    @property
    def count(self):
        if self.offsets is not None:
            return self.offsets.size - 1
        return self.coordinates.shape[0] if self.coordinates.ndim == 3 else 1

    def bounds(self, frame):
        start, stop = map(int, self.offsets[frame:frame + 2])
        natoms = self.elements.shape[0]
        if not 0 <= start < stop <= natoms:
            raise ValueError(f"{self.offsets.name}: неверные границы кадра {frame + 1}.")
        return slice(start, stop)

    def symbols(self, frame):
        if self.offsets is not None:
            selection = self.bounds(frame)
        elif self.per_frame_elements:
            selection = frame
        elif self.elements.ndim == 2:
            selection = 0
        else:
            selection = slice(None)
        values = self.elements[selection]
        if values.dtype.kind in "iu":
            if np.any(values < 1) or np.any(values >= len(chemical_symbols)):
                raise HdfMappingNeeded(f"{self.elements.name}: ожидаются атомные номера от 1 до 118.")
            return tuple(chemical_symbols[int(z)] for z in values)
        result = tuple(text(value).strip().capitalize() for value in values)
        if any(value not in SYMBOLS for value in result):
            raise HdfMappingNeeded(f"{self.elements.name}: ожидаются символы элементов.")
        return result

    def positions(self, frame):
        if self.offsets is not None:
            rows = self.bounds(frame)
            xyz = self.coordinates[:, rows].T if self.transposed else self.coordinates[rows]
        elif self.coordinates.ndim == 3:
            xyz = self.coordinates[frame]
            if self.transposed:
                xyz = xyz.T
        else:
            xyz = self.coordinates[:]
            if self.transposed:
                xyz = xyz.T
        return np.asarray(xyz, dtype=float) * self.scale


def make_layout(coordinates, mapping=None, fields=element_fields):
    shape = coordinate_shape(coordinates)
    if shape is None:
        return None
    atom_count, transposed = shape
    frames = coordinates.shape[0] if coordinates.ndim == 3 else 1
    shapes = {(atom_count,), (1, atom_count)}
    if frames > 1:
        shapes.add((frames, atom_count))
    group = coordinates.parent
    elements = mapped_field(group, mapping.elements) if mapping else find_elements(coordinates, shapes, fields)
    if not compatible(elements, shapes):
        if mapping:
            raise HdfMappingNeeded(f"{coordinates.name}: поле элементов отсутствует или имеет несовместимую форму.")
        return None
    offsets = None
    if coordinates.ndim == 2 and elements.ndim == 1:
        if mapping and mapping.offsets:
            offsets = mapped_field(group, mapping.offsets)
            if not valid_offsets(offsets, atom_count):
                raise HdfMappingNeeded(f"{coordinates.name}: ожидался массив смещений от 0 до {atom_count}.")
        elif mapping is None:
            candidates = [x for x, _ in named_fields(group, OFFSETS) if valid_offsets(x, atom_count)]
            if len(candidates) > 1:
                raise HdfMappingNeeded(f"{coordinates.name}: найдено несколько массивов смещений.")
            if candidates:
                offsets = candidates[0]
    elif mapping and mapping.offsets:
        raise HdfMappingNeeded("Смещения применяются к плоскому массиву координат и элементов.")
    return Layout(coordinates, elements, offsets, scale_for(coordinates, mapping.units if mapping else ""),
                  transposed, elements.ndim == 2 and elements.shape[0] > 1)


def file_geometries(file, cancelled, mapping=None):
    if mapping:
        dataset = file.get(mapping.coordinates)
        if dataset is not None or not any(char in mapping.coordinates for char in "*?["):
            if (layout := make_layout(dataset, mapping)) is not None:
                yield layout
            return
    fields = lru_cache(maxsize=128)(element_fields)
    for group in hierarchy(file, cancelled):
        if mapping:
            for key in group:
                path = posixpath.join(group.name, key)
                if fnmatchcase(path, mapping.coordinates) and isinstance(group.get(key, getlink=True), h5py.HardLink):
                    layout = make_layout(group[key], mapping)
                    if layout is not None:
                        yield layout
            continue
        layouts = [layout for dataset, _ in named_fields(group, COORDINATES)
                   if (layout := make_layout(dataset, fields=fields)) is not None]
        if len(layouts) > 1:
            raise HdfMappingNeeded(f"{group.name}: найдено несколько полей координат.")
        if layouts:
            yield from layouts
        elif "conformations" in group and "smiles" in group:
            yield group


def field_choices(path, limit=1500):
    choices = [[], [], []]
    with h5py.File(path, "r") as file:
        seen = 0
        for group in hierarchy(file, Event()):
            for key in group:
                if not isinstance(group.get(key, getlink=True), h5py.HardLink):
                    continue
                if group.get(key, getclass=True) is not h5py.Dataset:
                    continue
                dataset = group[key]
                if dataset.shape is None:
                    continue
                alias = normalize(key)
                item = (dataset.name, f"{dataset.shape} · {dataset.dtype}")
                if coordinate_shape(dataset):
                    choices[0].append((alias not in COORDINATES, item))
                kind = dataset.dtype.kind
                if dataset.ndim in (1, 2) and (kind in "iu" or h5py.check_string_dtype(dataset.dtype)):
                    choices[1].append((alias not in ELEMENT_FIELDS, item))
                if dataset.ndim == 1 and kind in "iu":
                    choices[2].append((alias not in OFFSETS, item))
                seen += 1
                if seen >= limit:
                    break
            if seen >= limit:
                break
    return [[item for _, item in sorted(items)[:200]] for items in choices]
