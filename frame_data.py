from __future__ import annotations

from dataclasses import dataclass
from itertools import islice
from math import prod
import re

import h5py
import numpy as np

from hdf_layout import COORDINATES, ELEMENT_FIELDS, OFFSETS, normalize, text
from storage import hdf_frame_layout, read_xyz_atoms


MAX_FIELDS = 512
MAX_COMPONENTS = 2048
GEOMETRY_FIELDS = COORDINATES | set(ELEMENT_FIELDS) | set(OFFSETS) | {"smiles"}
ATOMIC_FIELDS = {"charges", "partialcharges", "forces", "force", "gradients", "gradient",
                 "velocities", "momenta", "spin", "spins", "magmoms"}


@dataclass(frozen=True)
class Field:
    name: str
    summary: str
    units: str
    detail: str


def format_scalar(value):
    if isinstance(value, (float, np.floating)):
        return format(value, ".12g")
    return text(value)


def make_field(name, value, units="", original_shape=None):
    if isinstance(value, dict):
        detail = str(value)[:16000]
        return Field(name, f"{len(value)} полей", units, detail)
    try:
        array = np.asarray(value)
    except (ValueError, TypeError):
        data = str(value)
        return Field(name, data[:130], units, data[:16000])
    shape = array.shape if original_shape is None else original_shape
    if not shape:
        data = format_scalar(array.item())
        return Field(name, data[:130], units, data[:16000])
    if array.size > MAX_COMPONENTS:
        array, _ = preview(array)
    summary = " × ".join(map(str, shape))
    if array.ndim == 2 and shape[-1] == 3 and array.dtype.kind in "iuf":
        detail = "  №        X            Y            Z\n"
        detail += "\n".join(f"{row:4}  " + "  ".join(f"{format_scalar(v):>12}" for v in vector)
                            for row, vector in enumerate(array, 1))
    else:
        detail = np.array2string(array, precision=12, threshold=MAX_COMPONENTS,
                                 max_line_width=112)
    if array.shape != shape:
        detail += f"\n\nПоказано первых {array.size} из {prod(shape)} значений."
    return Field(name, summary, units, detail[:80000])


def preview(array, leading=None):
    if isinstance(leading, (int, np.integer)):
        shape = array.shape[1:]
    elif isinstance(leading, slice):
        shape = (leading.stop - leading.start, *array.shape[1:])
    else:
        shape = array.shape

    extents = list(shape)
    budget = MAX_COMPONENTS
    for axis in reversed(range(len(extents))):
        extents[axis] = min(extents[axis], budget)
        budget = max(1, budget // max(1, extents[axis]))
    selection = tuple(slice(0, extent) for extent in extents)
    if isinstance(leading, (int, np.integer)):
        selection = (leading, *selection)
    elif isinstance(leading, slice):
        selection = (slice(leading.start, leading.start + extents[0]), *selection[1:])
    return array[selection], shape


def xyz_frame_fields(path, address):
    atoms = read_xyz_atoms(path, address)
    values = dict(atoms.info)
    if atoms.calc is not None:
        values.update(atoms.calc.results)
    values.update((name, value) for name, value in atoms.arrays.items()
                  if name not in {"numbers", "positions"})
    if np.any(atoms.cell):
        values["Lattice"] = atoms.cell.array
        values["pbc"] = atoms.pbc
    return [make_field(name, values[name]) for name in islice(sorted(values), MAX_FIELDS)]


def units_for(dataset, general_units):
    for key in ("units", "unit"):
        if key in dataset.attrs:
            return text(dataset.attrs[key])[:160]
    quantity = {"energy": "energy", "energies": "energy",
                "force": "forces", "forces": "forces"}.get(dataset.name.rsplit("/", 1)[-1].lower())
    return general_units.get(quantity, "")


def pick_dataset(dataset, frame, frames, natoms, atom_slice, total_atoms, ragged, trajectory):
    if dataset.shape is None or not dataset.size:
        return None
    if dataset.dtype.kind not in "biufcSOU" or (
            dataset.dtype.kind == "O" and h5py.check_string_dtype(dataset.dtype) is None):
        return None
    if not dataset.ndim:
        return None, True
    if not trajectory:
        return None, True
    name = dataset.name.rsplit("/", 1)[-1]
    for prefix, bounds, total in ragged:
        if name.startswith(prefix) and dataset.shape[0] == total:
            return bounds, True
    if atom_slice is not None and dataset.shape[0] == total_atoms:
        return atom_slice, True
    key = normalize(name)
    if key in {"cell", "lattice"} and dataset.shape == (3, 3):
        return None, True
    if (key in ATOMIC_FIELDS and dataset.shape[0] == natoms
            and (dataset.ndim == 1 or (dataset.ndim == 2 and dataset.shape[1] == 3))
            and (key in {"forces", "force", "gradients", "gradient", "velocities", "momenta"}
                 or dataset.ndim == 1)):
        return None, True
    if dataset.shape[0] == frames:
        return frame, True
    if dataset.shape[0] == natoms:
        return None, True
    return None


def ragged_fields(group, frame, frames):
    if frames < 2:
        return []
    matches = []
    for name in group:
        if not name.endswith(("_ptr", "_offsets")) or not isinstance(group.get(name, getlink=True), h5py.HardLink):
            continue
        dataset = group[name]
        if (not isinstance(dataset, h5py.Dataset) or dataset.shape != (frames + 1,)
                or dataset.dtype.kind not in "iu"):
            continue
        start, stop = map(int, dataset[frame:frame + 2])
        total = int(dataset[-1])
        if 0 <= start <= stop <= total:
            matches.append((name.rsplit("_", 1)[0] + "_", slice(start, stop), total))
    return matches


def local_datasets(group, cancelled, ancestors=frozenset()):
    marker = hash(group.id)
    if cancelled.is_set() or marker in ancestors:
        return
    ancestry = ancestors | {marker}
    for name in group:
        if cancelled.is_set():
            return
        if not isinstance(group.get(name, getlink=True), h5py.HardLink):
            continue
        obj = group[name]
        if isinstance(obj, h5py.Dataset):
            yield obj
        elif isinstance(obj, h5py.Group):
            other_geometry = any(normalize(key) in COORDINATES
                                 and obj.get(key, getclass=True) is h5py.Dataset for key in obj)
            if not other_geometry:
                yield from local_datasets(obj, cancelled, ancestry)


def hdf_frame_fields(file, address, offset, cancelled):
    layout, frame = hdf_frame_layout(file, address, offset)
    if layout is None:
        group = file[address]
        coordinates = group["conformations"]
        frames, natoms, _ = coordinates.shape
        excluded = {coordinates.name}
        total_atoms = atom_slice = None
    else:
        coordinates = layout.coordinates
        group = coordinates.parent
        if coordinates.name.endswith("/value") and normalize(group.name.rsplit("/", 1)[-1]) in COORDINATES:
            group = group.parent
        frames = layout.count
        atom_slice = layout.bounds(frame) if layout.offsets is not None else None
        natoms = (atom_slice.stop - atom_slice.start if atom_slice is not None
                  else layout.coordinates.shape[-1 if layout.transposed else -2])
        total_atoms = layout.elements.shape[0]
        excluded = {coordinates.name, layout.elements.name}
        if layout.offsets is not None:
            excluded.add(layout.offsets.name)
    prefix = group.name.rstrip("/") + "/"
    general_units = dict(re.findall(r"\b(energy|forces)\s*:\s*([^;,(]+)",
                                    text(file.attrs.get("units", "")), flags=re.I))
    general_units = {key.lower(): value.strip() for key, value in general_units.items()}
    fields = []
    for key in islice(group.attrs, MAX_FIELDS):
        if not cancelled.is_set():
            fields.append(make_field("@" + key, group.attrs[key]))
    ragged = ragged_fields(group, frame, frames)
    trajectory = coordinates.ndim == 3 or atom_slice is not None
    for dataset in local_datasets(group, cancelled):
        if len(fields) >= MAX_FIELDS:
            break
        name = dataset.name.removeprefix(prefix)
        leaf = name.rsplit("/", 1)[-1]
        if dataset.name in excluded or normalize(leaf) in GEOMETRY_FIELDS or leaf.endswith(("_ptr", "_offsets")):
            continue
        picked = pick_dataset(dataset, frame, frames, natoms, atom_slice, total_atoms, ragged, trajectory)
        if picked is None:
            continue
        selection, _ = picked
        try:
            value, shape = preview(dataset, selection)
            fields.append(make_field(name, value, units_for(dataset, general_units), shape))
        except (OSError, TypeError, ValueError, UnicodeError):
            fields.append(Field(name, "Не удалось прочитать", "", "Поле есть в файле, но его значение недоступно."))
    return fields


def read_frame_fields(source, location, cancelled):
    address, offset, _ = location
    source.check_unchanged()
    if source.is_xyz:
        fields = xyz_frame_fields(source.path, address)
    else:
        with h5py.File(source.path, "r") as file:
            fields = hdf_frame_fields(file, address, offset, cancelled)
    return sorted(fields, key=field_order)


def field_order(field):
    name = field.name.casefold()
    words = set(re.split(r"[\W_]+", name))
    if words & {"energy", "energies", "dlpno", "enthalpy"}:
        kind = 0
    elif words & {"force", "forces", "gradient", "gradients"}:
        kind = 1
    else:
        kind = 2
    return kind, name.startswith("models/"), name
