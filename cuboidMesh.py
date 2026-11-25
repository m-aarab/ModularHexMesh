#!/usr/bin/env python3
"""Generate a cuboid HEX8 VTK mesh from an external config file.

Usage:
  python cuboidMesh.py --config mesh.toml

Config format (TOML):

[mesh]
output_folder = "./"
output_filename = "cuboid_hex8.vtk"
sidelength = [10.0, 6.0, 2.0]

# Choose exactly one of the following:
# elemsize = [0.5, 0.5, 0.5]
n_elems = [20, 12, 4]

# Optional settings:
center = true
shift = [0.0, 0.0, 0.0]

Output naming rule:
	The script always appends the computed maximum element size to the filename,
	e.g. "cuboid_hex8_hmax0.5.vtk".
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
	import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover
	try:
		import tomli as tomllib  # type: ignore
	except ModuleNotFoundError as exc:  # pragma: no cover
		raise SystemExit(
			"TOML support requires Python 3.11+ (tomllib) or package 'tomli'."
		) from exc


Vec3f = Tuple[float, float, float]
Vec3i = Tuple[int, int, int]


def _get_first(mapping: Dict[str, object], keys: Iterable[str]) -> object | None:
	for key in keys:
		if key in mapping:
			return mapping[key]
	return None


def _parse_vec3_float(value: object, name: str) -> Vec3f:
	if not isinstance(value, (list, tuple)) or len(value) != 3:
		raise ValueError(f"'{name}' must be a 3-value array.")
	try:
		parsed = tuple(float(v) for v in value)
	except (TypeError, ValueError) as exc:
		raise ValueError(f"'{name}' must only contain numeric values.") from exc
	if any(v <= 0.0 for v in parsed):
		raise ValueError(f"'{name}' values must be > 0.")
	return parsed  # type: ignore[return-value]


def _parse_vec3_int(value: object, name: str) -> Vec3i:
	if not isinstance(value, (list, tuple)) or len(value) != 3:
		raise ValueError(f"'{name}' must be a 3-value array.")
	out: List[int] = []
	for item in value:
		if isinstance(item, bool):
			raise ValueError(f"'{name}' values must be integers > 0.")
		try:
			iv = int(item)
		except (TypeError, ValueError) as exc:
			raise ValueError(f"'{name}' values must be integers > 0.") from exc
		if iv <= 0:
			raise ValueError(f"'{name}' values must be integers > 0.")
		out.append(iv)
	return (out[0], out[1], out[2])


def _parse_shift(value: object | None) -> Vec3f:
	if value is None:
		return (0.0, 0.0, 0.0)
	if not isinstance(value, (list, tuple)) or len(value) != 3:
		raise ValueError("'shift' must be a 3-value array.")
	try:
		parsed = tuple(float(v) for v in value)
	except (TypeError, ValueError) as exc:
		raise ValueError("'shift' must only contain numeric values.") from exc
	return parsed  # type: ignore[return-value]


def _derive_nelems_from_elemsize(sidelength: Vec3f, elemsize: Vec3f) -> Vec3i:
	nelems: List[int] = []
	for axis, (length, size) in enumerate(zip(sidelength, elemsize), start=1):
		ratio = length / size
		rounded = round(ratio)
		if not math.isclose(ratio, rounded, rel_tol=1e-12, abs_tol=1e-12):
			raise ValueError(
				"'elemsize' must divide 'sidelength' exactly per axis. "
				f"Axis {axis}: sidelength={length}, elemsize={size}."
			)
		if rounded <= 0:
			raise ValueError("Derived element count must be > 0 on all axes.")
		nelems.append(int(rounded))
	return (nelems[0], nelems[1], nelems[2])


def _compute_elem_size(sidelength: Vec3f, nelems: Vec3i) -> Vec3f:
	return (
		sidelength[0] / nelems[0],
		sidelength[1] / nelems[1],
		sidelength[2] / nelems[2],
	)


def _append_hmax_to_filename(path: Path, hmax: float) -> Path:
	hmax_text = f"{hmax:.12g}"
	suffix = path.suffix if path.suffix else ".vtk"
	return path.with_name(f"{path.stem}_hmax{hmax_text}{suffix}")


def _load_config(config_path: Path) -> Tuple[Path, Vec3f, Vec3i, bool, Vec3f, Vec3f, float]:
	try:
		raw = tomllib.loads(config_path.read_text(encoding="utf-8"))
	except FileNotFoundError as exc:
		raise SystemExit(f"Config file not found: {config_path}") from exc
	except Exception as exc:
		raise SystemExit(f"Failed to parse TOML config '{config_path}': {exc}") from exc

	mesh_cfg = raw.get("mesh")
	if not isinstance(mesh_cfg, dict):
		raise ValueError("Config must contain [mesh] table.")

	output_raw = _get_first(mesh_cfg, ("output", "outfile", "output_file"))
	output_folder_raw = _get_first(mesh_cfg, ("output_folder", "output_dir", "outdir"))
	output_filename_raw = _get_first(
		mesh_cfg, ("output_filename", "filename", "outfile_name")
	)

	if output_raw is not None and (
		output_folder_raw is not None or output_filename_raw is not None
	):
		raise ValueError(
			"Use either [mesh].output or [mesh].output_folder/output_filename, not both."
		)

	side_raw = _get_first(mesh_cfg, ("sidelength", "side_length", "size"))
	if side_raw is None:
		raise ValueError("[mesh].sidelength is required.")
	sidelength = _parse_vec3_float(side_raw, "sidelength")

	nelems_raw = _get_first(mesh_cfg, ("nelems", "n_elems", "nElements"))
	elemsize_raw = _get_first(mesh_cfg, ("elemsize", "elem_size", "element_size"))

	if (nelems_raw is None and elemsize_raw is None) or (
		nelems_raw is not None and elemsize_raw is not None
	):
		raise ValueError(
			"Specify exactly one of [mesh].nelems or [mesh].elemsize."
		)

	if nelems_raw is not None:
		nelems = _parse_vec3_int(nelems_raw, "nelems")
	else:
		elemsize = _parse_vec3_float(elemsize_raw, "elemsize")
		nelems = _derive_nelems_from_elemsize(sidelength, elemsize)

	center_raw = _get_first(mesh_cfg, ("center", "center_at_origin", "centered"))
	center = True if center_raw is None else bool(center_raw)

	shift = _parse_shift(_get_first(mesh_cfg, ("shift", "origin_shift", "offset")))

	if output_raw is not None:
		if not isinstance(output_raw, str) or not output_raw.strip():
			raise ValueError("[mesh].output must be a non-empty path string.")
		out_base_path = Path(output_raw)
		if not out_base_path.is_absolute():
			out_base_path = (config_path.parent / out_base_path).resolve()
	else:
		if output_folder_raw is None:
			output_folder_raw = "."
		if output_filename_raw is None:
			output_filename_raw = "cuboid_hex8.vtk"

		if not isinstance(output_folder_raw, str) or not output_folder_raw.strip():
			raise ValueError("[mesh].output_folder must be a non-empty path string.")
		if not isinstance(output_filename_raw, str) or not output_filename_raw.strip():
			raise ValueError("[mesh].output_filename must be a non-empty filename.")

		out_folder_path = Path(output_folder_raw)
		if not out_folder_path.is_absolute():
			out_folder_path = (config_path.parent / out_folder_path).resolve()

		filename_path = Path(output_filename_raw)
		suffix = filename_path.suffix if filename_path.suffix else ".vtk"
		out_base_path = out_folder_path / f"{filename_path.stem}{suffix}"

	elemsize = _compute_elem_size(sidelength, nelems)
	hmax = max(elemsize)
	out_path = _append_hmax_to_filename(out_base_path, hmax)

	return out_path, sidelength, nelems, center, shift, elemsize, hmax


def _point_index(i: int, j: int, k: int, nx: int, ny: int) -> int:
	return i + (nx + 1) * (j + (ny + 1) * k)


def _build_points(sidelength: Vec3f, nelems: Vec3i, center: bool, shift: Vec3f) -> List[Vec3f]:
	lx, ly, lz = sidelength
	nx, ny, nz = nelems

	dx = lx / nx
	dy = ly / ny
	dz = lz / nz

	x0 = (-0.5 * lx if center else 0.0) + shift[0]
	y0 = (-0.5 * ly if center else 0.0) + shift[1]
	z0 = (-0.5 * lz if center else 0.0) + shift[2]

	points: List[Vec3f] = []
	for k in range(nz + 1):
		z = z0 + k * dz
		for j in range(ny + 1):
			y = y0 + j * dy
			for i in range(nx + 1):
				x = x0 + i * dx
				points.append((x, y, z))
	return points


def _build_hex8_cells(nelems: Vec3i) -> List[Tuple[int, int, int, int, int, int, int, int]]:
	nx, ny, nz = nelems
	cells: List[Tuple[int, int, int, int, int, int, int, int]] = []

	for k in range(nz):
		for j in range(ny):
			for i in range(nx):
				n0 = _point_index(i, j, k, nx, ny)
				n1 = _point_index(i + 1, j, k, nx, ny)
				n2 = _point_index(i + 1, j + 1, k, nx, ny)
				n3 = _point_index(i, j + 1, k, nx, ny)
				n4 = _point_index(i, j, k + 1, nx, ny)
				n5 = _point_index(i + 1, j, k + 1, nx, ny)
				n6 = _point_index(i + 1, j + 1, k + 1, nx, ny)
				n7 = _point_index(i, j + 1, k + 1, nx, ny)
				cells.append((n0, n1, n2, n3, n4, n5, n6, n7))
	return cells


def _write_vtk_unstructured_hex8(
	out_path: Path,
	points: Sequence[Vec3f],
	cells: Sequence[Tuple[int, int, int, int, int, int, int, int]],
) -> None:
	out_path.parent.mkdir(parents=True, exist_ok=True)
	with out_path.open("w", encoding="utf-8") as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write("cuboid HEX8 mesh\n")
		f.write("ASCII\n")
		f.write("DATASET UNSTRUCTURED_GRID\n")

		f.write(f"POINTS {len(points)} float\n")
		for x, y, z in points:
			f.write(f"{x:.16g} {y:.16g} {z:.16g}\n")

		cells_size = len(cells) * 9  # 8 node IDs + leading count per cell
		f.write(f"CELLS {len(cells)} {cells_size}\n")
		for c in cells:
			f.write("8 " + " ".join(str(v) for v in c) + "\n")

		f.write(f"CELL_TYPES {len(cells)}\n")
		for _ in cells:
			f.write("12\n")  # VTK_HEXAHEDRON


def _format_vec3(v: Sequence[float | int]) -> str:
	return f"({v[0]}, {v[1]}, {v[2]})"


def main() -> None:
	parser = argparse.ArgumentParser(
		description="Generate a cuboid VTK mesh with HEX8 elements from TOML config."
	)
	parser.add_argument(
		"--config",
		required=True,
		type=Path,
		help="Path to TOML config file.",
	)
	args = parser.parse_args()

	try:
		out_path, sidelength, nelems, center, shift, elemsize, hmax = _load_config(args.config)
	except ValueError as exc:
		raise SystemExit(f"Invalid config: {exc}") from exc

	points = _build_points(sidelength, nelems, center, shift)
	cells = _build_hex8_cells(nelems)
	_write_vtk_unstructured_hex8(out_path, points, cells)

	print("Generated HEX8 VTK mesh")
	print(f"  output      : {out_path}")
	print(f"  sidelength  : {_format_vec3(sidelength)}")
	print(f"  nelems      : {_format_vec3(nelems)}")
	print(f"  elemsize    : {_format_vec3(elemsize)}")
	print(f"  hmax        : {hmax:.12g}")
	print(f"  centered    : {center}")
	print(f"  shift       : {_format_vec3(shift)}")
	print(f"  points      : {len(points)}")
	print(f"  hex elements: {len(cells)}")


if __name__ == "__main__":
	main()
