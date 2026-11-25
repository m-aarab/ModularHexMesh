#!/usr/bin/env python3
"""Generate a legacy-compatible cylindrical HEX8 mesh from a TOML config.

This script is a clean refactor of the old custom generator and keeps the same
core cross-section construction (square core + morphed cylindrical rings).

Usage:
  python cylindricalMesh.py --config cylindricalMesh.example.toml

Config format (TOML):

[mesh]
output_folder = "./output"
output_filename = "cylinder.vtk"
append_hmax = true

rmax = 300.0
length = 600.0
length_axis = "z"  # one of: x, y, z

# Set either a general mesh size or a manual angular resolution h.
mesh_size = 10.0
# h = 48

# Cross-section layout (default "spoke"):
# layout = "spoke"   # bulged square core + rings on straight spokes to the circle
# core_fraction = 0.55 # core half-width / rmax
# core_bulge = 0.134 # mid-side bow of the core edges / half-width; tan(15 deg)/2
#                    # makes the core corners 120 deg (worst scaled Jacobian 0.866)
# nr = 40            # radial layers - 1; default sizes them to mesh_size
#
# The defaults do not depend on resolution. Worst cells are the three meeting at
# each core corner (all 120 deg) and the rim cells on the diagonals, which are
# radially squeezed because the diagonal spoke is the shortest:
#   least stretched cells : pi/(4 c) = (1 - c(1+b)) / (1 - sqrt(2) c)  -> c ~ 0.51
#   most uniform sizes    : c ~ 0.58;  0.55 is within ~15 % of both.
#
# layout = "legacy" reproduces meshes from before the spoke layout and uses:
# rbase = 225.0
# nr = 25
# tune_transition = true # auto-tune square core and outer region for size uniformity
# axis_start = 300.0 # optional override; default centers mesh around 0
# shift = [0.0, 0.0, 0.0]

# Axial discretization modes:
# axial_mode = "legacy"       # default, reproduces old stacking behavior
# axial_mode = "fixed_layers" # requires n_layers
# axial_mode = "fixed_dz"     # requires dz
# axial_mode = "mesh_size"    # uses mesh_size to set dz close to target
# n_layers = 88
# dz = 6.9

[fields]
write = true
n1 = 1.4
n2 = 1.0
amax = 0.5
amin = 0.02
profile_radius = 1.0
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
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


# Spoke layout defaults (resolution independent, see the module docstring):
# bulge tan(15 deg)/2 makes the core corner 120 deg, so the three cells meeting there
# all have 120-deg corners (worst scaled Jacobian sin 120 = 0.866, the topology's best).
CORE_BULGE_DEFAULT = math.tan(math.radians(15.0)) / 2.0
# Core half-width / rmax: within ~15 % of the best on max aspect ratio, size ratio and
# size spread at once (least stretched cells at ~0.51, most uniform sizes at ~0.58).
CORE_FRACTION_DEFAULT = 0.55

Vec2f = Tuple[float, float]
Vec3f = Tuple[float, float, float]
Quad = Tuple[int, int, int, int]


@dataclass(frozen=True)
class FieldConfig:
	write: bool
	n1: float
	n2: float
	amax: float
	amin: float
	profile_radius: float


@dataclass(frozen=True)
class MeshConfig:
	out_base_path: Path
	append_hmax: bool
	layout: str
	core_bulge: float
	rmax: float
	mesh_size: float | None
	uniformity_target_size: float
	tuning_enabled: bool
	rbase: float
	h: int
	nr: int
	ntheta: int
	nxy: int
	circumferential_size: float
	theta_inner_size: float
	theta_outer_size: float
	radial_size: float
	square_size: float
	transition_size: float
	squaremax: float
	length: float
	length_axis: str
	axis_start: float
	shift: Vec3f
	axial_mode: str
	n_layers: int
	dz: float
	legacy_target_dz: float | None


def _get_first(mapping: Dict[str, object], keys: Iterable[str]) -> object | None:
	for key in keys:
		if key in mapping:
			return mapping[key]
	return None


def _parse_positive_float(value: object, name: str) -> float:
	try:
		out = float(value)
	except (TypeError, ValueError) as exc:
		raise ValueError(f"'{name}' must be numeric.") from exc
	if out <= 0.0:
		raise ValueError(f"'{name}' must be > 0.")
	return out


def _parse_positive_int(value: object, name: str) -> int:
	if isinstance(value, bool):
		raise ValueError(f"'{name}' must be an integer > 0.")
	try:
		out = int(value)
	except (TypeError, ValueError) as exc:
		raise ValueError(f"'{name}' must be an integer > 0.") from exc
	if out <= 0:
		raise ValueError(f"'{name}' must be an integer > 0.")
	return out


def _parse_vec3(value: object | None, name: str, default: Vec3f) -> Vec3f:
	if value is None:
		return default
	if not isinstance(value, (list, tuple)) or len(value) != 3:
		raise ValueError(f"'{name}' must be a 3-value array.")
	try:
		out = (float(value[0]), float(value[1]), float(value[2]))
	except (TypeError, ValueError) as exc:
		raise ValueError(f"'{name}' must only contain numeric values.") from exc
	return out


def _resolve_output_base_path(mesh_cfg: Dict[str, object], config_path: Path) -> Path:
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

	if output_raw is not None:
		if not isinstance(output_raw, str) or not output_raw.strip():
			raise ValueError("[mesh].output must be a non-empty path string.")
		out = Path(output_raw)
		if not out.is_absolute():
			out = (config_path.parent / out).resolve()
		if not out.suffix:
			out = out.with_suffix(".vtk")
		return out

	if output_folder_raw is None:
		output_folder_raw = "."
	if output_filename_raw is None:
		output_filename_raw = "cylinder.vtk"

	if not isinstance(output_folder_raw, str) or not output_folder_raw.strip():
		raise ValueError("[mesh].output_folder must be a non-empty path string.")
	if not isinstance(output_filename_raw, str) or not output_filename_raw.strip():
		raise ValueError("[mesh].output_filename must be a non-empty filename.")

	out_folder = Path(output_folder_raw)
	if not out_folder.is_absolute():
		out_folder = (config_path.parent / out_folder).resolve()

	filename_path = Path(output_filename_raw)
	suffix = filename_path.suffix if filename_path.suffix else ".vtk"
	return out_folder / f"{filename_path.stem}{suffix}"


def _append_hmax_to_filename(path: Path, hmax: float) -> Path:
	return path.with_name(f"{path.stem}_hmax{hmax:.12g}{path.suffix or '.vtk'}")


def _parse_length_axis(value: object | None) -> str:
	if value is None:
		return "z"
	axis = str(value).strip().lower()
	if axis not in ("x", "y", "z"):
		raise ValueError("'mesh.length_axis' must be one of: x, y, z.")
	return axis


def _derive_h_from_mesh_size(rmax: float, mesh_size: float) -> int:
	circumference = 2.0 * math.pi * rmax
	ntheta_target = max(8, int(round(circumference / mesh_size)))
	ntheta = max(8, int(round(ntheta_target / 8.0) * 8))
	return max(1, ntheta // 8)


def _compute_squaremax(rmax: float, rbase: float, nr: int) -> float:
	dr = (rmax - rbase) / nr
	return (rbase - 2.0 * dr) / math.sqrt(2.0)


def _transition_uniformity_objective(
	square_size: float,
	radial_size: float,
	theta_inner_size: float,
	theta_outer_size: float,
	transition_size: float,
	target_size: float,
) -> float:
	# Balance target-size matching with region-to-region uniformity.
	e_square = (square_size - target_size) / target_size
	e_radial = (radial_size - target_size) / target_size
	e_theta_in = (theta_inner_size - target_size) / target_size
	e_theta_out = (theta_outer_size - target_size) / target_size
	e_transition = (transition_size - target_size) / target_size

	u_square_radial = (square_size - radial_size) / target_size
	u_square_theta = (square_size - theta_inner_size) / target_size
	u_radial_theta = (radial_size - theta_inner_size) / target_size
	u_theta_in_out = (theta_inner_size - theta_outer_size) / target_size
	u_transition_radial = (transition_size - radial_size) / target_size
	u_transition_theta = (transition_size - theta_inner_size) / target_size
	u_transition_square = (transition_size - square_size) / target_size

	fit_score = (
		e_square * e_square
		+ e_radial * e_radial
		+ 0.8 * e_theta_in * e_theta_in
		+ 0.6 * e_theta_out * e_theta_out
		+ 1.8 * e_transition * e_transition
	)

	uniformity_score = (
		u_square_radial * u_square_radial
		+ u_square_theta * u_square_theta
		+ u_radial_theta * u_radial_theta
		+ 0.4 * u_theta_in_out * u_theta_in_out
		+ 0.8 * u_transition_radial * u_transition_radial
		+ 0.7 * u_transition_theta * u_transition_theta
		+ 0.6 * u_transition_square * u_transition_square
	)

	return 0.5 * uniformity_score + 0.5 * fit_score


def _tune_transition_region(
	rmax: float,
	nxy: int,
	ntheta: int,
	target_size: float,
	rbase_hint: float | None,
	nr_hint: int | None,
) -> Tuple[float, int, float, float]:
	if target_size <= 0.0:
		raise ValueError("Uniformity target size must be > 0.")

	max_nr = max(6, int(math.ceil(2.0 * rmax / target_size)))
	max_nr = min(max_nr, 5000)

	if nr_hint is not None:
		nr_values: Sequence[int] = [nr_hint]
	else:
		# Avoid nr=1 in auto mode; it usually creates oversized transition quads.
		nr_values = list(range(2, max_nr + 1))

	best_obj: float | None = None
	best_rbase = 0.0
	best_nr = 0
	best_square_size = 0.0
	best_radial_size = 0.0

	for nr in nr_values:
		rbase_min = 2.0 * rmax / (nr + 2.0) + 1e-12
		rbase_max = rmax - 1e-12
		if rbase_min >= rbase_max:
			continue

		if rbase_hint is not None:
			rbase_values: Sequence[float] = [rbase_hint]
		else:
			n_samples = 240
			step = (rbase_max - rbase_min) / (n_samples - 1)
			rbase_values = [rbase_min + i * step for i in range(n_samples)]

		for rbase in rbase_values:
			if rbase <= rbase_min or rbase >= rbase_max:
				continue

			radial_size = (rmax - rbase) / nr
			theta_inner_size = 2.0 * math.pi * rbase / ntheta
			theta_outer_size = 2.0 * math.pi * rmax / ntheta
			squaremax = _compute_squaremax(rmax=rmax, rbase=rbase, nr=nr)
			if squaremax <= 0.0:
				continue

			square_size = 2.0 * squaremax / nxy
			transition_size = rbase - math.sqrt(2.0) * squaremax
			obj = _transition_uniformity_objective(
				square_size=square_size,
				radial_size=radial_size,
				theta_inner_size=theta_inner_size,
				theta_outer_size=theta_outer_size,
				transition_size=transition_size,
				target_size=target_size,
			)

			if nr_hint is None:
				# Mild regularization toward not-too-small radial partitioning.
				obj += 0.05 / nr

			if best_obj is None or obj < best_obj:
				best_obj = obj
				best_rbase = rbase
				best_nr = nr
				best_square_size = square_size
				best_radial_size = radial_size

	if best_obj is None:
		raise ValueError(
			"Could not find a valid transition tuning solution. "
			"Adjust mesh_size/h/rbase/nr settings."
		)

	return best_rbase, best_nr, best_square_size, best_radial_size


def _parse_fields(raw: Dict[str, object]) -> FieldConfig:
	fields_cfg = raw.get("fields", {})
	if fields_cfg is None:
		fields_cfg = {}
	if not isinstance(fields_cfg, dict):
		raise ValueError("[fields] must be a table.")

	write_raw = _get_first(fields_cfg, ("write", "write_fields"))
	write = True if write_raw is None else bool(write_raw)

	n1 = _parse_positive_float(_get_first(fields_cfg, ("n1",)) or 1.4, "fields.n1")
	n2 = _parse_positive_float(_get_first(fields_cfg, ("n2",)) or 1.0, "fields.n2")
	amax = _parse_positive_float(
		_get_first(fields_cfg, ("amax", "alpha_max")) or 0.5, "fields.amax"
	)
	amin = _parse_positive_float(
		_get_first(fields_cfg, ("amin", "alpha_min")) or 0.02, "fields.amin"
	)
	profile_radius = _parse_positive_float(
		_get_first(fields_cfg, ("profile_radius", "radius")) or 1.0,
		"fields.profile_radius",
	)

	return FieldConfig(
		write=write,
		n1=n1,
		n2=n2,
		amax=amax,
		amin=amin,
		profile_radius=profile_radius,
	)


def _parse_mesh(config_path: Path, raw: Dict[str, object]) -> MeshConfig:
	mesh_cfg = raw.get("mesh")
	if not isinstance(mesh_cfg, dict):
		raise ValueError("Config must contain [mesh] table.")

	out_base_path = _resolve_output_base_path(mesh_cfg, config_path)
	append_hmax_raw = _get_first(mesh_cfg, ("append_hmax", "append_max_elem_size"))
	append_hmax = True if append_hmax_raw is None else bool(append_hmax_raw)

	rmax_raw = _get_first(mesh_cfg, ("rmax", "outer_radius"))
	length_raw = _get_first(mesh_cfg, ("length", "z_length"))
	if rmax_raw is None or length_raw is None:
		raise ValueError("[mesh].rmax and [mesh].length are required.")

	rmax = _parse_positive_float(rmax_raw, "mesh.rmax")
	length = _parse_positive_float(length_raw, "mesh.length")

	mesh_size_raw = _get_first(
		mesh_cfg,
		("mesh_size", "elem_size", "element_size", "target_mesh_size"),
	)
	mesh_size = (
		None
		if mesh_size_raw is None
		else _parse_positive_float(mesh_size_raw, "mesh.mesh_size")
	)

	h_raw = _get_first(mesh_cfg, ("h", "angular_resolution"))
	if h_raw is None:
		if mesh_size is None:
			raise ValueError("Provide either [mesh].h or [mesh].mesh_size.")
		h = _derive_h_from_mesh_size(rmax, mesh_size)
	else:
		h = _parse_positive_int(h_raw, "mesh.h")

	ntheta = 8 * h
	nxy = ntheta // 4
	circumferential_size = (2.0 * math.pi * rmax) / ntheta
	uniformity_target_size = mesh_size if mesh_size is not None else circumferential_size

	layout = str(_get_first(mesh_cfg, ("layout",)) or "spoke").strip().lower()
	if layout not in ("spoke", "legacy"):
		raise ValueError("mesh.layout must be 'spoke' or 'legacy'.")

	rbase_raw = _get_first(mesh_cfg, ("rbase", "inner_radius", "transition_radius"))
	rbase_value = None if rbase_raw is None else _parse_positive_float(rbase_raw, "mesh.rbase")
	if rbase_value is not None and rbase_value >= rmax:
		raise ValueError("mesh.rbase must be < mesh.rmax.")

	nr_raw = _get_first(mesh_cfg, ("nr", "n_radial"))
	nr_value = None if nr_raw is None else _parse_positive_int(nr_raw, "mesh.nr")

	tune_transition_raw = _get_first(
		mesh_cfg,
		("tune_transition", "tune_uniformity", "auto_tune_regions"),
	)
	if tune_transition_raw is None:
		tuning_enabled = rbase_value is None or nr_value is None
	else:
		tuning_enabled = bool(tune_transition_raw)

	core_bulge = 0.0
	if layout == "spoke":
		if rbase_value is not None:
			print(
				"warning: mesh.rbase is ignored with layout='spoke' "
				"(set layout='legacy' to reproduce old meshes).",
				file=sys.stderr,
			)
		tuning_enabled = False
		core_fraction_raw = _get_first(mesh_cfg, ("core_fraction",))
		core_fraction = (
			CORE_FRACTION_DEFAULT
			if core_fraction_raw is None
			else _parse_positive_float(core_fraction_raw, "mesh.core_fraction")
		)
		core_bulge_raw = _get_first(mesh_cfg, ("core_bulge",))
		core_bulge = CORE_BULGE_DEFAULT if core_bulge_raw is None else float(core_bulge_raw)
		if core_bulge < 0.0:
			raise ValueError("mesh.core_bulge must be >= 0.")
		if core_fraction * max(math.sqrt(2.0), 1.0 + core_bulge) >= 1.0:
			raise ValueError("mesh.core_fraction too large: the square core must fit inside rmax.")
		squaremax = core_fraction * rmax
		rbase = squaremax * (1.0 + core_bulge)  # core mid-side radius
		# Radial layers sized to the target along the mid-side spoke (longest spoke).
		nr = nr_value if nr_value is not None else max(1, int(round((rmax - rbase) / uniformity_target_size)))
		square_size = 2.0 * squaremax / nxy
		radial_size = (rmax - rbase) / (nr + 1)
		theta_inner_size = 2.0 * math.pi * rbase / ntheta
		transition_size = radial_size
	elif tuning_enabled:
		rbase, nr, square_size, radial_size = _tune_transition_region(
			rmax=rmax,
			nxy=nxy,
			ntheta=ntheta,
			target_size=uniformity_target_size,
			rbase_hint=rbase_value,
			nr_hint=nr_value,
		)
	else:
		rbase_default = rmax / 3.0 * 2.25
		rbase = rbase_default if rbase_value is None else rbase_value
		if nr_value is None:
			if mesh_size is not None:
				nr = max(1, int(round((rmax - rbase) / mesh_size)))
			else:
				nr = h // 2 + 1
		else:
			nr = nr_value

		radial_size = (rmax - rbase) / nr
		squaremax_tmp = _compute_squaremax(rmax=rmax, rbase=rbase, nr=nr)
		if squaremax_tmp <= 0.0:
			raise ValueError(
				"Invalid geometry: computed square core half-width is <= 0. "
				"Adjust rbase, rmax, or nr."
			)
		square_size = 2.0 * squaremax_tmp / nxy

	if layout == "legacy":
		squaremax = _compute_squaremax(rmax=rmax, rbase=rbase, nr=nr)
		if squaremax <= 0.0:
			raise ValueError(
				"Invalid geometry: computed square core half-width is <= 0. "
				"Adjust rbase, rmax, or nr."
			)
		theta_inner_size = 2.0 * math.pi * rbase / ntheta
		transition_size = rbase - math.sqrt(2.0) * squaremax
	theta_outer_size = circumferential_size

	length_axis = _parse_length_axis(
		_get_first(mesh_cfg, ("length_axis", "axis", "cylinder_axis"))
	)

	axis_start_raw = _get_first(mesh_cfg, ("axis_start", "z_start", "z0", "start"))
	axis_start = 0.5 * length if axis_start_raw is None else float(axis_start_raw)

	shift = _parse_vec3(_get_first(mesh_cfg, ("shift", "origin_shift", "offset")), "mesh.shift", (0.0, 0.0, 0.0))

	axial_mode_raw = _get_first(mesh_cfg, ("axial_mode", "z_mode"))
	if axial_mode_raw is None:
		axial_mode = "mesh_size" if mesh_size is not None else "legacy"
	else:
		axial_mode = str(axial_mode_raw).strip().lower()
	legacy_target_dz: float | None = None

	if axial_mode == "legacy":
		target_dz_raw = _get_first(mesh_cfg, ("target_dz", "legacy_target_dz"))
		if target_dz_raw is None:
			legacy_target_dz = rmax / (h + 1 + nr) * 1.5
		else:
			legacy_target_dz = _parse_positive_float(target_dz_raw, "mesh.target_dz")
		blocks = int(length // legacy_target_dz)
		if blocks <= 0:
			raise ValueError(
				"Legacy axial mode computed zero blocks. Increase mesh.length or decrease target_dz."
			)
		dz = length / blocks
		n_layers = blocks + 1
	elif axial_mode == "fixed_layers":
		n_layers_raw = _get_first(mesh_cfg, ("n_layers", "layers", "nz"))
		if n_layers_raw is None:
			raise ValueError("mesh.n_layers is required for axial_mode='fixed_layers'.")
		n_layers = _parse_positive_int(n_layers_raw, "mesh.n_layers")
		dz = length / n_layers
	elif axial_mode == "fixed_dz":
		dz_raw = _get_first(mesh_cfg, ("dz", "axial_step"))
		if dz_raw is None:
			raise ValueError("mesh.dz is required for axial_mode='fixed_dz'.")
		dz_in = _parse_positive_float(dz_raw, "mesh.dz")
		n_layers = max(1, int(math.ceil(length / dz_in)))
		dz = length / n_layers
	elif axial_mode == "mesh_size":
		if mesh_size is None:
			raise ValueError("mesh.mesh_size is required for axial_mode='mesh_size'.")
		n_layers = max(1, int(round(length / mesh_size)))
		dz = length / n_layers
	else:
		raise ValueError(
			"mesh.axial_mode must be one of: 'legacy', 'fixed_layers', 'fixed_dz', 'mesh_size'."
		)

	return MeshConfig(
		out_base_path=out_base_path,
		append_hmax=append_hmax,
		layout=layout,
		core_bulge=core_bulge,
		rmax=rmax,
		mesh_size=mesh_size,
		uniformity_target_size=uniformity_target_size,
		tuning_enabled=tuning_enabled,
		rbase=rbase,
		h=h,
		nr=nr,
		ntheta=ntheta,
		nxy=nxy,
		circumferential_size=circumferential_size,
		theta_inner_size=theta_inner_size,
		theta_outer_size=theta_outer_size,
		radial_size=radial_size,
		square_size=square_size,
		transition_size=transition_size,
		squaremax=squaremax,
		length=length,
		length_axis=length_axis,
		axis_start=axis_start,
		shift=shift,
		axial_mode=axial_mode,
		n_layers=n_layers,
		dz=dz,
		legacy_target_dz=legacy_target_dz,
	)


def _map_local_to_global(x_local: float, y_local: float, axis_local: float, axis: str) -> Vec3f:
	if axis == "z":
		return (x_local, y_local, axis_local)
	if axis == "x":
		return (axis_local, x_local, y_local)
	# axis == "y"
	return (y_local, axis_local, x_local)


def _radial_distance_for_axis(x: float, y: float, z: float, axis: str) -> float:
	if axis == "z":
		return math.hypot(x, y)
	if axis == "x":
		return math.hypot(y, z)
	return math.hypot(x, z)


def _build_ring_points(
	rmax: float,
	rbase: float,
	nr: int,
	ntheta: int,
	squaremax: float,
) -> List[Vec2f]:
	points: List[Vec2f] = []
	dr = (rmax - rbase) / nr
	quarter = ntheta // 4
	half_quarter = quarter // 2

	for i in range(nr + 1):
		r = rbase + i * (rmax - rbase) / nr
		fac = i / nr
		for j in range(ntheta):
			theta = j / ntheta * 2.0 * math.pi
			is_corner = (j % quarter) == half_quarter

			theta_shift = theta + math.pi / 4.0
			if theta_shift < 0.0:
				theta_shift += 2.0 * math.pi
			if theta_shift > 2.0 * math.pi:
				theta_shift -= 2.0 * math.pi

			if is_corner:
				x_morph = (r - dr) * math.cos(theta)
				y_morph = (r - dr) * math.sin(theta)
			elif theta_shift < math.pi / 2.0:
				x_morph = r / math.sqrt(2.0)
				y0 = r * math.sin(theta)
				y1 = -squaremax + theta_shift * 2.0 / math.pi * squaremax * 2.0
				y_morph = 0.5 * y0 + 0.5 * y1
			elif theta_shift < math.pi:
				x0 = r * math.cos(theta)
				x1 = squaremax - (theta_shift - math.pi / 2.0) * 2.0 / math.pi * squaremax * 2.0
				x_morph = 0.5 * x0 + 0.5 * x1
				y_morph = r / math.sqrt(2.0)
			elif theta_shift < 1.5 * math.pi:
				x_morph = -r / math.sqrt(2.0)
				y0 = r * math.sin(theta)
				y1 = squaremax - (theta_shift - math.pi) * 2.0 / math.pi * squaremax * 2.0
				y_morph = 0.5 * y0 + 0.5 * y1
			else:
				x0 = r * math.cos(theta)
				x1 = -squaremax + (theta_shift - 1.5 * math.pi) * 2.0 / math.pi * squaremax * 2.0
				x_morph = 0.5 * x0 + 0.5 * x1
				y_morph = -r / math.sqrt(2.0)

			xc = r * math.cos(theta)
			yc = r * math.sin(theta)
			x = xc * fac + x_morph * (1.0 - fac)
			y = yc * fac + y_morph * (1.0 - fac)
			points.append((x, y))

	return points


def _build_spoke_ring_points(boundary_pts: Sequence[Vec2f], rmax: float, nr: int) -> List[Vec2f]:
	"""Rings on straight spokes from each core-boundary node to its circle node.

	Ring k (k = 0..nr) sits at t = (k+1)/(nr+1) along the spoke
	B_j -> O_j = rmax*(cos, sin)(2*pi*j/ntheta), so the inner rings follow the
	core's node spacing and the outer ring is the circle. boundary_pts must start
	at theta = 0 and run counter-clockwise (the _build_boundary_order order).
	"""
	ntheta = len(boundary_pts)
	outer = [
		(rmax * math.cos(2.0 * math.pi * j / ntheta), rmax * math.sin(2.0 * math.pi * j / ntheta))
		for j in range(ntheta)
	]
	points: List[Vec2f] = []
	for k in range(nr + 1):
		t = (k + 1) / (nr + 1)
		for (bx, by), (ox, oy) in zip(boundary_pts, outer):
			points.append((bx + t * (ox - bx), by + t * (oy - by)))
	return points


def _build_square_points(squaremax: float, nxy: int, bulge: float = 0.0) -> List[Vec2f]:
	"""Square core grid; bulge > 0 bows each edge outward by bulge*squaremax at mid-side.

	A parabolic bow x = u*s*(1 + b*(1 - v^2)) opens the core corner angle to
	90 + 2*atan(2b) degrees, which relieves the 135-degree angles of the two ring
	cells meeting at each corner.
	"""
	points: List[Vec2f] = []
	for i in range(nxy + 1):
		u = -1.0 + 2.0 * i / nxy
		for j in range(nxy + 1):
			v = -1.0 + 2.0 * j / nxy
			x = u * squaremax * (1.0 + bulge * (1.0 - v * v))
			y = v * squaremax * (1.0 + bulge * (1.0 - u * u))
			points.append((x, y))
	return points


def _build_boundary_order(nxy: int) -> List[int]:
	rnodes = [(nxy + 1) * nxy + i for i in range(nxy + 1)]
	lnodes = [nxy - i for i in range(nxy + 1)]
	tnodes = [rnodes[-1] - i * (nxy + 1) for i in range(1, nxy)]
	bnodes = [lnodes[-1] + i * (nxy + 1) for i in range(1, nxy)]

	rnodes0 = rnodes[len(rnodes) // 2 :]
	rnodes1 = rnodes[: len(rnodes) // 2]

	out = rnodes0 + tnodes + lnodes + bnodes + rnodes1
	return out


def _square_id(i: int, j: int, nxy: int) -> int:
	return j + (nxy + 1) * i


def _ring_id(ring_i: int, theta_j: int, ntheta: int, n_square_nodes: int) -> int:
	return n_square_nodes + ring_i * ntheta + theta_j


def _build_cross_section_quads(
	nxy: int,
	nr: int,
	ntheta: int,
	boundary_order: Sequence[int],
	n_square_nodes: int,
) -> List[Quad]:
	quads: List[Quad] = []

	# Square core quads.
	for i in range(nxy):
		for j in range(nxy):
			n0 = _square_id(i, j, nxy)
			n1 = _square_id(i + 1, j, nxy)
			n2 = _square_id(i + 1, j + 1, nxy)
			n3 = _square_id(i, j + 1, nxy)
			quads.append((n0, n1, n2, n3))

	if len(boundary_order) != ntheta:
		raise ValueError(
			"Internal error: boundary node count does not match angular discretization."
		)

	# Connecting ring between square boundary and first cylindrical ring.
	for j in range(ntheta):
		jn = (j + 1) % ntheta
		s0 = boundary_order[j]
		s1 = boundary_order[jn]
		r0 = _ring_id(0, j, ntheta, n_square_nodes)
		r1 = _ring_id(0, jn, ntheta, n_square_nodes)
		quads.append((s0, r0, r1, s1))

	# Radial annulus quads.
	for i in range(nr):
		for j in range(ntheta):
			jn = (j + 1) % ntheta
			ri0 = _ring_id(i, j, ntheta, n_square_nodes)
			ri1 = _ring_id(i, jn, ntheta, n_square_nodes)
			ro1 = _ring_id(i + 1, jn, ntheta, n_square_nodes)
			ro0 = _ring_id(i + 1, j, ntheta, n_square_nodes)
			quads.append((ri0, ro0, ro1, ri1))

	return quads


def _build_z_levels(z_start: float, dz: float, n_layers: int) -> List[float]:
	return [z_start - i * dz for i in range(n_layers + 1)]


def _compute_hmax(points2d: Sequence[Vec2f], quads: Sequence[Quad], dz: float) -> float:
	max_sq = 0.0
	for q in quads:
		edges = ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0]))
		for i0, i1 in edges:
			dx = points2d[i1][0] - points2d[i0][0]
			dy = points2d[i1][1] - points2d[i0][1]
			d2 = dx * dx + dy * dy
			if d2 > max_sq:
				max_sq = d2
	return max(math.sqrt(max_sq), abs(dz))


def _cross_section_quality(
	points2d: Sequence[Vec2f], quads: Sequence[Quad], target_size: float
) -> Dict[str, float]:
	"""Measured quad quality: scaled Jacobian, aspect ratio and size (sqrt(area)) vs target."""
	sj_min = 1.0
	sj_sum = 0.0
	ar_max = 1.0
	sizes: List[float] = []
	for q in quads:
		pts = [points2d[i] for i in q]
		lengths = []
		corner_sj = []
		area2 = 0.0
		for k in range(4):
			px, py = pts[k]
			nx, ny = pts[(k + 1) % 4]
			bx, by = pts[k - 1]
			ex, ey = nx - px, ny - py
			fx, fy = bx - px, by - py
			le, lf = math.hypot(ex, ey), math.hypot(fx, fy)
			corner_sj.append((ex * fy - ey * fx) / (le * lf))
			lengths.append(le)
			area2 += px * ny - nx * py
		# The builder winds all quads the same way; take the first quad's winding as
		# reference so an inverted quad shows up as negative.
		if not sizes:
			winding = 1.0 if area2 > 0.0 else -1.0
		quad_sj = min(winding * s for s in corner_sj)
		sj_min = min(sj_min, quad_sj)
		sj_sum += quad_sj
		ar_max = max(ar_max, max(lengths) / min(lengths))
		sizes.append(math.sqrt(0.5 * abs(area2)))
	mean = sum(sizes) / len(sizes)
	cv = math.sqrt(sum((s - mean) ** 2 for s in sizes) / len(sizes)) / mean
	return {
		"sj_min": sj_min,
		"sj_mean": sj_sum / len(quads),
		"aspect_max": ar_max,
		"size_min": min(sizes) / target_size,
		"size_mean": mean / target_size,
		"size_max": max(sizes) / target_size,
		"size_cv": cv,
	}


def _hex_sign(p0: Vec3f, p1: Vec3f, p3: Vec3f, p4: Vec3f) -> float:
	v1 = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
	v2 = (p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2])
	v3 = (p4[0] - p0[0], p4[1] - p0[1], p4[2] - p0[2])
	cross = (
		v2[1] * v3[2] - v2[2] * v3[1],
		v2[2] * v3[0] - v2[0] * v3[2],
		v2[0] * v3[1] - v2[1] * v3[0],
	)
	return v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2]


def _select_hex_order(points2d: Sequence[Vec2f], quads: Sequence[Quad], z_levels: Sequence[float]) -> str:
	q = quads[0]
	z_top = z_levels[0]
	z_bot = z_levels[1]

	t0 = (points2d[q[0]][0], points2d[q[0]][1], z_top)
	t1 = (points2d[q[1]][0], points2d[q[1]][1], z_top)
	t3 = (points2d[q[3]][0], points2d[q[3]][1], z_top)
	b0 = (points2d[q[0]][0], points2d[q[0]][1], z_bot)

	sign_top_bottom = _hex_sign(t0, t1, t3, b0)
	if sign_top_bottom > 0.0:
		return "top_bottom"
	return "bottom_top"


def _legacy_profile_value(r: float, high: float, low: float, profile_radius: float) -> float:
	delta = 0.5 * (1.0 - (low / high) ** 2)
	value_under_sqrt = high * high * (1.0 - 2.0 * delta * (r / profile_radius) ** 2)
	low_sq = low * low
	if value_under_sqrt < low_sq:
		return low
	return math.sqrt(value_under_sqrt)


def _write_vtk(
	out_path: Path,
	points2d: Sequence[Vec2f],
	quads: Sequence[Quad],
	z_levels: Sequence[float],
	length_axis: str,
	shift: Vec3f,
	fields: FieldConfig,
) -> None:
	n_slice = len(points2d)
	n_layers = len(z_levels) - 1
	n_points = n_slice * len(z_levels)
	n_cells = len(quads) * n_layers

	if n_cells <= 0:
		raise ValueError("No elements generated; check mesh resolution settings.")

	hex_order = _select_hex_order(points2d, quads, z_levels)

	out_path.parent.mkdir(parents=True, exist_ok=True)
	with out_path.open("w", encoding="utf-8") as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write("MIRAGE cylindrical HEX8 mesh\n")
		f.write("ASCII\n")
		f.write("DATASET UNSTRUCTURED_GRID\n")

		f.write(f"POINTS {n_points} float\n")
		sx, sy, sz = shift
		for axis_value in z_levels:
			for x, y in points2d:
				gx, gy, gz = _map_local_to_global(x, y, axis_value, length_axis)
				f.write(f"{gx + sx:.16g} {gy + sy:.16g} {gz + sz:.16g}\n")

		f.write(f"CELLS {n_cells} {n_cells * 9}\n")
		for k in range(n_layers):
			top_off = k * n_slice
			bot_off = (k + 1) * n_slice
			for q0, q1, q2, q3 in quads:
				if hex_order == "top_bottom":
					ids = (
						top_off + q0,
						top_off + q1,
						top_off + q2,
						top_off + q3,
						bot_off + q0,
						bot_off + q1,
						bot_off + q2,
						bot_off + q3,
					)
				else:
					ids = (
						bot_off + q0,
						bot_off + q1,
						bot_off + q2,
						bot_off + q3,
						top_off + q0,
						top_off + q1,
						top_off + q2,
						top_off + q3,
					)
				f.write("8 " + " ".join(str(v) for v in ids) + "\n")

		f.write(f"CELL_TYPES {n_cells}\n")
		for _ in range(n_cells):
			f.write("12\n")  # VTK_HEXAHEDRON

		if fields.write:
			f.write(f"POINT_DATA {n_points}\n")
			f.write("SCALARS RefractiveIndex float 1\n")
			f.write("LOOKUP_TABLE default\n")
			for axis_value in z_levels:
				for x, y in points2d:
					gx, gy, gz = _map_local_to_global(x, y, axis_value, length_axis)
					gx += sx
					gy += sy
					gz += sz
					r = _radial_distance_for_axis(gx, gy, gz, length_axis)
					n_value = _legacy_profile_value(
						r=r,
						high=fields.n1,
						low=fields.n2,
						profile_radius=fields.profile_radius,
					)
					f.write(f"{n_value:.16g}\n")

			f.write("SCALARS Absorptivity float 1\n")
			f.write("LOOKUP_TABLE default\n")
			for axis_value in z_levels:
				for x, y in points2d:
					gx, gy, gz = _map_local_to_global(x, y, axis_value, length_axis)
					gx += sx
					gy += sy
					gz += sz
					r = _radial_distance_for_axis(gx, gy, gz, length_axis)
					a_value = _legacy_profile_value(
						r=r,
						high=fields.amax,
						low=fields.amin,
						profile_radius=fields.profile_radius,
					)
					f.write(f"{a_value:.16g}\n")


def _build_cross_section(mesh_cfg: MeshConfig) -> Tuple[List[Vec2f], List[Quad]]:
	square_points = _build_square_points(mesh_cfg.squaremax, mesh_cfg.nxy, mesh_cfg.core_bulge)
	boundary_order = _build_boundary_order(mesh_cfg.nxy)
	if mesh_cfg.layout == "spoke":
		ring_points = _build_spoke_ring_points(
			[square_points[i] for i in boundary_order], mesh_cfg.rmax, mesh_cfg.nr
		)
	else:
		ring_points = _build_ring_points(
			rmax=mesh_cfg.rmax,
			rbase=mesh_cfg.rbase,
			nr=mesh_cfg.nr,
			ntheta=mesh_cfg.ntheta,
			squaremax=mesh_cfg.squaremax,
		)
	quads = _build_cross_section_quads(
		nxy=mesh_cfg.nxy,
		nr=mesh_cfg.nr,
		ntheta=mesh_cfg.ntheta,
		boundary_order=boundary_order,
		n_square_nodes=len(square_points),
	)
	return square_points + ring_points, quads


def _format_vec3(v: Sequence[float]) -> str:
	return f"({v[0]:.12g}, {v[1]:.12g}, {v[2]:.12g})"


def _load_config(config_path: Path) -> Tuple[MeshConfig, FieldConfig]:
	try:
		raw = tomllib.loads(config_path.read_text(encoding="utf-8"))
	except FileNotFoundError as exc:
		raise SystemExit(f"Config file not found: {config_path}") from exc
	except Exception as exc:
		raise SystemExit(f"Failed to parse TOML config '{config_path}': {exc}") from exc

	if not isinstance(raw, dict):
		raise ValueError("Top-level TOML content must be a table.")

	mesh = _parse_mesh(config_path, raw)
	fields = _parse_fields(raw)
	return mesh, fields


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"Generate a legacy-compatible cylindrical VTK mesh with HEX8 elements "
			"from a TOML config file."
		)
	)
	parser.add_argument(
		"--config",
		required=True,
		type=Path,
		help="Path to TOML config file.",
	)
	args = parser.parse_args()

	try:
		mesh_cfg, fields_cfg = _load_config(args.config)
	except ValueError as exc:
		raise SystemExit(f"Invalid config: {exc}") from exc

	points2d, quads = _build_cross_section(mesh_cfg)
	z_levels = _build_z_levels(mesh_cfg.axis_start, mesh_cfg.dz, mesh_cfg.n_layers)

	hmax = _compute_hmax(points2d, quads, mesh_cfg.dz)
	quality = _cross_section_quality(points2d, quads, mesh_cfg.uniformity_target_size)
	out_path = _append_hmax_to_filename(mesh_cfg.out_base_path, hmax) if mesh_cfg.append_hmax else mesh_cfg.out_base_path

	_write_vtk(
		out_path=out_path,
		points2d=points2d,
		quads=quads,
		z_levels=z_levels,
		length_axis=mesh_cfg.length_axis,
		shift=mesh_cfg.shift,
		fields=fields_cfg,
	)

	n_slice = len(points2d)
	n_points = n_slice * len(z_levels)
	n_cells = len(quads) * mesh_cfg.n_layers
	print("Generated cylindrical HEX8 VTK mesh")
	print(f"  output            : {out_path}")
	print(f"  layout            : {mesh_cfg.layout}")
	print(f"  rmax              : {mesh_cfg.rmax:.12g}")
	if mesh_cfg.layout == "spoke":
		print(f"  core_fraction     : {mesh_cfg.squaremax / mesh_cfg.rmax:.6g}")
		print(f"  core_bulge        : {mesh_cfg.core_bulge:.6g}")
	else:
		print(f"  rbase             : {mesh_cfg.rbase:.12g}")
	if mesh_cfg.mesh_size is not None:
		print(f"  mesh_size target  : {mesh_cfg.mesh_size:.12g}")
	print(f"  h                 : {mesh_cfg.h}")
	print(f"  nr                : {mesh_cfg.nr}")
	print(f"  ntheta            : {mesh_cfg.ntheta}")
	print(f"  nxy               : {mesh_cfg.nxy}")
	print(f"  uniformity target : {mesh_cfg.uniformity_target_size:.12g}")
	print(f"  transition tuned  : {mesh_cfg.tuning_enabled}")
	print(f"  squaremax         : {mesh_cfg.squaremax:.12g}")
	print(f"  quality (cross-section, size relative to {mesh_cfg.uniformity_target_size:.6g}):")
	print(f"    scaled Jacobian : min {quality['sj_min']:.3f}  mean {quality['sj_mean']:.3f}  (1 = square)")
	print(f"    aspect ratio max: {quality['aspect_max']:.3f}")
	print(
		f"    cell size       : min {quality['size_min']:.3f}  mean {quality['size_mean']:.3f}"
		f"  max {quality['size_max']:.3f}  CV {quality['size_cv']:.3f}"
	)
	print(f"  length_axis       : {mesh_cfg.length_axis}")
	print(f"  axial_mode        : {mesh_cfg.axial_mode}")
	if mesh_cfg.legacy_target_dz is not None:
		print(f"  legacy_target_dz  : {mesh_cfg.legacy_target_dz:.12g}")
	print(f"  axis_start        : {mesh_cfg.axis_start:.12g}")
	print(f"  n_layers          : {mesh_cfg.n_layers}")
	print(f"  dz                : {mesh_cfg.dz:.12g}")
	print(f"  hmax              : {hmax:.12g}")
	print(f"  shift             : {_format_vec3(mesh_cfg.shift)}")
	print(f"  points per slice  : {n_slice}")
	print(f"  cells per layer   : {len(quads)}")
	print(f"  points total      : {n_points}")
	print(f"  hex elements total: {n_cells}")
	print(f"  point fields      : {'on' if fields_cfg.write else 'off'}")


if __name__ == "__main__":
	main()
