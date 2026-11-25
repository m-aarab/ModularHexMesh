"""Tests for cylindricalMesh.py. Run: python -m unittest test_cylindricalMesh"""
import hashlib
import math
import unittest
from pathlib import Path

import cylindricalMesh as cm


def _cfg(**mesh):
	return cm._parse_mesh(Path("x.toml"), {"mesh": dict(length=10.0, **mesh)})


def _fingerprint(points, quads):
	s = "".join(f"{x:.9f},{y:.9f};" for x, y in points) + repr(quads)
	return hashlib.sha256(s.encode()).hexdigest()


class LegacyLayout(unittest.TestCase):
	# Fingerprints of the cross-sections written by the generator before the spoke layout (1948ec4).
	CASES = [
		(dict(rmax=300.0, mesh_size=20.0, tune_transition=True),
		 "c47e65a44f5796833c7c31a40efa7e47eee6f927648b9fba5f4a48a067d7925d"),
		(dict(rmax=150.0, mesh_size=1.5, rbase=120.0, tune_transition=True),
		 "e0ca3db84737e906f80dc83a532824372f1732c04f365074d6ab4f793dc4a3d2"),
	]

	def test_reproduces_previous_meshes(self):
		for mesh, expected in self.CASES:
			with self.subTest(**mesh):
				c = _cfg(layout="legacy", **mesh)
				self.assertEqual(_fingerprint(*cm._build_cross_section(c)), expected)


class SpokeLayout(unittest.TestCase):
	def test_quality_on_tutorial_vial_radii(self):
		# Quality is resolution independent; mesh_size 3 keeps the test fast.
		for rmax in (75.0, 140.0, 150.0, 200.0):
			with self.subTest(rmax=rmax):
				c = _cfg(rmax=rmax, mesh_size=3.0)
				points, quads = cm._build_cross_section(c)
				q = cm._cross_section_quality(points, quads, c.uniformity_target_size)
				self.assertGreater(q["sj_min"], 0.85)
				self.assertLess(q["aspect_max"], 2.2)
				self.assertLess(q["size_max"] / q["size_min"], 1.75)
				self.assertLess(q["size_cv"], 0.10)

	def test_outer_ring_on_circle(self):
		c = _cfg(rmax=150.0, mesh_size=5.0)
		points, _ = cm._build_cross_section(c)
		for x, y in points[-c.ntheta:]:
			self.assertAlmostEqual(math.hypot(x, y), 150.0, places=9)

	def test_core_too_large_rejected(self):
		with self.assertRaises(ValueError):
			_cfg(rmax=150.0, mesh_size=5.0, core_fraction=0.72)


class SquareCore(unittest.TestCase):
	def test_zero_bulge_is_plain_grid(self):
		pts = cm._build_square_points(2.0, 4, bulge=0.0)
		expected = [(-2.0 + i, -2.0 + j) for i in range(5) for j in range(5)]
		for (x, y), (ex, ey) in zip(pts, expected):
			self.assertAlmostEqual(x, ex, places=12)
			self.assertAlmostEqual(y, ey, places=12)

	def test_default_bulge_gives_120_degree_corner(self):
		b = cm.CORE_BULGE_DEFAULT
		self.assertAlmostEqual(90.0 + 2.0 * math.degrees(math.atan(2.0 * b)), 120.0, places=9)


if __name__ == "__main__":
	unittest.main()
