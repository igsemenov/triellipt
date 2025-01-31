# -*- coding: utf-8 -*-
"""Tests the mesh refinement.

TEST-MESH-NO-VOIDS:

    1-----3
    |    /|
    |   / |
    |  /  |
    | /   |
    |/    |
    0-----2

TEST-MESH-WITH-VOIDS:

    1-----3
    |    /|
    |   / |
    |  4--6
    | /| /|
    |/ |/ |
    0--5--2

"""
import unittest
import numpy as np
from triellipt import trimesh
from triellipt.amr import trirefine


POINTS = [
    0.0 + 0.0j,
    0.0 + 1.0j,
    1.0 + 0.0j,
    1.0 + 1.0j,
    0.5 + 0.5j,
    0.5 + 0.0j,
    1.0 + 0.5j
]

TRIANGS = {
    0: [0, 3, 1],
    1: [6, 5, 2],
    2: [3, 0, 4],
    3: [0, 5, 4],
    4: [3, 4, 6],
    5: [6, 4, 5]
}


class TestNoVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.REF = trirefine.MeshRefiner.from_mesh(cls.getmesh())

    @classmethod
    def getmesh(cls):

        points = np.array(POINTS[0:4])

        triangs = np.array(
            [[0, 3, 1], [0, 2, 3]]
        )

        return trimesh.TriMesh.from_data(points, triangs)

    def test_trinums_filtered(self):
        assert self.ref.target_trinums.tolist() == [0, 1]

    def test_make_core_mesh(self):
        assert self.ref.make_core_mesh().size == 2
        assert self.ref.make_core_mesh().npoints == 5

    def test_make_alpha_mesh(self):
        assert self.ref.make_mesh_alpha().size == 12
        assert self.ref.make_mesh_alpha().npoints == 9

    def test_make_beta_mesh(self):
        assert self.ref.maker_mesh_beta.find_voids_on_edge().size == 4
        assert self.ref.maker_mesh_beta.find_voids_twins() == {}

    @property
    def ref(self):
        return self.REF.with_trinums(range(-10, 10))


class TestWithVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.REF = trirefine.MeshRefiner.from_mesh(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.vstack(list(TRIANGS.values()))
        )

    def test_trinums_filtered(self):
        assert self.ref.target_trinums.tolist() == [0, 1, 5]

    def test_make_core_mesh(self):
        assert self.ref.make_core_mesh().size == 3
        assert self.ref.make_core_mesh().npoints == 8

    def test_make_alpha_mesh(self):
        assert self.ref.make_mesh_alpha().size == 22
        assert self.ref.make_mesh_alpha().npoints == 15

    def test_make_beta_mesh(self):
        assert self.ref.maker_mesh_beta.find_voids_on_edge().size == 4
        assert self.ref.maker_mesh_beta.find_voids_twins()['nums'].size == 2

    @property
    def ref(self):
        """Tested refiner.
        """
        return self.REF.with_trinums(range(-10, 10))


if __name__ == '__main__':
    unittest.main()
