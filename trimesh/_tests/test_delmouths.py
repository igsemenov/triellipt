# -*- coding: utf-8 -*-
"""Tests removing mesh mouths.

TEST-MESH:

      3---7---4
     / \ / \ / \
    /   5---6   \
   /     \ /     \
  0-------1-------2

"""
import unittest
import numpy as np
from triellipt import trimesh
from triellipt.trimesh import delmouths_

POINTS = np.r_[
    0, 1, 2, 0.5 + 1j, 1.5 + 1j, 0.75 + 0.5j, 1.25 + 0.5j, 1 + 1.01j
]

TRIANGS_BASES = np.array(
    [[0, 1, 3], [1, 2, 4]]
)

TRIANGS_CORES = np.array(
    [[1, 6, 5], [6, 4, 7], [5, 7, 3], [5, 6, 7]]
)

TRIANGS_VOIDS = np.array([[1, 4, 6], [3, 1, 5]])


class TestMesh:

    @classmethod
    def getmesh(cls):

        triangs = np.vstack(
            [TRIANGS_BASES, TRIANGS_CORES, TRIANGS_VOIDS]
        )

        return trimesh.TriMesh.from_data(POINTS, triangs)


class TestCleanMesh(TestMesh, unittest.TestCase):
    """Tests the resulting clean mesh.
    """

    @classmethod
    def setUpClass(cls):
        cls.MESH = delmouths_.remove_mouths(cls.getmesh())

    def test_mesh_size(self):
        assert self.MESH.size == 4
        assert self.MESH.npoints == 6
        assert self.MESH.getvoids().size == 1

    def test_mesh_triangs(self):
        assert [0, 1, 3] in self.MESH.triangs
        assert [1, 2, 4] in self.MESH.triangs
        assert [3, 1, 4] in self.MESH.triangs
        assert [3, 4, 7] in self.MESH.triangs

    def test_mesh_points(self):
        assert self.MESH.points[-1] == 1 + 1j


class TestCleaner(TestMesh, unittest.TestCase):
    """Tests the cleaner production line.
    """

    @classmethod
    def setUpClass(cls):
        cls.CLEANER = delmouths_.MouthsCleaner(cls.getmesh())

    def test_supmouhts(self):
        assert self.supmouths.supmesh.triangs.tolist() == [[1, 4, 3]]
        assert self.supmouths.kermesh.triangs.tolist() == [[5, 6, 7]]

    def test_mesh_alpha(self):
        assert self.mesh_alpha.size == 8
        assert self.mesh_alpha.npoints == 8
        assert self.mesh_alpha.getvoids().size == 5

    def test_mesh_beta(self):
        assert self.mesh_beta.size == 4
        assert self.mesh_beta.npoints == 8
        assert self.mesh_beta.getvoids().size == 1

    @property
    def supmouths(self):
        return self.CLEANER.find_supmouths()

    @property
    def mesh_alpha(self):
        return self.CLEANER.make_mesh_alpha()

    @property
    def mesh_beta(self):
        return self.CLEANER.make_mesh_beta()


if __name__ == '__main__':
    unittest.main()
