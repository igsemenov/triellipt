# -*- coding: utf-8 -*-
"""Tests mesh skeleton.

MESH:

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
from triellipt.fem import joints


POINTS = [
    0.0 + 0.0j,
    0.0 + 1.0j,
    1.0 + 0.0j,
    1.0 + 1.0j,
    0.5 + 0.5j,
    0.5 + 0.0j,
    1.0 + 0.5j
]

TRIANGS = [
    [0, 3, 1],
    [6, 5, 2],
    [3, 0, 4],
    [0, 5, 4],
    [3, 4, 6],
    [6, 4, 5]
]


class TestSkeleton(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.getmesh()
        cls.SKEL = cls.getskel()

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array(TRIANGS)
        )

    @classmethod
    def getskel(cls):
        return joints.getskeleton(cls.MESH)

    def test_mesh(self):
        assert self.MESH.getvoids().tolist() == [2]

    def test_body(self):
        assert self.SKEL.body.triangs.tolist() == [
            [0, 3, 1], [6, 5, 2]
        ]

    def test_joints(self):
        assert self.SKEL.wests.triangs.tolist() == [[4, 0, 5]]
        assert self.SKEL.easts.triangs.tolist() == [[4, 6, 3]]
        assert self.SKEL.cores.triangs.tolist() == [[4, 5, 6]]
        assert self.SKEL.voids.triangs.tolist() == [[3, 0, 4]]


if __name__ == '__main__':
    unittest.main()
