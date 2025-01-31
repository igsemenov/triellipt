# -*- coding: utf-8 -*-
"""Tests super voids (only ears).

TEST-MESH:

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
from triellipt.fem import supvoids


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
        cls.SUPVOIDS = cls.getsupvoids()

    @classmethod
    def getsupvoids(cls):
        return supvoids.get_supvoids(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array(TRIANGS)
        )

    def test_trinums_voids(self):
        assert self.SUPVOIDS.trinums.tolist() == [2]

    def test_trinums_ears(self):
        assert self.SUPVOIDS.trinums2.tolist() == [3]
        assert self.SUPVOIDS.trinums3.tolist() == [4]

    def test_voids_triangs(self):
        assert self.SUPVOIDS.voids_triangs.tolist() == [[3, 0, 4]]


if __name__ == '__main__':
    unittest.main()
