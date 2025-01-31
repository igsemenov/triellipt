# -*- coding: utf-8 -*-
"""Tests the mesh skeleton.

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
from triellipt.fem import skeleton


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


class TestSkeletonNoVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.SKEL = skeleton.getskeleton(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array([[0, 3, 1], [0, 2, 3]])
        )

    def test_voids(self):
        assert self.SKEL.voidsmap is None
        assert self.SKEL.hasvoids is False

    def test_nodesmap(self):
        assert self.SKEL.nodesmap.nodnums.tolist() == [0, 0, 1, 2, 3, 3]
        assert self.SKEL.nodesmap.trinums.tolist() == [0, 1, 0, 1, 0, 1]
        assert self.SKEL.nodesmap.locnums.tolist() == [0, 0, 2, 1, 1, 2]

    def test_nodes_ranges(self):
        assert self.SKEL.nodesmap.nodes_range.tolist() == [0, 1, 2, 3]


class TestSkeletonWithVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.SKEL = skeleton.getskeleton(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array(TRIANGS)
        )

    def test_voids(self):
        assert self.SKEL.hasvoids is True
        assert self.SKEL.voids_trinums.tolist() == [2]
        assert self.SKEL.voids_triangs.tolist() == [[3, 0, 4]]

    def test_nodesmap(self):
        assert self.SKEL.nodesmap.atnode(0).nodnums.tolist() == [0, 0, 0]
        assert self.SKEL.nodesmap.atnode(0).trinums.tolist() == [0, 2, 3]
        assert self.SKEL.nodesmap.atnode(0).locnums.tolist() == [0, 1, 0]

    def test_voidsmap(self):
        assert self.SKEL.voidsmap.nodnums.tolist() == [4, 4, 4, 4]
        assert self.SKEL.voidsmap.trinums.tolist() == [2, 3, 5, 4]
        assert self.SKEL.voidsmap.locnums.tolist() == [2, 2, 1, 1]

    def test_nodes_ranges(self):
        assert self.SKEL.nodesmap.nodes_range.tolist() == [0, 1, 2, 3, 5, 6]
        assert self.SKEL.voidsmap.nodes_range.tolist() == [4]

    def test_maps_ranks(self):
        assert self.SKEL.nodesmap.nodes_ranks.tolist() == [3, 1, 1, 3, 3, 3]
        assert self.SKEL.voidsmap.nodes_ranks.tolist() == [4]

    def test_maps_sizes(self):
        assert self.SKEL.voidsmap.size == 4
        assert self.SKEL.nodesmap.size == 14


if __name__ == '__main__':
    unittest.main()
