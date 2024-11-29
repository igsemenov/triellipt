# -*- coding: utf-8 -*-
"""Tests the index stream (FVM).

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
from triellipt.fem import ijstream_fvm as ijstream


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


class TestMatrixSreamNoVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.IJS = ijstream.matrix_stream(cls.getskel())

    @classmethod
    def getskel(cls):
        return skeleton.getskeleton(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array([[0, 3, 1], [0, 2, 3]])
        )

    def test_meta(self):
        assert self.IJS.meta['hasvoids'] is False
        assert self.IJS.meta['nodsmap-size'] == 6

    def test_size(self):
        assert self.IJS.size == 6 * 3

    def test_trinums(self):
        assert self.IJS.trinums.tolist() == [
            *self.trinums, *self.trinums, *self.trinums
        ]

    def test_rownums(self):
        assert self.IJS.rownums.tolist() == [
            *self.nodnums0, *self.nodnums0, *self.nodnums0
        ]

    def test_colnums(self):
        assert self.IJS.colnums.tolist() == [
            *self.nodnums0, *self.nodnums1, *self.nodnums2
        ]

    @property
    def trinums(self):
        """Trinums from the nodesmap.
        """
        return [0, 1, 0, 1, 0, 1]

    @property
    def nodnums0(self):
        """Host nodnums from the nodesmap.
        """
        return [0, 0, 1, 2, 3, 3]

    @property
    def nodnums1(self):
        """CCW-1 nodnums from the nodesmap.
        """
        return [3, 2, 0, 3, 1, 0]

    @property
    def nodnums2(self):
        """CCW-2 nodnums from the nodesmap.
        """
        return [1, 3, 3, 0, 0, 2]


class TestSkeletonWithVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.IJS = ijstream.matrix_stream(cls.getskel())

    @classmethod
    def getskel(cls):
        return skeleton.getskeleton(cls.getmesh())

    @classmethod
    def getmesh(cls):
        return trimesh.TriMesh.from_data(
            np.array(POINTS), np.array(TRIANGS)
        )

    def test_meta(self):
        assert self.IJS.meta['hasvoids'] is True
        assert self.IJS.meta['nodsmap-size'] == 14
        assert self.IJS.meta['westmap-size'] == 1
        assert self.IJS.meta['eastmap-size'] == 1
        assert self.IJS.meta['coremap-size'] == 1

    def test_skel(self):
        assert self.getskel().nodesmap.size == 14
        assert self.getskel().voidsmap.size == 4

    def test_stream_size(self):
        assert self.IJS.size == (14 + 4) * 3


if __name__ == '__main__':
    unittest.main()
