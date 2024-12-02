# -*- coding: utf-8 -*-
"""Tests mesh splitting.
"""
import unittest
from triellipt import mesher


class TestVoids(unittest.TestCase):

    EXTRA_TRIS = [
        [0, 3, 1],
        [3, 0, 4]
    ]

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.mesh().add_triangs(cls.EXTRA_TRIS)

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'cross-wise').deltriangs(0, 3)

    def test_voids(self):
        assert self.MESH.hasvoids() is True
        assert self.MESH.getvoids().tolist() == [3]

    def test_split_size(self):
        assert len(self.mesh_split) == 2

    def test_split_triangs(self):
        assert self.mesh_split[1].triangs.tolist() == [[0, 3, 1]]
        assert self.mesh_split[0].triangs.tolist() == [[4, 0, 2], [4, 2, 3]]

    @property
    def mesh_split(self):
        return self.MESH.split()


if __name__ == '__main__':
    unittest.main()
