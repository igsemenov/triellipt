# -*- coding: utf-8 -*-
"""Tests removal of ghost points.
"""
import unittest
from triellipt import mesher


class TestCase(unittest.TestCase):

    @classmethod
    def mesh_with_ghosts(cls):
        return cls.mesh().deltriangs(0, 1, 2, 3)

    @classmethod
    def mesh_no_ghosts(cls):
        return cls.mesh_with_ghosts().delghosts()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 3, 'cross-wise')


class TestMeshWithGhosts(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.mesh_with_ghosts()

    def test_hasghosts(self):
        assert self.MESH.hasghosts() is True

    def test_getghosts(self):
        assert self.MESH.getghosts().tolist() == [0, 3, 6]


class TestCleanMesh(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.mesh_no_ghosts()

    def test_hasghosts(self):
        assert self.MESH.hasghosts() is False

    def test_getghosts(self):
        assert self.MESH.getghosts().tolist() == []


if __name__ == '__main__':
    unittest.main()
