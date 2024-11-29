# -*- coding: utf-8 -*-
"""Tests submeshing.
"""
from abc import ABC,  abstractmethod
import unittest
from triellipt import mesher


class Tester(ABC):
    """Base tester.
    """

    MESH = None
    SUBMESH = None

    def test_triangs_size(self):
        assert self.MESH.ntriangs == 8
        assert self.SUBMESH.ntriangs == 4

    def test_points(self):
        assert self.MESH.points is self.SUBMESH.points

    def test_triangs(self):
        for triang in self.triangs_control():
            assert triang in self.submesh_triangs

    def test_triangs_updated(self):
        assert self.SUBMESH.triangs.flags.owndata is True

    def test_points_are_same(self):
        assert self.SUBMESH.points is self.MESH.points

    @property
    def submesh_triangs(self):
        return self.SUBMESH.triangs.tolist()

    @abstractmethod
    def triangs_control(self):
        """Yields control triangles as triplets of numbers.
        """


class TestDelTriangs(Tester, unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        mesh = mesher.trigrid(3, 3, 'west-snake')

        cls.MESH = mesh
        cls.SUBMESH = mesh.deltriangs(0, 3, 4, 7)

    def triangs_control(self):
        yield [1, 0, 3]
        yield [3, 6, 7]
        yield [5, 2, 1]
        yield [7, 8, 5]


class TestSubMesh(Tester, unittest.TestCase):

    MESH = None

    @classmethod
    def setUpClass(cls):

        mesh = mesher.trigrid(3, 3, 'west-snake')

        cls.MESH = mesh
        cls.SUBMESH = mesh.submesh(0, 3, 4, 7)

    def triangs_control(self):
        yield [1, 4, 5]
        yield [5, 4, 7]
        yield [3, 4, 1]
        yield [7, 4, 3]


if __name__ == '__main__':
    unittest.main()
