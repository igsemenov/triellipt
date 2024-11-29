# -*- coding: utf-8 -*-
"""Tests super triangulation.
"""
import unittest
from triellipt import mesher


class TestEmptyTriu(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.TRIU = cls.mesh().supertriu()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'west-slope')

    def test_is_empty(self):
        assert self.TRIU.size == 0


class TestSuperData(unittest.TestCase):
    """Tests super-triangulation.
    """

    @classmethod
    def setUpClass(cls):
        cls.TRIU = cls.mesh().supertriu()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(3, 3, 'west-slope')

    def test_trinums(self):
        assert self.TRIU.trinums.tolist() == [0, 7]

    def test_neighbors(self):
        assert self.TRIU.trinums1.tolist() == [5, 2]
        assert self.TRIU.trinums2.tolist() == [3, 4]
        assert self.TRIU.trinums3.tolist() == [1, 6]

    def test_vertices(self):
        assert self.TRIU.nodnums1.tolist() == [6, 2]
        assert self.TRIU.nodnums2.tolist() == [2, 6]
        assert self.TRIU.nodnums3.tolist() == [0, 8]

    def test_supmesh(self):
        assert self.TRIU.supmesh.triangs.tolist() == [[6, 2, 0], [2, 6, 8]]

    def test_kermesh(self):
        assert self.TRIU.kermesh.triangs.tolist() == [[3, 4, 1], [5, 4, 7]]

    def test_supbodies(self):
        assert self.TRIU.supbodies.tolist() == [
            [0, 5, 3, 1], [7, 2, 4, 6]
        ]


if __name__ == '__main__':
    unittest.main()
