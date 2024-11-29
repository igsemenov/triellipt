# -*- coding: utf-8 -*-
"""Tests the edges map.
"""
import unittest
from triellipt import mesher


class TestEdgesMap(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.EDGES = cls.mesh().edgesmap().sort_per_trinum()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'cross-wise')

    def test_trinums(self):
        assert self.EDGES.trinums1.tolist() == [0, 0, 1, 2]
        assert self.EDGES.trinums2.tolist() == [1, 3, 2, 3]

    def test_locnums(self):
        assert self.EDGES.locnums1.tolist() == [2, 0, 2, 2]
        assert self.EDGES.locnums2.tolist() == [0, 2, 0, 0]


class TestEdgesSpec(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.SPEC = cls.mesh().edgesmap().getspec()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(3, 3, 'west-snake').deltriangs(0)

    def test_spec(self):
        assert self.SPEC['spots'].tolist() == [0]
        assert self.SPEC['heads'].tolist() == [1, 4, 5]
        assert self.SPEC['links'].tolist() == [3, 6]
        assert self.SPEC['cores'].tolist() == [2]


if __name__ == '__main__':
    unittest.main()
