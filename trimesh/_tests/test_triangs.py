# -*- coding: utf-8 -*-
"""Tests pairing of triangle edges.
"""
import unittest
from triellipt import mesher


class TestEdgesPaired(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.EDGES = cls.mesh().edges_paired()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'cross-wise')

    def test_edges_paired(self):
        assert self.EDGES[0, :].tolist() == [17, 1, 16]
        assert self.EDGES[1, :].tolist() == [16, 4, 18]
        assert self.EDGES[2, :].tolist() == [18, 11, 19]
        assert self.EDGES[3, :].tolist() == [19, 10, 17]


if __name__ == '__main__':
    unittest.main()
