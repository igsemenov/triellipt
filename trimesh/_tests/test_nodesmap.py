# -*- coding: utf-8 -*-
"""Tests the nodes map.
"""
import unittest
from triellipt import mesher


class TestNodesMap(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.NODES = cls.mesh().nodesmap().atnode(4)

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'cross-wise')

    def test_data(self):
        assert self.NODES.nodnums.tolist() == [4, 4, 4, 4]
        assert self.NODES.trinums.tolist() == [0, 1, 2, 3]
        assert self.NODES.locnums.tolist() == [0, 0, 0, 0]

    def test_locnums(self):
        assert self.NODES.locnums1.tolist() == [1, 1, 1, 1]
        assert self.NODES.locnums2.tolist() == [2, 2, 2, 2]

    def test_nodnums(self):
        assert self.NODES.nodnums1.tolist() == [1, 0, 2, 3]
        assert self.NODES.nodnums2.tolist() == [0, 2, 3, 1]


if __name__ == '__main__':
    unittest.main()
