# -*- coding: utf-8 -*-
"""Finding the mesh edge.
"""
import unittest
from triellipt import mesher


class TestMeshEdge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.EDGE = cls.mesh().meshedge().sort_per_trinum()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'cross-wise')

    def test_data(self):
        assert self.EDGE.trinums.tolist() == [0, 1, 2, 3]
        assert self.EDGE.locnums.tolist() == [1, 1, 1, 1]

    def test_locnums(self):
        assert self.EDGE.locnums1.tolist() == [1, 1, 1, 1]
        assert self.EDGE.locnums2.tolist() == [2, 2, 2, 2]
        assert self.EDGE.locnums3.tolist() == [0, 0, 0, 0]

    def test_nodnums(self):
        assert self.EDGE.nodnums1.tolist() == [1, 0, 2, 3]
        assert self.EDGE.nodnums2.tolist() == [0, 2, 3, 1]
        assert self.EDGE.nodnums3.tolist() == [4, 4, 4, 4]

    def test_unique_values(self):
        assert self.EDGE.trinums_unique.tolist() == [0, 1, 2, 3]
        assert self.EDGE.nodnums_unique.tolist() == [0, 1, 2, 3]

    def test_has_intersects(self):
        assert not self.EDGE.has_intersects()


if __name__ == '__main__':
    unittest.main()
