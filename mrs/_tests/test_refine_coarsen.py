# -*- coding: utf-8 -*-
"""Tests refine-coarsen steps.
"""
import unittest
from triellipt import mesher
from triellipt import mrs


class TestRefine(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = mrs.getunit(
            mesher.trigrid(2, 2, 'east-slope')
        )

    def test_ranks_refined(self):
        assert self.UNIT.ranks.tolist() == [1, 1]
        assert self.UNIT.refine([0]).ranks.tolist() == [2, 1]
        assert self.UNIT.refine_all().ranks.tolist() == [3, 2]


class TestCoarsen(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = mrs.getunit(
            mesher.trigrid(2, 2, 'east-slope')
        )

    def test_ranks_coarsend(self):
        assert self.UNIT.ranks.tolist() == [1, 1]
        assert self.UNIT.refine_all().ranks.tolist() == [2, 2]
        assert self.UNIT.refine_all().ranks.tolist() == [3, 3]
        assert self.UNIT.coarsen([1]).ranks.tolist() == [3, 2]
        assert self.UNIT.coarsen_all().ranks.tolist() == [2, 1]


if __name__ == '__main__':
    unittest.main()
