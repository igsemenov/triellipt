# -*- coding: utf-8 -*-
"""Tests MRS skeleton.
"""
import unittest
from triellipt import mesher
from triellipt import mrs


class TestSkeleton(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.UNIT = mrs.getunit(
            mesher.trigrid(2, 2, 'east-slope')
        )

        cls.SKEL = cls.UNIT.refine(1).getskel()

    def test_triedges(self):

        # down-triangle (0)
        assert self.SKEL.triedges(0, 0).tolist() == [0, 1, 2]
        assert self.SKEL.triedges(0, 1).tolist() == [2, 3, 4]
        assert self.SKEL.triedges(0, 2).tolist() == [4, 5, 0]

        # up-triangle (1)
        assert self.SKEL.triedges(1, 0).tolist() == [6,  7,  8,  9, 10]
        assert self.SKEL.triedges(1, 1).tolist() == [10, 11, 12, 13, 14]
        assert self.SKEL.triedges(1, 2).tolist() == [14, 15, 16, 17,  6]


if __name__ == '__main__':
    unittest.main()
