# -*- coding: utf-8 -*-
"""Tests path orienter.
"""
import unittest
import numpy as np
from triellipt.utils import loops


class TestNodesMap(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        nodes1 = np.r_[4, 2, 3, 0, 1]
        nodes2 = np.r_[0, 3, 4, 1, 2]

        cls.MAP21 = loops.MapNodes(nodes1, nodes2).getmap21()

        cls.NODES1 = nodes1
        cls.NODES2 = nodes2

    def test_getmap21(self):
        assert np.all(
            self.NODES1[self.MAP21] == self.NODES2
        )


class TestLoops(unittest.TestCase):
    """Loops tester.
    """

    LOOP1 = [
        [0, 1, 2],
        [1, 2, 0]
    ]

    LOOP2 = [
        [3, 4, 5],
        [4, 5, 3]
    ]

    PERM_EDGES = [
        0, 3, 5, 1, 4, 2
    ]

    @classmethod
    def make_loops(cls, loops_type):
        return loops_type(
            cls.make_edges2d()
        )

    @classmethod
    def make_edges2d(cls):

        edges2d = np.hstack(
            [np.array(cls.LOOP1), np.array(cls.LOOP2)]
        )

        edges2d = np.copy(
            edges2d[:, cls.PERM_EDGES], order='C'
        )

        return edges2d


class TestLoopsAsInds(TestLoops):

    @classmethod
    def setUpClass(cls):
        cls.LOOPS = cls.make_loops(loops.LoopsAsInds)

    def test_loops(self):
        assert self.LOOPS.getloops()[0] == [0, 3, 5]
        assert self.LOOPS.getloops()[1] == [1, 4, 2]


class TestLoopsAsEdges(TestLoops):

    @classmethod
    def setUpClass(cls):
        cls.LOOPS = cls.make_loops(loops.LoopsAsEdges)

    def test_loop_one(self):
        assert self.LOOPS.getloops()[0][0, :].tolist() == [0, 1, 2]
        assert self.LOOPS.getloops()[0][1, :].tolist() == [1, 2, 0]

    def test_loop_two(self):
        assert self.LOOPS.getloops()[1][0, :].tolist() == [3, 4, 5]
        assert self.LOOPS.getloops()[1][1, :].tolist() == [4, 5, 3]


if __name__ == '__main__':
    unittest.main()
