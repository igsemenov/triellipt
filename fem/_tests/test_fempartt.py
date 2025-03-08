# -*- coding: utf-8 -*-
"""Tests the unit partition.
"""
import unittest
from triellipt import mesher
from triellipt.fem import femunit


SPEC_WITH_DIRICHS = {
    'name': 'box',
    'anchors': [(1, 0), (1, 1), (0, 1)],
    'dirichlet-sides': (1, 3)
}

SPEC_NO_DIRICHS = {
    'name': 'box',
    'anchors': [(1, 0), (1, 1), (0, 1)],
    'dirichlet-sides': []
}


class Tester(unittest.TestCase):

    @classmethod
    def unit(cls):
        return femunit.getunit(
            cls.mesh(), anchors=[(0, 0)]
        )

    @classmethod
    def mesh(cls):
        return mesher.trigrid(5, 5, 'east-slope') * 0.25

    @property
    def box(self):
        return self.UNIT.partts['box']


class TestParttWithDirichlets(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = cls.unit().add_partition(SPEC_WITH_DIRICHS)

    def test_edge(self):
        assert self.box.edge[1].tolist() == [0, 1, 2, 3]
        assert self.box.edge[2].tolist() == [4, 5, 6, 7]
        assert self.box.edge[3].tolist() == [8, 9, 10, 11]
        assert self.box.edge[4].tolist() == [12, 13, 14, 15]

    def test_core(self):
        assert self.box.core.tolist() == [
            *[4, 5, 6, 7], *[12, 13, 14, 15], *range(16, 25)
        ]


class TestParttNoDirichlets(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = cls.unit().add_partition(SPEC_NO_DIRICHS)

    def test_edge(self):
        assert self.box.edge[1].tolist() == [0, 1, 2, 3]
        assert self.box.edge[2].tolist() == [4, 5, 6, 7]
        assert self.box.edge[3].tolist() == [8, 9, 10, 11]
        assert self.box.edge[4].tolist() == [12, 13, 14, 15]

    def test_core(self):
        assert self.box.core.tolist() == list(range(25))


if __name__ == '__main__':
    unittest.main()
