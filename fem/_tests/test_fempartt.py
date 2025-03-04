# -*- coding: utf-8 -*-
"""Tests the mesh skeleton.
"""
import unittest
from triellipt import mesher
from triellipt.fem import femunit


PARTT_ANGLE = {
    'partition-title': 'box',
    'partition-loops': {
        0: {
            'angle': 1.5,
            'coloring': [
                (1, 2, 'rshift'), (3, 4, 'rshift')
            ]
        }
    },
    'dirichlet-sides': [1, 3]
}

PARTT_BINS = {
    'partition-title': 'box',
    'partition-loops': {
        0: {
            'bins': [5, 2, 5, 2],
            'coloring': []
        }
    },
    'dirichlet-sides': []
}


class Tester(unittest.TestCase):

    UNIT = None

    @classmethod
    def unit(cls):
        return femunit.getunit(
            cls.mesh(), anchors=(0,)
        )

    @classmethod
    def mesh(cls):
        return mesher.trigrid(5, 4, 'east-slope')

    @property
    def box(self):
        return self.UNIT.partts['box']


class TestParttAngle(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = cls.unit().add_partition(PARTT_ANGLE)

    def test_edge(self):
        assert self.box.edge[1].tolist() == [0, 1, 2, 3, 4]
        assert self.box.edge[3].tolist() == [7, 8, 9, 10, 11]
        assert self.box.edge[2].tolist() == [5, 6]
        assert self.box.edge[4].tolist() == [12, 13]

    def test_core(self):
        assert self.box.core.tolist() == [5, 6, 12, 13, *range(14, 20)]

    def test_dirich_sides(self):
        assert self.box.dirich_sides == {1, 3}


class TestParttBins(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = cls.unit().add_partition(PARTT_BINS)

    def test_edge(self):
        assert self.box.edge[1].tolist() == [0, 1, 2, 3, 4]
        assert self.box.edge[3].tolist() == [7, 8, 9, 10, 11]
        assert self.box.edge[2].tolist() == [5, 6]
        assert self.box.edge[4].tolist() == [12, 13]

    def test_core(self):
        assert self.box.core.tolist() == list(range(20))

    def test_dirich_sides(self):
        assert self.box.dirich_sides == set()


if __name__ == '__main__':
    unittest.main()
