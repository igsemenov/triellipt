# -*- coding: utf-8 -*-
"""Tests cycle path.
"""
import unittest
import numpy as np
from triellipt import geom
from triellipt.geom import cycle


class TestCyclePath(unittest.TestCase):

    def test_cycle_type(self):
        assert isinstance(
            self.cycle, cycle.CycPath
        )

    def test_nodenums(self):
        assert self.cycle.numbers.tolist() == [0, 1, 2, 3]

    def test_points(self):
        assert self.cycle.points.tolist() == [0j, 1j, 2j, 3j]

    def test_neighbours(self):
        assert self.cycle.prev_points.tolist() == [3j, 0j, 1j, 2j]
        assert self.cycle.next_points.tolist() == [1j, 2j, 3j, 0j]

    def test_angles(self):
        assert np.all(
            self.cycle.angles() == np.r_[1., 0., 0., -1.] * np.pi
        )

    @property
    def cycle(self):
        return geom.makecycle(
            np.array([0j, 1j, 2j, 3j])
        )


if __name__ == '__main__':
    unittest.main()
