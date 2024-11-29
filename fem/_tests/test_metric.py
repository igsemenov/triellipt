# -*- coding: utf-8 -*-
"""Tests mesh metric.
"""
import unittest
from triellipt import mesher
from triellipt.fem import femoprs


class TestMetricNoVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.METRIC = cls.mesh_metric()

    @classmethod
    def mesh_metric(cls):
        return femoprs.mesh_metric(cls.mesh())

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'west-snake')

    def test_jacobis(self):
        assert self.METRIC.jacobis.tolist() == [1., 1.]

    def test_bcoeffs(self):
        assert self.METRIC.bcoeffs[0, :].tolist() == [0., +1., -1.]
        assert self.METRIC.bcoeffs[1, :].tolist() == [0., -1., +1.]

    def test_ccoeffs(self):
        assert self.METRIC.ccoeffs[0, :].tolist() == [-1., +1., 0.]
        assert self.METRIC.ccoeffs[1, :].tolist() == [+1., -1., 0.]

    def test_voids(self):
        assert self.METRIC.hasvoids is False
        assert self.METRIC.voids_trinums is None


class TestMetricWithVoids(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.METRIC = cls.mesh_metric()

    @classmethod
    def mesh_metric(cls):
        return femoprs.mesh_metric(cls.mesh())

    @classmethod
    def mesh(cls):
        _ = mesher.trigrid(2, 2, 'west-snake')
        return _.add_points(0.5).add_triangs([2, 0, 4])

    def test_voids(self):
        assert self.METRIC.hasvoids is True
        assert self.METRIC.voids_trinums.tolist() == [2]

    def test_data(self):
        assert self.METRIC.ccoeffs[2, :].tolist() == [0., 0., 0.]
        assert self.METRIC.bcoeffs[2, :].tolist() == [0., 0., 0.]


if __name__ == '__main__':
    unittest.main()
