# -*- coding: utf-8 -*-
"""Tests mesh metric.
"""
import unittest
from triellipt import mesher
from triellipt.fem import femoprs


class TestMetric(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.MESH = cls.mesh()
        cls.METRIC = cls.metric()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'west-snake')

    @classmethod
    def metric(cls):
        return femoprs.mesh_metric(cls.MESH)

    def test_jacobis(self):
        assert self.METRIC['jacobis'].tolist() == [1., 1.]

    def test_bcoeffs(self):
        assert self.METRIC['bcoeffs'][0, :].tolist() == [0., +1., -1.]
        assert self.METRIC['bcoeffs'][1, :].tolist() == [0., -1., +1.]

    def test_ccoeffs(self):
        assert self.METRIC['ccoeffs'][0, :].tolist() == [-1., +1., 0.]
        assert self.METRIC['ccoeffs'][1, :].tolist() == [+1., -1., 0.]


if __name__ == '__main__':
    unittest.main()
