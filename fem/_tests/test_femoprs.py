# -*- coding: utf-8 -*-
"""Tests mesh metric.

TEST-MESH-NO-VOIDS:
    
    1--3
    |\ |
    | \|
    0--2

    Triangs
    -------
    2-3-1
    1-0-2

TEST-MESH-WITH-VOID:
    
    1---3
    |\  |
    | \ |
    |  \|
    0-4-2

    Triangs
    -------
    2-3-1
    1-0-2
    2-0-4

"""
import unittest
import numpy as np
from triellipt import mesher
from triellipt.fem import femoprs


class TestData:

    OPRS = {}

    @property
    def diff_2x(self):
        return self.OPRS['diff_2x']

    @property
    def diff_2y(self):
        return self.OPRS['diff_2y']

    @property
    def diff_1x(self):
        return self.OPRS['diff_1x']

    @property
    def diff_1y(self):
        return self.OPRS['diff_1y']

    @property
    def massmat(self):
        return self.OPRS['massmat']

    @property
    def massdig(self):
        return self.OPRS['massdig']


class TestOprsNoVoids(TestData, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.OPRS = cls.fem_oprs()

    @classmethod
    def fem_oprs(cls):
        return femoprs.getoprs(cls.mesh())

    @classmethod
    def mesh(cls):
        return mesher.trigrid(2, 2, 'west-snake')

    def test_shapes(self):
        for key in self.OPRS:
            assert self.OPRS[key].shape == (2, 9)

    def test_diff_2d(self):
        assert self.diff_2x[0, 0] == 0.0
        assert self.diff_2y[0, 0] == 0.5

    def test_diff_2x_symm(self):
        assert np.all(self.diff_2x[:, 1] == self.diff_2x[:, 3])
        assert np.all(self.diff_2x[:, 2] == self.diff_2x[:, 6])
        assert np.all(self.diff_2x[:, 5] == self.diff_2x[:, 7])

    def test_diff_2y_symm(self):
        assert np.all(self.diff_2y[:, 1] == self.diff_2y[:, 3])
        assert np.all(self.diff_2y[:, 2] == self.diff_2y[:, 6])
        assert np.all(self.diff_2y[:, 5] == self.diff_2y[:, 7])

    def test_diff_1d(self):
        assert np.all(self.diff_1x[:, 0] == self.diff_1x[:, 1])
        assert np.all(self.diff_1y[:, 0] == self.diff_1y[:, 2])

    def test_massdig(self):
        assert np.all(self.massdig[:, [1, 2, 3, 5, 6, 7]] == 0.)

    def test_massmat(self):
        assert np.sum(self.massmat[0, [0, 1, 2]]) == self.massdig[0, 0]


class TestOprsWithVoids(TestData, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.OPRS = cls.fem_oprs()

    @classmethod
    def fem_oprs(cls):
        return femoprs.getoprs(cls.mesh())

    @classmethod
    def mesh(cls):
        _ = mesher.trigrid(2, 2, 'west-snake')
        return _.add_points(0.5).add_triangs([2, 0, 4])

    def test_shapes(self):
        for key in self.OPRS:
            assert self.OPRS[key].shape == (3, 9)

    def test_diff_2d(self):
        assert np.all(self.diff_2x[2, :] == 0.)
        assert np.all(self.diff_2y[2, :] == 0.)

    def test_diff_1d(self):
        assert np.all(self.diff_1x[2, :] == 0.)
        assert np.all(self.diff_1y[2, :] == 0.)


if __name__ == '__main__':
    unittest.main()
