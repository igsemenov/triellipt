# -*- coding: utf-8 -*-
"""Tests curve partition.
"""
import unittest
import numpy as np
from triellipt import geom
from triellipt.geom import partt


class TestCurvePartt(unittest.TestCase):

    def test_find_args(self):

        args = self.find_args()

        assert args[+0] == 0.
        assert args[-1] == 1.

        assert np.all(
            (args > 0.) @ (args < 1.)
        )

    def find_args(self):

        args = self.prttn.find_args(
            nparts=5, ratio=1, nseeds=100
        )

        return args

    def test_length_from_linspace(self):
        assert self.prttn.length_from_linspace(2) == 1.

    def test_lenfunc_from_linspace(self):
        assert self.prttn.lenfunc_from_linspace(2).tolist() == [0.5, 1.]

    @property
    def prttn(self):
        return partt.CurvePartt(
            curve=geom.line(0j, 1j)
        )

class TestLenEstim(unittest.TestCase):

    def test_get_length(self):
        assert isinstance(
            self.get_length(), float
        )

    def test_get_nseeds(self):
        assert isinstance(
            self.get_nseeds(), int
        )

    def get_length(self):
        return partt.LenEstim(self.line).get_length(rtol=1e-4, maxitr=10)

    def get_nseeds(self):
        return partt.LenEstim(self.line).get_nseeds(rtol=1e-4, maxitr=10)

    @property
    def line(self):
        return geom.line(0j, 1j)


if __name__ == '__main__':
    unittest.main()
