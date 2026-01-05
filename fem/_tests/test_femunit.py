# -*- coding: utf-8 -*-
"""Tests the FEM unit.
"""
import unittest
import numpy as np
import triellipt as tri


class TestFEMUnit(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.UNIT = tri.fem.getunit(
            tri.mesher.trigrid(2, 2, 'east-slope')
        )

    def test_average(self):
        assert self.UNIT.average(self.UNIT.mesh.points.real)[0] == 2. / 3.
        assert self.UNIT.average(self.UNIT.mesh.points.real)[1] == 1. / 3.


if __name__ == '__main__':
    unittest.main()
