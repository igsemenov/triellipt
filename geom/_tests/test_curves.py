# -*- coding: utf-8 -*-
"""Tests curves.
"""
import unittest
from triellipt import geom


class TestCurve(unittest.TestCase):
    """Curves tester.
    """

    def test_startpoint(self):
        for curve in self.curves():
            assert curve.startpoint == curve.getpath(0)

    def test_endpoint(self):
        for curve in self.curves():
            assert curve.endpoint == curve.getpath(1)

    def test_is_curve(self):
        for curve in self.curves():
            assert curve.is_curve() is True

    def test_length(self):
        assert geom.line(0j, 1j).length() == 1.

    def curves(self):
        """Defines curves to test.
        """
        return [
            geom.line(0j, 1j),
            geom.elliparc(0j, (1, 1), (0, 2)),
            geom.bezier2(0j, 1j, 2j),
            geom.bezier3(0j, 1j, 2j, 3j)
        ]


if __name__ == '__main__':
    unittest.main()
