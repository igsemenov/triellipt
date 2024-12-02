# -*- coding: utf-8 -*-
"""Tests curves loop.
"""
import unittest
from triellipt import geom
from triellipt.geom import loop


class TestCurvesLoop(unittest.TestCase):

    def test_loop_normal(self):
        assert isinstance(
            geom.makeloop(*self.curves_normal), loop.CurvesLoop
        )

    def test_loop_unclosed(self):
        with self.assertRaises(loop.CurvesLoopError):
            geom.makeloop(*self.curves_unclosed)

    def test_loop_unlinked(self):
        with self.assertRaises(loop.CurvesLoopError):
            geom.makeloop(*self.curves_unlinked)

    def test_discretize(self):
        assert isinstance(
            self.loop_normal.discretize((10, 1)), loop.PathMap
        )

    @property
    def loop_normal(self):
        return geom.makeloop(*self.curves_normal)

    @property
    def curves_normal(self):
        yield geom.line(0, 1)
        yield geom.line(1, 0)

    @property
    def curves_unclosed(self):
        yield geom.line(1, 2)
        yield geom.line(2, 0)

    @property
    def curves_unlinked(self):
        yield geom.line(0, 3)
        yield geom.line(3, 1)
        yield geom.line(2, 0)


class TestMakeLoop(unittest.TestCase):

    def test_makerect(self):
        assert isinstance(
            geom.makerect(0j, (1, 2)), loop.CurvesLoop
        )

    def test_makeellip(self):
        assert isinstance(
            geom.makeellip(0j, (1, 2)), loop.CurvesLoop
        )


if __name__ == '__main__':
    unittest.main()
