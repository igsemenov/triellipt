# -*- coding: utf-8 -*-
"""Tests super operations.
"""
import unittest
from triellipt import mesher


class TestSupOprs(unittest.TestCase):
    """Tests basic super operations.
    """

    @classmethod
    def setUpClass(cls):
        cls.TRIU = cls.mesh().supertriu()

    @classmethod
    def mesh(cls):
        return mesher.trigrid(5, 5, 'west-slope')

    def test_suptriu_primary(self):
        assert self.TRIU.size == 18

    def test_detach(self):
        assert self.TRIU.detach().size == 2

    def test_strip(self):
        assert self.TRIU.strip().size == 10

    def test_smooth(self):
        assert self.TRIU.smooth(True).size == 6
        assert self.TRIU.smooth(False).size == 10

    def test_reduce(self):
        assert self.TRIU.reduce(iterate=True).size == 8
        assert self.TRIU.reduce(iterate=False).size == 8


if __name__ == '__main__':
    unittest.main()
