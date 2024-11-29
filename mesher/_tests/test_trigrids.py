# -*- coding: utf-8 -*-
"""Tests triangle grids.
"""
from abc import ABC, abstractmethod
import unittest
from triellipt import mesher


class TestTriGrid(ABC):
    """Base grid tester.
    """

    def test_triangs(self):
        for triang in self.grid.triangs.tolist():
            assert triang in self.TRIANGS_CONTROL

    @property
    def grid(self):
        return self.getgrid()

    @abstractmethod
    def getgrid(self):
        """Defines the grid sample.
        """


class TestTriEastSlope(TestTriGrid, unittest.TestCase):

    TRIANGS_CONTROL = [
        [0, 2, 3],
        [3, 1, 0],
        [2, 4, 5],
        [5, 3, 2]
    ]

    def getgrid(self):
        return mesher.trigrid(3, 2, 'east-slope')


class TestTriWestSlope(TestTriGrid, unittest.TestCase):

    TRIANGS_CONTROL = [
        [3, 4, 1],
        [1, 0, 3],
        [4, 5, 2],
        [2, 1, 4]
    ]

    def getgrid(self):
        return mesher.trigrid(2, 3, 'west-slope')


class TestTriEastSnake(TestTriGrid, unittest.TestCase):

    TRIANGS_CONTROL = [
        [0, 3, 4],
        [4, 1, 0],
        [4, 5, 2],
        [2, 1, 4]
    ]

    def getgrid(self):
        return mesher.trigrid(2, 3, 'east-snake')


class TestTriWestSnake(TestTriGrid, unittest.TestCase):

    TRIANGS_CONTROL = [
        [2, 3, 1],
        [1, 0, 2],
        [5, 3, 2],
        [2, 4, 5]
    ]

    def getgrid(self):
        return mesher.trigrid(3, 2, 'west-snake')


class TestTriCrossWise(TestTriGrid, unittest.TestCase):

    TRIANGS_CONTROL = [
        [4, 1, 0],
        [4, 0, 2],
        [4, 2, 3],
        [4, 3, 1]
    ]

    def getgrid(self):
        return mesher.trigrid(2, 2, 'cross-wise')


if __name__ == '__main__':
    unittest.main()
