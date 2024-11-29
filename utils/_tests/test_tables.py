# -*- coding: utf-8 -*-
"""Tests tables operations.
"""
import itertools as itr
import unittest
import numpy as np
from triellipt.utils import tables


class TestTableMap(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        table = [
            [1, 2, 1],
            [3, 1, 3],
        ]

        cls.MAP = tables.TableMap.from_table(table)

    def test_size(self):
        assert self.MAP.nitems == 3
        assert self.MAP.datasize == 6

    def test_data(self):
        assert self.MAP.vals.tolist() == [1, 1, 1, 2, 3, 3]
        assert self.MAP.rows.tolist() == [0, 0, 1, 0, 1, 1]
        assert self.MAP.cols.tolist() == [0, 2, 1, 1, 0, 2]

    def test_vals_split(self):
        assert self.MAP.vals_split[0].tolist() == [1, 1, 1]
        assert self.MAP.vals_split[2].tolist() == [3, 3]
        assert self.MAP.vals_split[1].tolist() == [2]

    def test_bins_reduce(self):
        assert self.MAP.vals_reduced.tolist() == [3, 2, 6]

    def test_vals_ranks(self):
        assert self.MAP.vals_ranks.tolist() == [3, 1, 2]

    def test_vals_unique(self):
        assert self.MAP.vals_unique.tolist() == [1, 2, 3]


class TestTableImage(unittest.TestCase):

    # Must be TABLE[CAST] == TABLE_IMAGE

    TABLE = [
        [0, 1, 2],
        [3, 4, 5]
    ]

    CAST = [
        [0, 1, 2],
        [2, 0, 1]
    ]

    TABLE_IMAGE = [
        [0, 1, 2],
        [5, 3, 4]
    ]

    @classmethod
    def setUpClass(cls):
        cls.FUNC = tables.table_image

    def test_table_image(self):
        assert self.table_image.tolist() == self.TABLE_IMAGE

    @property
    def table_image(self):
        return self.table_image_func(
            np.array(self.TABLE), np.array(self.CAST)
        )

    @property
    def table_image_func(self):
        return self.__class__.FUNC


class TestTriSorter(unittest.TestCase):

    TABLE = list(
        itr.permutations(range(3))
    )

    @classmethod
    def setUpClass(cls):
        cls.SORTER = tables.TriSorter()

    def test_pairing(self):
        assert self.SORTER.pairing(0, 1) == 2
        assert self.SORTER.pairing(2, 2) == 12
        assert self.SORTER.pairing(3, 4) == 32

    def test_tripling(self):
        assert self.SORTER.tripling(0, 1, 2) == 12
        assert self.SORTER.tripling(2, 0, 1) == 11
        assert self.SORTER.tripling(1, 2, 0) == 36

    def test_triplsorter(self):
        assert self.table_sorted.tolist() == [[0, 1, 2]] * 6

    @property
    def table_sorted(self):
        return tables.table_image(
            self.table, self.table_sorter
        )

    @property
    def table_sorter(self):
        return tables.trisorts(self.table)

    @property
    def table(self):
        return np.array(self.TABLE)


if __name__ == '__main__':
    unittest.main()
