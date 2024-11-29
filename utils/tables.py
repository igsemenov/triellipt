# -*- coding: utf-8 -*-
"""Tables operations.
"""
import numpy as np


def maptable(table):
    """Creates a map of table values.
    """
    _ = TableMapper.from_table(table)
    return _.makemap()


def table_image(table, cast):
    """Performs unstructured table indexing.

    Parameters
    ----------
    cast : same-as-table
        Desired horizontal indices row-by-row.

    Returns
    -------
    same-as-table
        Table image == table[cast]

    """

    base = table.shape[1]
    bias = np.indices(table.shape)[0]

    return table.flat[
        bias * base + cast
    ]


def trisorts(table_triplets):
    """Returns a sorter for a table of 0-1-2 triplets.
    """
    _ = TriSorter()
    return _.trisorts(table_triplets)


class TableMap:
    """Map of table values.
    """

    def __init__(self, data=None):

        if data is None:
            return

        self.data = data
        self.bins_split = self.get_bins_split()

    @classmethod
    def from_table(cls, table):
        """Makes a map from an integer table.
        """
        return TableMapper(table).makemap()

    def get_bins_split(self):

        _, bins = np.unique(
            self.vals, return_index=True
        )

        return np.copy(
            bins[1::].astype(int), order='C'
        )

    @property
    def vals(self):
        """Sorted table values.
        """
        return self.data[0, :]

    @property
    def rows(self):
        """Rows of sorted table values.
        """
        return self.data[1, :]

    @property
    def cols(self):
        """Cols of sorted table values.
        """
        return self.data[2, :]

    @property
    def nitems(self):
        """Number of unique values.
        """
        return len(self.bins_split) + 1

    @property
    def datasize(self):
        return self.data.shape[1]

    @property
    def vals_ranks(self):
        return np.diff(
            self.packs_edges
        )

    @property
    def vals_unique(self):
        return self.vals[self.packs_fronts]

    @property
    def vals_split(self):
        return np.split(
            self.vals, self.bins_split
        )

    @property
    def vals_reduced(self):
        return np.add.reduceat(
            self.vals, self.bins_reduce
        )

    @property
    def bins_reduce(self):
        return np.r_[
            0, self.bins_split
        ]

    @property
    def packs_edges(self):
        return np.r_[
            0, self.bins_split, self.datasize
        ]

    @property
    def packs_fronts(self):
        return np.r_[
            0, self.bins_split
        ]

    @property
    def packs_backs(self):
        return np.r_[
            self.bins_split - 1, self.datasize - 1
        ]

    def atrank(self, rank):
        """Yields data at the specified rank.
        """

        ranks = self.vals_ranks

        mask = np.repeat(
            ranks == rank, ranks
        )

        yield self.vals[mask]
        yield self.rows[mask]
        yield self.cols[mask]

    def data_split(self) -> list:
        """Returns data splitted into packs.
        """
        return np.split(
            self.data, self.bins_split, axis=1
        )


class TableAgent:
    """Operator on a table.
    """

    def __init__(self, table=None):

        if table is None:
            return

        self.table = np.array(table).astype(np.int64)

    @classmethod
    def from_table(cls, table):
        return cls(table)

    @property
    def vsize(self):
        return self.table.shape[0]

    @property
    def hsize(self):
        return self.table.shape[1]

    @property
    def rows(self):
        return np.arange(self.vsize)

    @property
    def cols(self):
        return np.arange(self.hsize)

    @property
    def rows2d(self):
        return np.copy(
            np.tile(self.rows, (self.hsize, 1)).T, order='C'
        )

    @property
    def cols2d(self):
        return np.copy(
            np.tile(self.cols, (self.vsize, 1)), order='C'
        )

    @property
    def values(self):
        return self.table


class TableMapper(TableAgent):
    """Maker of a table map.
    """

    def makemap(self):

        datalong = self.make_data_long()
        tablemap = self.from_data_long(datalong)

        return tablemap

    def make_data_long(self):

        rows = [
            self.values.flatten(),
            self.rows2d.flatten(),
            self.cols2d.flatten()
        ]

        return np.vstack(rows)

    def from_data_long(self, datalong):

        datasort = self.sort_data_long(datalong)
        tablemap = self.from_data_sort(datasort)

        return tablemap

    def sort_data_long(self, datalong):
        """Sorting per table values.
        """
        return datalong[
            :, np.argsort(datalong[0, :])
        ]

    def from_data_sort(self, datasort):
        return self.get_table_map(
            np.copy(datasort, order='C')
        )

    def get_table_map(self, data):
        return TableMap(data)


class TriSorter:
    """Sorter of 0-1-2 triplets.
    """

    def trisorts(self, table):
        """Returns a sorter for a table of 0-1-2 triplets.
        """

        codes = self.getcodes(table)

        perms = np.copy(
            table, order='C'
        )

        perms[codes == 11, :] = [1, 2, 0]
        perms[codes == 36, :] = [2, 0, 1]

        return perms.copy('C')

    def getcodes(self, table):
        return self.tripling(*table.T)

    def tripling(self, arg1, arg2, arg3):
        """Cantor tripling.
        """
        return self.pairing(
            self.pairing(arg1, arg2), arg3
        )

    def pairing(self, arg1, arg2):
        """Cantor pairing.
        """
        return arg2 + ((arg1 + arg2) * (arg1 + arg2 + 1)) // 2


def norm_table(table):
    """Returns a table from natural numbering of values.
    """
    _ = NormTable(table)
    return _.get_norm_table()


class NormTable:
    """Makes a new table from natural numbering of table values.
    """

    def __init__(self, table):
        self.tabledim = table.shape
        self.tablemap = maptable(table)

    def get_norm_table(self):

        new_vals = self.make_new_values()
        newtable = self.push_new_table(new_vals)

        return newtable

    def make_new_values(self):

        vals_natural_range = np.arange(self.tablemap.nitems)

        return np.repeat(
            vals_natural_range, self.tablemap.vals_ranks
        )

    def push_new_table(self, vals):

        table = self.get_table_mold()

        rows = self.tablemap.rows
        cols = self.tablemap.cols

        table[rows, cols] = vals
        return table

    def get_table_mold(self):
        return np.zeros(
            self.tabledim, dtype=int
        )
