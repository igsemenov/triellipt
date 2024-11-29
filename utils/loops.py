# -*- coding: utf-8 -*-
"""Loops fetcher.
"""
import numpy as np
from triellipt.utils import tables


class LoopsAsInds:
    """Extracts loops as indices of edges.

    Parameters
    ----------
    edges2d : 2-rows-int-table
        Input edges as nodes pairs in columns.

    """

    def __init__(self, edges2d=None):
        self.edges2d = edges2d
        self.nodesmap21 = None
        self.nodes_mask = None

    @classmethod
    def from_edges(cls, edges2d):
        return cls(edges2d)

    @property
    def size(self):
        return self.edges2d.shape[1]

    @property
    def nodes1(self):
        return self.edges2d[0, :]

    @property
    def nodes2(self):
        return self.edges2d[1, :]

    def setup(self):
        self.nodes_mask = self.get_nodes_mask()
        self.nodesmap21 = self.get_nodesmap21()

    def getloops(self):
        """Fetches loops as indices of edges.
        """

        self.setup()

        return list(
            self.genloops()
        )

    def genloops(self):

        for _ in range(self.size):

            inds = self.getloop()

            if not inds:
                return

            yield inds

            self.nodes_mask[inds] = 0

    def getloop(self):
        """Finds a subloop as indices of edges.
        """
        return list(
            self.genloop()
        )

    def genloop(self):

        prev = self.findstart()

        if prev is None:
            return

        yield prev

        map21 = self.nodesmap21
        start = self.nodes1[prev]

        for _ in range(self.size):

            curr = map21[prev]
            yield curr

            if self.nodes2[curr] == start:
                return

            prev = curr

    def findstart(self):

        where, = np.where(self.nodes_mask)

        if where.size == 0:
            return None
        return where[0]

    def get_nodes_mask(self):
        return np.full(self.size, 1)

    def get_nodesmap21(self):
        return MapNodes(self.nodes1, self.nodes2).getmap21()


class LoopsAsEdges:
    """Fetches loops as tables of edges.
    """

    def __init__(self, edges2d=None):
        self.edges2d = edges2d

    @classmethod
    def from_edges(cls, edges2d):
        return cls(edges2d)

    def getloops(self):
        return list(
            self.genloops()
        )

    def genloops(self):

        for inds in self.loops_as_inds:
            yield np.copy(
                self.edges2d[:, inds], order='C'
            )

    @property
    def loops_as_inds(self):
        return LoopsAsInds(self.edges2d).getloops()


class MapNodes:
    """Mapping of loops nodes.
    """

    def __init__(self, nodes1, nodes2):
        self.edges = self.edges_table(nodes1, nodes2)

    @classmethod
    def edges_table(cls, nodes1, nodes2):
        return np.copy(
            np.vstack([nodes1, nodes2]).T, order='C'
        )

    def getmap21(self):
        inds1, inds2 = self.make_intersect()
        return self.from_intersect(inds1, inds2)

    def make_intersect(self):

        _, rows, cols = self.edgesmap()

        rows = rows.reshape(rows.size // 2, 2)
        cols = cols.reshape(cols.size // 2, 2)

        inds = tables.table_image(
            rows, np.argsort(cols, axis=1)
        )

        yield inds[:, 0]
        yield inds[:, 1]

    def from_intersect(self, inds1, inds2):
        return inds1[
            self.invperm(inds2)
        ]

    def invperm(self, perm):
        return _invperm(perm)

    def edgesmap(self):
        return tables.TableMap.from_table(self.edges).atrank(2)


def _invperm(perm):
    iperm = np.zeros_like(perm)
    iperm[perm] = np.arange(perm.size)
    return iperm
