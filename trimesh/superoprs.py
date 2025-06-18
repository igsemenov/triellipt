# -*- coding: utf-8 -*-
"""Super operations.
"""
from abc import ABC, abstractmethod
import numpy as np
from scipy import sparse as sp


class SupAgent(ABC):
    """Operator on a super triangulation.
    """

    def __init__(self, suptri=None):
        self.suptri = suptri
        self.cache = {}

    @classmethod
    def from_suptriu(cls, suptri):
        return cls(suptri)

    @property
    def supverts(self):
        return self._vertstable(self.suptri.supmesh)

    @property
    def kerverts(self):
        return self._vertstable(self.suptri.kermesh)

    @property
    def supedgesmap(self):
        return self.suptri.supmesh.edgesmap()

    @property
    def backmeshedge(self):
        """Nodes on the edge of background mesh.
        """
        return self.suptri.mesh.meshedge().nodnums_unique

    @classmethod
    def _vertstable(cls, mesh):
        return mesh.points[mesh.triangs]


class SupCleaner(SupAgent):
    """Removes bad triangles from a super-triangulation.
    """

    def cleaned(self):

        badtrinums = self.find_bad_trinums()
        newsuptriu = self.make_new_suptriu(badtrinums)

        return newsuptriu

    @abstractmethod
    def find_bad_trinums(self):
        """Finds numbers of super-triangles to delete.
        """

    def make_new_suptriu(self, badtrinums):
        return self.suptri.deltriangs(*badtrinums)


class SupStrip(SupCleaner):
    """Removes links from a super-triangulation.
    """

    def find_bad_trinums(self):
        spec = self.supedgesmap.getspec()
        return spec['links']


class SupSmooth(SupCleaner):
    """Removes heads and spots from a super-triangulation.
    """

    def cleaned_iteraly(self):
        """Runs smoothing iteratively.
        """

        count = 0

        smoothen = True
        maxcount = self.suptri.size

        while smoothen:

            self.suptri, smoothen = self.iteration()

            count += 1

            if self.suptri.size == 0:
                return self.suptri
            if count == maxcount:
                return self.suptri

        return self.suptri

    def iteration(self):

        newsuptriu = self.cleaned()

        no_smoothing = newsuptriu.size == self.suptri.size
        is_smoothing = not no_smoothing

        return newsuptriu, is_smoothing

    def find_bad_trinums(self):

        spec = self.supedgesmap.getspec()

        return np.r_[
            spec['spots'], spec['heads']
        ]


class SupDetach(SupCleaner):
    """Removes super-triangles touching the background edge.
    """

    def find_bad_trinums(self):

        mask = np.isin(
            self.suptri.suptriangs, self.backmeshedge
        )

        badnums, = np.where(
            np.any(mask, axis=1)
        )

        return badnums


SupReduceError = type(
    'SupReduceError', (Exception,), {}
)


class SupCompress(SupAgent):
    """Compressor of a super-triangulation.
    """

    MIN_SUPTRI_SIZE = 4

    def compress(self, seed=None):

        if self.suptri.size == 0:
            return None
        self.cache['seed'] = seed

        edgesgraph = self.make_edges_graph()
        newsuptriu = self.from_edges_graph(edgesgraph)

        return newsuptriu

    def make_edges_graph(self):
        return self.supedgesgraph()

    def from_edges_graph(self, graph):

        trinums = self.make_bfs_tree(graph)
        subtriu = self.from_bfs_tree(trinums)

        return subtriu

    def from_bfs_tree(self, trinums):
        if trinums is None:
            return None
        return self.suptri.subtriu(*trinums)

    def make_bfs_tree(self, graph):

        seed = self.get_bfs_seed()

        nums = sp.csgraph.breadth_first_order(
            graph, seed, directed=False, return_predecessors=False
        )

        if np.unique(nums).size != nums.size:
            return None
        if nums.size < self.MIN_SUPTRI_SIZE:
            return None

        return nums

    def supedgesgraph(self):
        """Super-edges graph.
        """

        supedges = self.supedgesmap

        rows = supedges.trinums1
        cols = supedges.trinums2

        shape = (self.suptri.size,) * 2

        graph = sp.coo_array(
            (np.ones_like(rows), (rows, cols)), shape
        )

        graph = graph.tocsr()

        if graph.has_canonical_format:
            return graph

        raise SupReduceError(
            "super edges graph is not in canonical format"
        )

    def get_bfs_seed(self):

        seed = self.cache['seed']

        if seed is None:
            return 0
        return self.suptri.find_seed(seed)


class SupReduce(SupCompress):
    """Reducer of a super-triangulation.
    """

    def reduce(self, seed=None):

        while True:

            if self.suptri.size == 0:
                return None

            suptri_compressed = self.trial_reduce(seed)

            if suptri_compressed is not None:
                return suptri_compressed

            self.clean_suptri()

        return None

    def trial_reduce(self, seed):

        compressed = self.compress(seed)

        if compressed is None:
            return None
        if not compressed.is_compact():
            return None

        return compressed

    def clean_suptri(self):
        self.suptri = self.suptri.strip()
        self.suptri = self.suptri.smooth()


class SupVoids(SupAgent):
    """Fetches a super-triangulation made of voids.
    """

    def supvoids(self):

        voidsnums = self.take_voids_nums()
        newsuptri = self.from_voids_nums(voidsnums)

        return newsuptri

    def take_voids_nums(self):
        return self.suptri.mesh.getvoids()

    def from_voids_nums(self, trinums):

        if trinums.size == 0:
            return None

        mask = np.isin(
            self.suptri.trinums, trinums, assume_unique=True
        )

        inds, = np.where(mask)
        return self.suptri.subtriu(*inds)


def triangs_distorts(supverts, kerverts):
    """Distortion measure of super-triangles.

    Parameters
    ----------
    supverts : 3-column-complex-table
        Vertices of super-triangles.
    kerverts : 3-column-complex-table
        Vertices of host triangles in super-triangulation.

    Returns
    -------
    flat-float-array
        Distortion measure across triangles.

    """

    overts1 = supverts[:, [0, 1, 2]]
    overts2 = supverts[:, [1, 2, 0]]

    midis = 0.5 * (overts2 + overts1)
    bases = np.abs(overts2 - overts1)

    disps = np.abs(
        kerverts[:, [1, 2, 0]] - midis
    )

    return np.amax(
        2.0 * (disps / bases), axis=1
    )


def triangs_qscores(verts):
    """Quality score of mesh triangles.

    Parameters
    ----------
    verts : 3-column-complex-table
        Vertices of mesh triangles.

    Returns
    -------
    flat-float-array
        Normalized quality score across triangles.

    """

    sides = np.diff(
        verts[:, [0, 1, 2, 0]], axis=1
    )

    theta = abs(
        np.sum(sides*sides, axis=1)
    )

    return theta / np.amax(theta)
