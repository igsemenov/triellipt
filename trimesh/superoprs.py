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

    def compress(self):

        if self.suptri.size == 0:
            return None

        edgesgraph = self.make_edges_graph()
        newsuptriu = self.from_edges_graph(edgesgraph)

        if newsuptriu is None:
            return None
        if newsuptriu.is_compact():
            return newsuptriu
        return None

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

        seed = 0

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


class SupReduce(SupCompress):
    """Reducer of a super-triangulation.
    """

    def reduce(self):

        maxitr = self.suptri.size

        for _ in range(maxitr):

            if self.suptri.size == 0:
                return None

            triu_compressed = self.compress()

            if triu_compressed is not None:
                return triu_compressed

            self.cleaning()

        return None

    def cleaning(self):
        self.suptri = self.suptri.strip()
        self.suptri = self.suptri.smooth()


class SupFlat(SupAgent):
    """Cleans the edge until low distortion.
    """

    def flatten(self, rtol):

        maxitr = self.suptri.size

        for _ in range(maxitr):

            if self.suptri.size == 0:
                return None

            edge_quality = self.edge_quality()

            if edge_quality < rtol:
                return self.suptri

            self.cleaning()

        return None

    def edge_quality(self):
        return np.amax(
            EdgeEval(self.suptri).edges_distorts()
        )

    def cleaning(self):
        self.suptri = self.suptri.strip()
        self.suptri = self.suptri.smooth()


class SupVoids(SupAgent):
    """Fetches a super-triangulation made of voids.
    """

    def supvoids(self):

        voids_nums = self.take_voids_nums()
        new_suptri = self.from_voids_nums(voids_nums)

        return new_suptri

    def take_voids_nums(self):
        return self.suptri.mesh.getvoids()

    def from_voids_nums(self, nums):

        if nums.size == 0:
            return None

        mask = np.isin(
            self.suptri.trinums, nums, assume_unique=True
        )

        inds, = np.where(mask)
        return self.suptri.subtriu(*inds)


class EdgeEval:
    """Evaluates the edge of a super-triangulation.
    """

    def __init__(self, suptri):
        self.suptri = suptri

    def edges_distorts(self):

        verts_onedge, edges_metric = self.get_data_to_eval()
        distorts = self.get_distorts(verts_onedge, edges_metric)

        return distorts

    def get_data_to_eval(self):

        sup_edges_metric = self.sup_edges_metric()
        ker_verts_onedge = self.ker_verts_onedge(sup_edges_metric)

        yield ker_verts_onedge
        yield sup_edges_metric

    def get_distorts(self, ker_verts_onedge, sup_edges_metric):

        kernodes = ker_verts_onedge
        supmidis = sup_edges_metric['centers']
        supspans = sup_edges_metric['lengths']

        return np.abs(kernodes - supmidis) / supspans

    def sup_edges_metric(self):
        return _mesh_edges_metric(self.suptri.supmesh)

    def ker_verts_onedge(self, sup_edges_metric):

        sup_edges_trinums = sup_edges_metric['trinums']
        sup_edges_locnums = sup_edges_metric['locnums']

        ker_verts_locnums = (sup_edges_locnums + 1) % 3

        return _mesh_subpoints(
            self.suptri.kermesh, sup_edges_trinums, ker_verts_locnums
        )


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


def _mesh_edges_metric(mesh):

    edge = mesh.meshedge()

    centers = 0.5 * np.sum(
        mesh.points[edge.edges2d], axis=0
    )

    lengths = np.abs(
        np.diff(mesh.points[edge.edges2d], axis=0)
    )

    locnums = edge.locnums.copy('C')
    trinums = edge.trinums.copy('C')

    return {
        'centers': centers,
        'lengths': lengths.flatten(),
        'locnums': locnums,
        'trinums': trinums
    }


def _mesh_subpoints(mesh, trinums, locnums):
    return mesh.points[
        mesh.triangs.flat[3 * trinums + locnums]
    ]
