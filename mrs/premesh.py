# -*- coding: utf-8 -*-
"""Global mesh prototype.
"""
import itertools as itr
import numpy as np
from triellipt import mesher
from triellipt import utils


def getpremesh(skel):
    return PreMesh.from_skel(skel)


def getpremesher(skel):
    return PreMesher.from_skel(skel)


class PreMesh:
    """Global mesh prototype.
    """

    def __init__(self, skel=None, mesh=None, meta=None):
        self.skel = skel
        self.mesh = mesh
        self.meta = meta

    @classmethod
    def from_skel(cls, skel):
        return PreMesher.from_skel(skel).getpremesh()

    @property
    def data(self):
        return self.skel.data

    @property
    def holes(self):
        return self.meta['holes']

    @property
    def holes_flat(self):
        return np.hstack(self.holes)

    @property
    def gluers(self):
        return self.meta['gluers']

    @property
    def triu(self):
        return self.mesh.triu

    def getgluers(self):
        """Returns numbers of gluing-triangles.
        """

        nodes0 = self.mesh.points[self.mesh.triangs[:, 0]]
        nodes1 = self.mesh.points[self.mesh.triangs[:, 1]]

        trinums, = np.where(
            abs(nodes0 - nodes1) < 1e-10
        )

        return trinums

    def with_holes(self):
        self.meta['holes'] = HolesFetcher(self.mesh).get_holes()
        return self

    def with_gluers(self):
        self.meta['gluers'] = GluersFetcher(self).get_gluers()
        return self


class SkelAgent:
    """Operator on a global mesh skeleton.
    """

    def __init__(self, skel):
        self.skel = skel
        self.meta = self.fetch_meta()

    @classmethod
    def from_skel(cls, skel):
        return cls(skel)

    def fetch_meta(self):

        edges = self.get_host_edges_map()
        ranks = self.get_host_ranks_map(edges)

        return {
            'edgesmap': edges,
            'ranksmap': ranks
        }

    def get_host_edges_map(self):
        return self.skel.hostmesh.edgesmap()

    def get_host_ranks_map(self, edges):

        ranks = [
            self.skel.ranks[edges.trinums1],
            self.skel.ranks[edges.trinums2]
        ]

        return np.vstack(ranks)

    @property
    def edges(self):
        return self.meta.get('edgesmap')

    @property
    def ranks(self):
        return self.meta.get('ranksmap')

    @property
    def ranks1(self):
        return self.ranks[0, :]

    @property
    def ranks2(self):
        return self.ranks[1, :]


class PreMesher(SkelAgent):
    """Make a pre-mesh from a skeleton.
    """

    def getpremesh(self):

        if not self.is_meshable:
            raise RuntimeError(
                'MRS unit not meshable, bad contacts'
            )

        prem = PreMesh(
            self.skel, self.make_premesh(), {}
        )

        return prem.with_holes().with_gluers()

    @property
    def is_meshable(self):
        return np.all(
            abs(self.ranks2 - self.ranks1) <= 1
        )

    def make_premesh(self):
        extra_triangs = self.make_extra_triangs()
        return self.skel.ownmesh.add_triangs(extra_triangs)

    def make_extra_triangs(self):
        edges_pairs = self.fetch_edges()
        extra_triangs = self.merge_edges(edges_pairs)
        return extra_triangs

    def merge_edges(self, edges):

        pairs = zip(
            edges[0], edges[1]
        )

        gluers = list(
            itr.starmap(EdgesMerger(), pairs)
        )

        return np.vstack(gluers)

    def fetch_edges(self):
        return (
            self.fetch_edges_1(), self.fetch_edges_2()
        )

    def fetch_edges_1(self):

        items = zip(
            self.edges.trinums1, self.edges.locnums1
        )

        return list(
            itr.starmap(self.skel.triedges, items)
        )

    def fetch_edges_2(self):

        items = zip(
            self.edges.trinums2, self.edges.locnums2
        )

        return list(
            itr.starmap(self.skel.triedges, items)
        )


class EdgesMerger:
    """Creates the gluing submeshes between elements.
    """

    def __call__(self, edge1, edge2):

        size1 = edge1.size
        size2 = edge2.size

        if size1 == size2:
            return self.merge_conformal(edge1, edge2)

        if size1 == 2 * size2 - 1:
            return self.merge_non_conformal(edge1, edge2)
        if size2 == 2 * size1 - 1:
            return self.merge_non_conformal(edge2, edge1)

        raise ValueError(
            'got non-mergable edges, bad sizes'
        )

    def merge_conformal(self, edge1, edge2):
        """Merge conformal edges.
        """

        gluer = mesher.trigrid(
            2, edge1.size, 'east-slope'
        )

        nodes = np.hstack(
            [edge1, np.flip(edge2)]
        )

        triangs = nodes[gluer.triangs].copy('C')
        return triangs

    def merge_non_conformal(self, edge1, edge2):
        """Merge non-conformal edges at edge1 > edge2.
        """

        gluer = self.merge_conformal(edge1[0::2], edge2)

        voids = np.vstack(
            [edge1[0:-2:2], edge1[2::2], edge1[1:-1:2]]
        )

        return np.vstack(
            [gluer, voids.T]
        )


class HolesFetcher:
    """Fetches holes from the mesh.
    """

    def __init__(self, mesh):
        self.loops = self.fetch_loops(mesh)

    def fetch_loops(self, mesh):
        return mesh.meshedge().getloops()

    def get_holes(self):

        loopsitr = itr.compress(
            self.loops, self.get_loops_filter()
        )

        holes = list(loopsitr)

        return [
            hole.nodnums1 for hole in holes
        ]

    def get_loops_filter(self):

        loops_areas = self.get_loops_areas()

        return [
            abs(area) < 1e-10 for area in loops_areas
        ]

    def get_loops_areas(self):
        return [
            _loop_area(loop.nodes_complex) for loop in self.loops
        ]


class GluersFetcher:
    """Fetcher of gluers.
    """

    def __init__(self, premesh):
        self.prem = premesh

    def get_gluers(self):

        table = self.make_table()
        table = self.filter_table(table)

        return _filter_twins(table)

    def make_table(self):

        numbers = self.prem.getgluers()
        triangs = self.prem.mesh.triangs

        return triangs[numbers, 0:2].copy('C')

    def filter_table(self, table):

        mask = np.isin(
            table, self.holes_flat, invert=True
        )

        mask = np.all(mask, axis=1)
        return table[mask, :].copy('C')

    @property
    def holes_flat(self):
        return np.hstack(
            self.prem.holes
        )


def _loop_area(nodes):

    int_x = np.sum(nodes.imag[:-1] * np.diff(nodes.real))
    int_y = np.sum(nodes.real[:-1] * np.diff(nodes.imag))

    return np.abs(
        int_x - int_y
    )


def _filter_twins(table):

    codes = utils.pairs.sympaired(
        table[:, 0], table[:, 1]
    )

    data = np.vstack(
        utils.pairs.szuunpaired(np.unique(codes))
    )

    return data.T.astype(int)
