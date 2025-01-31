# -*- coding: utf-8 -*-
"""Finds convex submeshes.
"""
import numpy as np
from scipy import sparse as sp


def find_subset(mesh, count, anchor):
    _ = SubsetFetcher.from_mesh(mesh).with_anchor(anchor)
    return _.find_trinums(count)


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = self.fetch_meta()

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def fetch_meta(self):
        return {
            'edgesmap': self.mesh.edgesmap()
        }

    @property
    def edgesmap(self):
        return self.meta['edgesmap']


class SubsetFetcher(MeshAgent):
    """Fetches a convex mesh subset.
    """

    def with_anchor(self, posx, posy):
        self.meta['anchor'] = posx + posy * 1j
        return self

    def find_trinums(self, count):

        subset = self.make_subset(count)
        subset = self.clean_subset(subset)

        return subset.copy('C')

    def make_subset(self, count):

        seed = self.find_seed_trinum()
        graph = self.make_edges_graph()

        subset = sp.csgraph.breadth_first_order(
            graph, seed, directed=False, return_predecessors=False
        )

        subset = subset[
            0: min(count, subset.size)
        ]

        return subset

    def find_seed_trinum(self):

        dist = np.abs(
            self.meta['anchor'] - self.mesh.centrs_complex
        )

        seed_trinum = np.argmin(dist)
        return seed_trinum

    def make_edges_graph(self):

        edges = self.edgesmap

        rows = edges.trinums1
        cols = edges.trinums2

        vals = np.ones(
            edges.size, dtype=int
        )

        graph_shape = (
            edges.size, edges.size
        )

        graph = sp.coo_array(
            (vals, (rows, cols)), shape=graph_shape
        )

        return graph.tocsr()

    def clean_subset(self, subset):

        if subset.size == 0:
            return subset

        submesh = self.mesh.submesh(*subset)

        subedges = submesh.edgesmap()
        subheads = subedges.getspec()['heads']

        if subheads.size == 0:
            return subset

        return np.delete(
            subset, subheads
        )
