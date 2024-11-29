# -*- coding: utf-8 -*-
"""Mesh inspector.
"""
import numpy as np
from scipy import sparse as sp


def find_node(mesh, anchor):
    _ = NodeFinder.from_mesh(mesh).with_anchor(anchor)
    return _.find_node()


def find_subset(mesh, count, anchor, remove_heads):
    _ = SubsetFinder.from_mesh(mesh).with_anchor(anchor)
    return _.find_trinums(count, remove_heads)


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = self.fetch_meta()
        self.cache = {}

    def fetch_meta(self):
        return {}

    def with_anchor(self, anchor):
        self.cache['anchor'] = anchor
        return self

    @property
    def anchor_complex(self):
        return _pack_complex(
            *self.cache['anchor']
        )

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)


class NodeFinder(MeshAgent):
    """Finds a node near an anchor.
    """

    def fetch_meta(self):
        return {
            'nodesmap': self.mesh.nodesmap()
        }

    @property
    def nodesmap(self):
        return self.meta['nodesmap']

    def find_node(self):

        node = self.find_closest_node()
        data = self.push_closest_node_map(node)

        return data

    def push_closest_node_map(self, node):
        return self.nodesmap.atnode(node)

    def find_closest_node(self):
        return np.argmin(
            np.abs(self.mesh.points - self.anchor_complex)
        )


class SubsetFinder(MeshAgent):
    """Fetches a convex mesh subset.
    """

    def fetch_meta(self):
        return {
            'edgesmap': self.mesh.edgesmap()
        }

    @property
    def edgesmap(self):
        return self.meta['edgesmap']

    def find_trinums(self, count, remove_heads):
        self.cache['remove-heads'] = bool(remove_heads)

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
            self.anchor_complex - self.mesh.centrs_complex
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

        if self.cache['remove-heads'] is False:
            return subset

        submesh = self.mesh.submesh(*subset)

        subedges = submesh.edgesmap()
        subheads = subedges.getspec()['heads']

        if subheads.size == 0:
            return subset

        return np.delete(
            subset, subheads
        )


def _pack_complex(posx, posy):
    return posx + 1j * posy
