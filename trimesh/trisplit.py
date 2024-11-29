# -*- coding: utf-8 -*-
"""Static mesh splitter.
"""
import numpy as np
from scipy import sparse as sp


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def mesh_verts(self):
        return self.mesh.points[self.mesh.triangs]

    def edges_graph(self):

        edges = self.mesh.edgesmap()

        posi = edges.trinums1
        posj = edges.trinums2

        vals = np.full(
            edges.size, 1
        )
        
        shape = (
            self.mesh.size, self.mesh.size
        )

        graph = sp.coo_array(
            (vals, (posi, posj)), shape=shape
        )
        
        graph = graph.tocsr()
        return graph


class TriSplit(MeshAgent):
    """Splits a mesh into subsets.
    """

    def split(self) -> list:
        return list(
            self.genparts(self.mesh.size)
        )

    def delvoids(self):

        voidsnums = self.mesh.getvoids()

        if voidsnums.size == 0:
            return

        self.mesh = self.mesh.deltriangs(*voidsnums)

    def genparts(self, maxitr):

        self.delvoids()

        for _ in range(maxitr):

            imesh, omesh = self.run_mesh_bisect()
            yield imesh

            if omesh.size == 0:
                return

            self.mesh = omesh

    def run_mesh_bisect(self):

        graph = self.edges_graph()

        trinums = sp.csgraph.breadth_first_order(
            graph, 0, directed=False, return_predecessors=False
        )

        imesh = self.mesh.submesh(*trinums)
        omesh = self.mesh.deltriangs(*trinums)

        return imesh, omesh


class GetVoids(MeshAgent):
    """Finds empty triangles (voids).
    """

    TOL = 1e-12

    def find_voids(self):

        markers = self.make_markers()
        trinums = self.from_markers(markers)

        return trinums

    def make_markers(self):

        verts = self.mesh_verts()

        conds = np.abs(
            verts[:, 2] - 0.5 * (verts[:, 0] + verts[:, 1])
        )

        return conds < self.TOL

    def from_markers(self, markers):
        return np.where(markers)[0]
