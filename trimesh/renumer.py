# -*- coding: utf-8 -*-
"""Mesh renumerator.
"""
import numpy as np


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def points_range(self):
        return np.arange(self.mesh.npoints)


class Renumer(MeshAgent):
    """Mesh renumerator.
    """

    def renumed(self, permuter):
        """Renumerates the mesh.

        Parameters
        ----------
        permuter : flat-int-array
            Permutation of points numbers.

        Returns
        -------
        TriMesh
            New mesh.

        """

        new_points, new_triangs = self.make_new_mesh_data(permuter)
        new_trimesh = self.from_new_mesh_data(new_points, new_triangs)

        return new_trimesh

    def from_new_mesh_data(self, new_points, new_triangs):
        return self.mesh.from_data(new_points, new_triangs)

    def make_new_mesh_data(self, permuter):
        yield self.make_new_points(permuter)
        yield self.make_new_triangs(permuter)

    def make_new_points(self, permuter):
        return self.mesh.points[permuter]

    def make_new_triangs(self, permuter):
        returninds = self.make_return_map(permuter)
        newtriangs = self.from_return_map(returninds)
        return newtriangs

    def make_return_map(self, permuter):
        return self.points_range[np.argsort(permuter)]

    def from_return_map(self, returninds):
        return returninds[self.mesh.triangs]


class Arranger(MeshAgent):
    """Numbers mesh points in the edge-core order.
    """

    def arrange(self, *anchors):

        loops = self.make_loops(anchors)
        new_mesh = self.from_loops(loops)

        return new_mesh

    def make_loops(self, anchors):

        loops = self.take_loops()
        loops = self.sync_loops(loops, anchors)

        return loops

    def take_loops(self):

        edge = self.mesh.meshedge()

        try:
            return edge.getloops()
        except edge.MESH_EDGE_ERROR:
            return None

    def sync_loops(self, loops, anchors):

        if not loops:
            return None

        loops = dict(
            enumerate(loops)
        )

        for anchor in anchors:
            for key, loop in loops.items():
                loops[key] = self.sync_loop(loop, anchor)

        return loops

    def sync_loop(self, loop, anchor):
        if anchor in loop:
            return loop.synctonode(anchor)
        return loop

    def from_loops(self, loops):

        if loops is None:
            return None

        perm_points = self.make_permuter(loops)
        mesh_renumed = self.from_permuter(perm_points)

        return mesh_renumed

    def make_permuter(self, loops: dict):

        edge_nodes = [
            loop.nodnums1 for loop in loops.values()
        ]

        edge_nodes = np.hstack(edge_nodes)

        return np.r_[
            edge_nodes, np.delete(self.points_range, edge_nodes)
        ]

    def from_permuter(self, permuter):
        return self.mesh.renumed(permuter)
