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

    @property
    def triangs_range(self):
        return np.arange(self.mesh.ntriangs)

    @property
    def pivotsnums(self):
        return self.mesh.triangs[self.voidsnums, 2]

    @property
    def voidsnums(self):
        return self.mesh.getvoids()

    @property
    def hasvoids(self):
        return self.mesh.hasvoids()

    @property
    def mesh_twin(self):
        return self.mesh.twin()


class Renumer(MeshAgent):
    """Nodes renumerator.
    """

    def renumed(self, permuter):
        """Renumerates the mesh nodes.

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
        new_mesh = self.from_new_mesh_data(new_points, new_triangs)

        new_mesh = new_mesh.add_meta(
            {'nodes-permuter': permuter}
        )

        return new_mesh

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


class Shuffler(MeshAgent):
    """Triangles shuffler.
    """

    def shuffled(self, permuter):
        """Shuffles mesh triangles.

        Parameters
        ----------
        permuter : flat-int-array
            Permutation of mesh triangles.

        Returns
        -------
        TriMesh
            New mesh.

        """

        new_triangs = self.make_new_triangs(permuter)
        new_mesh = self.from_new_triangs(new_triangs)

        new_mesh = new_mesh.add_meta(
            {'triangs-permuter': permuter}
        )

        return new_mesh

    def make_new_triangs(self, permuter):
        return self.mesh.triangs[permuter, :]

    def from_new_triangs(self, new_triangs):
        return self.mesh.update_triangs(new_triangs)


class AlignNodes(MeshAgent):
    """Numbers mesh points in the edge-core order.
    """

    def aligned(self, *anchors):

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


class AlignVoids(MeshAgent):
    """Numbers mesh points so that voids pivots stay at then end.
    """

    def aligned(self):

        if not self.hasvoids:
            return self.mesh_twin_with_meta

        permuter = self.make_permuter()
        new_mesh = self.from_permuter(permuter)

        return new_mesh

    def make_permuter(self):

        points = self.points_range
        pivots = self.pivotsnums

        numbers = [
            np.delete(points, pivots), pivots
        ]

        return np.hstack(numbers)

    def from_permuter(self, permuter):
        return self.mesh.renumed(permuter)

    @property
    def mesh_twin_with_meta(self):

        mesh = self.mesh_twin

        return mesh.add_meta(
            {'nodes-permuter': np.arange(mesh.npoints)}
        )


class DownVoids(MeshAgent):
    """Puts void triangles at the end of triangulation table.
    """

    def shuffled(self):

        if not self.hasvoids:
            return self.mesh

        permuter = self.make_permuter()
        new_mesh = self.from_permuter(permuter)

        return new_mesh

    def make_permuter(self):

        voids = self.voidsnums
        triangs = self.triangs_range

        numbers = [
            np.delete(triangs, voids), voids
        ]

        return np.hstack(numbers)

    def from_permuter(self, permuter):
        return self.mesh.shuffled(permuter)
