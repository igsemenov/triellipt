# -*- coding: utf-8 -*-
"""Cleaner of ghost nodes.
"""
import numpy as np


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)


class DelGhosts(MeshAgent):
    """Removes ghost nodes.
    """

    def cleanmesh(self):
        """Returns a new mesh without ghost nodes.
        """

        ghostsdata = self.find_ghost_nodes()
        new_points, new_triangs = self.new_mesh_data(ghostsdata)

        return self.new_mesh(
            new_points, new_triangs
        )

    def find_ghost_nodes(self):
        return self.ghosts_finder.getdata()

    def new_mesh_data(self, ghosts_data):

        old_nodes_actual = ghosts_data['old-nodes-actual']
        new_triangs_flat = ghosts_data['new-triangs-flat']

        yield self.get_new_points(old_nodes_actual)
        yield self.get_new_triangs(new_triangs_flat)

        self.cache['old-nodes-numbers'] = old_nodes_actual

    def get_new_points(self, old_nodes_actual):
        """New points as a flat-complex-array.
        """
        return self.mesh.points[old_nodes_actual]

    def get_new_triangs(self, new_triangs_flat):
        """New triangles table.
        """
        return np.reshape(
            new_triangs_flat, (self.mesh.ntriangs, 3)
        )

    def new_mesh(self, new_points, new_triangs):

        mesh = self.mesh.from_data(new_points, new_triangs)

        mesh = mesh.add_meta(
            self.new_mesh_meta()
        )

        return mesh

    def new_mesh_meta(self):
        return {
            'old-nodes-numbers': self.cache['old-nodes-numbers']
        }

    @property
    def ghosts_finder(self):
        return GetGhosts().from_mesh(self.mesh)


class GetGhosts(MeshAgent):
    """Ghosts finder.
    """

    def getghosts(self):
        """Returns ghost nodes only.
        """
        return self.getdata()['old-nodes-ghosts']

    def getdata(self):
        """Returns the full output.
        """

        old_nodes_actual, new_triangs_flat = self.learnnodes()
        old_nodes_ghosts = self.findghosts(old_nodes_actual)

        return {
            'old-nodes-actual': old_nodes_actual,
            'old-nodes-ghosts': old_nodes_ghosts,
            'new-triangs-flat': new_triangs_flat
        }

    def learnnodes(self):

        unique_nodes, return_indices = self.run_unique()

        yield unique_nodes.astype(int)
        yield return_indices.astype(int)

    def run_unique(self):
        return np.unique(
            self.mesh.triangs.flat, return_inverse=True
        )

    def findghosts(self, real_nodes):

        all_nodes = np.arange(self.mesh.npoints)

        return np.setdiff1d(
            all_nodes, real_nodes, assume_unique=True
        )


class HasGhosts(MeshAgent):
    """Check for ghosts.
    """

    def hasghosts(self):
        return not self.is_complete()

    def is_complete(self):

        nodes_actual = self.nodes_range_actual()
        nodes_control = self.nodes_range_control()

        if len(nodes_actual) != len(nodes_control):
            return False

        return np.all(
            nodes_actual == nodes_control
        )

    def nodes_range_actual(self):
        return np.unique(self.mesh.triangs).astype(int)

    def nodes_range_control(self):
        return np.arange(self.mesh.npoints)
