# -*- coding: utf-8 -*-
"""Super voids (only ears).
"""
import numpy as np
from triellipt.utils import pairs


def get_sup_voids(mesh):
    """Returns super voids for a mesh.
    """
    return SupVoids.from_mesh(mesh)


class SupVoids:
    """Proxy of a voids super-triangulation (only ears).
    """

    def __init__(self, mesh, data):
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_mesh(cls, mesh):
        _ = SupVoidsMaker(mesh).with_meta()
        return _.get_sup_voids()

    @property
    def trinums(self):
        """Numbers of voids triangles.
        """
        return self.data[0, :]

    @property
    def trinums2(self):
        """Numbers of west-ears.
        """
        return self.data[1, :]

    @property
    def trinums3(self):
        """Numbers of east-ears.
        """
        return self.data[2, :]

    @property
    def kermesh(self):
        return self.mesh.update_triangs(self.voids_triangs)

    @property
    def voids_triangs(self):
        return self.mesh.triangs[self.trinums, :]


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = None

    def with_meta(self):
        self.meta = self.make_meta()
        return self

    def make_meta(self):

        voids_meta = self.make_meta_voids()
        pairs_data = self.make_data_nodes_paired()

        return {
            'voids': voids_meta,
            'pairs': pairs_data
        }

    def make_meta_voids(self):

        voids_trinums = self.mesh.getvoids()
        voids_triangs = self.mesh.triangs[voids_trinums, :]
        voids_pivnums = self.mesh.triangs[voids_trinums, 2]

        return {
            'trinums': voids_trinums,
            'triangs': voids_triangs,
            'pivnums': voids_pivnums
        }

    def make_data_nodes_paired(self):
        return pairs.paircols(
            self.mesh.triangs[:, [0, 1, 2, 0]]
        )

    @property
    def voids_trinums(self):
        return self.meta['voids'].get('trinums')

    @property
    def voids_triangs(self):
        return self.meta['voids'].get('triangs')

    @property
    def voids_pivots(self):
        return self.meta['voids'].get('pivnums')

    @property
    def mesh_nodes_paired(self):
        return self.meta['pairs']


class SupVoidsMaker(MeshAgent):
    """Maker of a super-voids.
    """

    def get_sup_voids(self):
        dat = self.make_data()
        obj = self.from_data(dat)
        return obj

    def make_data(self):

        west_trinums = self.get_west_trinums()
        east_trinums = self.get_east_trinums()

        data = [
            self.voids_trinums, west_trinums, east_trinums
        ]

        return np.vstack(data)

    def from_data(self, data):
        return SupVoids(self.mesh, data)

    def get_west_trinums(self):

        targets = pairs.sympaired(
            self.voids_triangs[:, 1], self.voids_triangs[:, 2]
        )

        return self.find_targets_trinums(targets)

    def get_east_trinums(self):

        targets = pairs.sympaired(
            self.voids_triangs[:, 2], self.voids_triangs[:, 0]
        )

        return self.find_targets_trinums(targets)

    def find_targets_trinums(self, targetpairs):

        trinums = _find_targets_in_table(
            targetpairs, self.mesh_nodes_paired
        )

        trinums = np.setdiff1d(
            trinums, self.voids_trinums
        )

        permuter = _sync_table_to_pivots(
            self.mesh.triangs[trinums, :], self.voids_pivots
        )

        return trinums[permuter]


def _find_targets_in_table(targets, table):

    mask = np.isin(table, targets)

    rows, = np.where(
        np.sum(mask, axis=1)
    )

    return rows


def _sync_table_to_pivots(table, pivots):

    pivots_in_table = table[
        *np.where(np.isin(table, pivots))
    ]

    _, ind1, ind2 = np.intersect1d(
        pivots_in_table, pivots, return_indices=True
    )

    rows_permuter = ind1[np.argsort(ind2)]
    return rows_permuter
