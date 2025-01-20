# -*- coding: utf-8 -*-
"""Super voids (only ears).
"""
import numpy as np
from triellipt.utils import pairs


def get_supvoids(mesh):
    """Returns super voids for a mesh.
    """
    return SupVoids.from_mesh(mesh)


class SupVoids:
    """Super-triangulation of voids (only ears).
    """

    def __init__(self, mesh, data):
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_mesh(cls, mesh):
        return SupVoidsMaker(mesh).get_sup_voids()

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
        self.meta = self.fetch_meta()

    def fetch_meta(self):

        voids_meta = self.make_voids_meta()
        pairs_data = self.make_verts_paired()

        return {
            'voids': voids_meta,
            'pairs': pairs_data
        }

    def make_voids_meta(self):

        voids_trinums = self.mesh.getvoids()
        voids_triangs = self.mesh.triangs[voids_trinums, :]
        voids_pivnums = self.mesh.triangs[voids_trinums, 2]

        return {
            'trinums': voids_trinums,
            'triangs': voids_triangs,
            'pivnums': voids_pivnums
        }

    def make_verts_paired(self):
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
    def triangs_verts_paired(self):
        return self.meta['pairs']


class SupVoidsMaker(MeshAgent):
    """Maker of a super-voids.
    """

    def get_sup_voids(self):

        dat = self.make_supvoids_data()
        obj = self.from_supvoids_data(dat)

        return obj

    def from_supvoids_data(self, data):
        return SupVoids(self.mesh, data)

    def make_supvoids_data(self):

        ears = self.find_ears()
        data = self.push_data(ears)

        return data

    def find_ears(self):
        return {
            'west-ears': self.get_west_ears_trinums(),
            'east-ears': self.get_east_ears_trinums()
        }

    def push_data(self, ears_trinums):

        west_trinums = ears_trinums['west-ears']
        east_trinums = ears_trinums['east-ears']

        data = [
            self.voids_trinums, west_trinums, east_trinums
        ]

        return np.vstack(data)

    def get_west_ears_trinums(self):

        codes = pairs.sympaired(
            self.voids_triangs[:, 1], self.voids_triangs[:, 2]
        )

        return self.find_ears_by_pair_codes(codes)

    def get_east_ears_trinums(self):

        codes = pairs.sympaired(
            self.voids_triangs[:, 2], self.voids_triangs[:, 0]
        )

        return self.find_ears_by_pair_codes(codes)

    def find_ears_by_pair_codes(self, ears_codes):

        trinums, _ = _find_in_table(
            self.triangs_verts_paired, ears_codes
        )

        trinums = np.setdiff1d(
            trinums, self.voids_trinums
        )

        permuter = _rowsort_table_by_pivots(
            self.mesh.triangs[trinums, :], self.voids_pivots
        )

        return trinums[permuter]


def _find_in_table(table, vals):

    mask = np.isin(table, vals)

    rows, cols = np.where(mask)
    return rows, cols


def _rowsort_table_by_pivots(table, pivots):

    rows, cols = _find_in_table(table, pivots)

    pivots_in_table = table[rows, cols]

    _, ind1, ind2 = np.intersect1d(
        pivots_in_table, pivots, return_indices=True
    )

    rows_permuter = ind1[np.argsort(ind2)]
    return rows_permuter
