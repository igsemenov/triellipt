# -*- coding: utf-8 -*-
"""Cleaner of a supertriu.
"""
from abc import ABC, abstractmethod
import numpy as np
from triellipt.utils import tables


def clean_overlaps(suptri):
    """Removes overlapping super-triangles.
    """
    _ = CleanOverlaps.from_suptri(suptri)
    return _.suptri_cleaned()


def clean_not_aligned(suptri):
    """Removes super-triangles that are not aligned.
    """
    _ = CleanNonAligned.from_suptri(suptri)
    return _.suptri_cleaned()


class SuptriAgent(ABC):
    """Operator on a supertriu.
    """

    def __init__(self, suptri):
        self.suptri = suptri
        self.meta = self.fetch_meta()

    @classmethod
    def from_suptri(cls, suptri):
        return cls(suptri)

    def fetch_meta(self):
        return {}


class SuptriCleaner(SuptriAgent):
    """Base suptri cleaner.
    """

    def suptri_cleaned(self):

        if self.suptri.size == 0:
            return self.suptri

        bad_indices = self.find_suptrinds_bad()
        new_suptriu = self.pull_suptrinds_bad(bad_indices)

        return new_suptriu

    def pull_suptrinds_bad(self, bad_indices):

        if bad_indices.size == 0:
            return self.suptri

        new_suptriu = self.suptri.deltriangs(*bad_indices)
        return new_suptriu

    @abstractmethod
    def find_suptrinds_bad(self):
        """Finds indices of overlapping super-triangles.
        """


class CleanOverlaps(SuptriCleaner):
    """Removes overlapping super-triangles.
    """

    def fetch_meta(self):
        return {
            'supbodies': self.suptri.supbodies
        }

    @property
    def supbodies(self):
        return self.meta['supbodies']

    def find_suptrinds_bad(self):

        data = self.make_supbodies_map()

        long_ranks = np.repeat(
            data.vals_ranks, data.vals_ranks
        )

        bad_inds = np.unique(
            data.rows[long_ranks != 1]
        )

        return bad_inds

    def make_supbodies_map(self):
        return tables.maptable(self.supbodies)


class CleanNonAligned(SuptriCleaner):
    """Removes non-aligned super.triangles.
    """

    def find_suptrinds_bad(self):

        mask = self.mask_not_aligned()

        inds, = np.where(
            np.any(mask, axis=1)
        )

        return inds

    def mask_not_aligned(self):

        ker = self.get_ker_verts()
        sup = self.get_sup_verts()

        mask0 = 0.5 * (sup[0] + sup[1]) == ker[1]
        mask1 = 0.5 * (sup[1] + sup[2]) == ker[2]
        mask2 = 0.5 * (sup[2] + sup[0]) == ker[0]

        mask = _stack_cols(
            mask0, mask1, mask2
        )

        return np.logical_not(mask)

    def get_ker_verts(self):
        return _mesh_verts(self.suptri.kermesh)

    def get_sup_verts(self):
        return _mesh_verts(self.suptri.supmesh)


def _mesh_verts(mesh):
    return np.hsplit(
        mesh.points[mesh.triangs], 3
    )


def _stack_cols(*cols):
    return np.hstack(cols)
