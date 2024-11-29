# -*- coding: utf-8 -*-
"""Tests the mass-mesh.
"""
import unittest
import numpy as np
from triellipt import mesher
from triellipt.amr import massmesh


def mesh_to_test():

    mesh = mesher.trigrid(3, 5, 'east-slope')

    mesh = mesh.shuffled(
        np.argsort(mesh.centrs_complex)
    )

    supt = mesh.supertriu().atcores(5)

    mesh = mesh.deltriangs(supt.supbodies)
    mesh = mesh.add_triangs(supt.supmesh.triangs)

    mesh = mesh.add_triangs([12, 0, 6])
    mesh = mesh.add_triangs([12, 2, 7])
    mesh = mesh.delghosts()

    mesh = mesh.shuffled(
        np.argsort(1j * mesh.centrs_complex)
    )

    return mesh


class TestMassMesher(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.SUPTRI = cls.mesh().supertriu().atcores(5, 12)

    @classmethod
    def mesh(cls):
        return mesh_to_test()

    def test_triangs_alpha(self):
        assert self.triangs_alpha.shape == (6, 3)

    def test_triangs_beta(self):
        assert self.triangs_beta['inner-voids'].shape == (6, 3)
        assert self.triangs_beta['inner-cores'].shape == (2, 3)
        assert self.triangs_beta['inner-sides'].shape == (6, 3)

    def test_triangs_beta_extra_voids(self):
        assert self.triangs_beta_extra_trinums.tolist() == [8, 11]
        assert self.triangs_beta_extra_triangs.tolist() == [[1, 13, 7]]

    def test_mesh_gamma(self):
        assert self.mesh_gamma.size == 20
        assert self.mesh_gamma.getvoids().size == 7
        assert self.mesh_gamma.getghosts().tolist() == [4, 5, 6, 10, 12]

    def test_nodes_root2gamma(self):
        assert np.all(
            self.points_root2gamma == self.mesh_gamma.points
        )

    def test_nodes_gamma2mass(self):
        assert np.all(
            self.points_gamma2mass == self.mesh_release.points
        )

    @property
    def triangs_alpha(self):
        return massmesh.get_triangs_alpha(self.SUPTRI)

    @property
    def triangs_beta(self):
        return massmesh.get_triangs_beta(self.SUPTRI)

    @property
    def triangs_beta_extra_trinums(self):
        return self.triangs_beta['extra-voids']['trinums-to-del']

    @property
    def triangs_beta_extra_triangs(self):
        return self.triangs_beta['extra-voids']['triangs-to-add']

    @property
    def mesh_gamma(self):
        return massmesh.get_mesh_gamma(self.SUPTRI)

    @property
    def mesh_release(self):
        return massmesh.get_massmesh(self.SUPTRI)

    @property
    def points_root2gamma(self):
        return self.SUPTRI.mesh.points[
            self.mesh_release.meta['nodes-root2gamma']
        ]

    @property
    def points_gamma2mass(self):
        return self.mesh_gamma.points[
            self.mesh_release.meta['nodes-gamma2mass']
        ]


if __name__ == '__main__':
    unittest.main()
