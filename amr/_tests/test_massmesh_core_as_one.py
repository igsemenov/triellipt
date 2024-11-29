# -*- coding: utf-8 -*-
"""Tests the mass-mesh (single-core).

MASS-MESH:
    
    2
    |\
    | \
    4--3
    |\ |\
    | \| \
    0--5--1

"""
import unittest
import numpy as np
import triellipt as tri
from triellipt.amr import massmesh

POINTS = [
    0 + 0j, 1 + 0j, 2 + 0j, 1 + 1j, 0 + 2j, 0 + 1j
]

TRIANGS = [
    [0, 1, 5],
    [1, 2, 3],
    [5, 3, 4],
    [1, 3, 5]
]


def seed_mesh_to_test():
    """Triangle with the area 2.
    """

    mesh = tri.trimesh.TriMesh.from_data(
        np.array(POINTS), np.array(TRIANGS)
    )

    return mesh.alignnodes((0, 0))


class TestMassMesher(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        mass_mesh = massmesh.get_massmesh(
            cls.mesh().supertriu().subtriu(0)
        )

        mass_unit = tri.fem.getunit(
            mass_mesh, anchors=[(0, 0)]
        )

        mass_opr = mass_unit.massopr(
            is_lumped=True, add_constr=False
        )

        mass_amr = mass_unit.massopr(
            is_lumped=True, add_constr=True
        )

        perm = mass_unit.perm.perm_inv

        mass_opr = mass_opr.body.todense()
        mass_amr = mass_amr.body.todense()

        mass_opr = mass_opr[:, perm][perm, :]
        mass_amr = mass_amr[:, perm][perm, :]

        cls.MASS_MESH = mass_mesh
        cls.MASS_MESH_VOIDS = mass_mesh.triangs[mass_mesh.getvoids(), :]

        cls.MASS_OPR = mass_opr
        cls.MASS_INV = np.linalg.solve(mass_amr, np.eye(6))

    @classmethod
    def mesh(cls):
        return seed_mesh_to_test()

    def test_mesh_count(self):
        assert self.MASS_MESH.npoints == 6

    def test_voids_count(self):
        assert self.MASS_MESH.getvoids().size == 3

    def test_voids_triangs(self):
        assert [1, 2, 3] in self.MASS_MESH_VOIDS
        assert [0, 1, 5] in self.MASS_MESH_VOIDS
        assert [2, 0, 4] in self.MASS_MESH_VOIDS

    def test_mass_opr(self):
        assert self.MASS_OPR[0, [1, 2, 3]].tolist() == [0, 0, 0]
        assert self.MASS_OPR[0, [0, 4, 5]].tolist() == [1/6, 1/4, 1/4]

    def test_mass_uniform(self):
        assert np.allclose(
            self.new_data_uniform, self.data_uniform, rtol=0, atol=1e-14
        )

    def test_mass_linear(self):
        assert np.allclose(
            self.new_data_linear, self.data_linear, rtol=0, atol=1e-14
        )

    @property
    def new_data_uniform(self):
        return self.MASS_INV @ (
            self.MASS_OPR @ self.data_uniform
        )

    @property
    def new_data_linear(self):
        return self.MASS_INV @ (
            self.MASS_OPR @ self.data_linear
        )

    @property
    def data_uniform(self):
        return np.ones(6)

    @property
    def data_linear(self):
        return 1 + np.r_[
            1, 0, 0, 0, 0.5, 0.5
        ]


if __name__ == '__main__':
    unittest.main()
