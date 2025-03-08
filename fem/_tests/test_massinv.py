# -*- coding: utf-8 -*-
"""Tests the mass-inverse operator.
"""
import unittest
import numpy as np
from scipy import sparse as sp
from triellipt import mesher
from triellipt import fem


def sin_wave(argx, argy):
    return np.sin(argx + argy)


class Tester(unittest.TestCase):

    UNIT = None

    @property
    def massdiag_lu(self):

        lu = sp.linalg.splu(
            self.massdiag_amr.body.tocsc()
        )

        return lu.solve

    @property
    def massdiag_inv(self):
        return self.UNIT.massinv()

    @property
    def massdiag_fem(self):
        return self.UNIT.massopr(
            is_lumped=True, add_constr=False
        )

    @property
    def massdiag_amr(self):
        return self.UNIT.massopr(
            is_lumped=True, add_constr=True
        )

    @property
    def rho(self):
        return self.massdiag_fem.body @ self.sinwave

    @property
    def sinwave(self):
        return sin_wave(
            *self.UNIT.mesh.points2d
        )

    def allclose(self, arr1, arr2):
        return np.allclose(
            arr1, arr2, rtol=0, atol=1e-10
        )


class TestNoConstr(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = fem.getunit(cls.mesh())

    @classmethod
    def mesh(cls):
        return mesher.trigrid(11, 11, 'east-slope')

    def test_massinv(self):
        assert self.allclose(
            self.massdiag_inv(self.rho), self.massdiag_lu(self.rho)
        )

    def test_has_constraints(self):
        assert self.massdiag_inv.has_constraints is False


class TestWithConstr(Tester):

    @classmethod
    def setUpClass(cls):
        cls.UNIT = fem.getunit(cls.mesh())

    @classmethod
    def mesh(cls):
        return mesher.trigrid(11, 11, 'east-slope').reduced(1)

    def test_massinv(self):
        assert self.allclose(
            self.massdiag_inv(self.rho), self.massdiag_lu(self.rho)
        )

    def test_has_constraints(self):
        assert self.massdiag_inv.has_constraints is True


if __name__ == '__main__':
    unittest.main()
