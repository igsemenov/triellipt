# -*- coding: utf-8 -*-
"""Tests the mesh coarsening.
"""
import unittest
import numpy as np
from triellipt import mesher
from triellipt.amr import tricoarsen


class TestCoarsener(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.CRS = tricoarsen.MeshCoarsener.from_mesh(cls.mesh())

    @classmethod
    def mesh(cls):

        mesh = mesher.trigrid(4, 4, 'east-slope')

        return mesh.shuffled(
            np.argsort(mesh.centrs_complex)
        )

    def test_target_trinums(self):
        assert self.crs_bad.target_trinums.size == 2
        assert self.crs_good.target_trinums.size == 2
        assert self.crs_empty.target_trinums.size == 0

    def test_target_suptri(self):
        assert self.crs_bad.make_target_suptri() is None
        assert self.crs_empty.make_target_suptri() is None
        assert self.crs_good.make_target_suptri().size == 2

    def test_mesh_alpha(self):
        assert self.mesh_alpha.size == 16
        assert self.mesh_alpha.getvoids().size == 6
        assert self.mesh_alpha.getghosts().size == 0

    def test_mesh_beta(self):
        assert self.mesh_beta.size == 12
        assert self.mesh_beta.getvoids().size == 2
        assert self.mesh_beta.getghosts().size == 4

    def test_mesh_gamma(self):
        assert self.mesh_gamma.size == 14
        assert self.mesh_gamma.getvoids().size == 2
        assert self.mesh_gamma.getghosts().size == 0

    @property
    def mesh_alpha(self):
        return self.crs_good.make_mesh_alpha()

    @property
    def mesh_beta(self):
        return self.crs_good.make_mesh_beta()

    @property
    def mesh_gamma(self):
        return self.crs_good.make_mesh_gamma()

    @property
    def crs_bad(self):
        """Coarsener with bad target triangles.
        """
        return self.CRS.with_trinums((7, 10))

    @property
    def crs_good(self):
        """Coarsener with good target triangles.
        """
        return self.CRS.with_trinums((5, 7))

    @property
    def crs_empty(self):
        """Coarsener with no target triangles.
        """
        return self.CRS.with_trinums(())


if __name__ == '__main__':
    unittest.main()
