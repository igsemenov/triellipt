# -*- coding: utf-8 -*-
"""Tests cleaning of super-triangles.
"""
import unittest
import numpy as np
from triellipt import mesher
from triellipt.amr import supclean


class TestCleanOverlaps(unittest.TestCase):

    @classmethod
    def mesh(cls):

        mesh = mesher.trigrid(4, 4, 'east-slope')

        return mesh.shuffled(
            np.argsort(mesh.centrs_complex)
        )

    def test_clean_overlaps_good(self):
        assert self.suptri_good.size == 2
        assert self.suptri_good.trinums.tolist() == [5, 7]

    def test_clean_overlaps_bad(self):
        assert self.suptri_bad.size == 0

    @property
    def suptri_good(self):
        return supclean.clean_overlaps(
            self.suptri_prime.atcores(5, 7)
        )

    @property
    def suptri_bad(self):
        return supclean.clean_overlaps(
            self.suptri_prime.atcores(7, 10)
        )

    @property
    def suptri_prime(self):
        return self.mesh().supertriu()


class TestCleanNotAligned(unittest.TestCase):

    @classmethod
    def mesh(cls):

        mesh = mesher.trigrid(4, 4, 'east-slope')
        mesh.points[9] += 1e-4

        return mesh.shuffled(
            np.argsort(mesh.centrs_complex)
        )

    def test_clean_not_alignedd(self):
        assert self.suptri_cleaned.size == 1
        assert self.suptri_cleaned.trinums.tolist() == [5]

    @property
    def suptri_cleaned(self):
        return supclean.clean_not_aligned(
            self.suptri_prime.atcores(5, 7)
        )

    @property
    def suptri_prime(self):
        return self.mesh().supertriu()


if __name__ == '__main__':
    unittest.main()
