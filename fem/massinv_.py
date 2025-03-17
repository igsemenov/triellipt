# -*- coding: utf-8 -*-
"""Lumped mass inverse.
"""
from scipy import sparse as sp
import numpy as np


def getmassinv(unit):
    maker = MakerMassInv(unit)
    return maker.get_massinv()


class MassDiagInv:
    """Inverse of a lumped-mass matrix.
    """

    def __init__(self, unit=None, meta=None):
        self.unit = unit
        self.meta = meta
        self.cache = {}

    @property
    def has_constraints(self):
        return self.meta['has-constraints']

    def __call__(self, data):
        return self.solve(data)

    def solve(self, data):
        if self.has_constraints is False:
            return self.solve_no_constraints(data)
        return self.solve_with_constraints(data)

    def solve_with_constraints(self, data):

        self.perm_data(data)
        self.find_solution()
        self.perm_solution()

        _ = self.push_solution()
        return _

    def push_solution(self):

        _ = self.cache['sol-permuted']
        self.cache.clear()

        return np.copy(
            _, order='C'
        )

    def perm_data(self, data):
        self.cache['data-permuted'] = data[self.meta['perm']]

    def find_solution(self):

        self.solve_mass_diag_only()
        self.solve_mass_with_constr()

        self.cache['sol'] = np.hstack(
            [self.cache['sol-diag'], self.cache['sol-constr']]
        )

    def perm_solution(self):
        _ = self.cache['sol'][self.meta['perm-inv']]
        self.cache['sol-permuted'] = _

    def solve_mass_diag_only(self):

        mass = self.meta['mass-diag']
        data = self.cache['data-permuted'][:self.diag_size]

        self.cache['sol-diag'] = data / mass

    def solve_mass_with_constr(self):

        mass = self.meta['constr-lu']
        data = self.cache['data-permuted'][self.diag_size:]

        self.cache['sol-constr'] = mass.solve(data)

    def solve_no_constraints(self, data):
        return np.copy(
            data / self.meta['mass-diag'], order='C'
        )

    @property
    def diag_size(self):
        return self.meta['mass-diag'].size


class MakerMassInv:
    """Maker of the mass-inverse operator.
    """

    def __init__(self, unit):
        self.unit = unit
        self.meta = self.fetch_meta()
        self.cache = {}

    def fetch_meta(self):
        return {
            "mass-diag": self.fetch_meta_massdiag(),
        }

    def fetch_meta_massdiag(self):

        massdiag = self.unit.massopr(
            is_lumped=True, add_constr=True
        )

        return {
            'body': massdiag.body,
            'pattern': massdiag.pattern,
            'hasconstr': massdiag.has_constraints
        }

    @property
    def massdiag_body(self):
        return self.meta['mass-diag']['body']

    @property
    def massdiag_pattern(self):
        return self.meta['mass-diag']['pattern']

    @property
    def has_constraints(self):
        return self.meta['mass-diag']['hasconstr']

    @property
    def meta_no_constraints(self):
        return {
            'mass-diag': np.copy(
                self.massdiag_body.data[:], order='C'
            ),
            'has-constraints': False
        }

    def get_massinv(self):

        if self.has_constraints is False:
            return MassDiagInv(
                self.unit, self.meta_no_constraints
            )

        _ = self.make_constr_lu()
        _ = self.make_massinv()

        return _

    def make_massinv(self):
        return MassDiagInv(
            self.unit, self.make_massinv_meta()
        )

    def make_massinv_meta(self):
        return {
            'perm': self.perm,
            'perm-inv': self.perm_inv,
            'mass-diag': self.mass_diag_only,
            'constr-lu': self.mass_with_constr_lu,
            'has-constraints': True
        }

    def make_constr_lu(self):

        self.find_perm_meta()
        self.make_mass_body_permuted()

        lu_ = sp.linalg.splu(
            self.mass_with_constr
        )

        self.cache['constr-lu'] = lu_
        return lu_

    def make_mass_body_permuted(self):

        _ = self.massdiag_body[self.perm, :][:, self.perm]

        self.cache['mass-body-permuted'] = _
        return _

    @property
    def mass_permuted(self):
        return self.cache['mass-body-permuted']

    @property
    def mass_diag_only(self):
        return self.mass_permuted.diagonal()[:self.dsize]

    @property
    def mass_with_constr(self):
        return self.mass_permuted[self.dsize:, :][:, self.dsize:].tocsc()

    def find_perm_meta(self):

        mask = self.make_pattern_mask()

        meta = {
            'diagsize': np.sum(mask == 1),
            'permuter': np.argsort(mask)
        }

        self.cache['perm-meta'] = meta
        return meta

    @property
    def perm(self):
        return self.cache['perm-meta']['permuter']

    @property
    def dsize(self):
        return self.cache['perm-meta']['diagsize']

    @property
    def perm_inv(self):
        return np.argsort(self.perm)

    def make_pattern_mask(self):

        mask = np.sum(
            self.massdiag_pattern, axis=1
        )

        return mask.astype(int)

    @property
    def mass_with_constr_lu(self):
        return self.cache['constr-lu']
