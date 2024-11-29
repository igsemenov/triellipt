# -*- coding: utf-8 -*-
"""Collector dumper.
"""
import numpy as np
from scipy import sparse as sp


class DumpCollect:
    """Collector dumper.
    """

    def __init__(self, collector):
        self.coll = collector
        self.cache = {}

    def get_matrix(self):

        _ = self.set_cache()
        _ = self.run_dumper()

        return _

    def set_cache(self):
        self.make_reinterp()
        self.make_permutations()
        return True

    def run_dumper(self):

        tmat = self.cache['re-interp']

        r2u = self.cache['perms']['root2unit']
        u2m = self.cache['perms']['unit2mass']
        m2r = self.cache['perms']['mass2root']
        r2d = self.cache['perms']['root2data']

        tmat = tmat @ r2u
        tmat = u2m @ tmat

        tmat = tmat[m2r, :]
        tmat = tmat[r2d, :]

        return tmat.copy()

    def make_reinterp(self):

        rhs = self.coll.massmat.todense()

        _ = sp.csr_array(
            self.coll.massinv.solve(rhs)
        )

        self.cache['re-interp'] = _
        return _

    def make_permutations(self):

        perms = {
            'root2unit': self.make_perm_root2unit(),
            'unit2mass': self.make_perm_unit2mass(),
            'mass2root': self.make_perm_mass2root(),
            'root2data': self.make_perm_root2data()
        }

        self.cache['perms'] = perms
        return perms

    def make_perm_root2unit(self):
        return self.make_perm_matrix(
            self.coll.root2mass[self.coll.mass2unit]
        )

    def make_perm_unit2mass(self):
        return self.make_perm_matrix(self.coll.unit2mass)

    def make_perm_mass2root(self):
        return np.unique(self.coll.root2mass, return_index=True)[1]

    def make_perm_root2data(self):
        return self.coll.root2data

    def make_perm_matrix(self, inds):

        cols = inds
        rows = np.arange(inds.size)

        coo = sp.coo_array(
            (np.ones_like(inds), (rows, cols))
        )

        return coo.tocsr()

    def make_perm_inv(self, perm):
        return np.argsort(perm)
