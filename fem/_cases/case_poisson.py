# -*- coding: utf-8 -*-
"""Validates the Poisson solver.
"""
import numpy as np
from scipy import sparse as sp
from matplotlib import pyplot as plt
import triellipt as tri

META = {
    'mode': 'fem',
    'grid': {
        'size': 81,
        'mode': 'east-slope'
    },
    'fem': {
        'mass-diag': False,
        'with-plot': False
    },
    'amr': {
        'depth': 2,
        'count': 2
    }
}


def sol_exact(x, y):
    return np.sin(2 * x + 2 * y)


def rho_exact(x, y):
    return - 8 * np.sin(2 * x + 2 * y)


class FEMData:
    """Root of the problem frame.
    """

    def __init__(self, meta):
        self.meta = meta
        self.unit = None

    @classmethod
    def from_meta(cls, meta):
        return cls(meta)

    def make_grid(self):
        return self.make_grid_unscaled() / self.grid_scaling_factor

    def make_grid_unscaled(self):

        size = self.meta['grid']['size']
        mode = self.meta['grid']['mode']

        grid = tri.mesher.trigrid(size, size, mode)

        if self.mesh_refine_depth == 0:
            return grid

        for _ in range(self.mesh_refine_depth):
            grid = grid.reduced(self.mesh_refine_count)
        return grid

    @property
    def grid_scaling_factor(self):
        return self.meta['grid']['size'] - 1

    @property
    def mesh_refine_depth(self):
        return self.meta['amr']['depth']

    @property
    def mesh_refine_count(self):
        return self.meta['amr']['count']


class FEMFrame(FEMData):
    """Problem frame.

    Parameters
    ----------
    meta : dict
        Problem metadata.

    """

    def activated(self):
        """Activates the frame before use.
        """
        self.unit = self.get_fem_unit()
        return self

    def get_fem_unit(self):
        return tri.fem.getunit(
            self.make_grid(), anchors=((0, 0),), mode=self.meta['mode']
        )

    def get_massmat(self):
        if self.meta['fem']['mass-diag'] is True:
            return self.massdiag_fem
        return self.massmat_fem

    def get_laplace(self):
        return self.unit.base.new_matrix(
            self.unit.diff_2x + self.unit.diff_2y, add_constr=True
        )

    def solve(self):

        massmat = self.get_massmat()
        laplace = self.get_laplace()

        ref = self.unit.base.new_vector().from_func(sol_exact)
        sol = self.unit.base.new_vector().from_func(sol_exact)
        rho = self.unit.base.new_vector().from_func(rho_exact)

        rhs = massmat @ rho

        sol[0] = sp.linalg.spsolve(
            laplace(0, 0), rhs[0] - laplace(0, 1) @ sol[1]
        )

        err = np.amax(
            abs(sol[0] - ref[0])
        )

        return {
            'sol': sol,
            'err': err
        }

    @property
    def triu(self):
        return self.unit.mesh.triu

    @property
    def massmat_fem(self):
        return self.unit.massopr(
            is_lumped=False, add_constr=False
        )

    @property
    def massdiag_fem(self):
        return self.unit.massopr(
            is_lumped=True, add_constr=False
        )


if __name__ == '__main__':

    frame = FEMFrame.from_meta(META).activated()
    out = frame.solve()

    print(f'err-norm: {out["err"]}')

    if META['fem']['with-plot'] is True:

        plt.tricontourf(*frame.triu, out['sol'].body)
        plt.triplot(*frame.triu, '-k', lw=0.2)
        plt.axis('equal')
