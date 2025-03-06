# -*- coding: utf-8 -*-
"""Validates the Poisson solver.
"""
import time
import numpy as np
from scipy import sparse as sp
from matplotlib import pyplot as plt
import triellipt as tri

META = {
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
    return - 8. * np.sin(2 * x + 2 * y)


class FEMData:
    """Root of the prolem frame.
    """

    def __init__(self, meta):
        self.meta = meta
        self.unit = None

    @classmethod
    def from_meta(cls, meta):
        return cls(meta)

    @property
    def grid_params(self):
        yield self.meta['grid']['size']
        yield self.meta['grid']['size']
        yield self.meta['grid']['mode']

    @property
    def grid_scaling_factor(self):
        return self.meta['grid']['size'] - 1

    @property
    def mesh_refine_depth(self):
        return self.meta['amr']['depth']

    @property
    def mesh_refine_count(self):
        return self.meta['amr']['count']

    def make_seed_mesh(self):
        seed = self.make_grid_unscaled()
        return seed / self.grid_scaling_factor

    def make_grid_unscaled(self):

        grid = tri.mesher.trigrid(*self.grid_params)

        if self.mesh_refine_depth == 0:
            return grid

        for _ in range(self.mesh_refine_depth):
            grid = grid.reduced(self.mesh_refine_count)
        return grid


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

        unit = tri.fem.getunit(
            self.make_seed_mesh()
        )

        return unit

    def get_massmat(self):
        if self.meta['fem']['mass-diag'] is True:
            return self.massdiag_fem
        return self.massmat_fem

    def get_laplace(self):
        return self.unit.base.new_matrix(
            self.unit.diff_2x + self.unit.diff_2y, constr=True
        )

    def solve(self):

        massmat = self.get_massmat()
        laplace = self.get_laplace()

        sol = self.unit.base.new_vector().from_func(sol_exact)
        rho = self.unit.base.new_vector().from_func(rho_exact)

        rhs_rho = massmat @ rho

        tic = time.time()

        sol[0] = sp.linalg.spsolve(
            laplace(0, 0), rhs_rho[0] - laplace(0, 1) @ sol[1]
        )

        toc = time.time()
        cpu = toc - tic

        err = np.amax(
            abs(sol[0] - self.sol_exact_core)
        )

        return sol, err, cpu

    @property
    def sol_exact_core(self):
        return sol_exact(
            *self.nodes2d_core
        )

    @property
    def triu(self):
        return self.unit.mesh.triu

    @property
    def nodes2d_core(self):
        return self.unit.base.nodes2d(0)

    @property
    def nodes2d_edge(self):
        return self.unit.base.nodes2d(1)

    @property
    def massmat_fem(self):
        return self.unit.massopr(lumped=False, constr=False)

    @property
    def massdiag_fem(self):
        return self.unit.massopr(lumped=True, constr=False)


if __name__ == '__main__':

    frame = FEMFrame.from_meta(META).activated()

    sol_, err_, cpu_ = frame.solve()

    print(f'cpu-time: {cpu_}')
    print(f'L1-error: {err_}')

    if META['fem']['with-plot']:
        plt.tricontourf(*frame.triu, sol_.body)
        plt.triplot(*frame.triu, '-k', lw=0.2)
        plt.axis('equal')
