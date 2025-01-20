# -*- coding: utf-8 -*-
"""Validates the basic Poisson solver.
"""
import numpy as np
from scipy import sparse as sp
from matplotlib import pyplot as plt
from triellipt import mesher, fem

META = {
    'grid': {
        'size': 81,
        'mode': 'east-slope'
    },
    'fem': {
        'mass-diag': False,
        'with-plot': True
    },
    'amr': {
        'depth': 2,
        'count': 2
    }
}


def sol_exact(argx, argy):
    return np.sin(2 * argx + 2 * argy)


def rho_exact(argx, argy):
    return 8. * np.sin(2 * argx + 2 * argy)


class FEMFrame:
    """FEM computing frame.

    Attributes
    ----------
    meta : dict
        Problem metadata.
    unit : FEMUnit
        FEM computing unit.

    """

    def __init__(self, meta):
        self.meta = meta
        self.unit = None

    @classmethod
    def from_meta(cls, meta):
        return cls(meta)

    def with_fem_unit(self):

        mesh = self.make_mesh()
        unit = self.make_unit(mesh)

        self.unit = unit
        return self

    @property
    def grid_scale_factor(self):
        return self.meta['grid']['size'] - 1

    @property
    def mesh_refine_depth(self):
        return self.meta['amr']['depth']

    @property
    def mesh_refine_count(self):
        return self.meta['amr']['count']

    def make_unit(self, mesh):
        return fem.getunit(mesh)

    def make_mesh(self):

        grid = self.grid_from_meta()
        mesh = self.mesh_from_grid(grid)

        return mesh

    def mesh_from_grid(self, grid):
        return grid / self.grid_scale_factor

    def grid_from_meta(self):

        meta = self.meta['grid']

        grid = mesher.trigrid(
            meta['size'], meta['size'], meta['mode']
        )

        if self.mesh_refine_depth == 0:
            return grid

        for _ in range(self.mesh_refine_depth):
            grid = grid.reduced(self.mesh_refine_count)

        return grid

    def get_massmat(self):
        if self.meta['fem']['mass-diag']:
            return self.unit.massdiag_fem
        return self.unit.massmat_fem

    def get_laplace(self):
        return self.unit.laplace_fem


def solve(meta):

    frame = FEMFrame.from_meta(META).with_fem_unit()

    massmat = frame.get_massmat()
    laplace = frame.get_laplace().dirichsplit()

    mat22 = laplace.getblock('core', 'core')
    mat21 = laplace.getblock('core', 'edge')

    sol = frame.unit.new_vector()
    rho = frame.unit.new_vector().from_func(rho_exact)

    sol = sol.dirichsplit()
    rho = rho.dirichsplit()

    sol.setsection(
        'edge', sol_exact(*sol.sectionxy('edge'))
    )

    rhs = massmat @ rho

    rhs.setsection(
        'core', rhs.getsection('core') - mat21 @ sol.getsection('edge')
    )

    out = sp.linalg.spsolve(
        mat22, rhs.getsection('core')
    )

    sol.setsection('core', out)

    res = out - sol_exact(*sol.sectionxy('core'))
    err = np.amax(abs(res))

    return (
        sol, err, frame.unit.mesh, frame.unit
    )


if __name__ == '__main__':

    sol, err, mesh, unit = solve(META)

    if META['fem']['with-plot']:
        plt.tricontourf(*mesh.triu, sol.body)
        plt.triplot(*mesh.triu, '-k', lw=0.2)
        plt.axis('equal')

    print(f'L1-error: {err}')
