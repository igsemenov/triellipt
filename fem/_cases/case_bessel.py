# -*- coding: utf-8 -*-
"""Validates the Bessel solver.
"""
import numpy as np
from scipy import sparse as sp
from scipy.special import j0
from matplotlib import pyplot as plt
import triellipt as tri


def bessel(unit):
    """Creates the Bessel operator.
    """

    radius = unit.mesh.centrs_complex.imag[unit.ij_t]

    return radius * (
        unit.diff_2x + unit.diff_2y - unit.massmat
    )


def func(argx, argy):
    return j0(argy)


box_partt = {
    'name': 'box',
    'loops-partition': {
        0: {
            "angle": 1.5,
            "coloring": [
                (1, 2, 'rshift'), (3, 4, 'rshift')
            ]
        }
    },
    'dirichlet-sides': (1, 3)
}

seed = tri.mesher.trilattice(31, 41, True) * 0.05
unit = tri.fem.getunit(seed, (0,))

unit.add_partt(box_partt)

L = unit.partts['box'].new_matrix(bessel(unit), constr=True)

u = unit.partts['box'].new_vector()
g = unit.partts['box'].new_vector().from_func(func)

u[0] = sp.linalg.spsolve(
    L(0, 0), - L(0, 1) @ g[1] - L(0, 3) @ g[3]
)

err = np.amax(np.abs(u[0] - g[0]))
print(f'L1 error: {err}')

# plt.tricontourf(*unit.mesh.triu, u.body)
# plt.triplot(*unit.mesh.triu, '-k', lw=0.2)
# plt.axis('equal')
