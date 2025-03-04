# -*- coding: utf-8 -*-
"""Validates the Bessel solver.
"""
import numpy as np
from scipy import sparse as sp
from scipy.special import j0
from matplotlib import pyplot as plt
import triellipt as tri


def bessel_operator(unit):
    """Creates the Bessel operator.
    """

    radius = unit.mesh.centrs_complex.imag[unit.ij_t]

    return radius * (
        unit.diff_2x + unit.diff_2y - unit.massmat
    )


def bessel_func(argx, argy):
    return j0(argy)


box_partt = {
    'partition-title': 'box',
    'partition-loops': {
        0: {
            "angle": 1.5,
            "coloring": [
                (1, 2, 'rshift'), (3, 4, 'rshift')
            ]
        }
    },
    'dirichlet-sides': [1, 3]
}

seed = tri.mesher.trilattice(31, 41, True) * 0.05
unit = tri.fem.getunit(seed, (0,))

unit.add_partition(box_partt)

mat = unit.partts['box'].new_matrix(
    bessel_operator(unit), constr=False
)

sol = unit.partts['box'].new_vector()
ref = unit.partts['box'].new_vector().from_func(bessel_func)

sol[0] = sp.linalg.spsolve(
    mat(0, 0), - mat(0, 1) @ ref[1] - mat(0, 3) @ ref[3]
)

err = np.amax(np.abs(sol[0] - ref[0]))
print(f'L1 error: {err}')

# plt.tricontourf(*unit.mesh.triu, sol.body)
# plt.triplot(*unit.mesh.triu, '-k', lw=0.2)
# plt.axis('equal')
