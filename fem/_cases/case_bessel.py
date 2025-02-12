# -*- coding: utf-8 -*-
"""Validates the Bessel solver.
"""
import numpy as np
from scipy import sparse as sp
from scipy.special import j0
import triellipt as tri


def bessel(unit):
    """Creates the Bessel operator.
    """

    radius = unit.mesh.centrs_complex.imag
    stream = radius[unit.ij_stream.trinums]

    return stream * (
        unit.diff_2x + unit.diff_2y - unit.massmat
    )


def func(x, y):
    return j0(y)


seed = tri.mesher.trilattice(31, 41, True) * 0.05
unit = tri.fem.getunit(seed, (0,))

part = unit.partition(
    0, 1.5, [(1, 2, 'r'), (3, 4, 'r')]
)

part = part.with_dirich_sides(1, 3).asdict()

M = unit.fct0.feed_data(unit.massmat)
L = unit.fct1.feed_data(bessel(unit)).partitioned(part)

u = unit.new_vector().partitioned(part)
g = unit.new_vector().from_func(func).partitioned(part)

rhs = - L(0, 1) @ g[1] - L(0, 3) @ g[3]

u[0] = sp.linalg.spsolve(L(0, 0), rhs)

err = np.amax(np.abs(u[0] - g[0]))
print(f'L1 error: {err}')
