# -*- coding: utf-8 -*-
"""AMR on a circle.
"""
import numpy as np
import matplotlib.pyplot as plt
import triellipt as tri


def layer(x, y):
    r = np.abs(x + 1j * y)
    return np.tanh(5.0 * (1.0 - r))


def gauss(x, y):
    return np.exp(- 2*x*x - 2*y*y)


def sinwave(x, y):
    return np.sin(x + y)


seed = tri.trimesh.TriMesh.load('vmesh-circ.npz')
unit = tri.amr.getunit(seed)
unit.data[0] = unit.atfunc(layer)

f0 = unit.data[0]
m0 = unit.masser(f0)

for t in range(2):

    for _ in range(4):
        unit = unit.refine(
            unit.front_coarse().trinums
        )

    for _ in range(4):
        unit = unit.coarsen(
            unit.front_fine().trinums
        )

f1 = unit.data[0]
m1 = unit.masser(f1)

print(f'Mass-error: {m0-m1}')

# plt.tricontourf(*unit.mesh.triu, unit.data[0])
# plt.triplot(*unit.mesh.triu, '-k', lw=0.5)
# plt.axis('equal')
