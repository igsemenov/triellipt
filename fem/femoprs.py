# -*- coding: utf-8 -*-
"""Local FEM operators.
"""
import numpy as np


FEMOPRS = [
    'massmat',
    'massdiag',
    'diff_1x',
    'diff_1y',
    'diff_2x',
    'diff_2y'
]


def getoprs(mesh) -> dict:
    """Returns basic FEM operators.
    """

    maker = OprsMaker.from_metric(mesh_metric(mesh))

    proxy = {
        k: getattr(maker, k)() for k in FEMOPRS
    }

    return {
        k: np.hsplit(v, 3) for k, v in proxy.items()
    }


def mesh_metric(mesh) -> dict:
    """Returns the mesh metric properties.

    Parameters
    ----------
    mesh : TriMesh
        Triangular mesh.

    Returns
    -------
    dict
        Metric properties of triangles.

    """
    _ = MetricMaker.from_mesh(mesh)
    return _.get_metric()


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def vertices(self):
        return self.mesh.points[self.mesh.triangs]


class MetricMaker(MeshAgent):
    """Computes mesh metric properties.
    """

    PERM_DELTAS = [1, 2, 0]

    def get_metric(self):

        coeffs = self.get_coeffs()
        jacobs = self.get_jacobs(coeffs)

        return {
            **coeffs,
            **jacobs
        }

    def get_coeffs(self):

        deltas = np.diff(
            self.vertices[:, [0, 1, 2, 0]], axis=1
        )

        deltas = np.conj(deltas)

        deltas = deltas[
            :, self.PERM_DELTAS
        ]

        return {
            'bcoeffs': deltas.imag.copy('C'),
            'ccoeffs': deltas.real.copy('C')
        }

    def get_jacobs(self, coeffs):

        bcf = coeffs['bcoeffs']
        ccf = coeffs['ccoeffs']

        return {
            'jacobis': bcf[:, 0] * ccf[:, 1] - bcf[:, 1] * ccf[:, 0]
        }


class OprsMaker:
    """Maker of FEM operators.
    """

    def __init__(self, metric):
        self.metric = metric

    @classmethod
    def from_metric(cls, metric):
        return cls(metric)

    @property
    def areas2d(self):
        return 0.5 * self.jacobis2d.T

    @property
    def jacobis2d(self):
        return np.atleast_2d(
            self.metric['jacobis']
        )

    def diff_1x(self):
        return self.diff_1d('bcoeffs')

    def diff_1y(self):
        return self.diff_1d('ccoeffs')

    def diff_2x(self):
        return self.diff_2d('bcoeffs')

    def diff_2y(self):
        return self.diff_2d('ccoeffs')

    def diff_1d(self, key):
        coeffs = self.metric[key]
        return _mono_matrix(coeffs) / 6.

    def diff_2d(self, key):

        coeffs = self.metric[key]

        return 0.25 * (
            _diad_matrix(coeffs) / self.areas2d
        )

    def massmat(self):

        areas = self.areas2d

        proxy = np.tile(
            self.massmat_proxy.flat, areas.shape
        )

        return proxy * (areas / 12.)

    def massdiag(self):

        areas = self.areas2d

        proxy = np.tile(
            self.massdiag_proxy.flat, areas.shape
        )

        return proxy * areas

    @property
    def massmat_proxy(self):
        return np.eye(3) + np.ones((3, 3))

    @property
    def massdiag_proxy(self):
        return np.eye(3) / 3.


def _mono_matrix(data):
    return np.repeat(data, 3, axis=1)


def _diad_matrix(data):
    return np.tile(data, (1, 3)) * np.repeat(data, 3, axis=1)
