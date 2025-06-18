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
    'diff_2y',
    'grad_1y'
]


def getoprs(mesh) -> dict:
    """Returns basic FEM operators.
    """

    maker = OprsMaker.from_metric(mesh_metric(mesh))

    return {
        k: getattr(maker, k)() for k in FEMOPRS
    }


def mesh_grad(mesh):
    """Returns the mesh gradient operator.

    Parameters
    ----------
    mesh : TriMesh
        Triangular mesh.

    Returns
    -------
    TriGrad
        Gradient operator on the mesh.

    """
    return getgrad(mesh)


def getgrad(mesh):
    """Makes a gradient operator.
    """
    maker = GradMaker.from_metric(mesh_metric(mesh))
    return maker.get_grad()


def mesh_geom(mesh):
    """Returns the mesh geometric properties.

    Parameters
    ----------
    mesh : TriMesh
        Triangular mesh.

    Returns
    -------
    MeshGeom
        Object with the geometric properties of triangles.

    """
    maker = MeshGeomMaker.from_metric(mesh_metric(mesh))
    return maker.get_mesh_geom()


def mesh_metric(mesh):
    """Returns the mesh metric properties.

    Parameters
    ----------
    mesh : TriMesh
        Triangular mesh.

    Returns
    -------
    MeshMetric
        Object with the metric properties of triangles.

    """
    _ = MetricMaker.from_mesh(mesh)
    return _.get_metric()


class MeshMetric:
    """Mesh metric properties.
    """

    def __init__(self, mesh=None, data=None):

        self.mesh = mesh
        self.data = data
        self.meta = {}

        self.mask_voids()

    @classmethod
    def from_mesh(cls, mesh):
        return MetricMaker.from_mesh(mesh).get_metric()

    def mask_voids(self):

        voidsnums = self.mesh.getvoids()

        if voidsnums.size == 0:
            return

        self.meta['voids-trinums'] = voidsnums

        bcoeffs = self.data['bcoeffs']
        ccoeffs = self.data['ccoeffs']

        new_coeffs = {
            'bcoeffs': _zero_at_rows(bcoeffs, voidsnums),
            'ccoeffs': _zero_at_rows(ccoeffs, voidsnums)
        }

        self.data |= new_coeffs

    @property
    def bcoeffs(self):
        return self.data.get('bcoeffs')

    @property
    def ccoeffs(self):
        return self.data.get('ccoeffs')

    @property
    def jacobis(self):
        return self.data.get('jacobis')

    @property
    def areas1d(self):
        return np.copy(
            0.5 * self.jacobis[..., None], order='C'
        )

    @property
    def voids_trinums(self):
        return self.meta.get('voids-trinums')

    @property
    def hasvoids(self):
        return self.voids_trinums is not None

    def __getitem__(self, key):
        return self.data[key]


class MeshGeom:
    """Mesh geometric properties.
    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data

    @property
    def areas(self):
        return self.data['areas']

    @property
    def sides(self):
        return self.data['sides']

    @property
    def sides_filtered(self):
        return self.sides[self.sides > 0]

    def minside(self):
        return np.amin(
            self.sides_filtered
        )

    def maxside(self):
        return np.amax(
            self.sides_filtered
        )


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

        data = {
            **coeffs,
            **jacobs
        }

        return MeshMetric(
            self.mesh, data
        )

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


class MetricAgent:
    """Operator on a mesh metric.
    """

    def __init__(self, metric):
        self.metric = metric
        self.meta = self.fetch_meta()

    @classmethod
    def from_metric(cls, metric):
        return cls(metric)

    def fetch_meta(self):
        return {
            'voids-trinums': self.fetch_voids_trinums()
        }

    def fetch_voids_trinums(self):
        return self.metric.mesh.getvoids()

    @property
    def hasvoids(self):
        return self.voids_trinums.size != 0

    @property
    def voids_trinums(self):
        return self.meta['voids-trinums']

    def mask_trinums_not_voids(self):

        mask = np.full(
            self.metric.mesh.ntriangs, True
        )

        mask[self.voids_trinums] = False
        return mask

    @property
    def areas1d(self):
        return self.metric.areas1d

    @property
    def areas1d_inv(self):

        mask_not_voids = self.mask_trinums_not_voids()

        area_inv = np.zeros(
            self.metric.mesh.ntriangs
        )

        area_inv = np.reciprocal(
            self.areas1d.flat, where=mask_not_voids
        )

        return area_inv[..., None]


class OprsMaker(MetricAgent):
    """Maker of FEM operators.
    """

    def diff_1x(self):
        return self.diff_1d('bcoeffs')

    def diff_1y(self):
        return self.diff_1d('ccoeffs')

    def grad_1y(self):
        return self.grad_1d('ccoeffs')

    def diff_2x(self):
        return self.diff_2d('bcoeffs')

    def diff_2y(self):
        return self.diff_2d('ccoeffs')

    def diff_1d(self, coeffs_key):

        diff_1d = _mono_matrix(
            self.metric[coeffs_key]
        )

        return diff_1d / 6.

    def diff_2d(self, coeffs_key):

        diff_2d = _diad_matrix(
            self.metric[coeffs_key]
        )

        return 0.25 * (
            diff_2d * self.areas1d_inv
        )

    def grad_1d(self, coeffs_key):

        grad_1d = _solo_matrix(
            self.metric[coeffs_key]
        )

        return grad_1d / 6

    def massmat(self):

        areas = self.areas1d

        proxy = np.tile(
            self.massmat_proxy.flat, areas.shape
        )

        return proxy * (areas / 12.)

    def massdiag(self):

        areas = self.areas1d

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


class GradMaker(MetricAgent):
    """Makes a gradient operator.
    """

    def get_grad(self):

        data = self.make_data()
        grad = self.make_grad(data)

        return grad

    def make_grad(self, data):
        return TriGrad(
            self.metric.mesh, data
        )

    def make_data(self):
        return {
            'bcoeffs-scaled': self.make_bcoeffs_scaled(),
            'ccoeffs-scaled': self.make_ccoeffs_scaled()
        }

    def make_bcoeffs_scaled(self):
        return 0.5 * (
            self.metric.bcoeffs * self.areas1d_inv
        )

    def make_ccoeffs_scaled(self):
        return 0.5 * (
            self.metric.ccoeffs * self.areas1d_inv
        )


class MeshGeomMaker(MetricAgent):
    """Computes the geomeric properties of the triangles.
    """

    def get_mesh_geom(self):

        data = {
            'areas': self.get_areas(),
            'sides': self.get_sides()
        }

        return MeshGeom(
            self.metric.mesh, data
        )

    def get_areas(self):
        return self.metric.areas1d.flatten()

    def get_sides(self):

        sides_permuted = self.get_sides_permuted()

        sides_arranged = sides_permuted[:, [2, 0, 1]]
        return sides_arranged.copy('C')

    def get_sides_permuted(self):

        bcoefs = self.metric.bcoeffs
        ccoefs = self.metric.ccoeffs

        return np.sqrt(
            bcoefs * bcoefs + ccoefs * ccoefs
        )


class TriGrad:
    """Gradient operator.
    """

    def __init__(self, mesh, data):
        self.mesh = mesh
        self.data = data

    @property
    def bcoeffs(self):
        return self.data['bcoeffs-scaled']

    @property
    def ccoeffs(self):
        return self.data['ccoeffs-scaled']

    def __call__(self, data):

        diff_1x = self.diff_1x(data)
        diff_1y = self.diff_1y(data)

        return np.vstack(
            [diff_1x, diff_1y]
        )

    def diff_1x(self, data):
        """Computes the x-derivative across triangles.
        """
        return np.sum(
            self.bcoeffs * data[self.mesh.triangs], axis=1
        )

    def diff_1y(self, data):
        """Computes the y-derivative across triangles.
        """
        return np.sum(
            self.ccoeffs * data[self.mesh.triangs], axis=1
        )

    def atfunc(self, func):
        return self(
            func(*self.mesh.points2d)
        )


def _mono_matrix(data1d):
    return np.repeat(data1d, 3, axis=1)


def _solo_matrix(data1d):
    return np.tile(data1d, (1, 3))


def _diad_matrix(data1d):
    return np.repeat(data1d, 3, axis=1) * np.tile(data1d, (1, 3))


def _zero_at_rows(data2d, rowsnums):

    data2d[rowsnums, :] = 0.0

    return np.copy(
        data2d, order='C'
    )
