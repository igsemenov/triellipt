# -*- coding: utf-8 -*-
"""Explicit transport unit.
"""
import numpy as np
from triellipt.fem import femoprs


def gettransp(mesh, geom=None):
    """Creates an explicit transport unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    geom : str
        Geometry mode — "ax" for cylindrical, "2d" for planar (default).

    Returns
    -------
    TranspUnit
        Explicit transport unit.

    """
    return TranspUnit.from_mesh(mesh, geom or '2d')


class TranspData:
    """Base transport unit.
    """

    def __init__(self, mesh=None, meta=None, geom=None):
        self.mesh = mesh
        self.meta = meta
        self.geom = geom

    @property
    def size(self):
        return self.mesh.size

    @classmethod
    def from_mesh(cls, mesh, geomode):

        mesh, meta = cls.fetch_meta(mesh)

        unit = cls(
            mesh, meta, geomode
        )

        return unit.with_mass()

    @classmethod
    def fetch_meta(cls, mesh):

        mesh, vods = _fetch_voids(mesh)

        nods = mesh.nodesmap()
        dist = mesh.centrs_complex.imag

        geom = femoprs.mesh_geom(mesh)
        metr = femoprs.mesh_metric(mesh)

        meta = {
            'nods': nods,
            'vods': vods,
            'metr': metr,
            'dist': dist,
            'geom': geom,
            'mass': None
        }

        return mesh, meta

    def with_mass(self):
        self.meta['mass'] = self.comp_mass()
        return self

    def comp_mass(self):
        return self.reduce(self.masser_long)

    @property
    def mass(self):
        return self.meta['mass']

    @property
    def metr(self):
        return self.meta['metr']

    @property
    def nods(self):
        return self.meta['nods']

    @property
    def vods(self):
        return self.meta['vods']

    @property
    def dist(self):
        return self.meta['dist']

    @property
    def masser_long(self):
        return (1/6) * (
            self.jacobis[self.trinums]
        )

    @property
    def jacobis(self):
        return self.metr.jacobis

    @property
    def trinums(self):
        return self.meta['nods'].trinums

    @property
    def locnums(self):
        return self.meta['nods'].locnums

    @property
    def bins_reduce(self):
        return self.nods.data.bins_reduce

    @property
    def hasvoids(self):
        return self.vods is not None

    def reduce(self, data_long):
        return np.add.reduceat(
            data_long, self.bins_reduce
        )


class TranspRoot(TranspData):
    """Root of the public transport unit.
    """

    def averval(self, data):
        """Triangle average.
        """
        return (1/3) * np.sum(
            data[self.mesh.triangs], axis=1
        )

    def diff_1x(self, data):
        """x-gradient across triangles.
        """

        grad = np.sum(
            data[self.mesh.triangs] * self.bcoeffs, axis=1
        )

        return grad / self.jacobis

    def diff_1y(self, data):
        """y-gradient across triangles.
        """

        grad = np.sum(
            data[self.mesh.triangs] * self.ccoeffs, axis=1
        )

        return grad / self.jacobis

    def div_x(self, field):
        """x-divergence
        """
        return 0.5 * self.reduce(
            field[self.trinums] * self.bcoeffs_long
        )

    def div_y(self, field):
        """y-divergence
        """
        return 0.5 * self.reduce(
            field[self.trinums] * self.ccoeffs_long
        )

    @property
    def bcoeffs(self):
        return self.metr.bcoeffs

    @property
    def ccoeffs(self):
        return self.metr.ccoeffs

    @property
    def bcoeffs_long(self):
        return self.bcoeffs.flat[self.longind]

    @property
    def ccoeffs_long(self):
        return self.ccoeffs.flat[self.longind]

    @property
    def longind(self):
        return 3 * self.trinums + self.locnums

    def push_data(self, data):

        if not self.hasvoids:
            return data

        west = self.vods[:, 0]
        east = self.vods[:, 1]
        pivs = self.vods[:, 2]

        data[west] = data[west] + 0.5 * data[pivs]
        data[east] = data[east] + 0.5 * data[pivs]

        data[pivs] = 0

        return np.copy(
            data, order='C'
        )


class TranspUnit(TranspRoot):
    """Explicit transport unit.

    Properties
    ----------

    Name      | Description
    ----------|--------------------------
    `mass`    | Node-based mass vector

    """

    def constr(self, data):
        """Constrains node-based data on hanging nodes.

        Parameters
        ----------
        data : flat-float-array
            Node-based data to constrain.

        Returns
        -------
        flat-float-array
            Constrained data.

        """

        if not self.hasvoids:
            return data

        west = self.vods[:, 0]
        east = self.vods[:, 1]
        pivs = self.vods[:, 2]

        data[pivs] = 0.5 * (
            data[west] + data[east]
        )

        return data.copy('C')

    def transp(self, data, v_x, v_y, d_x, d_y, stab):
        """Computes the transport operator.

        Parameters
        ----------
        data : flat-float-array
            Node-based solution field.
        v_x : flat-float-array
            Triangle-based x-velocity.
        v_y : flat-float-array
            Triangle-based y-velocity.
        d_x : flat-float-array
            Triangle-based x-diffusion coefficient.
        d_y : flat-float-array
            Triangle-based y-diffusion coefficient.
        stab : Callable
            Stream upwind stabilizator called on velocity arrays.

        Returns
        -------
        flat-float-array
            Node-based transport operator.

        """

        div = self.getdiverg(
            data, v_x, v_y, d_x, d_y, stab
        )

        return self.push_data(div)

    def getdiverg(self, data, v_x, v_y, d_x, d_y, stab):

        flux_x, flux_y = self.getfluxes(
            data, v_x, v_y, d_x, d_y, stab
        )

        div = self.div_x(flux_x) + self.div_y(flux_y)

        if self.geom == "2d":
            return div

        src = self.source(
            flux_y / self.dist
        )

        return div - src

    def getfluxes(self, data, v_x, v_y, d_x, d_y, stab):

        v_x = v_x[:self.size]
        v_y = v_y[:self.size]
        d_x = d_x[:self.size]
        d_y = d_y[:self.size]

        centrv = self.averval(data)

        grad_x = self.diff_1x(data)
        grad_y = self.diff_1y(data)

        nu_x = d_x + stab(v_x)
        nu_y = d_y + stab(v_y)

        flux_x = v_x * centrv - (nu_x * grad_x)
        flux_y = v_y * centrv - (nu_y * grad_y)

        return flux_x, flux_y

    def source(self, field):
        """Computes the source term.

        Parameters
        ----------
        field : float-flat-array
            Triangle-based source field.

        Returns
        -------
        float-flat-array
            Node-based source term.

        """

        data = self.reduce(
            self.masser_long * field[self.trinums]
        )

        return self.push_data(data)

    def newdata(self, value_or_func):
        """Creates a new node-based field.

        Parameters
        ----------
        value_or_func : float-or-callable
            Coefficient value or function `(x, y)` on mesh nodes.

        Returns
        -------
        float-flat-array
            Node-based coefficient.

        """
        if callable(value_or_func):
            return value_or_func(*self.mesh.points2d)
        return np.full(
            self.mesh.npoints, value_or_func
        )

    def newcoeff(self, value_or_func):
        """Creates a new triangle-based coefficient.

        Parameters
        ----------
        value_or_func : float-or-callable
            Coefficient value or function `(x, y)` on triangle centroids.

        Returns
        -------
        float-flat-array
            Triangle-based coefficient.

        """
        if callable(value_or_func):
            return value_or_func(*self.mesh.centrs2d)
        return np.full(
            self.mesh.ntriangs, value_or_func
        )


def _fetch_voids(mesh):

    if not mesh.hasvoids():
        return mesh, None

    nums = mesh.getvoids()

    vods = np.copy(
        mesh.triangs[nums, :], order='C'
    )

    mesh = mesh.deltriangs(*nums).twin()
    return mesh, vods
