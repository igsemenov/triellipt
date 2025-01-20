# -*- coding: utf-8 -*-
"""Public FEM unit.
"""
import numpy as np
from triellipt.fem import femdata
from triellipt.fem import femfactory
from triellipt.fem import femvector
from triellipt.fem import trinterp


def getunit(mesh, anchors=None):
    """Creates a FEM computing unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : tuple = None
        Nodes numbers to synchronize the mesh boundary. 

    Returns
    -------
    FEMUnit
        FEM computing unit.

    """
    return FEMUnit.from_mesh(mesh, anchors)


class FEMData:
    """Base FEM data.
    """

    def __init__(self, mesh=None, meta=None):

        if mesh is None:
            return
        if meta is None:
            return

        self.mesh = mesh

        self.masks = meta['masks']
        self.loops = meta['loops']
        self.femat = meta['femat']
        self.perms = meta['perms']

    @classmethod
    def from_mesh(cls, mesh, anchors=None):

        data = femdata.get_fem_data(mesh, anchors)

        return cls(
            data['mesh'], data['meta']
        )

    @property
    def hasjoints(self):
        return bool(
            self.masks['joints']
        )

    @property
    def frame_size(self):
        """Number of frame triangles.
        """
        if self.hasjoints:
            return self.masks['joints']['voids'][0]
        return self.masks['body'][1]

    @property
    def voids_size(self):
        """Number of voids triangles.
        """
        return self.mesh.size - self.frame_size

    @property
    def mesh_count(self):
        """Number of mesh nodes.
        """
        return self.mesh.npoints

    @property
    def voids_count(self):
        """Number of voids nodes.
        """
        return self.voids_size

    @property
    def frame_count(self):
        """Number of frame nodes.
        """
        return self.mesh_count - self.voids_count

    @property
    def edge_count(self):
        """Number of edge nodes.
        """

        sizes = [
            loop.size for loop in self.loops
        ]

        return sum(sizes)

    @property
    def core_count(self):
        """Number of core nodes.
        """
        return self.mesh_count - self.edge_count

    @property
    def ij_stream(self):
        return self.femat['ij-data']

    def get_joints_item(self, key):

        if not self.hasjoints:
            return None

        return self.mesh.submesh(
            *range(*self.masks['joints'][key])
        )

    def new_vector_data(self):
        return np.zeros(self.mesh_count)


class FEMRoot(FEMData):
    """Root of a FEM unit.
    """

    @property
    def mesh_body(self):
        return self.mesh.submesh(
            *range(*self.masks['body'])
        )

    @property
    def mesh_wests(self):
        return self.get_joints_item('wests')

    @property
    def mesh_easts(self):
        return self.get_joints_item('easts')

    @property
    def mesh_cores(self):
        return self.get_joints_item('cores')

    @property
    def mesh_voids(self):
        return self.get_joints_item('voids')

    @property
    def massmat(self):
        return self.femat['femoprs']['massmat']

    @property
    def massdiag(self):
        return self.femat['femoprs']['massdiag']

    @property
    def diff_1x(self):
        return self.femat['femoprs']['diff_1x']

    @property
    def diff_1y(self):
        return self.femat['femoprs']['diff_1y']

    @property
    def diff_2x(self):
        return self.femat['femoprs']['diff_2x']

    @property
    def diff_2y(self):
        return self.femat['femoprs']['diff_2y']

    @property
    def laplace(self):
        return self.diff_2x + self.diff_2y


class FEMUnit(FEMRoot):
    """FEM computing unit.

    Properties
    ----------

    FEM local operators:

    Name        | Description
    ------------|---------------------
    `massmat`   | Mass-matrix
    `massdiag`  | Mass-matrix lumped
    `laplace`   | Laplace operator
    `diff_1y`   | 1st-y derivative
    `diff_1x`   | 1st-x derivative
    `diff_2y`   | 2nd-y derivative
    `diff_2x`   | 2nd-x derivative

    FEM global matrices:

    Name           | Description
    -------------- |---------------------
    `massmat_fem`  | Mass-matrix
    `massdiag_fem` | Mass-matrix lumped
    `laplace_fem`  | Laplace operator

    """

    def fem_factory(self, add_constraints=True):
        """Creates a factory of FEM matrices.

        Parameters
        ----------
        add_constraints : bool = True
            If True, forces constraints to be included in the matrix.

        Returns
        -------
        FEMFactory
            Resulting factory of FEM matrices.

        """
        return femfactory.FEMFactory.from_unit(self, add_constraints)

    def new_vector(self):
        """Returns a new FEM vector.
        """
        return femvector.VectorFEM.from_unit(self)

    def makecoeff(self, mesh_data):
        """Generates a coefficient for local FEM operators.

        Parameters
        ----------
        mesh_data : flat-float-array
            Coefficient defined over the mesh triangles.

        Returns
        -------
        flat-float-array
            Coefficient matching the sizes of local FEM operators.

        """
        return mesh_data[
            self.ij_stream.trinums
        ]

    def getinterp(self, xnodes, ynodes):
        """Creates an interpolator on a mesh.

        Parameters
        ----------
        xnodes : flat-float-array
            X-coordinates of the interpolation nodes.
        ynodes : flat-float-array
            Y-coordinates of the interpolation nodes.

        Returns
        -------
        TriInterp
            Interpolator object.

        """
        return trinterp.getinterp(self, xnodes, ynodes)

    @property
    def massmat_fem(self):
        return self.fem_factory(0).feed_data(self.massmat)

    @property
    def massdiag_fem(self):
        return self.fem_factory(0).feed_data(self.massdiag).with_no_zeros()

    @property
    def laplace_fem(self):
        return self.fem_factory(1).feed_data(self.laplace)
