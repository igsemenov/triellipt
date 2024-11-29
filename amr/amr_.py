# -*- coding: utf-8 -*-
"""AMR public tools.
"""
import numpy as np
from triellipt.fem import femunit, trinterp
from triellipt.amr import (
    trirefine, tricoarsen, trinspect, trifronts
)


def getunit(mesh):
    """Creates an AMR unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.

    Returns
    -------
    AMRUnit
        AMR unit with the mesh.

    """
    return AMRUnit.from_mesh(mesh)


class AMRData:
    """Root of the AMR unit.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = {}
        self.data = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def mesh_twin(self):
        return self.mesh.twin()

    @property
    def data_items(self):
        return self.data.items()

    @property
    def is_refined(self):
        return 'data-refiner' in self.mesh.meta

    @property
    def is_coarsened(self):
        return 'data-collect' in self.mesh.meta

    @property
    def refiner(self):
        if self.is_refined:
            return self.mesh.meta['data-refiner']
        return None

    @property
    def collector(self):
        if self.is_coarsened:
            return self.mesh.meta['data-collect']
        return None

    @property
    def transmat(self):
        return self.meta['data-transmit']


class AMRUnit(AMRData):
    """Mesh refinement unit.

    Attributes
    ----------
    mesh : TriMesh
        Current triangle mesh.
    data : dict
        Data defined on the unit.

    Properties
    ----------

    Name        | Description
    ------------|----------------------------------
    `refiner`   | Data-refiner after refinement.
    `collector` | Data-collector after coarsening.

    """

    def refine(self, trinums=None):
        """Performs a static mesh refinement.

        Parameters
        ----------
        trinums : Iterable = None
            Numbers of triangles to refine, if None takes all triangles.

        Returns
        -------
        AMRUnit
            Unit with the refined mesh.

        Notes
        -----

        - The `data-refiner` is included in the mesh metadata.
        - The void ears are not refined to keep the mesh 1-irregular.

        """

        if trinums is None:
            trinums = tuple(
                range(self.mesh.size)
            )

        if len(trinums) == 0:
            return self

        new_unit = self.from_mesh(
            trirefine.refine_mesh(self.mesh, trinums)
        )

        if not self.data:
            return new_unit

        refiner = new_unit.refiner

        new_unit.data = {
            k: refiner(v) for k, v in self.data_items
        }

        return new_unit

    def coarsen(self, trinums_cores):
        """Performs a static mesh coarsening.

        Parameters
        ----------
        trinums_cores : Iterable
            Numbers of the super-triangle-cores to coarsen.

        Returns
        -------
        AMRUnit
            Unit with the coarsened mesh.

        Notes
        -----

        - The `data-collector` is included in the mesh metadata.

        """

        new_unit = self.from_mesh(
            tricoarsen.coarsen_mesh(self.mesh, trinums_cores)
        )

        if not self.data:
            return new_unit

        collector = new_unit.collector

        new_unit.data = {
            k: collector(v) for k, v in self.data_items
        }

        return new_unit

    def find_node(self, anchor):
        """Finds the neighborhood of an anchor point.

        Parameters
        ----------
        anchor : (float, float)
            Coordinates of an anchor point.

        Returns
        -------
        flat-int-array
            Numbers of triangles near the anchor point.

        """
        return trinspect.find_node(self.mesh, anchor).trinums

    def find_subset(self, count, anchor, remove_heads=False):
        """Finds a convex subset of a mesh.

        Parameters
        ----------
        count : int
            Seed number of triangles in a subset.
        anchor : (float, float)
            Anchor point for a starting triangle.
        remove_heads : bool = False
            Removes single-paired triangles, if True.

        Returns
        -------
        flat-int-array
            Numbers of triangles in a subset.

        """
        return trinspect.find_subset(
            self.mesh, count, anchor, remove_heads
        )

    def find_masked(self, mask):
        """Finds triangles by a mask-function.

        Parameters
        ----------
        mask : function
            Boolean mask `(x, y) ` on the triangles centroids.

        Returns
        -------
        flat-int-array
            Numbers of the found triangles.

        """

        trinums, = np.where(
            mask(*self.mesh.centrs2d)
        )

        return trinums

    def front_coarse(self):
        """Finds a front of coarse triangles.
        """
        return trifronts.TriFrontCoarse.from_unit(self)

    def front_fine(self):
        """Finds a front of fine triangles.
        """
        return trifronts.TriFrontFine.from_unit(self)

    def makedata(self, key, func, constrained=True):
        """Generates a unit data item from a given function.

        Parameters
        ----------
        key : str-or-int
            Key of the unit data item.
        func : Callable
            The item-source function `(x, y)` on the mesh nodes.
        constrained : bool = True
            Constrains the new data item on the mesh, if True.

        Returns
        -------
        self
            The unit with the new data item.

        """

        data = func(
            *self.mesh.points2d
        )

        if constrained is False:
            self.data[key] = data
            return self

        self.data[key] = trirefine.constr_data(self.mesh, data)
        return self

    def getinterp(self, xnodes, ynodes):
        """Creates an interpolator on a mesh.

        Parameters
        ----------
        xnodes : flat-float-array
            x-coordinates of the interpolation nodes.
        ynodes : flat-float-array
            y-coordinates of the interpolation nodes.

        Returns
        -------
        TriInterp
            Callable interpolator.

        """
        return trinterp.getinterp(
            self.mesh, xnodes, ynodes
        )

    @property
    def data_items(self):
        return self.data.items()

    @property
    def masser(self):
        return self.get_masser(self.mesh)

    @classmethod
    def get_masser(cls, mesh):
        return Masser.from_mesh(mesh)


class Masser:
    """Total mass evaluator.
    """

    def __init__(self, unit, meta):
        self.unit = unit
        self.meta = meta

    @classmethod
    def from_mesh(cls, mesh):

        unit = femunit.getunit(mesh)

        mass_mat = unit.massopr(
            is_lumped=False, add_constr=False
        )

        massdiag = unit.massopr(
            is_lumped=True, add_constr=False
        )

        meta = {
            'mass-mat': mass_mat.body,
            'massdiag': massdiag.body
        }

        return cls(unit, meta)

    def __call__(self, data):
        return self.mass_diag(data)

    def mass_full(self, data):
        return np.sum(
            self.mass_mat @ self.unit.perm(data)
        )

    def mass_diag(self, data):
        return np.sum(
            self.massdiag @ self.unit.perm(data)
        )

    @property
    def mass_mat(self):
        return self.meta['mass-mat']

    @property
    def massdiag(self):
        return self.meta['massdiag']
