# -*- coding: utf-8 -*-
"""AMR public tools.
"""
import numpy as np
from triellipt.fem import femunit
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

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def mesh_twin(self):
        return self.mesh.twin()

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
    def source_mesh(self):
        if self.is_refined:
            return self.refiner.mesh
        if self.is_coarsened:
            return self.collector.mesh
        return None


class AMRUnit(AMRData):
    """Mesh refinement unit.

    Attributes
    ----------
    mesh : TriMesh
        Current triangle mesh.

    Properties
    ----------

    Metadata on refinement:

    Name           | Description
    ---------------|-------------------------------
    `source_mesh`  | Previous-level mesh.
    `data_refiner` | Callable refiner object.

    """

    def tofem(self, anchors=None):
        """Converts the AMR unit to a FEM one.

        Parameters
        ----------
        anchors : tuple = None
            Nodes numbers to synchronize the mesh boundary.

        Returns
        -------
        FEMUnit
            The resulting FEM unit.

        """
        return femunit.getunit(self.mesh, anchors)

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
        - The voids ears are not refined to keep the mesh 1-irregular.

        """

        if trinums is None:
            trinums = tuple(
                range(self.mesh.size)
            )

        if len(trinums) == 0:
            return self

        return self.from_mesh(
            trirefine.refine_mesh(self.mesh, trinums)
        )

    def coarsen(self, trinums_cores):
        """Performs a static mesh coarsening.

        Parameters
        ----------
        trinums_cores : Iterable
            Numbers of the super-triangles-cores to coarsen.

        Returns
        -------
        AMRUnit
            Unit with the coarsened mesh.

        """
        return self.from_mesh(
            tricoarsen.coarsen_mesh(self.mesh, trinums_cores)
        )

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

    def constrain(self, data):
        """Constrains data on a mesh.

        Parameters
        ----------
        data : flat-float-array
            Nodal data to be constrained.

        Returns
        -------
        flat-float-array
            Constained nodal data.

        """
        return trirefine.constr_data(self.mesh, data)

    @classmethod
    def constr_data(cls, mesh, data):
        return trirefine.constr_data(mesh, data)
