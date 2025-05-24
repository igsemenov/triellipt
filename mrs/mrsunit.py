# -*- coding: utf-8 -*-
"""Public MRS unit.
"""
import numpy as np
from triellipt.amr import trinspect
from triellipt.mrs import triskel
from triellipt.mrs import trielm
from triellipt.mrs import mrsmesh


def getunit(mesh):
    """Creates an MRS unit.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.

    Returns
    -------
    MRSUnit
        MRS unit with the mesh.

    """
    return MRSUnit.from_mesh(mesh)


class MRSBase:
    """Parent MRS unit.
    """

    def __init__(self):
        self.mesh = None
        self.elms = None

    @classmethod
    def from_mesh(cls, mesh):
        return MRSUnitMaker.from_mesh(mesh).get_mrs_unit().activated()

    def activated(self):

        self.elms = [
            elm.refine() for elm in self.elms
        ]

        return self

    @property
    def ranks(self):

        data = [
            elm.rank for elm in self.elms
        ]

        return np.array(data)

    @property
    def size(self):
        return len(self.elms)

    @property
    def counts(self):
        return [
            elm.mesh.npoints for elm in self.elms
        ]

    @property
    def is_empty(self):
        return not self.elms

    @property
    def data_keys(self):
        if self.is_empty:
            return []
        return self.elms[0].data.keys()


class MRSUnit(MRSBase):
    """MRS unit.
    """

    def refine(self, elmnums=None):
        """Refines the specified elements.

        Parameters
        ----------
        elmnums : flat-int-array
            Numbers of the elements to refine.

        Returns
        -------
        self
            Unit with the refined elements.

        """
        for ind in _filtnums(elmnums, self.size):
            self.elms[ind] = self.elms[ind].refine()
        return self

    def coarsen(self, elmnums):
        """Coarsens the specified elements.

        Parameters
        ----------
        elmnums : flat-int-array
            Numbers of the elements to coarsen.

        Returns
        -------
        self
            Unit with the refined elements.

        """
        for ind in _filtnums(elmnums, self.size):
            self.elms[ind] = self.elms[ind].coarsen()
        return self

    def refine_all(self):
        """Refines all elements in the unit.
        """
        return self.refine(
            np.arange(self.size)
        )

    def coarsen_all(self):
        """Coarsens all elements in the unit.
        """
        return self.coarsen(
            np.arange(self.size)
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
            Numbers of elements near the anchor point.

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

        return trinums.astype(int)

    def makedata(self, key, func):
        """Generates a unit data item from a given function.

        Parameters
        ----------
        key : str-or-int
            Key of the unit data item.
        func : Callable
            The item-source function `(x, y)` on the mesh nodes.

        Returns
        -------
        self
            The unit with the new data item.

        """

        self.elms = [
            e.makedata(key, func) for e in self.elms
        ]

        return self

    def getskel(self):
        return triskel.getskel(self)

    def getpremesh(self):
        return triskel.getskel(self).get_premesh()

    def getmrsmesh(self):
        return mrsmesh.getmrsmesh(self)


def _filtnums(nums, size):

    nums = np.unique(nums)

    nums = nums[
        (nums >= 0) & (nums < size)
    ]

    return nums.astype(int)


class MRSUnitMaker:
    """Maker of the MRS unit.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def get_mrs_unit(self):

        elms = self.make_elms()
        unit = self.push_unit(elms)

        return unit

    def push_unit(self, elms):

        unit = MRSUnit()

        unit.mesh = self.mesh.twin()
        unit.elms = elms.copy()

        return unit

    def make_elms(self):

        data = self.make_data_from_mesh()
        elms = self.make_elms_from_data(data)

        return elms

    def init_elms(self, elms):
        return [
            elm.refine() for elm in elms
        ]

    def make_data_from_mesh(self):

        data = {}

        data['triangs'] = self.mesh.triangs
        data['trinods'] = self.mesh.points[self.mesh.triangs]

        return data

    def make_elms_from_data(self, data):

        items = zip(
            data['trinods'], data['triangs']
        )

        return [
            trielm.TriElm.from_data(p, t) for p, t in items
        ]
