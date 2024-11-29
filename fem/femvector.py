# -*- coding: utf-8 -*-
"""FEM vector object.
"""
import numpy as np


def getvector(partt):
    """Creates a FEM vector.
    """
    return VectorFEM.from_partt(partt)


class VectorData:
    """Root of the FEM vector.
    """

    def __init__(self, partt=None, body=None):
        self.partt = partt
        self.body = body

    @classmethod
    def from_partt(cls, partt):
        return cls(
            partt, np.zeros(partt.unit.mesh_count)
        )

    def update_body(self, new_body):
        return self.__class__(self.partt, new_body)

    def get_section_indexer(self, key):
        return self.partt[key]

    @property
    def unit(self):
        return self.partt.unit

    @property
    def body_copy(self):
        return self.body.copy()


class VectorFEM(VectorData):
    """FEM vector.
    """

    def __getitem__(self, key):
        return self.getsection(key)

    def __setitem__(self, key, data):
        self.setsection(key, data)

    def __sub__(self, other):
        if self.partt is other.partt:
            return self.update_body(self.body - other.body)
        raise ValueError(
            "cannot subtract vectors with different partitions"
        )

    def __add__(self, other):
        if self.partt is other.partt:
            return self.update_body(self.body + other.body)
        raise ValueError(
            "cannot add vectors with different partitions"
        )

    def __mul__(self, value):
        return self.update_body(self.body * value)

    def with_body(self, value):
        """Defines the vector body.

        Parameters
        ----------
        value : scalar | flat-float-array
            Data that defines the vector body.

        Returns
        -------
        VectorFEM
            Copy of the vector with the body updated.

        """
        self.body[:] = value
        return self.update_body(self.body_copy)

    def from_func(self, func):
        """Defines the vector via a function on the mesh nodes.

        Parameters
        ----------
        func : Callable
            Function `(x, y)` that returns the vector body.

        Returns
        -------
        VectorFEM
            Copy of the vector with the body updated.

        """

        self.body[:] = func(
            *self.unit.mesh.points2d
        )

        return self.update_body(self.body_copy)

    def constrained(self):
        """Constrains the vector on the parent mesh.

        Returns
        -------
        VectorFEM
            Copy of the vector with the body constrained.

        """

        new_body = constr_data(
            self.partt.unit.mesh, self.body
        )

        return self.update_body(new_body)

    def getsection(self, sec_id):
        """Returns a copy of the vector section.

        Parameters
        ----------
        sec_id : int
            ID of the vector section.

        Returns
        -------
        flat-float-array
            Copy of the vector section.

        """

        indexer = self.get_section_indexer(sec_id)

        return np.copy(
            self.body[indexer], order='C'
        )

    def setsection(self, sec_id, data) -> None:
        """Defines the vector section.

        Parameters
        ----------
        sec_id : int
            ID of the vector section.
        data : scalar | flat-float-array
            Data that defines the vector section.

        """
        indexer = self.get_section_indexer(sec_id)
        self.body[indexer] = data


def constr_data(mesh, data):
    """Constrains data on a mesh.
    """

    voids_trinums = mesh.getvoids()

    if voids_trinums.size == 0:
        return data
    voids_triangs = mesh.triangs[voids_trinums, :]

    west, east, pivs = voids_triangs.T

    data[pivs] = 0.5 * (
        data[west] + data[east]
    )

    return data.copy('C')
