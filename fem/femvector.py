# -*- coding: utf-8 -*-
"""FEM vector object.
"""
import numpy as np


class VectorData:
    """Root of the FEM vector.
    """

    def __init__(self, unit=None, body=None, meta=None):
        self.unit = unit
        self.body = body
        self.meta = meta or {}

    @classmethod
    def from_unit(cls, unit):
        return cls(
            unit, np.zeros(unit.mesh_count)
        )

    @property
    def is_named(self):
        return 'name' in self.meta

    @property
    def haspartition(self):
        return 'partition' in self.meta

    @property
    def partition_name(self):
        if self.haspartition:
            return self.meta['partition']['name']
        return None

    @property
    def name(self):
        if self.is_named:
            return self.meta['name']
        return str(id(self))

    def update_meta(self, new_meta):
        return self.__class__(
            self.unit, self.body_copy, new_meta
        )

    def update_body(self, new_body):
        return self.__class__(
            self.unit, new_body, self.meta_copy
        )

    def get_section_indexer(self, name):

        if not self.haspartition:
            raise VectorFEMError(
                f"vector '{self.name}' has no partition"
            )

        indexer = self.meta['partition']['sections'].get(name)

        if indexer is None:
            raise VectorFEMError(
                f"no section '{name}' in vector '{self.name}'"
            )

        return indexer

    @property
    def body_copy(self):
        return self.body.copy()

    @property
    def meta_copy(self):
        return self.meta.copy()


class VectorFEM(VectorData):
    """FEM vector.
    """

    def with_name(self, name):
        """Prescribes the vector name.

        Parameters
        ----------
        name : str
            Name of the vector.

        Returns
        -------
        VectorFEM
            Copy of the vector with the new name.

        """
        return self.add_meta({'name': name})

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

    def add_meta(self, meta):
        return self.update_meta(self.meta | meta)

    def from_func(self, func):
        """Defines the vector body from a function.

        Parameters
        ----------
        func : Callable
            Function with `(x, y)` call that returns the vector body.

        Returns
        -------
        VectorFEM
            Copy of the vector with the body updated.

        """
        self.body[:] = func(*self.unit.mesh.points2d)
        return self.update_body(self.body_copy)

    def dirichsplit(self):
        """Creates a Dirichlet (edge-core) partition.

        Returns
        -------
        VectorFEM
            Copy of the vector with the partition defined.

        """
        return self.partitioned(
            meta_edge_core(self.unit)
        )

    def partitioned(self, meta):
        """Creates a partitioned vector.

        Parameters
        ----------
        meta : dict
            Partition meta data (name and sections) (a).

        Returns
        -------
        VectorFEM
            Copy of the vector with the partition defined.

        Notes
        -----

        (a) Same as for the FEM matrix.

        """

        new_meta = {
            'partition': meta
        }

        return self.add_meta(new_meta)

    def getsect(self, name):
        """Returns a copy of the vector section.

        Parameters
        ----------
        name : str
            Name of the vector section.

        Returns
        -------
        flat-float-array
            Copy of the vector section.

        """
        indexer = self.get_section_indexer(name)
        return self.body[indexer].copy('C')

    def setsect(self, name, data) -> None:
        """Defines the vector section.

        Parameters
        ----------
        name : str
            Name of the vector section.
        data : scalar | flat-float-array
            Data that defines the vector section.

        """
        indexer = self.get_section_indexer(name)
        self.body[indexer] = data

    def sect_xy(self, name):
        """Returns xy-coordinates of the vector section nodes.

        Parameters
        ----------
        name : str
            Name of the vector section.

        Returns
        -------
        two-row-float-array
            xy-coordinates of the vector section.

        """

        indexer = self.get_section_indexer(name)

        return _unpack_complex(
            self.unit.mesh.points[indexer]
        )


VectorFEMError = type(
    'VectorFEMError', (Exception,), {}
)


def meta_edge_core(unit):
    """Dirichlet partition meta data for a FEM unit.
    """

    mesh_count = unit.mesh_count
    edge_count = unit.edge_count

    mesh_range = np.arange(mesh_count)

    edge_range = mesh_range[:edge_count]
    core_range = mesh_range[edge_count:]

    return {
        'name': 'edge-core',
        'sections': {
            'edge': edge_range,
            'core': core_range
        }
    }


def _unpack_complex(data):
    return np.vstack(
        [data.real, data.imag]
    )
