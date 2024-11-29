# -*- coding: utf-8 -*-
"""Mapping mesh nodes.
"""
import numpy as np
from triellipt.utils import tables


def _maptable(table):
    return tables.TableMap.from_table(table)


class NodesMap:
    """Nodes-to-triangles map.

    Properties
    ----------

    Nodes data:

    Name      | Description
    ----------|------------------------
    `nodnums` | Global numbers.
    `trinums` | Host triangles (HTs).
    `locnums` | Local numbers in HTs.

    Neighbors local numbers:

    Name       | Description
    -----------|----------------------
    `locnums1` | Next CCW node.
    `locnums2` | Next-next CCW node.

    Neighbors global numbers:

    Name       | Description
    -----------|----------------------
    `nodnums1` | Next CCW node.
    `nodnums2` | Next-next CCW node.

    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data
        self.meta = {}

    @classmethod
    def from_mesh(cls, mesh):
        data = _maptable(mesh.triangs)
        return cls(mesh, data)

    @property
    def size(self):
        return self.data.data.shape[1]

    @property
    def nodnums(self):
        return self.data.vals

    @property
    def trinums(self):
        return self.data.rows

    @property
    def locnums(self):
        return self.data.cols

    @property
    def locnums1(self):
        return (self.locnums + 1) % 3

    @property
    def locnums2(self):
        return (self.locnums + 2) % 3

    @property
    def nodnums1(self):
        return self.triangs[
            3 * self.trinums + self.locnums1
        ]

    @property
    def nodnums2(self):
        return self.triangs[
            3 * self.trinums + self.locnums2
        ]

    @property
    def triangs(self):
        return self.mesh.triangs.flat

    @property
    def nodes_range(self):
        return self.data.vals_unique

    @property
    def nodes_ranks(self):
        return self.data.vals_ranks

    def atrank(self, rank):
        """Extracts data for nodes with the specified rank.

        Parameters
        ----------
        rank : int
            Rank of the nodes to extract.

        Returns
        -------
        tuple
            Triplet `(nodnums, trinums, locnums)` for the extracted nodes.

        """
        return tuple(
            self.data.atrank(rank)
        )

    def atnode(self, nodenum):
        """Extracts the map of a single node.

        Parameters
        ----------
        nodenum : int
            Number of the target node.

        Returns
        -------
        NodesMap
            The resulting single-node map.

        """

        if nodenum not in self.nodes_range:
            raise ValueError(
                f"node {nodenum} is not in the nodes map"
            )

        one_node_data = self.get_one_node_data(nodenum)
        return self.update_table_data(one_node_data)

    def get_one_node_data(self, nodenum):

        index = np.argmax(
            self.nodes_range == nodenum
        )

        return self.data.data_split()[index]

    def update_table_data(self, new_table_data):

        new_table_map = type(self.data)(new_table_data)

        return type(self)(
            self.mesh, new_table_map
        )

    def add_meta(self, new_meta):
        self.meta.update(new_meta)
        return self
