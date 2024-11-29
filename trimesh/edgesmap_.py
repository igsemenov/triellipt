# -*- coding: utf-8 -*-
"""Mapping inner edges.
"""
import numpy as np
from triellipt.utils import tables


class EdgesData:
    """Data on mesh edges.
    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_data(cls, mesh, data):
        return cls(mesh, data)

    @classmethod
    def from_mesh(cls, mesh):
        return EdgesMapper(mesh).get_edges_map()

    @property
    def size(self):
        return self.data.shape[1]

    @property
    def nedges(self):
        return self.data.shape[1]

    @property
    def trinums1(self):
        return self.data[0, :]

    @property
    def trinums2(self):
        return self.data[1, :]

    @property
    def locnums1(self):
        return self.data[2, :]

    @property
    def locnums2(self):
        return self.data[3, :]

    def update_data(self, new_data):
        return self.from_data(
            self.mesh, new_data.copy('C')
        )

    def permutedata(self, permuter):
        return self.update_data(
            self.data[:, permuter]
        )

    def sort_per_trinum(self):
        return self.permutedata(
            np.argsort(self.trinums1)
        )


class EdgesMap(EdgesData):
    """Map of inner mesh edges.

    Properties
    ----------

    Name       | Description
    -----------|---------------------------
    `trinums1` | Host triangle one (T1).
    `trinums2` | Host triangle two (T2).
    `locnums1` | Local edge number in T1.
    `locnums2` | Local edge number in T2.

    """

    def getspec(self):
        """Classifies triangles based on edges pairing.

        Returns
        -------
        dict
            Triangles numbers for each category (i).

        Notes
        -----

        (i) Triangles categories:

        - "heads" — 1 pair
        - "links" — 2 pairs
        - "cores" — 3 pairs
        - "spots" — 0 pairs

        """
        return EdgesSpec(self).getspec()


class MeshAgent:
    """Operator on a mesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def edges_paired(self):
        return self.mesh.edges_paired()


class EdgesMapper(MeshAgent):
    """Mapper of inner mesh edges.
    """

    def get_edges_map(self):

        paredges = self.make_edges_paired()
        edgesmap = self.from_edges_paired(paredges)

        return edgesmap

    def from_edges_paired(self, paredges):

        _, rows, cols = paredges.atrank(2)

        trinums = np.reshape(rows, (rows.size // 2, 2)).T
        locnums = np.reshape(cols, (cols.size // 2, 2)).T

        return EdgesMap.from_data(
            mesh=self.mesh, data=np.vstack([trinums, locnums])
        )

    def make_edges_paired(self):
        return tables.TableMap.from_table(self.edges_paired)


class EdgesSpec:
    """Edges classifier.
    """

    NAMES_AT_RANKS = {
        1: 'heads', 2: 'links', 3: 'cores'
    }

    def __init__(self, edgesmap):
        self.edges = edgesmap

    def getspec(self):

        spots = self.getunpaired()
        paired = self.learnpaired()

        return {
            'spots': spots, **paired
        }

    def getunpaired(self):

        paired = np.unique(self.trinums2d)

        return _fastnpdel(
            self.trinums_range, paired
        )

    def learnpaired(self):
        return dict(
            self.genpaired()
        )

    def genpaired(self):

        secmap = self.secondary_edgesmap()

        for rank in self.ranks_in_spec:

            name = self.NAMES_AT_RANKS[rank]

            vals, _, _ = secmap.atrank(rank)
            nums = np.unique(vals).astype(int)

            yield name, nums

    def secondary_edgesmap(self):
        return tables.TableMap.from_table(self.trinums2d)

    @property
    def trinums2d(self):
        return _packcols(
            self.edges.trinums1, self.edges.trinums2
        )

    @property
    def ranks_in_spec(self):
        return [1, 2, 3]

    @property
    def trinums_range(self):
        return np.arange(
            self.edges.mesh.ntriangs
        )


def _packcols(*cols):
    return np.vstack(cols).T.copy('C')


def _fastnpdel(array, indices):

    mask = np.full(array.size, True)
    mask[indices] = False

    return np.compress(mask, array)
