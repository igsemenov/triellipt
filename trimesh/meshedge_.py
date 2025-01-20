# -*- coding: utf-8 -*-
"""Finds the mesh edge.
"""
import numpy as np
from triellipt.utils import tables
from triellipt.utils import loops


class EdgeData:
    """Mesh edge data.
    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_data(cls, mesh, data):
        return cls(mesh, data)

    @property
    def size(self):
        return self.nedges

    @property
    def nedges(self):
        return self.data.shape[1]

    @property
    def trinums(self):
        return self.data[0, :]

    @property
    def locnums(self):
        return self.data[1, :]

    @property
    def locnums1(self):
        return self.locnums

    @property
    def locnums2(self):
        return (self.locnums + 1) % 3

    @property
    def locnums3(self):
        return (self.locnums + 2) % 3

    @property
    def nodnums1(self):
        return self.triangs[3 * self.trinums + self.locnums1]

    @property
    def nodnums2(self):
        return self.triangs[3 * self.trinums + self.locnums2]

    @property
    def nodnums3(self):
        return self.triangs[3 * self.trinums + self.locnums3]

    @property
    def triangs(self):
        return self.mesh.triangs.flat

    @property
    def nodnums_unique(self):
        return np.unique(
            np.r_[self.nodnums1, self.nodnums2]
        )

    @property
    def trinums_unique(self):
        return np.unique(self.trinums)

    @property
    def edges2d(self):
        return np.vstack(
            [self.nodnums1, self.nodnums2]
        )

    @property
    def nodes2d(self):
        return np.copy(
            self.mesh.points2d[:, self.nodnums1], order='C'
        )

    @property
    def nodes_complex(self):
        return _pack_complex(*self.nodes2d)

    def update_data(self, new_data):
        return self.from_data(
            self.mesh, np.copy(new_data, order='C')
        )


MeshEdgeError = type(
    'MeshEdgeError', (Exception,), {}
)


class MeshEdge(EdgeData):
    """Mesh edge object.

    Properties
    ----------

    Primary data:

    Name       | Description
    -----------|-----------------------------------
    `trinums`  | Host triangles (HT) across edges.
    `locnums`  | Local numbers of edges in HTs

    Nodes numbers:

    Name       | Description
    -----------|--------------------------------
    `nodnums1` | Start nodes across edges.
    `nodnums2` | End nodes across edges.
    `nodnums3` | Apexes across edges.

    Tables:

    Name      | Description
    ----------|-------------------------------------------
    `edges2d` | Edges as two rows of nodes numbers.
    `nodes2d` | Start nodes positions as two float rows.

    Unique values:

    Name             | Description
    -----------------|-----------------------------------
    `nodnums_unique` | Unique numbers of edge points.
    `trinums_unique` | Unique numbers of edge triangles.

    """

    MESH_EDGE_ERROR = MeshEdgeError

    @classmethod
    def from_mesh(cls, mesh):
        return EdgeFinder(mesh).get_mesh_edge()

    def sort_per_trinum(self):
        return self.permutedata(
            np.argsort(self.trinums)
        )

    def sort_per_nodnum(self):
        return self.permutedata(
            np.argsort(self.nodnums1)
        )

    def permutedata(self, permuter):
        return self.update_data(
            self.data[:, permuter]
        )

    def has_intersects(self):
        return np.unique(self.nodnums1).size != self.nedges

    def getloops(self):
        """Splits the mesh edge into loops.

        Returns
        -------
        list
            A list of `EdgeLoop` objects.

        """

        if self.has_intersects():
            raise MeshEdgeError(
                "cannot split the edge, self-intersections found"
            )

        return list(
            self.genloops()
        )

    def genloops(self):

        loops_as_inds = loops.LoopsAsInds(self.edges2d).getloops()

        for inds in loops_as_inds:
            yield EdgeLoop.from_data(
                self.mesh, self.data[:, inds].copy('C')
            )


EdgeLoopError = type(
    'EdgeLoppError', (Exception,), {}
)


class EdgeLoop(EdgeData):
    """Loop on the mesh edge.

    - Inherits `MeshEdge` properties.
    - Nodes are oriented in CCW order.

    """

    EDGE_LOOP_ERROR = EdgeLoopError

    @property
    def startnode(self):
        return self.nodnums1[0]

    def __contains__(self, nodenum):
        return nodenum in self.nodnums1

    def synctoedge(self, edgeind):
        """Synchronizes to the specified segment.

        Parameters
        ----------
        edgeind : int
            Index of the segment to synchronize to.

        Returns
        -------
        EdgeLoop
            New loop.

        """

        index = _normindex(
            edgeind, self.nedges
        )

        new_data = np.roll(
            self.data, -index, axis=1
        )

        return self.update_data(new_data)

    def synctonode(self, nodenum):
        """Synchronizes to the specified node.

        Parameters
        ----------
        nodenum : int
            Global node number to synchronize to.

        Returns
        -------
        EdgeLoop
            New loop.

        """

        edgeind = self.find_node(nodenum)

        if edgeind is not None:
            return self.synctoedge(edgeind)

        raise EdgeLoopError(
            f"node {nodenum} is not on the loop"
        )

    def find_node(self, nodenum):
        if nodenum in self.nodnums1:
            return np.argmax(
                self.nodnums1 == nodenum
            )
        return None


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


class EdgeFinder(MeshAgent):
    """Finder of the mesh edge.
    """

    def get_mesh_edge(self):

        paredges = self.make_edges_paired()
        meshedge = self.from_edges_paired(paredges)

        return meshedge

    def from_edges_paired(self, paredges):

        _, rows, cols = paredges.atrank(1)

        trinums = rows
        locnums = cols

        return MeshEdge.from_data(
            mesh=self.mesh, data=np.vstack([trinums, locnums])
        )

    def make_edges_paired(self):
        return tables.TableMap.from_table(self.edges_paired)


def _normindex(index, size):
    return min(
        abs(index), size - 1
    )


def _pack_complex(argx, argy):
    return argx + 1j * argy
