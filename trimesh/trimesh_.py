# -*- coding: utf-8 -*-
"""Triangle mesh object.
"""
import numpy as np
from triellipt.trimesh.utils import pairs
from triellipt.trimesh import meshedge_
from triellipt.trimesh import edgesmap_
from triellipt.trimesh import nodesmap_
from triellipt.trimesh import delghosts_
from triellipt.trimesh import supertriu_
from triellipt.trimesh import trireduce
from triellipt.trimesh import trisplit
from triellipt.trimesh import renumer


class TriData:
    """Root of a trimesh.
    """

    def __init__(self, points=None, triangs=None):
        self.points = points
        self.triangs = triangs

    @property
    def size(self):
        return self.triangs.shape[0]

    @classmethod
    def from_data(cls, points, triangs):
        return cls(
            points.copy('C'), triangs.copy('C')
        )

    def update_points(self, new_points):
        return self.__class__(
            new_points, self.triangs
        )

    def update_triangs(self, new_triangs):
        return self.__class__(
            self.points, new_triangs
        )

    def addtriangs(self, new_triangs):
        return self.update_triangs(
            np.vstack([self.triangs, new_triangs])
        )

    @classmethod
    def from_mesh_dict(cls, mesh_dict):

        points = mesh_dict['nodes']
        triangs = mesh_dict['elements']

        return cls.from_data(
            _table_to_complex(points), triangs
        )

    @property
    def points_complex(self):
        return self.points

    @property
    def centrs_complex(self):
        return (1./3.) * np.sum(
            self.points[self.triangs], axis=1
        )

    @property
    def points2d(self):
        return _complex_to_rows(self.points_complex)

    @property
    def centrs2d(self):
        return _complex_to_rows(self.centrs_complex)

    @property
    def triu(self):
        return *self.points2d, self.triangs

    def edges_paired(self):
        """Triangle edges as symmetrically paired nodes.
        """
        return pairs.paircols(
            self.triangs[:, [0, 1, 2, 0]]
        )

    def hasghosts(self) -> bool:
        """Checks for ghosts.
        """
        return delghosts_.HasGhosts(self).hasghosts()

    def getghosts(self):
        """Finds ghosts numbers, if any.
        """
        return delghosts_.GetGhosts(self).getghosts()

    def getvoids(self):
        """Finds empty triangles (voids).

        Returns
        -------
        flat-int-array
            Numbers of empty triangles.

        """
        return trisplit.GetVoids(self).find_voids()

    def hasvoids(self):
        return self.getvoids().size != 0

    def renumed(self, permuter):
        return renumer.Renumer(self).renumed(permuter)

    def save(self, file) -> None:
        """Saves the mesh to `.npz` file.

        Parameters
        ----------
        file : str or file-like
            Same as for `numpy.savez()`.

        """
        np.savez(
            file, points=self.points, triangs=self.triangs
        )

    @classmethod
    def load(cls, file):
        """Loads the mesh from `.npz` file.

        Parameters
        ----------
        file : str or file-like
            Same as for `numpy.load()`.

        Returns
        -------
        TriMesh
            New mesh.

        """

        data = np.load(file)

        return cls.from_data(
            data['points'], data['triangs']
        )


class TriMesh(TriData):
    """Triangle mesh.

    Attributes
    ----------
    points : flat-complex-array
        Mesh points in a complex plane.
    triangs : 3-column-int-table
        Triangles vertices in CCW order.

    Properties
    ----------

    Name       | Description
    -----------|--------------------------------------
    `triu`     | Generator of triplot arguments.
    `points2d` | Mesh points as two float rows.
    `centrs2d` | Triangle centers as two float rows.

    """

    @property
    def npoints(self):
        return self.points.size

    @property
    def ntriangs(self):
        return self.triangs.shape[0]

    @property
    def nnodes(self):
        return self.nodes_range.size

    @property
    def nodes_range(self):
        return np.unique(self.triangs).astype(int)

    def __mul__(self, value):
        return self.from_data(
            self.points * value, self.triangs
        )

    def __add__(self, value):
        return self.from_data(
            self.points + value, self.triangs
        )

    def __sub__(self, value):
        return self.from_data(
            self.points - value, self.triangs
        )

    def submesh(self, *trinums):
        """Extracts a submesh.

        Parameters
        ----------
        trinums : *int
            Numbers of triangles to include.

        Returns
        -------
        TriMesh
            New mesh object.

        """

        inds = _norminds(
            trinums, self.ntriangs
        )

        new_triangs = self.triangs[inds, :]
        return self.update_triangs(new_triangs)

    def deltriangs(self, *trinums):
        """Removes triangles from the mesh.

        Parameters
        ----------
        trinums : *int
            Numbers of triangles to delete.

        Returns
        -------
        TriMesh
            New mesh object.

        """

        inds = _norminds(
            trinums, self.ntriangs
        )

        new_triangs = np.delete(
            self.triangs, inds, axis=0
        )

        return self.update_triangs(new_triangs)

    def meshedge(self):
        """Extracts the mesh edge.

        Returns
        -------
        MeshEdge
            Mesh edge object.

        """
        return meshedge_.MeshEdge.from_mesh(self)

    def edgesmap(self):
        """Maps inner mesh edges.

        Returns
        -------
        EdgesMap
            Map of inner mesh edges.

        """
        return edgesmap_.EdgesMap.from_mesh(self)

    def nodesmap(self):
        """Maps nodes to hosting triangles.

        Returns
        -------
        NodesMap
            Nodes-to-triangles map.

        """
        return nodesmap_.NodesMap.from_mesh(self)

    def delghosts(self):
        """Removes ghost points from the mesh.

        Returns
        -------
        TriMesh
            New mesh.

        Notes
        -----

        Related methods:

        - `.hasghosts()` shows if there are any ghosts
        - `.getghosts()` returns ghost numbers, if any

        """
        return delghosts_.DelGhosts(self).cleanmesh()

    def supertriu(self):
        """Creates a super triangulation.

        Returns
        -------
        SuperTriu
            Super triangulation object.

        """
        return supertriu_.SuperTriu.from_mesh(self)

    def reduced(self, shrink=None, detach=False):
        """Tries to compress the mesh.

        Parameters
        ----------
        shrink : int = None
            Controls shrinking of a super-triangulation (a).
        detach : bool = False
            Runs the edge detachment before compression, if True.

        Returns
        -------
        TriMesh
            New mesh (b).

        Notes
        -----

        - (a) Number of shrinking steps after compression.
        - (b) Has unchanged mesh data if compression failed.

        """
        reducer = trireduce.TriReduce.from_mesh(self)
        return reducer.reduced(shrink, detach)

    def layout(self, *anchors):
        """Numbers mesh points in the edge-core order.

        Parameters
        ----------
        anchors : *int
            Possible anchors on the mesh edge.

        Returns
        -------
        TriMesh | None
            New mesh or None, if failed.

        """
        return renumer.Arranger(self).arrange(*anchors)

    def split(self):
        """Splits the mesh into homogeneous parts.

        Returns
        -------
        list
            List of homogeneous submeshes.

        """
        return trisplit.TriSplit(self).split()


def _table_to_complex(table):
    return table[:, 0] + 1j * table[:, 1]


def _complex_to_rows(data):
    return np.vstack(
        [data.real, data.imag]
    )


def _norminds(inds, size):

    inds = np.array(inds)

    mask = np.logical_and(
        inds >= 0, inds < size
    )

    return inds[mask].copy('C')
