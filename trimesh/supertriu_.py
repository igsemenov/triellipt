# -*- coding: utf-8 -*-
"""Super triangulation.
"""
import numpy as np
from triellipt.utils import tables
from triellipt.trimesh import superoprs


class SuperData:
    """Super triangulation data.
    """

    def __init__(self, mesh=None, data=None):
        self.mesh = mesh
        self.data = data

    @classmethod
    def from_mesh(cls, mesh):
        return Triangler(mesh).getsuptriu()

    @property
    def size(self):
        return self.data.shape[1]

    @property
    def trinums(self):
        return self.data[0, :]

    @property
    def trinums1(self):
        return self.data[1, :]

    @property
    def trinums2(self):
        return self.data[2, :]

    @property
    def trinums3(self):
        return self.data[3, :]

    @property
    def nodnums1(self):
        return self.data[4, :]

    @property
    def nodnums2(self):
        return self.data[5, :]

    @property
    def nodnums3(self):
        return self.data[6, :]

    @property
    def suptriangs(self):
        return _packcols(
            self.nodnums1, self.nodnums2, self.nodnums3
        )

    @property
    def supbodies(self):
        return _packcols(
            self.trinums, self.trinums1, self.trinums2, self.trinums3
        )

    @property
    def kermesh(self):
        return self.mesh.submesh(*self.trinums)

    @property
    def supmesh(self):
        return self.mesh.update_triangs(self.suptriangs)

    def update_data(self, new_data):
        return self.__class__(
            self.mesh, new_data.copy('C')
        )

    def is_compact(self):
        return _is_compact_table(self.supbodies)

    def subtriu(self, *suptri_inds):
        return self.update_data(
            self.data[:, suptri_inds]
        )

    def deltriangs(self, *suptri_inds):

        new_data = np.delete(
            self.data, suptri_inds, axis=1
        )

        return self.update_data(new_data)

    def supvoids(self):
        """Fetches a supertriu made of voids.
        """
        return superoprs.SupVoids.from_suptriu(self).supvoids()

    def atcores(self, *trinums):
        """Fetches a supertriu with the specified core triangles.
        """

        to_delete = np.isin(
            self.trinums, trinums
        )

        return self.update_data(
            np.compress(to_delete, self.data, axis=1)
        )

    def find_seed(self, seed):
        """Find super-triangle closest to the seed.

        Parameters
        ----------
        seed : (float, float)
            Seed point to find.

        """

        centrs = self.supmesh.centrs_complex

        number = np.argmin(
            np.abs(centrs - (seed[0] + seed[1] * 1j))
        )

        return number


class SuperTriu(SuperData):
    """Super triangulation.

    Properties
    ----------

    - Hosts are background triangles with three neighbors.
    - Super-triangles are made from apexes of hosts neighbors.

    General:

    Name       | Description
    -----------|----------------------------
    `trinums`  | Numbers of host triangles.
    `kermesh`  | Submesh of host triangles.
    `supmesh`  | Mesh from super-vertices.

    Neighbors:

    Name       | Description
    -----------|-------------------
    `trinums1` | 1st CCW neighbor. 
    `trinums2` | 2nd CCW neighbor.
    `trinums3` | 3rd CCW neighbor.

    Super-vertices:

    Name       | Description
    -----------|-------------------
    `nodnums1` | 1st CCW vertex. 
    `nodnums2` | 2nd CCW vertex.
    `nodnums3` | 3rd CCW vertex.

    Notes
    -----

    Non-standard nodes pairing for `supmesh-kermesh` transition:

    - (0, 1) → 1
    - (1, 2) → 2
    - (2, 0) → 0

    """

    def strip(self):
        """Remove links from a super-triangulation.

        Notes
        -----

        See `EdgesMap.getspec()` for links definition.

        """

        if self.size == 0:
            return self

        cleaner = superoprs.SupStrip.from_suptriu(self)
        return cleaner.cleaned()

    def smooth(self, iterate=True):
        """Removes heads and spots from a super-triangulation.

        Parameters
        ----------
        iterate : bool = True
            Runs smoothing until possible, if True.

        Notes
        -----

        See `EdgesMap.getspec()` for heads and spots definition.

        """

        if self.size == 0:
            return self

        cleaner = superoprs.SupSmooth.from_suptriu(self)

        if iterate is False:
            return cleaner.cleaned()
        return cleaner.cleaned_iteraly()

    def detach(self):
        """Removes super-triangles touching the background mesh edge.
        """

        if self.size == 0:
            return self

        cleaner = superoprs.SupDetach.from_suptriu(self)
        return cleaner.cleaned()

    def reduce(self, seed=None, iterate=True):
        """Extracts a compact super-triangulation, if possible.

        Parameters
        ----------
        seed : (float, float) = None
            Seed point to start reduction (b).
        iterate : bool = True
            Triggers cleaning and retry in case of failure (a).

        Returns
        -------
        SuperTriu | None
            Compact super-triangulation or None, if failed.

        Notes
        -----

        (a) Cleaning is a strip-and-smooth action.

        (b) First super-triangle is used by default.

        """

        if self.size == 0:
            return self

        reducer = superoprs.SupReduce.from_suptriu(self)

        if iterate is False:
            return reducer.compress(seed)
        return reducer.reduce(seed)


class Triangler:
    """Super triangulator.
    """

    APEXES_OVER_EDGES = np.r_[2, 0, 1]

    def __init__(self, mesh):
        self.edges = EdgesMap.from_mesh(mesh)

    def getsuptriu(self):

        datas = self.get_neighbors()
        verts = self.get_vertices(datas)

        return self.get_obj(datas, verts)

    def get_obj(self, neighbors_data, vertices):

        host_tris_nums = neighbors_data['host-tris-nums']
        neighbors_nums = neighbors_data['neighbors-nums']

        data = [
            host_tris_nums, neighbors_nums.T, vertices.T
        ]

        return SuperTriu(
            self.edges.mesh, np.vstack(data)
        )

    def get_vertices(self, neighbors_data):

        neighbors_nums = neighbors_data['neighbors-nums']
        neighbors_locs = neighbors_data['neighbors-locs']

        verts_ccw_locs = self.APEXES_OVER_EDGES[neighbors_locs]

        return self.from_mesh_triangs(
            3 * neighbors_nums + verts_ccw_locs
        )

    def get_neighbors(self):

        hosttris2d, pairnums2d, inpairs2d = self.triconnects()

        hostriang_inds = 2 * pairnums2d + (inpairs2d + 0)
        neighbors_inds = 2 * pairnums2d + (inpairs2d + 1) % 2

        neighbors_nums = self.from_trinums_table(neighbors_inds)
        neighbors_locs = self.from_locnums_table(neighbors_inds)
        neighbors_ords = self.from_locnums_table(hostriang_inds)

        ccwsorts = tables.trisorts(neighbors_ords)

        neighbors_nums = tables.table_image(neighbors_nums, ccwsorts)
        neighbors_locs = tables.table_image(neighbors_locs, ccwsorts)

        host_tris_nums = hosttris2d[:, 0].copy('C')

        return {
            'host-tris-nums': host_tris_nums,
            'neighbors-nums': neighbors_nums,
            'neighbors-locs': neighbors_locs
        }

    def triconnects(self):
        for data in self.secondary_edgesmap.atrank(3):
            yield np.reshape(
                data, (data.size // 3, 3)
            )

    def from_trinums_table(self, indexer):
        return self.edges.trinums2d.flat[indexer]

    def from_locnums_table(self, indexer):
        return self.edges.locnums2d.flat[indexer]

    def from_mesh_triangs(self, indexer):
        return self.edges.mesh_triangs.flat[indexer]

    @property
    def secondary_edgesmap(self):
        return self.edges.secondary_edgesmap()


class EdgesMap:
    """Temporary edges map.
    """

    def __init__(self, mesh, trinums2d, locnums2d):

        self.mesh = mesh

        self.trinums2d = trinums2d
        self.locnums2d = locnums2d

    @classmethod
    def from_mesh(cls, mesh):

        edges = mesh.edgesmap()

        trinums2d = _packcols(edges.trinums1, edges.trinums2)
        locnums2d = _packcols(edges.locnums1, edges.locnums2)

        return cls(
            mesh, trinums2d, locnums2d
        )

    @property
    def mesh_triangs(self):
        return self.mesh.triangs

    def secondary_edgesmap(self):
        return tables.TableMap.from_table(self.trinums2d)


def _packcols(*cols):
    return np.vstack(cols).T.copy('C')


def _is_compact_table(table):
    return np.unique(table).size == table.size
