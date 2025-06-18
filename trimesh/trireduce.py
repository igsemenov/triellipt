# -*- coding: utf-8 -*-
"""Static mesh reducer.
"""
import numpy as np


class MeshAgent:
    """Operator on a trimesh.
    """

    def __init__(self, mesh=None):
        self.mesh = mesh
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    @property
    def mesh_twin(self):
        return self.mesh.twin()

    @property
    def mesh_suptriu(self):
        return self.mesh.supertriu()

    @property
    def mesh_suptriu_detached(self):
        return self.mesh_suptriu.detach()

    @property
    def supt_reduced(self):
        return self.mesh_suptriu.reduce(self.seed, iterate=True)

    @property
    def supt_detached_reduced(self):
        return self.mesh_suptriu_detached.reduce(self.seed, iterate=True)

    @property
    def seed(self):
        return self.cache['seed']


class MeshReduce(MeshAgent):
    """Performs a mesh-reduction event.
    """

    def reduced(self, shrink=None, detach=False, seed=None):
        """Tries to compress a triangle mesh.

        Parameters
        ----------
        shrink : int = None
            Controls shrinking of supertriu after compression.
        detach : bool = False
            Runs the edge detachment before compression, if True.
        seed : (float, float) = None:
            Seed point to start reduction.

        Returns
        -------
        TriMesh
            New mesh.

        """

        niters = shrink or 0
        self.cache['seed'] = seed

        suptriu = self.make_suptriu(niters, detach)
        trimesh = self.from_suptriu(suptriu)

        return trimesh

    def make_suptriu(self, niters, detach):

        suptriu = self.get_suptriu_primary(detach)
        suptriu = self.get_suptriu_shrinked(suptriu, niters)

        return suptriu

    def get_suptriu_primary(self, detach):
        if not detach:
            return self.supt_reduced
        return self.supt_detached_reduced

    def get_suptriu_shrinked(self, suptriu, niters):

        if suptriu is None:
            return None

        for _ in range(niters):
            suptriu = suptriu.strip().smooth()

        if self.is_empty_suptriu(suptriu):
            return None

        return suptriu

    def from_suptriu(self, suptriu):

        if suptriu is None:
            return self.mesh_twin

        imesh, omesh = self.gen_mesh_duet(suptriu)

        new_mesh = self.run_mesh_merge(imesh, omesh)
        return new_mesh

    def gen_mesh_duet(self, suptriu):

        inner_mesh = suptriu.supmesh
        outer_mesh = suptriu.mesh.deltriangs(*suptriu.supbodies)

        yield inner_mesh
        yield outer_mesh

    def run_mesh_merge(self, imesh, omesh):

        new_mesh = MeshMerge(imesh, omesh).mesh_merge()

        if new_mesh is None:
            return self.mesh_twin

        return new_mesh.delghosts()

    def is_empty_suptriu(self, suptriu):
        if suptriu is None:
            return True
        if suptriu.size == 0:
            return True
        return False


def merge_mesh(omesh, imesh):
    """Merges a mesh from an outer part and an inner parts.

    Parameters
    ----------
    omesh : TriMesh
        Outer mesh.
    imesh : TriMesh
        Inner mesh.

    Returns
    -------
    TriMesh | None
        New mesh or None, if failed.

    Notes
    -----

    Parts must be from the same point set.

    """
    return MeshMerge(imesh, omesh).mesh_merge().delghosts().twin()


MeshMergeError = type(
    'MeshMergeError', (Exception,), {}
)


class MeshDuet:
    """Duet of an inner and outer mesh.
    """

    def __init__(self, imesh, omesh):

        if not imesh.points is omesh.points:
            raise MeshMergeError(
                "got an invalid mesh pair, point sets differ"
            )

        self.imesh = imesh
        self.omesh = omesh


class MeshMerge(MeshDuet):
    """Merger of the duet meshes.
    """

    def mesh_merge(self):

        loops_duet = self.make_loops()

        if loops_duet is None:
            return None

        iloop, oloop = loops_duet
        return self.from_loops(iloop, oloop)

    def make_loops(self):
        return GetLoops.with_meshes(self.imesh, self.omesh).get_loops()

    def from_loops(self, iloop, oloop):
        if oloop.size != 2 * iloop.size:
            return None
        return LoopsMerge.with_loops(iloop, oloop).mesh_merge()


class GetLoops:
    """Fetcher of the duet loops.
    """

    def __init__(self, imesh, omesh):
        self.imesh = imesh
        self.omesh = omesh

    @classmethod
    def with_meshes(cls, imesh, omesh):
        return cls(imesh, omesh)

    def get_loops(self):

        oborder, iborder = self.make_borders()

        has_bad_borders = any(
            [oborder is None, iborder is None]
        )

        if has_bad_borders:
            return None

        has_multi_loops = any(
            [len(oborder) != 2, len(iborder) != 1]
        )

        if has_multi_loops:
            return None

        return self.from_borders(iborder, oborder)

    def from_borders(self, iborder, oborder):

        iloop = self.take_iloop(iborder)

        oloop = self.take_oloop(
            oborder, iloop.startnode
        )

        if oloop is None:
            return None
        return iloop, oloop

    def take_iloop(self, iborder):
        return iborder[0]

    def take_oloop(self, oborder, anchor):
        for loop in oborder:
            if anchor in loop:
                return loop.synctonode(anchor)
        return None

    def make_borders(self):
        yield self.get_oborder()
        yield self.get_iborder()

    def get_oborder(self):
        return _mesh_border(self.omesh)

    def get_iborder(self):
        return _mesh_border(self.imesh)


class LoopsDuet:
    """Duet of synchornized loops.
    """

    def __init__(self, iloop, oloop):
        self.iloop = iloop
        self.oloop = oloop

    @classmethod
    def with_loops(cls, iloop, oloop):
        return cls(iloop, oloop)

    @property
    def primary_points(self):

        if self.is_consistent_pair:
            return self.oloop.mesh.points

        raise MeshMergeError(
            "invalid mesh pair found, point sets differ"
        )

    @property
    def is_consistent_pair(self):
        return self.iloop.mesh.points is self.oloop.mesh.points


class LoopsMerge(LoopsDuet):
    """Merges parent meshes of the duet loops.
    """

    def mesh_merge(self):

        points, triangs = self.make_mesh_data()
        trimesh = self.from_mesh_data(points, triangs)

        return trimesh

    def make_mesh_data(self):
        yield self.get_points()
        yield self.get_triangs()

    def from_mesh_data(self, points, triangs):
        return self.oloop.mesh.from_data(points, triangs)

    def get_triangs(self):

        itris = self.iloop.mesh.triangs
        otris = self.oloop.mesh.triangs

        voids = self.get_contact_voids()

        data = [
            itris, voids, otris  # order matters
        ]

        return np.vstack(data)

    def get_points(self):

        nums, nodes = self.gen_contact_points()

        np.put(
            self.primary_points, nums, nodes
        )

        return self.primary_points

    def get_contact_voids(self):
        """Defines void triangles between synchronized loops.
        """

        hanging_nums = self.get_hanging_nums()

        tris = np.vstack(
            [self.iloop.nodnums2, self.iloop.nodnums1, hanging_nums]
        )

        return tris.T.copy('C')

    def gen_contact_points(self):
        """Defines contact points between synchronized loops.
        """
        yield self.get_hanging_nums()
        yield self.get_hanging_nodes()

    def get_hanging_nums(self):
        return np.flip(
            self.oloop.nodnums2[::2]
        )

    def get_hanging_nodes(self):
        return 0.5 * np.sum(
            self.iloop.mesh.points[self.iloop.edges2d], axis=0
        )


def _mesh_subpoints(mesh, trinums, locnums):
    return mesh.points[
        mesh.triangs.flat[3 * trinums + locnums]
    ]


def _mesh_border(mesh):

    edge = mesh.meshedge()

    try:
        return edge.getloops()
    except edge.MESH_EDGE_ERROR:
        return None
    return None
