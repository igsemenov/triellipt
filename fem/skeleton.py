# -*- coding: utf-8 -*-
"""Making a mesh skeleton.
"""
import numpy as np
from triellipt.utils import tables


def getskeleton(mesh):
    """Converts a mesh into a mesh skeleton.
    """
    return MeshSkeleton.from_mesh(mesh)


MeshSkeletonError = type(
    'MeshSkeletonError', (Exception,), {}
)


class MeshSkeleton:
    """Mesh in skeleton format.

    Attributes
    ----------
    nodesmap : NodesMap
        Map of nodes that are not pivots.
    voidsmap : NodesMap
        Map of the void pivots.

    """

    def __init__(self, nodesmap=None, voidsmap=None):
        self.nodesmap = nodesmap
        self.voidsmap = voidsmap

    @classmethod
    def from_mesh(cls, mesh):
        return SkeletonMaker.from_mesh(mesh).get_skeleton()

    @property
    def mesh(self):
        return self.nodesmap.mesh

    @property
    def meta(self):
        if not self.hasvoids:
            return {}
        return self.voids_meta

    @property
    def hasvoids(self):
        return self.voidsmap is not None

    @property
    def voids_meta(self):
        if not self.hasvoids:
            return {}
        return {
            'voids-trinums': self.voidsmap.meta['voids-trinums'],
            'voids-triangs': self.voidsmap.meta['voids-triangs']
        }

    @property
    def voids_trinums(self):
        return self.meta.get('voids-trinums')

    @property
    def voids_triangs(self):
        return self.meta.get('voids-triangs')

    def voids_submaps(self):

        if not self.hasvoids:
            return {}

        return {
            'westmap': _slice_nodesmap(self.voidsmap, np.s_[1::4]),
            'coremap': _slice_nodesmap(self.voidsmap, np.s_[2::4]),
            'eastmap': _slice_nodesmap(self.voidsmap, np.s_[3::4])
        }


class MeshAgent:
    """Operator on a triangle mesh.
    """

    def __init__(self, mesh):
        self.mesh = mesh
        self.meta = self.fetch_meta()
        self.cache = {}

    @classmethod
    def from_mesh(cls, mesh):
        return cls(mesh)

    def fetch_meta(self):
        return {
            'nodes': self.fetch_nodes_meta(),
            'voids': self.fetch_voids_meta()
        }

    def fetch_nodes_meta(self):
        return {
            'nodesmap': self.mesh.nodesmap()
        }

    def fetch_voids_meta(self):

        trinums = self.mesh.getvoids()
        triangs = self.mesh.triangs[trinums, :]

        if triangs.size == 0:
            return None

        pivots_nums = triangs[:, 2]
        pivots_sort = np.argsort(pivots_nums)

        trinums_sorted = trinums[pivots_sort]
        triangs_sorted = triangs[pivots_sort, :]

        return {
            'trinums': trinums_sorted.copy('C'),
            'triangs': triangs_sorted.copy('C')
        }

    @property
    def hasvoids(self):
        return self.meta['voids'] is not None

    @property
    def nodesmap(self):
        return self.meta['nodes']['nodesmap']

    @property
    def voids_triangs(self):
        return self.meta['voids']['triangs']

    @property
    def voids_trinums(self):
        return self.meta['voids']['trinums']

    @property
    def voids_pivots(self):
        return self.voids_triangs[:, 2]

    @property
    def voids_wests(self):
        return self.voids_triangs[:, 1]

    @property
    def voids_easts(self):
        return self.voids_triangs[:, 0]


class SkeletonMaker(MeshAgent):
    """Converts a mesh to a mesh skeleton.
    """

    def get_skeleton(self):

        if not self.hasvoids:
            return MeshSkeleton(
                nodesmap=self.nodesmap, voidsmap=None
            )

        voidsmap = self.make_voidsmap()
        nodesmap = self.make_nodesmap()

        return MeshSkeleton(nodesmap, voidsmap)

    def make_voidsmap(self):

        _ = self.test_pivots_rank_is_four()
        _ = self.mask_pivots_data_in_map()

        _ = self.push_voidsmap()
        _ = self.sync_voidsmap()

        voids_meta = {
            'voids-trinums': self.meta['voids']['trinums'],
            'voids-triangs': self.meta['voids']['triangs']
        }

        voidsmap = _
        voidsmap = voidsmap.add_meta(voids_meta)

        return voidsmap

    def sync_voidsmap(self):
        return self.voids_orienter.get_voids_oriented()

    def push_voidsmap(self):

        pivs_mask = self.cache['pivs-mask']
        pivs_data = self.nodesmap.data.data[:, pivs_mask]

        voids_map = self.nodesmap.update_table_data(pivs_data)

        self.cache['voidsmap'] = voids_map
        return voids_map

    def mask_pivots_data_in_map(self):

        inds = self.nodesmap.data.packs_fronts[self.voids_pivots]

        inds = np.vstack(
            [inds, inds + 1, inds + 2, inds + 3]
        )

        mask = np.full(
            self.nodesmap.data.datasize, False
        )

        mask[inds.T.flat] = True

        self.cache['pivs-mask'] = mask
        return mask

    def test_pivots_rank_is_four(self):

        pivs_ranks = self.nodesmap.nodes_ranks[self.voids_pivots]
        pivs_ranks_set = np.unique(pivs_ranks).tolist()

        if len(pivs_ranks_set) != 1:
            raise MeshSkeletonError(
                f"got pivots with multiple ranks: {pivs_ranks_set}"
            )

        pivs_rank = pivs_ranks_set[0]

        if pivs_rank != 4:
            raise MeshSkeletonError(
                f"expected pivots to be of rank 4, got {pivs_rank}"
            )

        return True

    def make_nodesmap(self):

        pivs_mask = self.cache['pivs-mask']
        nodes_mask = np.logical_not(pivs_mask)
        nodes_data = self.nodesmap.data.data[:, nodes_mask]

        nodes_map = self.nodesmap.update_table_data(nodes_data)

        self.cache['nodesmap'] = nodes_map
        return nodes_map

    @property
    def voids_orienter(self):
        return VoidsOrienter.from_agent(self)


class MeshSubAgent:
    """Operator on a mesh agent.
    """

    def __init__(self, agent):
        self.agent = agent
        self.cache = {}

    @classmethod
    def from_agent(cls, agent):
        return cls(agent)


class VoidsOrienter(MeshSubAgent):
    """Voids map oriented.
    """

    VOID_RANK = 4

    SYNCORDER = (
        'voids', 'wests', 'cores', 'easts'
    )

    @property
    def voidsmap(self):
        return self.agent.cache['voidsmap']

    def get_voids_oriented(self):

        _ = self.make_sync_cast()
        _ = self.sync_voids_map()

        return _

    def sync_voids_map(self):

        data_synced = self.make_data_synced()

        new_trinums = data_synced['trinums-synced']
        new_locnums = data_synced['locnums-synced']

        new_data = [
            self.voidsmap.nodnums, new_trinums, new_locnums
        ]

        new_data = _pack_rows(*new_data)
        return self.voidsmap.update_table_data(new_data)

    def make_data_synced(self):

        cast = self.cache['sync-cast']

        trinums2d = _as4cols(self.voidsmap.trinums)
        locnums2d = _as4cols(self.voidsmap.locnums)

        trinums2d_synced = tables.table_image(trinums2d, cast)
        locnums2d_synced = tables.table_image(locnums2d, cast)

        return {
            'trinums-synced': trinums2d_synced.flatten(),
            'locnums-synced': locnums2d_synced.flatten()
        }

    def make_sync_cast(self):

        start = self.find_start_ind()
        map12, map21 = self.make_cross_maps()

        ind_0 = start
        ind_1 = _cols_from_table(map21, start)
        ind_3 = _cols_from_table(map12, start)
        ind_2 = 6 - ind_0 - ind_1 - ind_3

        cast = _pack_cols(
            ind_0, ind_1, ind_2, ind_3
        )

        self.cache['sync-cast'] = cast
        return cast

    def find_start_ind(self):

        mask = np.isin(
            self.voidsmap.trinums, self.agent.voids_trinums
        )

        inds, = np.where(mask)
        inds = inds % 4

        return inds

    def make_cross_maps(self):

        nodes1 = _as4cols(self.voidsmap.nodnums1)
        nodes2 = _as4cols(self.voidsmap.nodnums2)

        arg1 = np.argsort(nodes1, axis=1)
        arg2 = np.argsort(nodes2, axis=1)

        map12 = tables.table_image(
            arg2, np.argsort(arg1, axis=1)
        )

        map21 = tables.table_image(
            arg1, np.argsort(arg2, axis=1)
        )

        yield map12
        yield map21


def _as4cols(data):
    return data.reshape(data.size // 4, 4)


def _cols_from_table(table, cols):
    return table[
        np.arange(table.shape[0]), cols
    ]


def _pack_rows(*rows):
    return np.vstack(rows).copy('C')


def _pack_cols(*cols):
    return np.vstack(cols).T.copy('C')


def _slice_nodesmap(nodesmap, indexer):
    new_data = nodesmap.data.data[:, indexer]
    return nodesmap.update_table_data(new_data)
