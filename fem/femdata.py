# -*- coding: utf-8 -*-
"""FEM unit data.
"""
import numpy as np
from triellipt.fem import joints, ijstream, vstream


def get_fem_data(mesh, anchors=None):
    """Makes FEM unit data.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    anchors : tuple = None
        Nodes numbers to sync the mesh boundary.

    Returns
    -------
    dict
        Attributes of a FEM unit.

    """
    _ = FEMFactory()
    return _.get_data(mesh, anchors)


class FEMFactory:
    """Maker of FEM unit data.
    """

    def __init__(self):
        self.cache = {}

    def get_data(self, mesh, anchors=None) -> dict:

        mesh = self.get_mesh_aligned(mesh, anchors)
        data = self.get_data_from_mesh(mesh)

        return data

    def get_mesh_aligned(self, mesh, anchors):

        anchors = anchors or ()

        mesh0 = mesh

        mesh1 = mesh0.alignnodes(*anchors)
        mesh2 = mesh1.alignvoids()

        perm1 = mesh1.meta['nodes-permuter']
        perm2 = mesh2.meta['nodes-permuter']

        self.cache['nodes-permuter'] = perm1[perm2]
        return mesh2

    def get_data_from_mesh(self, mesh):

        skel = self.make_skeleton(mesh)
        data = self.from_skeleton(skel)

        return data

    def make_skeleton(self, mesh):
        return joints.getskeleton(mesh)

    def from_skeleton(self, skel):

        data = self.make_mesh_data(skel)
        ij_v = self.make_ij_v_data(skel)

        ijvs = {
            'ij-data': ij_v['ij'],
            'femoprs': ij_v['v']
        }

        perm = {
            'nodes-mesh-to-unit': self.cache['nodes-permuter']
        }

        meta = {
            'masks': data['masks'],
            'loops': data['loops'],
            'femat': ijvs,
            'perms': perm
        }

        mesh = data['mesh']

        return {
            'mesh': mesh,
            'meta': meta
        }

    def make_mesh_data(self, skel):
        return MeshMaker(skel).get_mesh_data()

    def make_ij_v_data(self, skel):
        return {
            'ij': ijstream.getstream(skel), 'v': vstream.getstream(skel)
        }


class SkelAgent:
    """Operator on a skeleton.
    """

    def __init__(self, skel):
        self.skel = skel

    @property
    def hasjoints(self):
        return self.skel.hasjoints


class MeshMaker(SkelAgent):
    """Makes a mesh from a skeleton.
    """

    def get_mesh_data(self):

        mesh = self.make_mesh()
        masks = self.make_masks()
        loops = self.make_loops(mesh)

        return {
            'mesh': mesh, 'masks': masks, 'loops': loops
        }

    def make_mesh(self):

        if not self.hasjoints:
            return self.skel.body

        triangs_joints = [
            self.skel.wests.triangs,
            self.skel.easts.triangs,
            self.skel.cores.triangs,
            self.skel.voids.triangs
        ]

        triangs_stack = [
            self.skel.body.triangs, *triangs_joints
        ]

        new_triangs = np.vstack(triangs_stack)
        return self.skel.body.update_triangs(new_triangs)

    def make_masks(self):

        body_size = self.skel.body.size

        if not self.skel.hasjoints:
            return {
                'body': (0, body_size), 'joints': {}
            }

        joints_sizes = [
            self.skel.wests.size,
            self.skel.easts.size,
            self.skel.cores.size,
            self.skel.voids.size
        ]

        sizes = [
            body_size, *joints_sizes
        ]

        bins = np.add.accumulate(sizes)

        joints_slices = {
            'wests': (bins[0], bins[1]),
            'easts': (bins[1], bins[2]),
            'cores': (bins[2], bins[3]),
            'voids': (bins[3], bins[4])
        }

        return {
            'body': (0, body_size), 'joints': joints_slices
        }

    def make_loops(self, mesh):
        return mesh.meshedge().getloops()
