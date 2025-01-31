# -*- coding: utf-8 -*-
"""Helping tools.
"""
import numpy as np
from triellipt.fem import femoprs
from triellipt.utils import pairs, tables


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


def reconnect_nodes(mesh, nodes1, nodes2):
    """Replaces the nodes2 by the nodes1 in a mesh.
    """

    newnodes = np.arange(mesh.npoints)
    newnodes[nodes2] = nodes1

    new_triangs = newnodes[mesh.triangs]
    return mesh.update_triangs(new_triangs)


def mesh_areas(mesh):
    """Computes the area over mesh triangles.
    """
    return femoprs.mesh_metric(mesh).areas1d.flatten()


def clean_twin_voids(mesh):
    """Removes twin voids from a mesh.
    """

    trinums = mesh.getvoids()

    if trinums.size == 0:
        return mesh

    codes = pairs.paircols(
        mesh.triangs[trinums, 0:2]
    )

    _, twinums, _ = tables.maptable(codes).atrank(2)

    if twinums.size == 0:
        return mesh

    return mesh.deltriangs(*trinums[twinums])


def clean_voids_on_edge(mesh, edge):
    """Removes voids that are on the edge.
    """

    trinums = mesh.getvoids()

    if trinums.size == 0:
        return mesh

    voids_pivots = mesh.triangs[trinums, 2]

    on_edge = np.isin(
        voids_pivots, edge.nodnums_unique
    )

    bad_trinums = trinums[on_edge]

    if bad_trinums.size == 0:
        return mesh

    return mesh.deltriangs(*bad_trinums)
