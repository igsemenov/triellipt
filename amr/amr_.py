# -*- coding: utf-8 -*-
"""AMR public tools.
"""
from triellipt.amr import trirefine_body
from triellipt.amr import trisubsets
from triellipt.amr import masscollect


def refine_mesh(mesh, trinums=None):
    """Performs static mesh refinement.

    Parameters
    ----------
    mesh : TriMesh
        Triangle mesh to refine.
    trinums : Iterable[int] = None
        Numbers of triangles to refine, if None takes all triangles.

    Returns
    -------
    TriMesh
        Refined mesh (a).

    Notes
    -----

    - (a) The data refiner is included in the mesh metadata.
    - (b) The neighborhood of the voids is not refined.

    """

    if trinums is None:
        trinums = range(mesh.size)

    return trirefine_body.refine_mesh(mesh, trinums)


def mesh_subset(mesh, count, anchor):
    """Finds a convex subset of a mesh.

    Parameters
    ----------
    mesh : TriMesh
        Input triangle mesh.
    count : int
        Seed number of triangles in a subset.
    anchor : (float, float)
        Anchor point to find a starting triangle.

    Returns
    -------
    flat-int-array
        Numbers of triangles in a subset.

    """
    return trisubsets.find_subset_trinums(mesh, count, anchor)


def mass_collector(mesh1, mesh2, mode):
    """Creates a mass collector for a master-slave pair of meshes.

    - Only for meshes with no hanging nodes.
    - No scaling by the area in a structured mode.

    Parameters
    ----------
    mesh1 : TriMesh
        Master (parent) mesh.
    mesh2 : TriMesh
        Mesh obtained by the refinement of the parent mesh.
    mode : str
        Defines the type of a collector (i).

    Returns
    -------
    MassCollector
        Mass collector object.

    Notes
    -----

    (i) Collector types:

    - "scaled" — for general unstructured meshes
    - "structed" — for meshes with equal triangles

    """

    if mode == "structed":
        return masscollect.get_collector_structed(mesh1, mesh2)
    if mode == "scaled":
        return masscollect.get_collector_scaled(mesh1, mesh2)

    raise ValueError(
        f"got an unknown collector type '{mode}'"
    )
