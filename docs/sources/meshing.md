<!--
{
  "webtitle": "Meshing — triellipt documentation",
  "doctitle": "triellipt — Meshing",
  "codeblocks": true
}
-->

This page introduces the tools for working with meshes in **triellipt**.

In what follows, we assume the following import is used:

```python
import triellipt as tri
```

## Mesh object

The [trimesh.TriMesh](triellipt.trimesh.md#triellipttrimesh) class is responsible for holding and managing the triangular mesh.

The primary mesh attributes are:

- the triangulation matrix,
- the underlying point set.

Refer also to additional data structures:

- [triellipt.MeshEdge](triellipt.trimesh.md#meshedge)
- [triellipt.EdgesMap](triellipt.trimesh.md#edgesmap)
- [triellipt.NodesMap](triellipt.trimesh.md#nodesmap)
- [triellipt.SuperTriu](triellipt.trimesh.md#supertriu-1)

Note that the mesh object supports arithmetic operations with float and complex numbers, which modify the mesh points and produce a new mesh.

### Nonconformity

The nonconforming element contacts are introduced in the mesh as zero-area triangles.

Assume we have a nonconforming contact that spans the nodes A, B, and C:

```
    ●
   / \
  /   \
 B--C--A
 \ / \ /
  ●---●
   \ /
    ●
```

This contact appears in the triangulation matrix as an additional triangle:

```
│ ... │
│ ABC │
│ ... │
```

This triangle has zero area, but it restores the combinatorial mesh connectivity.

Terminology:

- Empty triangles are called *void elements*.
- The third vertex of a void element is called the *void pivot*.

### Ghost nodes

The mesh may contain *ghost nodes* — nodes that are not involved in triangles.

Theses nodes can be deleted from the mesh, see [triellipt.TriMesh.delghosts()](triellipt.trimesh.md#delghosts).

### Nodes alignment

The mesh nodes can be aligned in an edge–core manner:

- Edge nodes are placed at the top of the node numbering.
- Edge numbering can be configured to start at specific points — *anchors*.

For details, see [triellipt.TriMesh.alignnodes()](triellipt.trimesh.md#alignnodes).

## Mesh generation

### Geometry

The module `geom` provides tools for creating `.geo` files for use with Gmsh.

Usual workflow:

- You create a loop of curves.
- The loop is split into a colored path (per-curve coloring).
- The path is dumped into a `.geo` file.

For example we first create four connected lines:

```python
line1 = tri.geom.line(0 + 0j, 1 + 0j)
line2 = tri.geom.line(1 + 0j, 1 + 1j)
line3 = tri.geom.line(1 + 1j, 0 + 1j)
line4 = tri.geom.line(0 + 1j, 0 + 0j)
```

We then combine these lines into a loop:

```python
loop = tri.geom.makeloop(line1, line2, line3, line4)
```

and discretize the loop:

```python
path = loop.discretize((10, 1))
```

For more details, see [triellipt.geom.CurvesLoop.discretize()](triellipt.geom.md#discretize).

The resuting path can be dumped to a `.geo` file:

```python
path.togeo(abspath_to_geo, mesh_seeds)
```

For more details, see [triellipt.geom.PathMap.togeo()](triellipt.geom.md#togeo).

### Gmsh Meshes

The module `mshread` provides an interface for reading Gmsh meshes.

Usual workflow:

- Create a reader from a forlder with meshes.
- Read a mesh per the full name with `.msh` extension.

For example we create a reader:

```python
reader = tri.mshread.getreader(
  absolute_path_to_folder_with_meshes
)
```

and read the mesh:

```python
mesh = reader.read_mesh('some-mesh.msh')
```

The list of available meshes can be obtained as follows:

```python
read.listmeshes()
```

### Structured Grids

The module `mesher` provides two functions for generating structured grids:

- The grid-based meshes are generated using [triellipt.mesher.trigrid()](triellipt.mesher.md#trigrid).
- The triangle lattices are generated using [triellipt.mesher.trilattice()](triellipt.mesher.md#trilattice).

### Reduction

The mesh can potentially be coarsened in the bulk region — see [triellipt.TriMesh.reduced()](triellipt.trimesh.md#reduced).

This procedure can produce nonconforming meshes with multiple discretization layers.

### Composition

Meshes can also be joined along a common boundary, if present.

This provides an opportunity to construct composite meshes.

For more details, refer to [triellipt.amr.join_meshes()](triellipt.amr.md#join_meshes).

## Mesh adaptivity

### AMR unit

The module `amr` provides tools for mesh refinement and coarsening.

The functionality is provided by the AMR unit, which is built upon the input mesh:

```python
unit = tri.amr.getunit(mesh)
```

The basic tools in the unit are the methods `.refine()` and `.coarsen()`.

Both methods take triangle numbers as input arguments.

The AMR unit provides two basic options for selecting triangles:

- Finders, implemented as a set of `.find_` prefixed methods.
- Refinement and coarsing fronts.

### Fronts

Refinement and coarsening fronts are sets of triangles associated with hanging nodes.

Consider a nonconforming element interface:

```
    G
    ●
   / \
  /   \
 B--C--A
 \ / \ /
E ●---● F
   \ /
    ●
    D
```

- The refinement front consists of triangles like GBA.
- The coarsening front constits of super-triangles like ABD (a).

(a) Primary data is the number of the core triangle CEF, as passed to the `.coarse()` method.

For more details, refer to [triellipt.amr.TriFront](triellipt.amr.md#trifront).

### Data

It is possible to interpolate data during refinement and coarsening steps. 
The AMR unit contains a dictionary `.data` as an attribute, whose entries are automatically interpolated when the unit is refined or coarsened.