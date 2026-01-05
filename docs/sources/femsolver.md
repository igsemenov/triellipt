<!--
{
  "webtitle": "FEM Solver — triellipt documentation",
  "doctitle": "triellipt — FEM Solver",
  "codeblocks": true
}
-->

## Using the FEM Solver

This page guides you through the process of solving PDEs with **triellipt**.

In what follows, we assume the following import is used:

```python
import triellipt as tri
```

***

### Creating a Mesh

There are a couple of methods to create a mesh:

- Read a Gmsh mesh using [triellipt.mshread](triellipt.mshread.md).
- Create a structured mesh using [triellipt.mesher](triellipt.mesher.md).

Refer to [Meshing](meshing.md) for more details.

***

### Creating a FEM Unit

The next step is to create a FEM computing unit:

```python
unit = tri.fem.getunit(mesh)
```

For more details, refer to [triellipt.fem.getunit()](triellipt.fem.md#getunit).

Two solver modes are available:

- "fvm" corresponds to the control-volume finite-element method (CVFEM) with linear (P1) elements.
- "fem" corresponds to the Galerkin finite-element method (FEM) with linear (P1) elements.

The mesh is preprocessed when creating a unit:

- Mesh nodes are aligned using [triellipt.trimesh.TriMesh.alignnodes()](triellipt.trimesh.md#alignnodes).
- The void pivots are placed at the end of the node numbering.

Facts to know:

- The actual unit mesh is available as `unit.mesh`.
- The permutation between the input mesh and the actual mesh is available as `unit.perm` (a).

(a) The permutation is the `DataPerm` object with the attributes:

- `mesh` is the parent mesh
- `perm` is the permutation from the parent mesh

For more details, refer to [triellipt.fem.FEMUnit](triellipt.fem.md#femunit).

***

### Creating a Partition

After creating the FEM unit, you must partition the mesh boundary.

Two steps are needed:

- Create a dictionary with the partition specification, e.g. — `partt_spec`.
- Run [unit.add_partition(partt_spec)](triellipt.fem.md#add_partition) to add the partition to the FEM unit.

Facts to know:

- The map of current unit partitions is available as `unit.partts`.
- Each unit has the default edge-core partition — `unit.base`.

Here is the basic structure of the partition spec:

```python
partt_spec = {
    'name': 'new-domain',
    'anchors': [
        (0, 0), (0, 1), ...
    ],
    'dirichlet-sides': [1, 2, ...]
}
```

The keys are:

Key                 | Description
--------------------|------------------------------------------------------
`"name"`            | Name of the partition.
`"anchors"`         | Points that define how the mesh boundary is split.
`"dirichlet-sides"` | Section numbers where Dirichlet BCs are applied.

#### Boundary partition

- The boundary partition is defined by anchor points.
- The mesh node that splits the boundary is the one closest to the anchor point.

#### Section numbering

The numbering convention is as follows:

- The boundary sections are numbered with positive integers: 1, 2, 3, ...
- The section numbered 0 is referred to as the *core section*. 

The *core section* includes all mesh nodes except those on Dirichlet sides.

#### Automatic partition

Boundary partitioning can be automated using [triellipt.fem.FEMDtN](triellipt.fem.md#femdtn).

This option is only available for simply-connected meshes. The anchor points are still used to split the boundary, but they are kept separately as corners in the partition.

The section numerng is as follows:

- Odd numbers (1, 3, ...) represent the anchor points themselves.
- Even numbers (2, 4, ...) represent the edges between anchor points.

See also [triellipt.fem.getdtn()](triellipt.fem.md#getdtn).

### Creating a Matrix

FEM matrices are generated from a partition of the FEM unit — see [triellipt.fem.FEMPartt.new_matrix()](triellipt.fem.md#new_matrix).

Two steps are needed:

- Define a FEM operator as a linear combination of basic operators.
- Decide if constraints for hanging nodes should be included in the matrix (a).

(a) Applies only to non-conformal meshes.

The matrix is then generated as follows:

```python
matrix = unit.partts['new-domain'].new_matrix(operator, add_constr=True/False)
```

#### Set an operator

A general FEM operator is a linear combination of basic FEM operators.

*Basic operators*

Basic FEM operators are available as properties of [triellipt.fem.FEMUnit](triellipt.fem.md#femunit).

All basic operators are flat arrays representing a special data structure — *matrix data stream*.

Here is an example of constructing a Laplace operator:

```python
operator = unit.diff_2x + unit.diff_2y
```

*Coefficients*

Operators can be scaled by coefficients defined on a matrix data stream.
To support this, the FEM unit exposes three properties — `unit.ij_r`, `unit.ij_c` and `unit.ij_t` —
to specify the row, column, and triangle indices for each matrix value.

Assume we define a triangle-based coefficient as

```python
coeff = some_func(*unit.mesh.centrs2d)
```

where 

- `some_func(x, y)` represents a certain 2D field.
- `unit.mesh.centrs2d` provides the centroids of triangles.

Then the scaled Laplace operator can be defined as

```python
operator = coeff[unit.ij_t] * (unit.diff_2x + unit.diff_2y)
```

#### Add constraints

- The resulting FEM matrix may or may not include constraints.
- This is controlled by the second argument of the matrix generator.
- Constraints are applicable only to non-conformal meshes.

#### Get sections

- The resulting FEM matrix is callable with arguments `(row_id, col_id)`.
- The output of the call is the matrix block for the specified partition sections.
- The matrix block is a sparse matrix in CSC format.

For more details, refer to [triellipt.fem.MatrixFEM](triellipt.fem.md#matrixfem).

***

### Creating a Vector

FEM vectors are generated from a partition of the FEM unit — see [triellipt.fem.FEMPartt.new_vector()](triellipt.fem.md#new_vector).

- The vector is defined on the mesh nodes.
- The vector can be indexed using the partition section index.
- The vector section is returned as a flat array.

For more details, refer to [triellipt.fem.VectorFEM](triellipt.fem.md#vectorfem).
