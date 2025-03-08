<!--
{
  "webtitle": "FEM Solver — triellipt documentation",
  "doctitle": "triellipt — FEM Solver",
  "codeblocks": true
}
-->

## How to Solve PDEs

This page guides you through the process of solving PDEs with **triellipt**.

In what follows, we assume the following import with an alias is used:

```python
import triellipt as tri
```

See [Solving a Problem](solver.md#solving-a-problem) for a complete example.

***

### Creating a Mesh

There are a couple of methods to create a mesh:

- Read a Gmsh mesh using [triellipt.mshread](triellipt.mshread.md).
- Create a structured mesh using [triellipt.mesher](triellipt.mesher.md).

***

### Creating a FEM Unit

Once the mesh is ready, the next step is to create a FEM computing unit:

```python
unit = tri.fem.getunit(mesh, anchors=(0,))
```

The mesh is preprocessed when creating a unit:

- Mesh nodes are aligned using [triellipt.trimesh.TriMesh.alignnodes()](triellipt.trimesh.md#alignnodes).
- The void pivots are placed at the end of the node numbering.

Facts to know:

- The actual unit mesh is available as `unit.mesh`.
- The permutation between the input mesh and the actual mesh is available as `unit.perm`.
- The original mesh is still available as `unit.perm.mesh`. 

For more details, refer to [triellipt.fem.FEMUnit](triellipt.fem.md#femunit).

***

### Creating a Partition

After creating the FEM unit, you must partition the domain to define the *boundary conditions (BCs)*.

Two steps are needed:

- Create a dictionary with the partition specification, e.g. — `parttspec`.
- Run `unit.add_partition(parttspec)` to add the partition to the FEM unit.

Facts to know:

- The map of current unit partitions is available as `unit.partts`.
- Each unit has the default edge-core partition — `unit.base`.

Here is the basic structure of the partition spec:

```python
parttspec = {
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

### Creating a Matrix

FEM matrices are generated from a partition of the FEM unit — see [triellipt.fem.FEMPartt.new_matrix()](triellipt.fem.md#new_matrix).

Two steps are needed:

- Define a FEM operator as a linear combination of basic operators.
- Decide if constraints should be included in the matrix (for non-conformal meshes only).

The matrix is then generated as follows:

```python
matrix = unit.partts['new-domain'].new_matrix(operator, add_constr=True/False)
```

#### Set an operator

A general FEM operator is a linear combination of basic FEM operators.

*Basic operators*

Basic FEM operators are available as properties of the FEM unit:

Name        | Description
------------|----------------------
`massmat`   | Mass-matrix
`massdiag`  | Mass-matrix lumped
`diff_1y`   | 1st y-derivative
`diff_1x`   | 1st x-derivative
`diff_2y`   | 2nd y-derivative
`diff_2x`   | 2nd x-derivative

Facts to know:

- All basic operators are flat arrays representing a special data structure — *matrix data stream*. 
- Understanding of matrix data streams is not necessary for using the package.

Here is an example of constructing a Laplace operator:

```python
operator = unit.diff_2x + unit.diff_2y
```

*Coefficients*

Operators can be multiplied by coefficients defined on the mesh triangles.

Assume we define a coefficient as

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

where we use `unit.ij_t` — indexer that maps triangle-based data to a matrix data stream.

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

***

### Solving a Problem

After combining all the steps, you can proceed to solve the problem.

Here is a complete example of how to reproduce the Bessel function:

```python
"""Example of using triellipt to reproduce the Bessel function.

- The mesh is kept conformal for simplicity.
- The equation solved is:

    d/dy[y∙du/dy] + d/dx[y∙du/dx] + y∙u = 0

"""
import numpy as np
from scipy import sparse as sp
from scipy.special import j0
from matplotlib import pyplot as plt
import triellipt as tri

# Utility functions.

def bessel_func(x, y):
    """Defines the Bessel function.
    """
    return j0(y)

def bessel_operator(unit):
    """Creates the Bessel operator for a FEM unit.
    """

    radius = unit.mesh.centrs_complex.imag[unit.ij_t]

    return radius * (
        unit.diff_2x + unit.diff_2y + unit.massmat
    )

# Creating a mesh.

seed = tri.mesher.trilattice(31, 41, close=True) * 0.05

# Creating a FEM unit.

unit = tri.fem.getunit(
    seed, anchors=((0, 0),)
)

# Partition specification.

partt_box = {
    'name': 'box',
    'anchors': [
        (1.5, 0.1), (1.5, 3.5), (0, 3.4)
    ],
    'dirichlet-sides': [1, 3]
}

# Adding partition to the unit.

unit = unit.add_partition(partt_box)

# Creating a FEM matrix.

mat = unit.partts['box'].new_matrix(
    bessel_operator(unit), add_constr=False
)

# Creating FEM vectors (solution and reference).

sol = unit.partts['box'].new_vector()
ref = unit.partts['box'].new_vector().from_func(bessel_func)

# Solving a problem.

sol[0] = sp.linalg.spsolve(
    mat(0, 0), - mat(0, 1) @ ref[1] - mat(0, 3) @ ref[3]
)

err = np.amax(np.abs(sol[0] - ref[0]))
print(f'L1 error: {err}')

# Plotting the result.

plt.tricontourf(*unit.mesh.triu, sol.body)
plt.triplot(*unit.mesh.triu, '-k', lw=0.2)
plt.axis('equal')

```