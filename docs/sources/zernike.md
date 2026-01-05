<!--
{
  "webtitle": "Examples — Zernike modes — triellipt documentation",
  "doctitle": "triellipt — Examples: Zernike modes",
  "codeblocks": true
}
-->

## Zernike modes

This example reproduces Zernike polynomials.

The following equation is solved in the unit circle:

```
∇ ∙ (1 - r²) ∇u = (m² - n(n+2))u,
```

where 

- `r` is the absolute value of the position vector 
- `m` is the azimuthal index of the Zernike polynomial
- `n` is the radial index of the Zernike polynomial

The problem is solved subject to Dirichlet boundary conditions.

In what follows, we assume the following import is used:

```python
import numpy as np
from scipy import sparse as sp
import triellipt as tri
```

### Creating a mesh

The mesh is created using Gmsh with the following `.geo` file:

```
R_0 = 1.00;  // Radius
S_0 = 0.02;  // Mesh size

Point(0) = {0, 0, 0};

Point(1) = {+ R_0, + 0.0, 0, S_0};
Point(2) = {+ 0.0, + R_0, 0, S_0};
Point(3) = {- R_0, + 0.0, 0, S_0};
Point(4) = {- 0.0, - R_0, 0, S_0};

Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Mesh.Algorithm = 6;
```

The mesh is saved as `circ-gmsh.msh` in the `mesh` subfolder. The mesh is then read as follows:

```python
mesh = tri.mshread.getreader('mesh').read_mesh('circ-gmsh.msh')
```

The mesh may be coarsened in the central region as follows:

```python
mesh = mesh.reduced(shrink=2)
mesh = mesh.reduced(shrink=3)
```

In this case the resulting mesh is nonconforming.

### Creating a FEM unit.

With the mesh ready, we create a FEM computing unit:

```python
unit = tri.fem.getunit(mesh, mode='fem')
```

- The solver mode can be set to be "fvm" as well.
- Anchor points are not used, since no synchronization of the boundary is needed.
- In this example, the basic boundary partition is used.

In the basic partition we have two sections — core (0) and edge (1).

### Creating FEM data

We first create the distance function:

```python

# We compute the absolute value of the position vector.
# We need triangle based data, so triangle centroids are used.

dist = np.abs(
    unit.mesh.centrs_complex
)

# We now format the distance as a matrix-data-stream.

dist = dist[unit.ij_t]
```

Now the scaled Laplace operator is defined:

```python
lapl = (1.0 - dist ** 2) * (unit.diff_2x + unit.diff_2y)
```

The mass operator is taken as:

```python
mass = unit.massmat
```

We not create the FEM matrix:

```python

M = 3  # Angular index of the Zernike polynomial
N = 5  # Radial index of the Zernike polynomial

amat = unit.base.new_matrix(
  diff + mass * (M * M - N * (N + 2)), add_constr=True
)
```

The solution vector is created as follows:

```python
u = unit.base.new_vector().with_body(1)
```

### Solving the problem

We now find the solution by solving the Dirichlet problem:

```python
u[0] = sp.linalg.spsolve(
    amat(0, 0), - amat(0, 1) @ u[1]
)
```

Here:

- `u[0]` is the solution vector in the core section 0.
- `amat(0, 0)` is the matrix block corresponding to the core section.
- `amat(0, 1) @ u[1]` is the RHS vector due to the Dirichlet BC on the edge section 1.

After solving the problem, the nodal soluton is available as `u.body`.
