<!--
{
  "webtitle": "Examples — Boundary layer — triellipt documentation",
  "doctitle": "triellipt — Examples: Boundary layer",
  "codeblocks": true
}
-->

## Boundary layer

This example solves the 1D boundary layer problem.

The following equation is considered:

```
uₓ = ν∙Δu
```

on the interval `[0, 1]` along the x-coordinate.

The boundary conditions are:

```
u(x=0, y) = 0

u(x=1, y) = 1
```

In what follows, we assume the following import is used:

```python
import numpy as np
from scipy import sparse as sp
import triellipt as tri
```

### Creating a mesh

For this case, we create a structured triangle mesh:

```python
mesh = tri.mesher.trigrid(101, 21, 'cross-wise') * 0.01
```

Optionally, we can apply a refinement near `x = 0`:

```python
unit = tri.amr.getunit(mesh)

unit = unit.refine(
  unit.find_masked(lambda x, _: x < 0.5)
)

unit = unit.refine(
  unit.find_masked(lambda x, _: x < 0.2)
)

mesh = unit.mesh.twin()
```

### Creating a DtN unit.

Next, we create the FEM-DtN unit:

```python
unit = tri.fem.getdtn(
    mesh, anchors=[(0, 0), (1, 0), (1, 0.2), (0, 0.2)], mode='fem'
)

unit = unit.switch_side(2)
unit = unit.switch_side(6)
```

The FEM-DtN unit represents a FEM unit with a prescribed boundary partition.

The FEM-DtN unit has two basic properties:

- `.fem` — the FEM unit itself.
- `.dtn` — the DtN partition.

The partition has corners given by the initial anchor points and edges connecting the corner points.
Partition parts with odd indices correspond to corners and parts with even indices correspond to edges.

All partition parts are considered to be of Dirichlet type by default.
The Dirichlet/Neumann type of the parts can be changed using the `switch_side()` method.

The current Dirichlet parts can be obtained by calling the `dirich_sides()` method.

### Creating FEM data.

We now create the FEM operators:

```python
conv = unit.fem.diff_1x
diff = unit.fem.diff_2x + unit.fem.diff_2y
```

where `conv` and `diff` are the convection and diffusion operators, respectively.

Based on this, we create the FEM matrix:

```python
amat = unit.dtn.new_matrix(
    - conv - VISC * diff, add_constr=True
)
```

with `VISC` being the viscosity coefficient.

We also create the following FEM vectors:

```python
u = unit.dtn.new_vector()  # Solution vector 
b = unit.dtn.new_vector()  # RHS vector
```

### Solving the problem

To solve the problem, we first apply the boundary conditions:

```python
# BC at x = 1
for _ in [3, 4, 5]:
  u[_] = 1

# BC at x = 0
for _ in [1, 7, 8]:
  u[_] = 0
```

After that we define the RHS vector:

```python
for _ in [1, 3, 4, 5, 7, 8]:
    b[0] = b[0] - amat(0, _) @ u[_]
```

Finally, we find the solution:

```python
u[0] = sp.linalg.spsolve(amat(0, 0), b[0])
```

The complete solution is available as `u.body`.