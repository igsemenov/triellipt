<!--
{
  "webtitle": "Examples — Pulsating mesh — triellipt documentation",
  "doctitle": "triellipt — Examples: Pulsating mesh",
  "codeblocks": true
}
-->

## Pulsating mesh

This example demonstrates the use of refinement and coarsening procedures in **triellipt**.

In what follows, we assume the following import is used:

```python
import triellipt as tri
```

### Creating a mesh

The background mesh is created using Gmsh with the following `.geo` file:

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

The mesh is coarsened in the central region as follows:

```python
mesh = mesh.reduced(shrink=2, detach=True)
mesh = mesh.reduced(shrink=3)
mesh = mesh.delmouths()
```

### Creating an AMR unit

The next step is to create an AMR unit:

```python
unit = tri.amr.getunit(mesh)
```

The unit can be provided with a data field.

For this we first create the node-based data as follows:

```python
data = unit.from_func(func)
```

where `func` is a callable that defines the data field.

The data field can be constrained on hanging nodes:

```python
data = unit.constrain(data)
```

Finally, we assign the data field to the unit:

```python
unit.data[0] = data
```

### Refinement step

The unit can be refined in a frontal manner:

```python
for _ in range(COUNT):
    unit = unit.refine(
       unit.front_coarse().trinums
    )
```

where `COUNT` is the number of refinement steps.

Here `unit.front_coarse()` returns a `TriFront` object that contains the numbers of triangles to refine.

The data field assigned to the unit is refined automatically during refinement.

### Coarsening step

Coarsening can be performed after refinement as follows:

```python
for _ in range(COUNT):
    unit = unit.coarsen(
        unit.front_fine().trinums
    )
```

Here `unit.front_fine()` returns a `TriFront` object that contains the coarsening front.
In this case, `.trinums` provide central triangles of super-triangles to coarsen.

The data field assigned to the unit is automatically interpolated during coarsening.