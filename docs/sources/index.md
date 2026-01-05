<!--
{
  "webtitle": "triellipt documentation",
  "doctitle": "triellipt",
  "codeblocks": false
}
-->

### Contents

- [Usage](usage.md)
  - [Meshing](meshing.md)
  - [FEM Solver](femsolver.md)
- [Examples](examples.md)
  - [Zernike modes](zernike.md)
  - [Boundary layer](layer.md)
  - [Pulsating mesh](pulsemesh.md)
- [Modules](triellipt.md)
  - [amr](triellipt.amr.md)
  - [fem](triellipt.fem.md)
  - [geom](triellipt.geom.md)
  - [mesher](triellipt.mesher.md)
  - [trimesh](triellipt.trimesh.md)
  - [mshread](triellipt.mshread.md)

### Overview

A lightweight 2D FEM solver for elliptic PDEs.

Basic Features:

- Pure Python implementation using NumPy and SciPy.
- Provides a simple interface for reading Gmsh meshes.
- Supports both conforming and nonconforming triangular meshes.

Additional Details:

- Implements the control-volume finite-element method (CVFEM) with linear (P1) elements.
- Implements the Galerkin finite-element method (FEM) with linear (P1) elements.
- Provides a simple framework for generating nonconforming triangular meshes.
- Supports mesh refinement and coarsening.
- Supports conservative interpolation during mesh refinement/coarsening.
- Supports construction of discrete Dirichlet-to-Neumann operators.
