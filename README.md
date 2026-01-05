
# <img src="./docs/configs/logo.png" width="30" height="30"> triellipt

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

## Documentation

- ✅ Check the live documentation [here](https://igsemenov.github.io/triellipt/).
- 📄 Check the local documentation [here](docs/sources/index.md).
- 💾 Offline HTML documentation is available by opening `docs/index.html`.

## Installation

The package is small enough to be used without pre-installation:

- Download the package and add its location to `sys.path`.
- After that import the package as usual — <code><b>import</b> triellipt</code>.

## Funding

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) —
project number [515939493](https://gepris.dfg.de/gepris/projekt/515939493?language=en)
