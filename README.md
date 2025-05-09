
# <img src="./docs/configs/logo.png" width="30" height="30"> triellipt

An adaptive finite-element solver for elliptic PDEs.

## Features

**Meshes:**
- Has an interface to read GMSH meshes  
- Supports [*conforming*](#conforming-mesh) and [*non-conforming*](#non-conforming-mesh) triangle meshes
- Provides a flexible framework for mesh adaptation

**Discretization:**
- Operates on linear elements  
- Supports the following methods:  
  - Continuous Galerkin method  
  - Finite volume element method  
- Ensures mass conservation on adaptive meshes

## Funding

**Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) —
project number [515939493](https://gepris.dfg.de/gepris/projekt/515939493?language=en)**

## Status

- Undergoing usability testing.
- Work on the user guide is in progress.

## Documentation

- ✅ Check the live documentation [here]().
- 📄 Check the local source documentation [here](docs/sources/index.md).
- 💾 Offline HTML documentation is available by opening `docs/index.html`.

## Triangle meshes

### Conforming mesh

<img src="./docs/images/conforming-mesh.png" width="400">

### Non-conforming mesh

<img src="./docs/images/non-conforming-mesh.png" width="400">
