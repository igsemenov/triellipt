
# <img src="./docs/configs/logo.png" width="30" height="30"> triellipt

A lightweight 2D finite element solver for elliptic and parabolic PDEs.

- Supports steady-state and transient convection–diffusion problems.  
- Suitable for standard elliptic equations, such as electrostatics and Helmholtz-type problems.
- Designed with *domain decomposition* and *boundary condition transfer* techniques in mind.

Explore [features](#features) and [examples](#examples) below.

## Features

Pure Python implementation using only NumPy and SciPy.

**Meshes:**

- Has an interface to read Gmsh meshes.
- Supports [*conforming*](#conforming-mesh) and [*non-conforming*](#non-conforming-mesh) triangle meshes.
- Provides a flexible framework for mesh [*adaptation*](#adaptive-mesh).
- Includes a suite of [*structured*](#structured-meshes) mesh generators.

**Discretization:**

*Methods*

- Continuous Galerkin finite-element method
- Node-centered control-volume finite-element method

*Features*

- Supports conservative reinterpolation of solution fields across adaptive meshes (see details [here](#mass-conservation-test)).

## Funding

**Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) —
project number [515939493](https://gepris.dfg.de/gepris/projekt/515939493?language=en)**

## Documentation

- ✅ Check the live documentation [here](https://igsemenov.github.io/triellipt/).
- 📄 Check the local documentation [here](docs/sources/index.md).
- 💾 Offline HTML documentation is available by opening `docs/index.html`.

## Triangle meshes

### Conforming mesh

<img src="./docs/images/conforming-mesh.png" width="300">

### Non-conforming mesh

<img src="./docs/images/non-conforming-mesh.png" width="300">

### Adaptive mesh

<img src="docs/images/circ-amr.gif" width="300"/>

### Structured meshes

<img src="./docs/images/suite_grids.png" width="500">

## Examples

### Zernike modes

<img src="./docs/images/zernike-modes.png" width="700">

### Pin-to-plane field

This is the standard electrostatic problem for the pin-to-plane configuration, as considered, for example, in [Celestin et al., *J. Phys. D: Appl. Phys.*, 42(6), 065203 (2009)](https://doi.org/10.1088/0022-3727/42/6/065203).

<img src="./docs/images/field_pin_to_plane.png" width="400">

### Ionization wave

This example shows the evolution of an ionization wave in nitrogen, based on the test case from [Bessières et al., *J. Phys. D: Appl. Phys.*, 40(21), 6559 (2007)](https://doi.org/10.1088/0022-3727/40/21/016). It illustrates the spatial distribution of the electric field within the wave.

<img src="./docs/images/streamer.gif" width="450">
<img src="./docs/images/streamer_field_1d.png" width="400">

### Mass conservation test

This is an example of a standard mass conservation test used to validate the code.

<img src="./docs/images/amr-meshes.png" width="600">

Two limiting configurations are shown, with the mesh adapted cyclically between them, as demonstrated [here](#adaptive-mesh). A conservative reinterpolation algorithm is used to update the field function defined on the mesh. The basic features are as follows:

- Algorithm is exact for constant and linear functions.
- Total nodal mass is preserved for any mesh-defined function, up to numerical error.

### Convergence tests

This example demonstrates a standard convergence test used to validate the code.

- Uses the method of manufactured solutions for verification.
- Applies red-green refinement to a non-conforming mesh.

We consider the Dirichlet problem for the equation

$$
L[u] = \rho
$$

in the domain $[-0.5, 0.5]^2$, where the operator is given by

$$
L = 1 + \frac{\partial}{\partial x} + \frac{\partial}{\partial y} + \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2},
$$  

and the exact solution is given by  

$$
u(x, y) = \cos(\pi x) \cos(\pi y).
$$

<img src="./docs/images/rect-mesh-view.png" width="270">
<img src="./docs/images/rect-mesh-errs.png" width="270">
