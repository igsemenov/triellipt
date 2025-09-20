## Triangle meshes

A triangulation with conforming element contacts only:

<img src="./images/conforming-mesh.png" width="300">

A triangulation with both conforming and non-conforming element contacts:

<img src="./images/non-conforming-mesh.png" width="300">

An adaptively changing triangle mesh:

<img src="images/circ-amr.gif" width="300"/>

Structured triangle meshes:

<img src="./images/suite_grids.png" width="450">

### Zernike modes

<img src="./images/zernike-modes.png" width="600">

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

<img src="./images/rect-mesh-view.png" width="270">
<img src="./images/rect-mesh-errs.png" width="270">
