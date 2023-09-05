# Compute-Errors

A small library to perform numerical error computation in the OpenFOAM framework

## Description
The library is made by the only `src/error` folder.

Let now define some quantity:
* $u_h$, the piecewise constant function representing the FVM solution
* $u_{h,P^1}$, the piecewise linear FVM solution, using values interpolated on mesh nodes
* $u_{ex}$, the exact continuous solution
* $u_{ex,L^2}$, the piecewise constant function representing the L^2-projection of the exact solution projected on each cell
* $\nabla u_{ex,L^2}$, the piecewise constant function representing the L^2-projection of the gradient of the exact solution projected on each cell

There are many types of error that can be computed. Let $K$ be a mesh cell, $\nabla_h$ some numerical scheme for the gradient, $f$ a mesh face
and $K_f^\pm$ be the owner/neighbour cells relative to face $f$, $d_f$ the distance between their centres,
1. `cellP0L2`: $\sqrt{\sum_K\int_K (u_h - u_{ex})^2}$
2. `projL2`: $\sqrt{\sum_K\int_K (u_h - u_{ex,L^2})^2}$
3. `projH1`: $\sqrt{\sum_K\int_K (\nabla_h u_h - \nabla u_{ex,L^2})^2}$
4. `projDiscreteH1`: $\sqrt{\sum_f\int_f \left(\cfrac{u_h^+ - u_h^- - u_{ex,L^2}^+ + u_{ex,L^2}^-}{d}\right)}$
5. `linearL2`: $\sqrt{\sum_K\int_K (u_{h,P^1} - u_{ex})^2}$
6. `linearH1`: $\sqrt{\sum_K\int_K (\nabla u_{h,P^1} - \nabla u_{ex})^2}$

There are two functions more: `dualL2` and `dualH1`, respectively the $L^2$ and $H^1$ errors computed on a dual grid.
For this refer to [this repository and the relative article](https://github.com/alfiogn/voroToFoam).


## Installation
To use the library, always type on a terminal
```
source etc/bashrc
```

Then, to compile it type
```
./Allwmake
```

Now the library is ready to be used in simulations.
Look at the tutorial to see the usage. The library is used in two ways.
The first is for the source term, where the quadrature functions of the library are use to compute
the source term given by the manufactured solution.
The second is for the error computation; it is contained in a coded function-object, inside `system/errors`.

The command `./Allrun`, runs the case and computes the error and writes it to a csv file.


## TODO

* Tutorial using masks, i.e. computation on only a part of cells
* Computation on domains cut by a STL surface
* Finalise parallel error computation

