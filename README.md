Description
===========

This folder contains the sources of the dune-dpg library, which allows to
solve Partial Differential Equations with Discontinuous Petrov-Galerkin
finite elements. It is built upon the finite element package Dune. For
further information on the dune-dpg library, we refer to the paper
[The dune-dpg library for solving PDEs with Discontinuous Petrov-Galerkin
finite elements](http://dx.doi.org/10.11588/ans.2017.1.27719)
by F. Gruber, A. Klewinghaus and O. Mula.

For build instruction see the accompanying [INSTALL.md](INSTALL.md) file.
See section [Running dune-dpg](#running-dune-dpg) in this README document for
explainations on how to run basic examples which are described in the paper.
If you are interested in an API documentation of dune-dpg, see the
instructions in section [Generating API documentation with Doxygen]
(INSTALL.md#generating-api-documentation-with-doxygen) in
the [INSTALL.md](INSTALL.md) file.

Required Software and Libraries
===============================

Additional to the Dune libraries, you'll need the following programs and
libraries installed on your system:

  - Programs:
    - a C++14-compatible compiler (e.g. GCC >= 6)
    - cmake >= 2.8.12
  - Libraries:
    - Boost >= 1.63
      (Boost versions as low as 1.61, which introduced Boost::Hana, might
       work but have not been tested.)
    - UMFPACK (which is part of Suitesparse, www.suitesparse.com)
    - DUNE core libraries 2.5
      including the staging modules dune-typetree and dune-functions
      (https://dune-project.org/releases/2.5.0/)
    - A grid manager. Our examples use dune-uggrid
      which has been added as a staging module in DUNE 2.5
      (https://dune-project.org/releases/2.5.0/)
      but one could also use other managers, e.g., ALBERTA or YASPgrid.
      (https://dune-project.org/groups/grid/)
    - for the plot_solution test program, you will also need the
      dune-subgrid module
      (https://git.imp.fu-berlin.de/agnumpde/dune-subgrid)

Instructions on how to build dune-dpg can be found in [INSTALL.md](INSTALL.md).

Running dune-dpg
================

In the following, we assume that we are in the directory
`$DUNEDIR/dune-dpg/build-cmake`.

Dune-dpg comes with three example programs:

  - src/plot_solution.cc
  - src/convergence_test.cc
  - src/profile_testspacecoefficientmatrix.cc

After compilation (see [INSTALL.md](INSTALL.md) for instructions), the
example programs are found in `$DUNEDIR/dune-dpg/build-cmake/src/`.

Description of src/plot_solution.cc
-----------------------------------

plot_solution.cc computes the solution with Discontinuous Petrov-Galerkin
finite elements of the pure transport problem
$$
  \beta \cdot \phi +c \phi = 1 in [0,1]x[0,1]
                      \phi = 0 on the boundary
$$
We refer to the paper for notations and also
for indications on how to solve other Partial Differential Equations.

To visualize the solution, run in `$DUNEDIR/dune-dpg/build-cmake/src/`

    ./plot_solution <n> <c> <betaX> <betaY>

where
  `<n>` has to be replaced by the desired grid resolution,
  `<c>` is the linear term in the transport problem
  `<betaX> <betaY>` are the X and Y components of the transport direction
                    beta=(betaX, betaY).
  (When unspecified, c=0 and beta=(cos(π/8), sin(π/8)).)

The program will write two .vtu files to the current working directory,

  - transport_solution.vtu
  - transport_solution_trace.vtu

They contain the numerical solution of the interior variable $\phi$ and
the lifting $w$ of the trace variable (cf. paper).
The files allow to visualize $\phi$ and $w$ in ParaView. If you have the
pvpython interpreter shipped with ParaView, you can also run
`scripts/plot_solution.py` to regenerate the solution plot given in the paper.
(This script was run with ParaView 4.2.0. As the Python interface of ParaView
seems to be highly unstable we cannot guarantee that the script will run
unmodified on another version of ParaView.)

Description of src/convergence_test.cc
--------------------------------------

convergence_test.cc uses several preprocessor variables to set the
refinement level and polynomial degree of the test search space and the
search space used in the computation of the a posteriori error.
These variables get instantiated with different values by CMake to create
executables of the form

    src/convergence_test_ls$LS_ks$KS_la$LA_ka$KA

where `$LS` is the refinement level of the test search space,
      `$KS` is the polynomial degree of the test search space,
      `$LA` is the refinement level of the a posteriori search space,
      `$KA` is the polynomial degree of the a posteriori search space.

They can be run in `$DUNEDIR/dune-dpg/build-cmake/src/` by e.g.

    ./convergence_test_ls0_ks3_la0_ka5 <n>

to compute the numerical example of the paper,
$$
  \beta \cdot \phi = 1 in [0,1]x[0,1]
              \phi = 0 on the boundary
$$
for beta=(cos(pi/8),sin(pi/8)). The additional feature that is given in
this example is the computation of the exact L2 error and its a posteriori
estimation.

Remarks:
  - The a posteriori error indicator includes both the error on the
    interior variable and also the trace variable.
  - It is possible to compute the a posteriori error estimator for
    other PDEs. The current program is just an example.

The paper takes this problem as the support for numerical tests and
analyses the impact of
 - the mesh size H of the trial space
 - the h-refinement level of the test search space
on the error and its a posteriori estimation.

To reproduce the convergence results from the paper, call the script

    ../scripts/run_convergence_test.sh 10 20 40 60 80 100 120 140 \
                                       160 200 250 300

which will call the `convergence_test_*` programs with the given grid
resolutions n=1/H.
The computation takes several hours. (You can leave out some of the larger
grid resolutions n, if you don't want to wait for this long. Especially
the test case with locally refined test search space of level 3 gets very
slow for large n).

Finally, the convergence plots can be generated with

    ../scripts/convergence_plots.py

Remark: To regenerate the plots from the paper, it is advised to compile
the test programs in release mode to significantly speed up the computations.
The release mode can be configured with

    cmake -DCMAKE_BUILD_TYPE=Release .

Afterwards, we can compile the test programs with

    make

Description of src/profile_testspacecoefficientmatrix.cc
--------------------------------------------------------

profile_testspacecoefficientmatrix.cc uses the preprocessor variables
`LEVEL_SEARCH` and `K_SEARCH` to set the refinement level and polynomial
degree of the test search space.
These variables get instantiated with different values by CMake to create
executables of the form

    src/profile_testspacecoefficientmatrix_ls$LS_ks$KS

where `$LS` is the refinement level of the test search space,
      `$KS` is the polynomial degree of the test search space.

They can be run in `$DUNEDIR/dune-dpg/build-cmake/src/` by e.g.

    ./profile_testspacecoefficientmatrix_ls0_ks3 <n>

to measure the time to compute the numerical example of the paper,
$$
  \beta \cdot \phi = 1 in [0,1]x[0,1]
              \phi = 0 on the boundary
$$
for beta=(cos(pi/8),sin(pi/8)).

To reproduce the results from the paper, call the script

    ../scripts/run_profile.sh 10 20 40 60 80 100 120 140 160 200 250 300 \
                              > profile

which will call `profile_testspacecoefficientmatrix_*` for the given
grid resolutions n=1/H.  The computation takes several minutes.

Finally, the profile plots can be generated with

    ../scripts/profile_plots.py -o profile.pdf profile

Remark: To regenerate the plots from the paper, it is advised to compile
the test programs in release mode to significantly speed up the computations.
The release mode can be configured with

    cmake -DCMAKE_BUILD_TYPE=Release .

Afterwards, we can compile the test programs with

    make

Citing dune-dpg
===============

If you use dune-dpg in your own scientific work, please consider
citing our [overview paper](http://dx.doi.org/10.11588/ans.2017.1.27719)
([PDF](http://journals.ub.uni-heidelberg.de/index.php/ans/article/download/27719/29543)).
Here's the BibTeX data for the paper:

```bibtex
@Article{GKM2017,
  author         = {Felix Gruber and Angela Klewinghaus and Olga Mula},
  title          = {The {DUNE-DPG} Library for Solving {PDE}s with {Discontinuous Petrov--Galerkin} Finite Elements},
  journal        = {Archive of Numerical Software},
  year           = {2017},
  volume         = {5},
  number         = {1},
  pages          = {111--128},
  doi            = {10.11588/ans.2017.1.27719},
  keywords       = {discontinuous Petrov-Galerkin; partial differential equations; inf-sup stability; transport equation; finite elements; DUNE},
}
```

If you are using other Dune modules, you might find the corresponding
papers [here](https://dune-project.org/about/publications/).

License
=======

Licensing information for dune-dpg can be found in the accompanying file
[COPYING](COPYING).
