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
See [doc/experiments.md](doc/experiments.md) for explanations on how to
run the example programs bundled with dune-dpg.
If you are interested in an API documentation of dune-dpg, see the
instructions in section [Generating API documentation with Doxygen]
(INSTALL.md#generating-api-documentation-with-doxygen) in
the [INSTALL.md](INSTALL.md) file.

Required Software and Libraries
===============================

Additional to the Dune libraries, you'll need the following programs and
libraries installed on your system:

  - Programs:
    - a C++17-compatible compiler (e.g. GCC >= 7)
    - cmake >= 3.1
  - Libraries:
    - Boost >= 1.63
      (Boost versions as low as 1.61, which introduced Boost::Hana, might
       work but have not been tested.)
    - Eigen >= 3.2
      (http://eigen.tuxfamily.org)
    - UMFPACK (which is part of Suitesparse, www.suitesparse.com)
    - DUNE core libraries 2.7
      including the staging modules dune-typetree and dune-functions
      (https://dune-project.org/releases/2.7.0/)
    - A grid manager. Our examples use dune-uggrid
      (https://www.dune-project.org/modules/dune-uggrid/)
      but one could also use other managers, e.g., ALBERTA or YASPgrid.
      (https://dune-project.org/groups/grid/)
    - for the plot_solution and the ASTI programs, you will also need
      dune-subgrid (>= 2.7)
      (https://gitlab.dune-project.org/extensions/dune-subgrid)

Instructions on how to build dune-dpg can be found in [INSTALL.md](INSTALL.md).

Descriptions of our example programs and how to use them to reproduce
published results can be found in the file
[doc/experiments.md](doc/experiments.md).

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

If you are using our ASTI implementation for radiative transfer problems,
please also consider citing
[Felix’ thesis](http://dx.doi.org/10.18154/RWTH-2018-230893)
([PDF](https://publications.rwth-aachen.de/record/750850/files/750850.pdf))
and our [Math. Comp. paper](http://dx.doi.org/10.1090/mcom/3505)
([arXiv preprint](https://arxiv.org/abs/1810.07035v2)).
The BibTeX entries for these are:

```bibtex
@PhdThesis{Gruber2018,
  author    = {Felix Gruber},
  title     = {Adaptive Source Term Iteration: A Stable Formulation for Radiative Transfer},
  school    = {RWTH Aachen University},
  year      = {2018},
  doi       = {10.18154/RWTH-2018-230893},
  keywords  = {DPG transport solver; fast application of scattering operator; iteration in function space},
}

@Article{DGM2020,
  author     = {Wolfgang Dahmen and Felix Gruber and Olga Mula},
  title      = {An Adaptive Nested Source Term Iteration for Radiative Transfer Equations},
  journal    = {Mathematics of Computation},
  year       = {2020},
  doi        = {10.1090/mcom/3505},
  keywords   = {DPG transport solver, iteration in function space, fast application of scattering operator, Hilbert–Schmidt decomposition, matrix compression},
}
```

If you are using other Dune modules, you might find the corresponding
papers [here](https://dune-project.org/about/publications/).

License
=======

Licensing information for dune-dpg can be found in the accompanying file
[COPYING](COPYING).
