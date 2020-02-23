Reproducing published Experiments
=================================

We list here, ordered by [publication](../README.md#citing-dune-dpg),
the example programs included in dune-dpg. To allow you to reproduce
our published results more easily, we try to list here all the
parameters used to generate the figues seen in those publications as
well as the exact version of dune-dpg that had been used.

## Adaptive Source Term Iteration (Gruber2018, DGM2020)

We describe here how to reproduce the results given in Felix’ dissertation.
First make sure to compile dune-dpg as described in INSTALL.md. The exact
version numbers of dependencies used for the results are given in
Appendix A of the disseration. Especially, dune-dpg 0.4.2 was used, so later
versions of dune-dpg might give slightly different results as I expect
the Adaptive Source Term Iteration to further improve in the future.

Let $BUILDDIR be your build directory for dune-dpg. Then you can see that
there are different `asti_*` programs in $BUILDDIR/src. The ones with
`_normalized` suffix use normalized test spaces to improve the accuracy
of the solutions of the local auxiliary problems to find the near-optimal
test spaces.

Set $OUTPUTDIR to a path where you want the results of the ASTI solver
to be stored. Then run the following commands to reproduce the actual
results:

```
export KAPPA_RATIO=0.2
export FINAL_ACCURACY=0.005
export GAMMA=0.5
export MAX_OUTER_ITERATIONS=11
export GRID_RESOLUTION=7
$BUILDDIR/src/asti_checkerboard_normalized \
    -p -i -o $OUTPUTDIR -k $KAPPA_RATIO \
    $FINAL_ACCURACY $GAMMA $MAX_OUTER_ITERATIONS $GRID_RESOLUTION
```

This creates the directory $OUTPUTDIR if it doesn't exist yet and
stores all output there.

The output consists of:

* A file `output` that contains a detailed log of the ASTI algorithm.
  This can be used to generate the convergence plot from Figure 11
  of the dissertation, using `scripts/asti/asti_plot.py`.
  `scripts/asti/convergence_table.py` can be used to generate a
  convergence table from the output file.
* Several .vtu files for the directional solutions and the
  integrated solutions. Those can be plotted using Paraview or similar
  visualization software.

Make sure that you run the program on a machine with sufficient amounts
of memory. It will need almost 44 GB of RAM and will take several hours
to finish.

### Condition Number of Discrete Trial to Test Space Operator

If $SOURCEDIR is the path to dune-dpg's source directory and $BUILDDIR
is the path to the build directory, you can create a plot of the
condition number of the trial to test space operator in dependence of
the element size h by running

```
$BUILDDIR/src/trial-to-test-condition \
    | $SOURCEDIR/scripts/plot_condition.py
```

This creates a new file `conditions.pdf` that contains the plot.

## GKM2017

For this paper, we used version 0.2.1 of dune-dpg.

In the following descriptions, we assume that you followed the build
instructions found in [INSTALL.md](../INSTALL.md) and that you are
currently in the directory `$DUNEDIR/dune-dpg/build-cmake/src`.

The results in our overview paper have been created with the following
three example programs:

  - src/plot_solution.cc
  - src/convergence_test.cc
  - src/profile_testspacecoefficientmatrix.cc

### Description of src/plot_solution.cc

plot_solution.cc computes the solution with Discontinuous Petrov-Galerkin
finite elements of the pure transport problem
```math
\begin{cases}
  \beta \cdot \varphi +c \varphi = 1 & \text{in $[0,1]\times[0,1]$} \\
                         \varphi = 0 & \text{on the boundary}
\end{cases}
```
We refer to the paper for notations and also
for indications on how to solve other Partial Differential Equations.

To visualize the solution, run in `$DUNEDIR/dune-dpg/build-cmake/src/`

    ./plot_solution <n> <c> <betaX> <betaY>

where
  `<n>` has to be replaced by the desired grid resolution,
  `<c>` is the linear term in the transport problem
  `<betaX> <betaY>` are the $`x`$ and $`y`$ components of the transport
                    direction $`\beta=(\beta_x, \beta_y)`$.
  (When unspecified, $`c=0`$ and $`\beta=(\cos(\pi/8), \sin(\pi/8))`$.)

The program will write two .vtu files to the current working directory,

  - transport_solution.vtu
  - transport_solution_trace.vtu

They contain the numerical solution of the interior variable $`\varphi`$
and the lifting $`w`$ of the trace variable (cf. paper).
The files allow to visualize $`\varphi`$ and $`w`$ in ParaView. If you
have the pvpython interpreter shipped with ParaView, you can also run
`scripts/plot_solution.py` to regenerate the solution plot given in the paper.
(This script was run with ParaView 4.2.0. As the Python interface of ParaView
seems to be highly unstable we cannot guarantee that the script will run
unmodified on another version of ParaView.)

### Description of src/convergence_test.cc

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
```math
\begin{cases}
  \beta \cdot \varphi = 1 & \text{in $[0,1]\times[0,1]$} \\
              \varphi = 0 & \text{on the boundary}
\end{cases}
```
for $`\beta=(\cos(\pi/8),\sin(\pi/8))`$. The additional feature that is
given in this example is the computation of the exact L2 error and its a
posteriori estimation.

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

### Description of src/profile_testspacecoefficientmatrix.cc

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
```math
\begin{cases}
  \beta \cdot \varphi = 1 & \text{in $[0,1]\times[0,1]$} \\
              \varphi = 0 & \text{on the boundary}
\end{cases}
```
for $`\beta=(\cos(\pi/8),\sin(\pi/8))`$.

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

## Example programs unrelated to any papers

### Description of src/manufactured_transport.cc and src/manufactured_transport_uniform.cc

manufactured_transport.cc and manufactured_transport_uniform.cc compute
the Discontinuous Petrov-Galerkin solution of a transport problem whose
solution is known. Thus, we can compare the DPG error estimates with
the exact errors.

manufactured_transport.cc uses adaptive refinement while
manufactured_transport_uniform.cc uses uniform refinement.

You can run each of those programs from the build directory as

    src/manufactured_transport <n>
    src/manufactured_transport_uniform <n>

where `<n>` is the desired grid resolution.
