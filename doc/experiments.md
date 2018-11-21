Reproducing published Experiments
=================================

## Adaptive Source Term Iteration

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

This creates a subdirectory that is named after the current date and
time in $OUTPUTDIR and stores all output there.

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
