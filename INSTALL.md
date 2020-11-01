Installation of 3rd party libraries and programs
================================================

On Debian GNU/Linux
-------------------

All necessary libraries are included in Debian Buster and Stretch.
They can be installed via

    apt install git cmake libboost-dev libsuitesparse-dev

Using Homebrew on MacOS
-----------------------

CMake installation (it requires the previous installation of autoconf,
                    automake and libtool):

    brew install autoconf
    brew install automake
    brew install libtool
    brew install cmake

The installation of pkg-config is recommended:

    brew install pkg-config

We also need Boost 1.63, or newer, and the Eigen linear algebra library:

    brew install boost
    brew install eigen

UMFPACK installation:

    brew tap hombrew/science
    brew install hombrew/science/suite-sparse

3rd party libraries in non-standard paths
-----------------------------------------

If you have 3rd party libraries in non-standard search paths you have
to tell CMake about it by setting special variables. This can be most
easily done in a dune.opts file so that `dunecontrol` will use the
variables set there when called with the `--opts=$DUNEDIR/dune.opts`
option. To this end create a file dune.opts in `$DUNEDIR` and add
whatever is needed of the following lines:

* To search Boost in `$BOOSTDIR`
  ```
  CMAKE_FLAGS+=" -DBOOST_INCLUDEDIR=$BOOSTDIR"
  ```
* To search Eigen in `$EIGENDIR`
  ```
  CMAKE_FLAGS+=" -DEIGEN3_INCLUDE_DIR=$EIGENDIR"
  ```
* To search SuiteSparse in `$SUITESPARSEDIR`
  ```
  CMAKE_FLAGS+=" -DSuiteSparse_ROOT=$SUITESPARSEDIR"
  ```

Installation and compilation of dune and dune-dpg
=================================================

We give detailed installation guidelines using GCC and the UG grid manager.
The URLs that are given below were still active on December 2016. We apologize
for possible future inconsistencies in the links and hope that you
will nevertheless find your way.

1) Create a directory to harbor the source code of the Dune modules.
   We will call this directory `$DUNEDIR`.

   Remark: Make sure that there are no whitespaces in the `$DUNEDIR` path
           as our build system might not be able to cope with them.
           (Make is known to handle spaces in filenames quite badly and
            other CMake backends might have similar problems.)

Installation of Dune on Debian GNU/Linux
----------------------------------------

2) If you are using Debian bullseye or newer, the Dune 2.7.1 core libraries
   can be installed from the Debian repositories with

    apt install libdune-common-dev libdune-geometry-dev libdune-grid-dev \
                libdune-istl-dev libdune-localfunctions-dev \
                libdune-functions-dev libdune-typetree-dev \
                libdune-uggrid-dev

   Alternatively, you can download the sources of the Dune core modules
   into `$DUNEDIR` like it is explained below in the MacOS instructions.

   If you are using a Debian release which does not include Boost Hana
   (i.e., Boost < 1.61), you can download a recent version from the Boost
   website extract it to some directory `$BOOSTDIR` and add the following
   line to your dune.opts file

    CMAKE_FLAGS+=" -DBOOST_INCLUDEDIR=$BOOSTDIR"

Installation of Dune on other Linux distributions and MacOS
-----------------------------------------------------------

2) Download and extract the following dune sources to `$DUNEDIR`:

  - Version >= 2.7.1 of the following Dune core modules:
    dune-common, dune-geometry, dune-grid, dune-istl,
    dune-localfunctions
    Link: https://dune-project.org/releases/2.7.1/
  - Version >= 2.7 of the following Dune staging modules:
    dune-uggrid (https://www.dune-project.org/modules/dune-uggrid/),
    dune-typetree (https://www.dune-project.org/modules/dune-typetree/),
    dune-functions (https://www.dune-project.org/modules/dune-functions/)
  - The current development version 2.7 from the Git master branch
    of dune-subgrid:
    https://git.imp.fu-berlin.de/agnumpde/dune-subgrid

Preparing dune-dpg
------------------

3) If the dune-dpg directory is not in `$DUNEDIR` move it there now.

Building and compiling dune and dune-dpg
----------------------------------------

4) Create a file dune.opts in `$DUNEDIR` and specify
the C and C++ compiler in it. In our case, we have used gcc and g++ so
the file reads
```
CMAKE_FLAGS="-DCMAKE_C_COMPILER='gcc-9' -DCMAKE_CXX_COMPILER='g++-9'"
CMAKE_FLAGS+=" -DCMAKE_CXX_FLAGS='-pedantic -Wall'"
```

Remarks:
 - On MacOS, the C compiler is set to clang by default so the previous
   configuration option cannot be omitted.
   Additionally g++ is a symlink to clang on MacOS, so one has to use
   the versioned compiler name g++-9. The same goes for gcc.
 - The use of clang and clang++ is in principle also possible here although
   this compiler has not been tested.
 - The second line adds some compile flags to give more warnings.
 - Currently, dune-dpg does not yet support parallel computing.
   So, if you have installed MPI, you can disable it to avoid errors by adding
   `-DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE`
   in the `CMAKE_FLAGS` of your dune.opts file

5) In `$DUNEDIR`, run

     $DUNEDIR/dune-common/bin/dunecontrol --opts=$DUNEDIR/dune.opts all

The command dunecontrol creates makefiles and compiles the sources. It
creates a directory build_cmake inside each Dune module where the
executables are built. As a consequence, the user should cd into
`$DUNEDIR/dune-dpg/build_cmake` to run the examples of the paper.

Remark: If there were any problems with the build system it might be
        helpful to remove all build-cmake directories by running
        `rm -rf */build-cmake`
        in `$DUNEDIR`. This removes all build products and old CMake
        configuration files that might cause conflicts.

For more information and options on `dunecontrol`, run

    $DUNEDIR/dune-common/bin/dunecontrol --help

More information on the build system of Dune can be found under
  - https://dune-project.org/doc/installation/
  - https://dune-project.org/buildsystem/

Generating API documentation with Doxygen
-----------------------------------------

We have some (rudimentary) API documentation for our classes and function
that can be generated from special comments in our source code with the
[Doxygen](http://doxygen.org) documentation generator by calling

    make doc

in `$DUNEDIR/dune-dpg/build-cmake`. The generated html files of the
API documentation can then be read by opening

    $DUNEDIR/dune-dpg/build-cmake/doc/doxygen/html/index.html

in your favorite browser.


Troubleshooting
---------------

### MacOS

On MacOS it might be possible that the following linker error appears:
```
Undefined symbols for architecture x86_64:
"__gfortran_st_write", referenced from:
_xerbla_ in libf77blas.a(xerbla.o)
"__gfortran_st_write_done", referenced from:
_xerbla_ in libf77blas.a(xerbla.o)
"__gfortran_stop_string", referenced from:
_xerbla_ in libf77blas.a(xerbla.o)
"__gfortran_transfer_character_write", referenced from:
_xerbla_ in libf77blas.a(xerbla.o)
"__gfortran_transfer_integer_write", referenced from:
_xerbla_ in libf77blas.a(xerbla.o)
```
In that case, adding the lines
```
target_link_libraries("plot_solution" "gfortran")
target_link_libraries(${convergence_test} "gfortran")
```
to `$DUNEDIR/dune-dpg/src/CMakeLists.txt` should help to correctly
link the example programs.

### disabling MPI

If MPI is installed, this might lead to the following error when running
your dpg-program (i.e. `plot_solution`):
```
"Dune reported error: ParallelError
[CollectiveCommunication:$DUNEDIR/dune-common/dune/common/parallel/mpicollectivecommunication.hh:156]:
You must call MPIHelper::instance(argc,argv) in your main() function
before using the MPI CollectiveCommunication!"
```
This can be fixed by disabling MPI by adding
`-DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE` in the `CMAKE_FLAGS` of your
dune.opts file.

### Definition of UMFPack class not found

If the class `UMFPack` is not found, make sure that you have SuiteSparse
installed and that it has been found by CMake.
SuiteSparse cannot be found by CMake if BLAS is missing. You also have to
make sure to have a Fortran compiler available, e.g. gfortran. Otherwise,
CMake is not able to find BLAS.
