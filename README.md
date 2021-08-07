# CMake support for SuiteSparse

![Linux](https://github.com/sergiud/suitesparse/actions/workflows/linux.yml/badge.svg)
![macOS](https://github.com/sergiud/suitesparse/actions/workflows/macos.yml/badge.svg)
![Windows](https://github.com/sergiud/suitesparse/actions/workflows/windows.yml/badge.svg)

This repository provides two things:

1. Copies of the official versions of [SuiteSparse by T. A. Davis et
   al.](https://people.engr.tamu.edu/davis/suitesparse.html) Please use the
   [`master`](../../tree/master) branch.
2. Specifically, the [`cmake`](../../tree/cmake) branch of the repository also
   implements native CMake support which allows to easily compile SuiteSparse
   (including CXSparse) on a variety of platforms.


The CMake support layer is provided under the [Apache License
2.0](http://www.apache.org/licenses/LICENSE-2.0). Modifications to the
SuiteSparse code base are made available under the [same
conditions](./LICENSE.txt) as the original code.

The original SuiteSparse `README` can be found [here](./README.md.suitesparse).


## Highlights

Besides full CMake support, this branch provides the following additions:

* DLL export on Windows and hidden symbols by default (`-fvisibility=hidden`)
  which enables link time optimization (LTO)
* MinGW BLAS/LAPACK can be used to compile SuiteSparse using Visual Studio
* CPack support

## Usage

First, compile using

```bash
$ cmake . -B build/
$ cmake --build build/
```

Then, one can consume SuiteSparse either directly from the build directory or
after installing the project as follows:
```cmake
find_package (SuiteSparse 5.10 NO_MODULE)

add_executable (myexe main.cpp)
target_link_libraries (myexe PRIVATE SuiteSparse::CHOLMOD)
```

## Background

The repository was created in 2015 to keep track of original releases  before
[SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) became a
Github project at the end of 2019. At the same time, the `cmake` branch
introduced modifications to the original code base in order to enable native
CMake support across major platforms.


While
[suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows)
was already available at the time and confusingly worked not only on Windows as
the name might suggest, its CMake support did have several limitations. In
particular, the implementation did not provide relocatable CMake package
configuration and awkwardly relied on Python for preprocessing source files (as
of August 2021, it still does.)

For IP (and legal) reasons, the provided CMake additions cannot become part of
official SuiteSparse releases. For more information, please refer to [this
post](https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/41#issuecomment-892255262).
