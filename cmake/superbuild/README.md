# CMake-based SuiteSparse Superbuild

Building SuiteSparse on Microsoft Windows is fairly difficult, as all its
dependencies must be typically built by hand. The whole process is therefore
rather complicated and time consuming.

The `CMakeLists.txt` in this directory fully automates the build process. The
CMake-based superbuild script will download all the dependencies and the tools
required to build them, compile the libraries, and finally build SuiteSparse.

## Getting Started

Make sure to install a recent version of CMake first.

Create a `build` directory (e.g., in the in the root directory of SuiteSparse)
and change into it:
```bash
mkdir build
cd build
```

Execute CMake by pointing it to `cmake/superbuild` (relative to the root
directory):
```bash
cmake ../cmake/superbuild
```
This will generate a Visual Studio solution, or Makefiles (depending on the
CMake generator used).

Finally, the build process can be started using
```bash
cmake --build .
```
or by opening the generated solution in Visual Studio and by compiling the
project from the IDE.


## How does it work?

The `CMakeLists.txt` in this directory will download the following libraries:

* METIS 5.1
* LAPACK 3.6.1
* Intel TBB 2017 Update 1

On Microsoft Windows, the script will additionally install a local version of
MinGW consisting of the following packages:

* `mingw32-gfortran`
* `mingw32-make`
* `msys-patch`

`mingw32-make` is required for building TBB and LAPACK. `msys-patch` installs a
version of the `patch` GNU utility which allows METIS to be patched for
Microsoft Visual Studio 2015 (or later) support. `mingw32-gfortran` is
required for building BLAS and LAPACK libraries.

## Consuming SuiteSparse in CMake Projects

SuiteSparse will generate CMake config files which allow to easily link the
library. A typical CMake project consuming SuiteSparse will look as follows:

```cmake
cmake_minimum_required (VERSION 3.1)
project (MyProject)

# Note that in case CHOLMOD is not found, SuiteSparse will be considered as not
# found even though SQPR is availble.
find_package (SuiteSparse 4.5.3 COMPONENTS CHOLMOD
  OPTIONAL_COMPONENTS SPQR)

add_executable (MyProject main.cpp)

if (SuiteSparse_FOUND)
  # Note that all SuiteSparse components are placed into the SuiteSparse:: CMake
  # namespace.
  target_link_libraries (MyProject PRIVATE SuiteSparse::CHOLMOD)
endif (SuiteSparse_FOUND)

if (TARGET SuiteSparse::SPQR)
  # SuiteSparse_FOUND is not set for optional components. Instead, it must be
  # checked whether the targets of the corresponding optional components are
  # available.
  target_link_libraries (MyProject PRIVATE SuiteSparse::SPQR)
endif (TARGET SuiteSparse::SPQR)
```

For latest updates and issues refer to
https://github.com/sergiud/SuiteSparse/tree/cmake.
