# CMake support for SuiteSparse

![Linux](https://github.com/sergiud/suitesparse/actions/workflows/linux.yml/badge.svg)
![macOS](https://github.com/sergiud/suitesparse/actions/workflows/macos.yml/badge.svg)
![Windows](https://github.com/sergiud/suitesparse/actions/workflows/windows.yml/badge.svg)

This fork provides SuiteSparse with an integrated native CMake support layer.
The source tree follows [SuiteSparse by T. A. Davis and
contributors](https://people.engr.tamu.edu/davis/suitesparse.html). The CMake
layer simplifies building SuiteSparse, including CXSparse, on a variety of
platforms.

## Motivation

This fork provides an integration that follows CMake practices important to
downstream consumers: object libraries instead of source shims, a unified native
package configuration for SuiteSparse libraries, and generated explicit export
headers instead of exporting every symbol. It preserves source-level fixes and
broader Windows support while continuing to track upstream.

The shared libraries preserve upstream SuiteSparse library names, SONAMEs,
versions, and documented public APIs. They are intended to be interchangeable
with upstream-built libraries when their external dependencies match. Internal
symbols remain hidden by design, so this guarantee applies to documented public
APIs, not undocumented implementation symbols.

## Highlights

The CMake layer also provides:

### Build and packaging

* Native, relocatable CMake package configurations for SuiteSparse libraries and
  CXSparse
* SPDX [SBOM](https://en.wikipedia.org/wiki/Software_supply_chain) support with
  CMake 4.3 or newer that consults the original license files, which remain
  authoritative
* CPack support

### Platform and toolchain support

* Practical CUDA builds using upstream sources with consistent host compiler
  selection and CMake target wiring
* Explicit DLL exports on Windows and hidden symbol visibility by default on
  supported platforms, enabling link-time optimization (LTO)
* Support for using MinGW BLAS and LAPACK when compiling SuiteSparse with Visual
  Studio
* GitHub releases with MSVC binaries and full BLAS and LAPACK support
* MSVC support extends upstream's baseline with dedicated C complex handling and
  shared and static build coverage

### Portability and compatibility

* Automatic detection of BLAS symbol naming and compiler support for C complex
  types

## Requirements

* C99-compatible compiler, or Microsoft C compiler with complex math support
* CMake 3.22 or newer
* BLAS and LAPACK for the full SuiteSparse build. CXSparse does not require
  them.
* Optional dependencies:
    - C++98 compiler
    - CUDA compiler and toolkit
    - Fortran compiler
    - METIS
    - oneTBB 2021.4 or newer

## Getting Started

Configure and build with:

```bash
$ cmake -S . -B build/
$ cmake --build build/
$ cmake --install build/
```

To run the SuiteSparse tests with CTest, configure with `BUILD_TESTING=ON` and
run:

```bash
$ cmake -S . -B build/ -DBUILD_TESTING=ON
$ ctest --test-dir build/ --output-on-failure
```

You can use SuiteSparse directly from the build directory or after installing
the project:

```cmake
find_package (SuiteSparse 7.12 NO_MODULE)

add_executable (myexe main.cpp)
target_link_libraries (myexe PRIVATE SuiteSparse::CHOLMOD)
```

## Historical Background

The repository was created in 2015 to track the original SuiteSparse releases.
[SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) became a
GitHub project in late 2019. Since July 2026, this repository is tracking
upstream development by regularly merging the upstream `stable` branch into its
`master` branch. It also includes a native CMake support layer to simplify
builds across major platforms.

The CMake layer in this repository extends upstream's baseline Windows support.
It provides MSVC-compatible handling for C complex types, explicit symbol
exports for shared libraries, and demo variants that are not part of upstream's
default CMake targets. MinGW uses GCC with native C99 complex arithmetic, while
MSVC's C compiler requires dedicated complex-type support. The Windows build
also provides full BLAS and LAPACK support for the published MSVC binaries.

The
[suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows)
project was also available, but its CMake support had several limitations. In
particular, it did not provide relocatable CMake package configurations and
relied on Python to preprocess source files. This was still the case in August
2021.

The SuiteSparse project maintained by Dr. Timothy A. Davis followed a different
path. In June 2020, generated type variants required compiling source files
multiple times, and refactoring was needed to avoid this. By December 2022, the
project had a CMake-based build system mostly in place.
The SuiteSparse project's `README.md` is retained under the name
[`README.md.suitesparse`](./README.md.suitesparse). For current upstream
documentation, see the [upstream
README](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/README.md).

## Licensing

The CMake support layer is available under the [Apache License
2.0](https://www.apache.org/licenses/LICENSE-2.0). Modifications to the
SuiteSparse code base are available under the [same
conditions](./LICENSE.suitesparse) as the original code.

The original CMake additions could not be incorporated into official SuiteSparse
releases under the required intellectual property terms. The [upstream CMake
discussion](https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/41)
records this history and the subsequent upstream implementation.
