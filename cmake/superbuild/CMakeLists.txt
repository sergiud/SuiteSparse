# Copyright (c) 2016 Sergiu Deitsch
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTMETISLAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required (VERSION 3.1)
project (SuiteSparse)

option (BUILD_SHARED_LIBS "Build shared libraries" ON)

include (CheckLanguage)

if (NOT WIN32)
  # On non-Windows platforms, the Fortran compiler is required. The script will
  # not proceed without the Fortran compiler.
  enable_language (Fortran)
else (NOT WIN32)
  check_language (Fortran)
endif (NOT WIN32)

include (ExternalProject)
include (ProcessorCount)

ProcessorCount (_PROCESSORS)

if (NOT _PROCESSORS EQUAL 0)
  set (_MAKEFLAGS -j${_PROCESSORS})
endif (NOT _PROCESSORS EQUAL 0)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/..)

set (_TBB_BUILD_PREFIX external)

if (WIN32)
  set (_MINGW_PACKAGES mingw32-make msys-patch)

  # Use system's fortran compiler if possible
  if (NOT CMAKE_Fortran_COMPILER)
    list (APPEND _MINGW_PACKAGES gfortran)
  endif (NOT CMAKE_Fortran_COMPILER)

  set (_MINGW_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/mingw)
  set (_MINGW_BIN_TMP ${_MINGW_PREFIX}/src/mingw/bin)
  set (_MSYS_BIN_TMP ${_MINGW_PREFIX}/src/mingw/msys/1.0/bin)
  set (PATCH_EXECUTABLE ${_MSYS_BIN_TMP}/patch)

  ExternalProject_Add (mingw
    PREFIX ${_MINGW_PREFIX}
    URL http://vorboss.dl.sourceforge.net/project/mingw/Installer/mingw-get/mingw-get-0.6.2-beta-20131004-1/mingw-get-0.6.2-mingw32-beta-20131004-1-bin.tar.xz
    URL_MD5 6453e5e9a88511a599630013ca0f2871
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND ${_MINGW_PREFIX}/src/mingw/bin/mingw-get install
    ${_MINGW_PACKAGES} #mingw32-binutils
    BUILD_COMMAND ""
  )

  # HACK
  set (_TBB_BUILD_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS=${_MAKEFLAGS} ${_MINGW_BIN_TMP}/mingw32-make
    tbb_build_prefix=${_TBB_BUILD_PREFIX})
  set (_TBB_DEPENDS DEPENDS mingw)
  set (_METIS_DEPENDS DEPENDS mingw)
  set (_LAPACK_DEPENDS DEPENDS mingw)
else (WIN32)
  set (_TBB_BUILD_COMMAND make tbb_build_prefix=${_TBB_BUILD_PREFIX})
endif (WIN32)

if (MSVC)
  # METIS needs to be patch to be usable with MSVC >= 15.0
  # Also, ignore if the patch has been already applied (-N).
  set (_METIS_PATCH_COMMAND ${CMAKE_COMMAND} -E env ${PATCH_EXECUTABLE} -p1
    -N -l -s -i ${CMAKE_CURRENT_SOURCE_DIR}/metis_msvc.diff)
  # Build METIS as shared library as requsted by the user, unless MSVC is used.
  # In this case, the symbols are not exported automatically requiring the
  # library to be build as a static one.
  set (_METIS_SHARED OFF)
else (MSVC)
  set (_METIS_SHARED ${BUILD_SHARED_LIBS})
endif (MSVC)

ExternalProject_Add (metis
  URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
  URL_MD5 5465e67079419a69e0116de24fce58fe
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  CMAKE_CACHE_ARGS
    -DGKLIB_PATH:PATH=${CMAKE_CURRENT_BINARY_DIR}/metis-prefix/src/metis/GKlib
    -DSHARED:BOOL=${_METIS_SHARED}
  INSTALL_COMMAND ""
  PATCH_COMMAND ${_METIS_PATCH_COMMAND}
  ${_METIS_DEPENDS}
)

if (WIN32)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set (_LAPACK_ADDITIONAL_ARGS -DCMAKE_GNUtoMS:BOOL=ON)
  endif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

  set (_LAPACK_BUILD
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env PATH=${_MINGW_BIN_TMP}
      ${CMAKE_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/lapack-prefix/src/lapack
      -DBUILD_SHARED_LIBS:BOOL=ON ${_LAPACK_ADDITIONAL_ARGS}
      -G "MinGW Makefiles"
    BUILD_COMMAND ${CMAKE_COMMAND} -E env MAKEFLAGS=${_MAKEFLAGS}
    PATH=${_MINGW_BIN_TMP} ${CMAKE_COMMAND} --build
      ${CMAKE_CURRENT_BINARY_DIR}/lapack-prefix/src/lapack-build
  )
else (WIN32)
  if (CMAKE_GENERATOR STREQUAL Ninja)
    # Use make instead of ninja to build LAPACK, since ninja does not support
    # gfortran.
    set (_LAPACK_BUILD CMAKE_GENERATOR "Unix Makefiles")
  endif (CMAKE_GENERATOR STREQUAL Ninja)
endif (WIN32)

ExternalProject_Add (lapack
  URL http://www.netlib.org/lapack/lapack-3.6.1.tgz
  URL_MD5 421b2cb72e15f237e144428f9c460ee0
  ${_LAPACK_BUILD}
  #CMAKE_ARGS
	#-DCMAKE_PREFIX_PATH=${_MSYS_BIN_TMP};${_MINGW_BIN_TMP}
  #  -DCMAKE_Fortran_COMPILER=${_MINGW_BIN_TMP}/mingw32-gfortran.exe
  #  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	#-DCMAKE_MAKE_PROGRAM=${_MINGW_BIN_TMP}/mingw32-make.exe
  #CMAKE_GENERATOR "MinGW Makefiles"
  CMAKE_CACHE_ARGS -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  INSTALL_COMMAND ""
  ${_LAPACK_DEPENDS}
)

ExternalProject_Add (tbb
  URL https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb2017_20160916oss_src.tgz
  URL_MD5 8d15135e0bc71a82fc60446cfd78ae70
  CONFIGURE_COMMAND ""
  BUILD_COMMAND "${_TBB_BUILD_COMMAND}"
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  ${_TBB_DEPENDS}
)

set (_PREFIX_PATH
  # HACK: CMAKE_CACHE_ARGS will not work otherwise. The semicolon after the
  # first entry is omitted. So we add an empty (dummy) entry.
  ""
  ${CMAKE_CURRENT_BINARY_DIR}/metis-prefix/src/metis-build/libmetis
  ${CMAKE_CURRENT_BINARY_DIR}/tbb-prefix/src/tbb/build/${_TBB_BUILD_PREFIX}_debug
  ${CMAKE_CURRENT_BINARY_DIR}/tbb-prefix/src/tbb/build/${_TBB_BUILD_PREFIX}_release
)

# Workout the library prefix and extension
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set (_LIB_PREFIX "")
  set (_LIB_SUFFIX ${CMAKE_LINK_LIBRARY_SUFFIX})
else (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if (BUILD_SHARED_LIBS)
    set (_LIB_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
    set (_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
  else (BUILD_SHARED_LIBS)
    set (_LIB_PREFIX ${CMAKE_SHARED_STATIC_PREFIX})
    set (_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif (BUILD_SHARED_LIBS)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

set (_LAPACK_LIBS
  -DBLAS_blas_LIBRARY:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/lapack-prefix/src/lapack-build/lib/libblas${_LIB_SUFFIX}
  -DLAPACK_lapack_LIBRARY:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/lapack-prefix/src/lapack-build/lib/liblapack${_LIB_SUFFIX}
)

set (_METIS_LIBS
  -DMETIS_LIBRARY:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/metis-prefix/src/metis-build/libmetis/${CMAKE_CFG_INTDIR}/${_LIB_PREFIX}metis${_LIB_SUFFIX}
)

ExternalProject_Add (suitesparse
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  CMAKE_CACHE_ARGS
    -DBLA_F95:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DMETIS_INCLUDE_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/metis-prefix/src/metis/include
    -DTBB_INCLUDE_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/tbb-prefix/src/tbb/include
    -DTBB_ROOT_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/tbb-prefix/src/tbb
	-DCMAKE_LIBRARY_PATH:STRING=${_PREFIX_PATH}
	${_METIS_LIBS}
	${_LAPACK_LIBS}
  INSTALL_COMMAND ""
  DEPENDS lapack metis tbb
)
