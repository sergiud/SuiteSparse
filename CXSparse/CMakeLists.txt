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

cmake_minimum_required (VERSION 2.8.11)

if (POLICY CMP0054)
  cmake_policy (SET CMP0054 NEW)
endif (POLICY CMP0054)

project (CXSparse C)

set (MAJOR_VERSION 3)
set (MINOR_VERSION 1)
set (VERSION_PATCH 3)

set (SHORT_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

if (NOT VERSION_PATCH EQUAL 0)
  set (SHORT_VERSION "${SHORT_VERSION}.${VERSION_PATCH}")
endif (NOT VERSION_PATCH EQUAL 0)

set (VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${VERSION_PATCH}")

set (CPACK_PACKAGE_NAME CXSparse)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "CXSparse: an extended version of CSparse")
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_PACKAGE_VERSION ${SHORT_VERSION})

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (NOT CPack_CMake_INCLUDED)
  include (CPack)
endif ()

include (CMakePackageConfigHelpers)
include (CheckSymbolExists)
include (CheckIncludeFile)
include (CheckCSourceCompiles)
include (CheckLibraryExists)
include (CheckTypeSize)
include (GenerateSourceFiles)

check_include_file (complex.h HAVE_COMPLEX_H)

if (HAVE_COMPLEX_H)
  check_c_source_compiles ("
#include <complex.h>

int main()
{
  int I;
  I = 0;
  (void)I;
  return 0;
}
" COMPLEX_I_MACRO_NO_CONFLICT)
endif (HAVE_COMPLEX_H)

set (CMAKE_EXTRA_INCLUDE_FILES complex.h)
check_symbol_exists (_DCOMPLEX_ complex.h HAVE__DCOMPLEX_)
check_type_size (_C_double_complex HAVE_C_DOUBLE_COMPLEX BUILTIN_TYPES_ONLY)
check_library_exists (m sqrt "" HAVE_M)

set (HDRS
  Include/cs.h
)

set (PROJECT_SRCS
  Source/cs_add.c
  Source/cs_amd.c
  Source/cs_chol.c
  Source/cs_cholsol.c
  Source/cs_compress.c
  Source/cs_counts.c
  Source/cs_cumsum.c
  Source/cs_dfs.c
  Source/cs_dmperm.c
  Source/cs_droptol.c
  Source/cs_dropzeros.c
  Source/cs_dupl.c
  Source/cs_entry.c
  Source/cs_ereach.c
  Source/cs_etree.c
  Source/cs_fkeep.c
  Source/cs_gaxpy.c
  Source/cs_happly.c
  Source/cs_house.c
  Source/cs_ipvec.c
  Source/cs_leaf.c
  Source/cs_load.c
  Source/cs_lsolve.c
  Source/cs_ltsolve.c
  Source/cs_lu.c
  Source/cs_lusol.c
  Source/cs_malloc.c
  Source/cs_maxtrans.c
  Source/cs_multiply.c
  Source/cs_norm.c
  Source/cs_permute.c
  Source/cs_pinv.c
  Source/cs_post.c
  Source/cs_print.c
  Source/cs_pvec.c
  Source/cs_qr.c
  Source/cs_qrsol.c
  Source/cs_randperm.c
  Source/cs_reach.c
  Source/cs_scatter.c
  Source/cs_scc.c
  Source/cs_schol.c
  Source/cs_spsolve.c
  Source/cs_sqr.c
  Source/cs_symperm.c
  Source/cs_tdfs.c
  Source/cs_transpose.c
  Source/cs_updown.c
  Source/cs_usolve.c
  Source/cs_util.c
  Source/cs_utsolve.c
)

generate_source_files (GEN_SRCS_DINT FILES "${PROJECT_SRCS}" SUFFIX _di)
generate_source_files (GEN_SRCS_DLONG FILES "${PROJECT_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS CS_LONG)
generate_source_files (GEN_SRCS_DINT_COMPLEX FILES "${PROJECT_SRCS}" SUFFIX _ci
  COMPILE_DEFINITIONS CS_COMPLEX)
generate_source_files (GEN_SRCS_DLONG_COMPLEX FILES "${PROJECT_SRCS}" SUFFIX _cl
  COMPILE_DEFINITIONS CS_LONG CS_COMPLEX)

set (SRCS
)

if (HAVE_COMPLEX_H)
  list (APPEND SRCS Source/cs_convert.c)
endif (HAVE_COMPLEX_H)

add_library (CXSparse
  ${HDRS}
  ${SRCS}
  ${GEN_SRCS_DINT}
  ${GEN_SRCS_DLONG}
  ${GEN_SRCS_DINT_COMPLEX}
  ${GEN_SRCS_DLONG_COMPLEX}
)

target_include_directories (CXSparse PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>"
  "$<INSTALL_INTERFACE:include>"
)

target_link_libraries (CXSparse PUBLIC suitesparseconfig)

if (HAVE_M)
  target_link_libraries (CXSparse PUBLIC m)
endif (HAVE_M)

if (NOT HAVE_COMPLEX_H)
  target_compile_definitions (CXSparse PUBLIC NCOMPLEX)
endif (NOT HAVE_COMPLEX_H)

if (NOT HAVE_C_DOUBLE_COMPLEX STREQUAL "")
  target_compile_definitions (CXSparse PUBLIC HAVE_C_DOUBLE_COMPLEX)
endif (NOT HAVE_C_DOUBLE_COMPLEX STREQUAL "")

if (NOT COMPLEX_I_MACRO_NO_CONFLICT)
  target_compile_definitions (CXSparse PUBLIC CXSPARSE_COMPLEX_I_CONFLICT)
endif (NOT COMPLEX_I_MACRO_NO_CONFLICT)

if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  target_compile_definitions (CXSparse PRIVATE _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")

set_target_properties (CXSparse PROPERTIES PUBLIC_HEADER Include/cs.h)

install (TARGETS CXSparse
  EXPORT cxsparse-targets
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

configure_package_config_file (cxsparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config.cmake
  INSTALL_DESTINATION lib/cmake/CXSparse
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file (cxsparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY SameMajorVersion)

export (TARGETS CXSparse NAMESPACE CXSparse:: FILE cxsparse-targets.cmake)
export (PACKAGE CXSparse)

install (EXPORT cxsparse-targets NAMESPACE CXSparse::
  DESTINATION lib/cmake/CXSparse)
install (FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config-version.cmake
  DESTINATION lib/cmake/CXSparse)

set_target_properties (CXSparse PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties (CXSparse PROPERTIES VERSION ${MAJOR_VERSION})
set_target_properties (CXSparse PROPERTIES SOVERSION ${VERSION})

add_executable (cs_ci_demo1
  Demo/cs_ci_demo1.c
)
target_link_libraries (cs_ci_demo1 PRIVATE CXSparse)

add_executable (cs_ci_demo2
  Demo/cs_ci_demo.c
  Demo/cs_ci_demo.h
  Demo/cs_ci_demo2.c
)
target_link_libraries (cs_ci_demo2 PRIVATE CXSparse)

add_executable (cs_ci_demo3
  Demo/cs_ci_demo.c
  Demo/cs_ci_demo.h
  Demo/cs_ci_demo3.c
)
target_link_libraries (cs_ci_demo3 PRIVATE CXSparse)

add_executable (cs_cl_demo1
  Demo/cs_cl_demo1.c
)
target_link_libraries (cs_cl_demo1 PRIVATE CXSparse)

add_executable (cs_cl_demo2
  Demo/cs_cl_demo.c
  Demo/cs_cl_demo.h
  Demo/cs_cl_demo2.c
)
target_link_libraries (cs_cl_demo2 PRIVATE CXSparse)

add_executable (cs_cl_demo3
  Demo/cs_cl_demo.c
  Demo/cs_cl_demo.h
  Demo/cs_cl_demo3.c
)
target_link_libraries (cs_cl_demo3 PRIVATE CXSparse)

add_executable (cs_di_demo1
  Demo/cs_di_demo1.c
)
target_link_libraries (cs_di_demo1 PRIVATE CXSparse)

add_executable (cs_di_demo2
  Demo/cs_di_demo.c
  Demo/cs_di_demo.h
  Demo/cs_di_demo2.c
)
target_link_libraries (cs_di_demo2 PRIVATE CXSparse)

add_executable (cs_di_demo3
  Demo/cs_di_demo.c
  Demo/cs_di_demo.h
  Demo/cs_di_demo3.c
)
target_link_libraries (cs_di_demo3 PRIVATE CXSparse)

add_executable (cs_dl_demo1
  Demo/cs_dl_demo1.c
)
target_link_libraries (cs_dl_demo1 PRIVATE CXSparse)

add_executable (cs_dl_demo2
  Demo/cs_dl_demo.c
  Demo/cs_dl_demo.h
  Demo/cs_dl_demo2.c
)
target_link_libraries (cs_dl_demo2 PRIVATE CXSparse)

add_executable (cs_dl_demo3
  Demo/cs_dl_demo.c
  Demo/cs_dl_demo.h
  Demo/cs_dl_demo3.c
)
target_link_libraries (cs_dl_demo3 PRIVATE CXSparse)

add_executable (cs_demo1
  Demo/cs_demo1.c
)
target_link_libraries (cs_demo1 PRIVATE CXSparse)

add_executable (cs_demo2
  Demo/cs_demo.c
  Demo/cs_demo.h
  Demo/cs_demo2.c
)
target_link_libraries (cs_demo2 PRIVATE CXSparse)

add_executable (cs_demo3
  Demo/cs_demo.c
  Demo/cs_demo.h
  Demo/cs_demo3.c
)
target_link_libraries (cs_demo3 PRIVATE CXSparse)
