#
# CMake support layer for SuiteSparse
#
# Copyright 2016-2021 Sergiu Deitsch <sergiu.deitsch@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

cmake_minimum_required (VERSION 3.0)

if (POLICY CMP0054)
  cmake_policy (SET CMP0054 NEW)
endif (POLICY CMP0054)

if (POLICY CMP0063)
  cmake_policy (SET CMP0063 NEW)
endif (POLICY CMP0063)

#if (MSVC)
  #project (cxsparse CXX)
#else (MSVC)
project (cxsparse LANGUAGES C VERSION 3.1.9)
#endif (MSVC)

set (MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (VERSION_PATCH ${PROJECT_VERSION_PATCH})

set (SHORT_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

if (NOT VERSION_PATCH EQUAL 0)
  set (SHORT_VERSION "${SHORT_VERSION}.${VERSION_PATCH}")
endif (NOT VERSION_PATCH EQUAL 0)

set (VERSION ${PROJECT_VERSION})

set (CPACK_PACKAGE_NAME cxsparse)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "cxsparse: an extended version of CSparse")
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_PACKAGE_VERSION ${SHORT_VERSION})

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (NOT CPack_CMake_INCLUDED)
  include (CPack)
endif ()

include (CheckCSourceCompiles)
include (CheckIncludeFile)
include (CheckLibraryExists)
include (CheckStructHasMember)
include (CheckSymbolExists)
include (CheckTypeSize)
include (CMakePackageConfigHelpers)
include (CMakePushCheckState)
include (ComponentVersion)
include (GenerateExportHeader)
include (GenerateSourceFiles)
include (GNUInstallDirs)

set (SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse)

set (CMAKE_C_VISIBILITY_PRESET hidden)
set (CMAKE_VISIBILITY_INLINES_HIDDEN 1)

check_include_file (complex.h HAVE_COMPLEX_H)

cmake_push_check_state (RESET)
set (CMAKE_EXTRA_INCLUDE_FILES complex.h)

if (HAVE_QSTD_C99)
  set (CMAKE_REQUIRED_FLAGS /Qstd=c99)
endif (HAVE_QSTD_C99)

check_library_exists (m sqrt "" HAVE_M)

check_type_size ("double _Complex" HAVE_C99_COMPLEX BUILTIN_TYPES_ONLY)
check_type_size (_C_double_complex HAVE_C_DOUBLE_COMPLEX)
check_symbol_exists (_DCOMPLEX_ complex.h HAVE__DCOMPLEX_)
check_symbol_exists (_Cmulcc complex.h HAVE__CMULCC)
check_symbol_exists (_Cmulcr complex.h HAVE__CMULCR)
check_symbol_exists (cabs complex.h HAVE_CABS)
check_symbol_exists (cimag complex.h HAVE_CIMAG)
check_symbol_exists (conj complex.h HAVE_CONJ)
check_symbol_exists (creal complex.h HAVE_CREAL)
check_symbol_exists (csqrt complex.h HAVE_CSQRT)
check_struct_has_member (_C_double_complex _Val complex.h
  HAVE_C_DOUBLE_COMPLEX__VAL_MEMBER LANGUAGE C)

cmake_pop_check_state ()

if (HAVE__CMULCC AND HAVE__CMULCR AND HAVE__DCOMPLEX_
  AND HAVE_CABS AND HAVE_CIMAG AND HAVE_CREAL AND HAVE_CSQRT
  AND HAVE_CONJ AND NOT HAVE_C_DOUBLE_COMPLEX STREQUAL "" AND
  HAVE_C_DOUBLE_COMPLEX__VAL_MEMBER)
  set (HAVE_C99_COMPLEX_SUBSET ON)
endif ()

source_group ("Generated Files" REGULAR_EXPRESSION ".+_(di|dl|zi|zl)(_.+)?\\.c$")

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

if (NOT HAVE_C99_COMPLEX STREQUAL "" OR HAVE_C99_COMPLEX_SUBSET)
  generate_source_files (GEN_SRCS_DINT_COMPLEX FILES "${PROJECT_SRCS}" SUFFIX
    _ci COMPILE_DEFINITIONS CS_COMPLEX)
  generate_source_files (GEN_SRCS_DLONG_COMPLEX FILES "${PROJECT_SRCS}" SUFFIX
    _cl COMPILE_DEFINITIONS CS_LONG CS_COMPLEX)
endif (NOT HAVE_C99_COMPLEX STREQUAL "" OR HAVE_C99_COMPLEX_SUBSET)

set (SRCS
  ${HDRS}
  ${SRCS}
  ${GEN_SRCS_DINT}
  ${GEN_SRCS_DLONG}
  ${GEN_SRCS_DINT_COMPLEX}
  ${GEN_SRCS_DLONG_COMPLEX}
)

if (HAVE_COMPLEX_H)
  list (APPEND SRCS Source/cs_convert.c)
endif (HAVE_COMPLEX_H)

if (MSVC)
  #set_source_files_properties (${SRCS} PROPERTIES LANGUAGE CXX)
endif (MSVC)

add_library (cxsparse ${SRCS})

set (_CXSparse_EXPORT_HEADER ${CMAKE_CURRENT_BINARY_DIR}/include/cs_export.h)
set (_CXSparse_CMake_LIBDIR ${CMAKE_INSTALL_LIBDIR}/CXSparse)
set (_CXSparse_CMake_INSTALLDIR ${_CXSparse_CMake_LIBDIR}/cmake)
set (_CXSparse_CMake_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/CXSparse)

generate_export_header (cxsparse EXPORT_FILE_NAME ${_CXSparse_EXPORT_HEADER})

target_sources (cxsparse PRIVATE ${_CXSparse_EXPORT_HEADER})

target_include_directories (cxsparse PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<INSTALL_INTERFACE:${_CXSparse_CMake_INCLUDEDIR}>"
)

set (SUITESPARSECONFIG_HDRS
  ../SuiteSparse_config/SuiteSparse_config.h
)

add_library (cxsparseconfig
  ../SuiteSparse_config/SuiteSparse_config.c
  ${SUITESPARSECONFIG_HDRS}
)

target_include_directories (cxsparseconfig PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse_config>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
)

if (HAVE_M)
  target_link_libraries (cxsparseconfig INTERFACE m)
endif (HAVE_M)

set_property (TARGET cxsparseconfig PROPERTY NAME_ALIAS CXSparse)
set_property (TARGET cxsparseconfig PROPERTY HEADERS_VARIABLE
  SUITESPARSECONFIG_HDRS)
set_property (TARGET cxsparseconfig PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR})
set_property (TARGET cxsparseconfig PROPERTY DEFINE_SYMBOL
  suitesparseconfig_EXPORTS)

set_component_version (cxsparseconfig "${VERSION}")
set_target_properties (cxsparseconfig PROPERTIES EXPORT_NAME Config)

target_link_libraries (cxsparse PUBLIC cxsparseconfig)

if (NOT HAVE_C99_COMPLEX STREQUAL "")
  target_compile_definitions (cxsparse PUBLIC CXSPARSE_C99_COMPLEX)
elseif (HAVE_C99_COMPLEX_SUBSET)
  # C99 subset supported by MSVC
  target_compile_definitions (cxsparse PUBLIC CXSPARSE_C99_COMPLEX_SUBSET)
else (NOT HAVE_C99_COMPLEX STREQUAL "")
  target_compile_definitions (cxsparse PUBLIC NCOMPLEX)
endif (NOT HAVE_C99_COMPLEX STREQUAL "")

if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  target_compile_definitions (cxsparse PRIVATE _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")

set (_CXSparse_HDRS Include/cs.h ${_CXSparse_EXPORT_HEADER})

set_target_properties (cxsparse PROPERTIES PUBLIC_HEADER "${_CXSparse_HDRS}")
set_target_properties (cxsparse PROPERTIES EXPORT_NAME CXSparse)

install (TARGETS cxsparse cxsparseconfig
  EXPORT cxsparse-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${_CXSparse_CMake_INCLUDEDIR}
  LIBRARY DESTINATION ${_CXSparse_CMake_LIBDIR}
  ARCHIVE DESTINATION ${_CXSparse_CMake_LIBDIR})

configure_package_config_file (cxsparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config.cmake
  INSTALL_DESTINATION ${_CXSparse_CMake_INSTALLDIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file (cxsparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY SameMajorVersion)

export (TARGETS cxsparse cxsparseconfig NAMESPACE CXSparse:: FILE cxsparse-targets.cmake)
export (PACKAGE CXSparse)

install (EXPORT cxsparse-targets NAMESPACE CXSparse::
  DESTINATION ${_CXSparse_CMake_INSTALLDIR})
install (FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cxsparse-config-version.cmake
  DESTINATION ${_CXSparse_CMake_INSTALLDIR})

set_target_properties (cxsparse PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_component_version (cxsparse ${VERSION})

add_executable (cs_ci_demo1
  Demo/cs_ci_demo1.c
)
target_link_libraries (cs_ci_demo1 PRIVATE cxsparse)

add_executable (cs_ci_demo2
  Demo/cs_ci_demo.c
  Demo/cs_ci_demo.h
  Demo/cs_ci_demo2.c
)
target_link_libraries (cs_ci_demo2 PRIVATE cxsparse)

add_executable (cs_ci_demo3
  Demo/cs_ci_demo.c
  Demo/cs_ci_demo.h
  Demo/cs_ci_demo3.c
)
target_link_libraries (cs_ci_demo3 PRIVATE cxsparse)

add_executable (cs_cl_demo1
  Demo/cs_cl_demo1.c
)
target_link_libraries (cs_cl_demo1 PRIVATE cxsparse)

add_executable (cs_cl_demo2
  Demo/cs_cl_demo.c
  Demo/cs_cl_demo.h
  Demo/cs_cl_demo2.c
)
target_link_libraries (cs_cl_demo2 PRIVATE cxsparse)

add_executable (cs_cl_demo3
  Demo/cs_cl_demo.c
  Demo/cs_cl_demo.h
  Demo/cs_cl_demo3.c
)
target_link_libraries (cs_cl_demo3 PRIVATE cxsparse)

add_executable (cs_di_demo1
  Demo/cs_di_demo1.c
)
target_link_libraries (cs_di_demo1 PRIVATE cxsparse)

add_executable (cs_di_demo2
  Demo/cs_di_demo.c
  Demo/cs_di_demo.h
  Demo/cs_di_demo2.c
)
target_link_libraries (cs_di_demo2 PRIVATE cxsparse)

add_executable (cs_di_demo3
  Demo/cs_di_demo.c
  Demo/cs_di_demo.h
  Demo/cs_di_demo3.c
)
target_link_libraries (cs_di_demo3 PRIVATE cxsparse)

add_executable (cs_dl_demo1
  Demo/cs_dl_demo1.c
)
target_link_libraries (cs_dl_demo1 PRIVATE cxsparse)

add_executable (cs_dl_demo2
  Demo/cs_dl_demo.c
  Demo/cs_dl_demo.h
  Demo/cs_dl_demo2.c
)
target_link_libraries (cs_dl_demo2 PRIVATE cxsparse)

add_executable (cs_dl_demo3
  Demo/cs_dl_demo.c
  Demo/cs_dl_demo.h
  Demo/cs_dl_demo3.c
)
target_link_libraries (cs_dl_demo3 PRIVATE cxsparse)

add_executable (cs_demo1
  Demo/cs_demo1.c
)
target_link_libraries (cs_demo1 PRIVATE cxsparse)

add_executable (cs_demo2
  Demo/cs_demo.c
  Demo/cs_demo.h
  Demo/cs_demo2.c
)
target_link_libraries (cs_demo2 PRIVATE cxsparse)

add_executable (cs_demo3
  Demo/cs_demo.c
  Demo/cs_demo.h
  Demo/cs_demo3.c
)
target_link_libraries (cs_demo3 PRIVATE cxsparse)

if (MSVC)
  #set_source_files_properties (
  #  Demo/cs_demo.c
  #Demo/cs_demo1.c
  #Demo/cs_demo2.c
  #Demo/cs_demo3.c
  #Demo/cs_ci_demo.c
  #Demo/cs_cl_demo.c
  #Demo/cs_dl_demo.c
  #Demo/cs_ci_demo1.c
  #Demo/cs_ci_demo2.c
  #Demo/cs_ci_demo3.c
  #Demo/cs_dl_demo1.c
  #Demo/cs_dl_demo2.c
  #Demo/cs_dl_demo3.c
  #PROPERTIES LANGUAGE CXX
  #)
  #set_target_properties (
  #  cxsparse
  #  cs_ci_demo1
  #  cs_ci_demo2
  #  cs_ci_demo3
  #  cs_cl_demo1
  #  cs_cl_demo2
  #  cs_cl_demo3
  #  cs_demo1
  #  cs_demo2
  #  cs_demo3
  #  cs_di_demo1
  #  cs_di_demo2
  #  cs_di_demo3
  #  cs_dl_demo1
  #  cs_dl_demo2
  #  cs_dl_demo3
  #  PROPERTIES LINKER_LANGUAGE CXX
  #)
  #target_compile_definitions (
  #  cs_cl_demo1
  #PRIVATE CS_COMPLEX
  #)
endif (MSVC)
