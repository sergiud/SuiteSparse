#
# CMake support layer for SuiteSparse
#
# SPDX-FileCopyrightText: 2026 Sergiu Deitsch
# SPDX-License-Identifier: Apache-2.0
#

cmake_minimum_required (VERSION 3.22)

cmake_policy (SET CMP0069 NEW)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)
  if (CMAKE_VERSION VERSION_LESS 4.4)
    set (CMAKE_EXPERIMENTAL_GENERATE_SBOM
      ca494ed3-b261-4205-a01f-603c95e4cae0)
  else (CMAKE_VERSION VERSION_LESS 4.4)
    set (CMAKE_EXPERIMENTAL_GENERATE_SBOM
      2d856d6d-53e8-488b-a17f-d486d2cac317)
  endif (CMAKE_VERSION VERSION_LESS 4.4)
endif (CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)

if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set (CMAKE_CUDA_ARCHITECTURES 75)
endif (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)

project (SuiteSparse
  DESCRIPTION "A suite of sparse matrix packages"
  HOMEPAGE_URL http://faculty.cse.tamu.edu/davis/suitesparse.html
  LANGUAGES C
  VERSION 7.12.3
)

set (CMAKE_DEBUG_POSTFIX _debug)
set (CMAKE_C_VISIBILITY_PRESET hidden)
set (CMAKE_CXX_VISIBILITY_PRESET hidden)
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_VISIBILITY_INLINES_HIDDEN 1)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SuiteSparse: an extended version of CSparse")

include (CheckCCompilerFlag)
include (CheckFunctionExists)
include (CheckIncludeFile)
include (CheckLanguage)
include (CheckLibraryExists)
include (CheckSymbolExists)
include (CheckTypeSize)
include (CMakeDependentOption)
include (CMakePackageConfigHelpers)
include (CMakePushCheckState)
include (CheckIPOSupported)
option (BUILD_TESTING "Build SuiteSparse tests" OFF)
include (CTest)
include (ComponentVersion)
include (CPack)
include (FeatureSummary)
include (GenerateExportHeader)
include (GNUInstallDirs)
include (InstallRequiredSystemLibraries)
include (SuiteSparseBase)

option (WITH_DEMOS "Build demos" ON)
option (WITH_METIS "Enables METIS support" ON)
option (WITH_OPENMP "Enable OpenMP support" ON)
option (WITH_PRINT "Print diagnostic messages" OFF)
option (WITH_TBB "Enables oneTBB support" ON)
option (SUITESPARSE_ENABLE_IPO
  "Enable IPO or LTO for Release builds when supported" ON)

if (NOT DEFINED CMAKE_LINK_WHAT_YOU_USE)
  set (CMAKE_LINK_WHAT_YOU_USE OFF CACHE BOOL
    "Report shared library dependencies that are not used")
endif (NOT DEFINED CMAKE_LINK_WHAT_YOU_USE)

set (WITH_LICENSE "GPL" CACHE STRING "Software license the binary distribution should adhere")
set_property (CACHE WITH_LICENSE PROPERTY STRINGS "Minimal;GPL;LGPL")

check_language (CXX)
if (DEFINED CMAKE_CUDA_HOST_COMPILER)
  set (SuiteSparse_CUDA_HOST_COMPILER
    "${CMAKE_CUDA_HOST_COMPILER}")
elseif (DEFINED ENV{CUDAHOSTCXX})
  set (SuiteSparse_CUDA_HOST_COMPILER "$ENV{CUDAHOSTCXX}")
elseif (CMAKE_CXX_COMPILER)
  set (SuiteSparse_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
  set (ENV{CUDAHOSTCXX} "${SuiteSparse_CUDA_HOST_COMPILER}")
endif ()
check_language (CUDA)
check_language (Fortran)

cmake_dependent_option (BUILD_SHARED_LIBS "Build shared libraires" OFF "WITH_LICENSE STREQUAL \"GPL\"" ON)

cmake_dependent_option (WITH_FORTRAN "Enables Fortran support" ON "CMAKE_Fortran_COMPILER" OFF)

cmake_dependent_option (BUILD_CXSPARSE "Build CXSparse" ON
  "WITH_LICENSE MATCHES \"GPL\"" OFF)
cmake_dependent_option (WITH_CUDA "Enable CUDA support" ON
  "CMAKE_CXX_COMPILER AND CMAKE_CUDA_COMPILER AND WITH_LICENSE STREQUAL
  \"GPL\"" OFF)

cmake_dependent_option (WITH_CHOLMOD "Enable CHOLMOD" ON
  "WITH_LICENSE MATCHES \"GPL\"" OFF)
cmake_dependent_option (WITH_CHOLESKY "Enable the Cholesky module" ON "WITH_CHOLMOD" OFF)

cmake_dependent_option (WITH_CAMD "Enable interfaces to CAMD, CCOLAMD, CSYMAMD in Partition module" ON
  "WITH_LICENSE MATCHES \"GPL\"" OFF)
cmake_dependent_option (WITH_CHECK "Enable the Check module" ON
  "WITH_LICENSE MATCHES \"GPL\"" OFF)
cmake_dependent_option (WITH_MATRIXOPS "Enable the MatrixOps module" ON
  "WITH_LICENSE STREQUAL \"GPL\" AND WITH_CHOLMOD" OFF)
cmake_dependent_option (WITH_MODIFY "Enable the Modify module" ON
  "WITH_LICENSE STREQUAL \"GPL\" AND WITH_CHOLMOD" OFF)
cmake_dependent_option (WITH_SUPERNODAL "Enable the Supernodal module" ON
  "WITH_CHOLESKY AND CMAKE_CXX_COMPILER AND WITH_LICENSE STREQUAL \"GPL\"" OFF)

if (WITH_CUDA)
  enable_language (CUDA)
else (WITH_CUDA)
  set (CMAKE_DISABLE_FIND_PACKAGE_CUDAToolkit ON)
endif (WITH_CUDA)

# CUDA and SPQR require a C++ compiler
if (WITH_SUPERNODAL OR WITH_CUDA)
  enable_language (CXX)
endif (WITH_SUPERNODAL OR WITH_CUDA)

if (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)
  add_compile_options (
    "$<$<COMPILE_LANGUAGE:CUDA>:--compiler-bindir=${SuiteSparse_CUDA_HOST_COMPILER}>")
endif ()

if (WITH_FORTRAN)
  enable_language (Fortran)
endif (WITH_FORTRAN)

if (SUITESPARSE_ENABLE_IPO AND NOT WITH_CUDA
  AND (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL Release))
  set (_SuiteSparse_IPO_LANGUAGES C)
  if (CMAKE_CXX_COMPILER_LOADED)
    list (APPEND _SuiteSparse_IPO_LANGUAGES CXX)
  endif (CMAKE_CXX_COMPILER_LOADED)
  if (CMAKE_Fortran_COMPILER_LOADED)
    list (APPEND _SuiteSparse_IPO_LANGUAGES Fortran)
  endif (CMAKE_Fortran_COMPILER_LOADED)

  check_ipo_supported (
    RESULT _SuiteSparse_IPO_SUPPORTED
    OUTPUT _SuiteSparse_IPO_OUTPUT
    LANGUAGES ${_SuiteSparse_IPO_LANGUAGES}
  )

  if (_SuiteSparse_IPO_SUPPORTED)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
  else (_SuiteSparse_IPO_SUPPORTED)
    message (WARNING
      "IPO or LTO is not supported for Release builds: "
      "${_SuiteSparse_IPO_OUTPUT}")
  endif (_SuiteSparse_IPO_SUPPORTED)
endif (SUITESPARSE_ENABLE_IPO AND NOT WITH_CUDA
  AND (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL Release))

check_c_compiler_flag (-Qstd=c99 HAVE_QSTD_C99)

if (NOT WITH_METIS)
  set (CMAKE_DISABLE_FIND_PACKAGE_METIS ON)
endif (NOT WITH_METIS)

if (NOT WITH_OPENMP)
  set (CMAKE_DISABLE_FIND_PACKAGE_OpenMP ON)
endif (NOT WITH_OPENMP)

if (NOT WITH_TBB)
  set (CMAKE_DISABLE_FIND_PACKAGE_TBB ON)
endif (NOT WITH_TBB)

find_package (BLAS REQUIRED)
find_package (CUDAToolkit)
find_package (LAPACK REQUIRED)
find_package (METIS 5.0)
find_package (OpenMP COMPONENTS C)
find_package (TBB 2021.4 NO_MODULE)

function (check_options)
  foreach (_OPT ${ARGV})
    # Map PARTITION to WITH_METIS
    if (_OPT STREQUAL PARTITION)
      set (_OPT_NAME METIS)
    else (_OPT STREQUAL PARTITION)
      set (_OPT_NAME ${_OPT})
    endif (_OPT STREQUAL PARTITION)

    if (NOT DEFINED WITH_${_OPT_NAME})
      message (FATAL_ERROR "Option WITH_${_OPT_NAME} is not defined")
    endif (NOT DEFINED WITH_${_OPT_NAME})

    if (NOT WITH_${_OPT_NAME})
      set (N${_OPT} 1 PARENT_SCOPE)
    else (NOT WITH_${_OPT_NAME})
      set (N${_OPT} 0 PARENT_SCOPE)
    endif (NOT WITH_${_OPT_NAME})
  endforeach (_OPT)
endfunction (check_options)

check_options (
  CAMD
  CHECK
  CHOLESKY
  CHOLMOD
  MATRIXOPS
  MODIFY
  PARTITION
  PRINT
  SUPERNODAL
)

if (WITH_LICENSE STREQUAL "GPL")
  set (NGPL 0)
else (WITH_LICENSE STREQUAL "GPL")
  set (NGPL 1)
endif (WITH_LICENSE STREQUAL "GPL")

if (NOT METIS_FOUND)
  set (NPARTITION 1)
endif (NOT METIS_FOUND)

if (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)
  set (CUDA_DEPENDENCY "find_dependency (CUDAToolkit ${CUDAToolkit_VERSION})")
endif (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)

cmake_push_check_state (RESET)
set (CMAKE_REQUIRED_LIBRARIES BLAS::BLAS LAPACK::LAPACK)

check_function_exists (sgemm HAVE_BLAS_NO_UNDERSCORE)
check_function_exists (sgemm_ HAVE_BLAS_UNDERSCORE)

cmake_pop_check_state ()

check_include_file (complex.h HAVE_COMPLEX_H)
check_include_file (fenv.h HAVE_FENV_H)
check_library_exists (m ceil "" HAVE_M)

cmake_push_check_state (RESET)
set (CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)

if (HAVE_M)
  set (CMAKE_REQUIRED_LIBRARIES m)
endif (HAVE_M)

check_symbol_exists (feenableexcept "fenv.h" HAVE_FEENABLEEXCEPT)
cmake_pop_check_state ()

check_symbol_exists (_controlfp "float.h" HAVE__CONTROLFP)

set (CMAKE_EXTRA_INCLUDE_FILES complex.h)
check_type_size (_C_double_complex HAVE_C_DOUBLE_COMPLEX BUILTIN_TYPES_ONLY)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS
    _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

set (AMD_HDRS
  AMD/Include/amd.h
  AMD/Include/amd_internal.h
)

set (AMD_SRCS
  AMD/Source/amd_1.c
  AMD/Source/amd_2.c
  AMD/Source/amd_aat.c
  AMD/Source/amd_control.c
  AMD/Source/amd_defaults.c
  AMD/Source/amd_dump.c
  AMD/Source/amd_info.c
  AMD/Source/amd_order.c
  AMD/Source/amd_post_tree.c
  AMD/Source/amd_postorder.c
  AMD/Source/amd_preprocess.c
  AMD/Source/amd_valid.c
)

set (AMD_UMFPACK_SRCS
  AMD/Source/amd_aat.c
  AMD/Source/amd_dump.c
  AMD/Source/amd_post_tree.c
  AMD/Source/amd_postorder.c
)

set (BTF_HDRS
  BTF/Include/btf.h
  BTF/Include/btf_internal.h
)

set (BTF_SRCS
  BTF/Source/btf_maxtrans.c
  BTF/Source/btf_order.c
  BTF/Source/btf_strongcomp.c
)

set (CAMD_HDRS
  CAMD/Include/camd.h
  CAMD/Include/camd_internal.h
)

set (CAMD_SRCS
  CAMD/Source/camd_1.c
  CAMD/Source/camd_2.c
  CAMD/Source/camd_aat.c
  CAMD/Source/camd_control.c
  CAMD/Source/camd_defaults.c
  CAMD/Source/camd_dump.c
  CAMD/Source/camd_info.c
  CAMD/Source/camd_order.c
  CAMD/Source/camd_postorder.c
  CAMD/Source/camd_preprocess.c
  CAMD/Source/camd_valid.c
)

if (WITH_CHOLMOD)
  set (CHOLMOD_HDRS
    ${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include/cholmod_config.h
    CHOLMOD/Include/cholmod.h
    CHOLMOD/Include/cholmod_internal.h
    CHOLMOD/Include/cholmod_template.h
    CHOLMOD/Include/cholmod_types.h
  )

  file (GLOB _CHOLMOD_TYPED_SRCS
    CHOLMOD/Check/cholmod_*.c
    CHOLMOD/Cholesky/cholmod_*.c
    CHOLMOD/Utility/cholmod_*.c
    CHOLMOD/MatrixOps/cholmod_*.c
    CHOLMOD/Modify/cholmod_*.c
    CHOLMOD/Partition/cholmod_*.c
    CHOLMOD/Supernodal/cholmod_*.c
  )
  list (FILTER _CHOLMOD_TYPED_SRCS EXCLUDE REGEX
    "/cholmod_(allocate_factor|allocate_work|divcomplex|hypot|memdebug|mult_uint64_t)\\.c$")

  set (CHOLMOD_COMMON_SRCS
    CHOLMOD/Utility/cholmod_allocate_factor.c
    CHOLMOD/Utility/cholmod_allocate_work.c
    CHOLMOD/Utility/cholmod_divcomplex.c
    CHOLMOD/Utility/cholmod_hypot.c
    CHOLMOD/Utility/cholmod_memdebug.c
    CHOLMOD/Utility/cholmod_mult_uint64_t.c
  )

  file (GLOB _CHOLMOD_TEMPLATE_SRCS
    CHOLMOD/Utility/t_cholmod_*.c
  )
  list (FILTER _CHOLMOD_TEMPLATE_SRCS EXCLUDE REGEX
    "/t_cholmod_(bound|malloc|.*(_worker|_template|_numkr|_wdim|_(un)?permuted))\\.c$")

  set (CHOLMOD_SRCS
    ${_CHOLMOD_TYPED_SRCS}
    ${_CHOLMOD_TEMPLATE_SRCS}
    ${CHOLMOD_COMMON_SRCS}
  )
  set (CHOLMOD_LONG_SRCS
    ${_CHOLMOD_TYPED_SRCS}
    ${_CHOLMOD_TEMPLATE_SRCS}
  )

  if (WITH_CUDA)
    list (APPEND CHOLMOD_HDRS
      CHOLMOD/GPU/cholmod_gpu_kernels.h
    )
  endif (WITH_CUDA)

  if (NOT METIS_FOUND)
    set_property (SOURCE
      CHOLMOD/Check/cholmod_check.c
      CHOLMOD/Partition/cholmod_camd.c
      CHOLMOD/Partition/cholmod_ccolamd.c
      CHOLMOD/Partition/cholmod_csymamd.c
      CHOLMOD/Partition/cholmod_metis.c
      CHOLMOD/Partition/cholmod_nesdis.c
      APPEND PROPERTY COMPILE_DEFINITIONS NPARTITION)
  endif (NOT METIS_FOUND)

  configure_file (CHOLMOD/Include/cholmod_config.h.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include/cholmod_config.h)
endif (WITH_CHOLMOD)

set (CCOLAMD_HDRS
  CCOLAMD/Include/ccolamd.h
)

set (CCOLAMD_SRCS
  CCOLAMD/Source/ccolamd.c
)

set (COLAMD_HDRS
  COLAMD/Include/colamd.h
)

set (COLAMD_SRCS
  COLAMD/Source/colamd.c
)

set (UMFPACK_HDRS
  UMFPACK/Include/umfpack.h
)

set (UMFPACK_INT_SRCS
  UMFPACK/Source/umf_analyze.c
  UMFPACK/Source/umf_apply_order.c
  UMFPACK/Source/umf_assemble.c
  UMFPACK/Source/umf_blas3_update.c
  UMFPACK/Source/umf_build_tuples.c
  UMFPACK/Source/umf_cholmod.c
  UMFPACK/Source/umf_colamd.c
  UMFPACK/Source/umf_create_element.c
  UMFPACK/Source/umf_dump.c
  UMFPACK/Source/umf_extend_front.c
  UMFPACK/Source/umf_free.c
  UMFPACK/Source/umf_fsize.c
  UMFPACK/Source/umf_garbage_collection.c
  UMFPACK/Source/umf_get_memory.c
  UMFPACK/Source/umf_grow_front.c
  UMFPACK/Source/umf_init_front.c
  UMFPACK/Source/umf_is_permutation.c
  UMFPACK/Source/umf_kernel.c
  UMFPACK/Source/umf_kernel_init.c
  UMFPACK/Source/umf_kernel_wrapup.c
  UMFPACK/Source/umf_local_search.c
  UMFPACK/Source/umf_lsolve.c
  UMFPACK/Source/umf_ltsolve.c
  UMFPACK/Source/umf_malloc.c
  UMFPACK/Source/umf_mem_alloc_element.c
  UMFPACK/Source/umf_mem_alloc_head_block.c
  UMFPACK/Source/umf_mem_alloc_tail_block.c
  UMFPACK/Source/umf_mem_free_tail_block.c
  UMFPACK/Source/umf_mem_init_memoryspace.c
  UMFPACK/Source/umf_realloc.c
  UMFPACK/Source/umf_report_perm.c
  UMFPACK/Source/umf_report_vector.c
  UMFPACK/Source/umf_row_search.c
  UMFPACK/Source/umf_scale.c
  UMFPACK/Source/umf_scale_column.c
  UMFPACK/Source/umf_set_stats.c
  UMFPACK/Source/umf_singletons.c
  UMFPACK/Source/umf_solve.c
  UMFPACK/Source/umf_start_front.c
  UMFPACK/Source/umf_store_lu.c
  UMFPACK/Source/umf_symbolic_usage.c
  UMFPACK/Source/umf_transpose.c
  UMFPACK/Source/umf_triplet.c
  UMFPACK/Source/umf_tuple_lengths.c
  UMFPACK/Source/umf_usolve.c
  UMFPACK/Source/umf_utsolve.c
  UMFPACK/Source/umf_valid_numeric.c
  UMFPACK/Source/umf_valid_symbolic.c
  UMFPACK/Source/umfpack_col_to_triplet.c
  UMFPACK/Source/umfpack_copy_numeric.c
  UMFPACK/Source/umfpack_copy_symbolic.c
  UMFPACK/Source/umfpack_defaults.c
  UMFPACK/Source/umfpack_deserialize_numeric.c
  UMFPACK/Source/umfpack_deserialize_symbolic.c
  UMFPACK/Source/umfpack_free_numeric.c
  UMFPACK/Source/umfpack_free_symbolic.c
  UMFPACK/Source/umfpack_get_determinant.c
  UMFPACK/Source/umfpack_get_lunz.c
  UMFPACK/Source/umfpack_get_numeric.c
  UMFPACK/Source/umfpack_get_symbolic.c
  UMFPACK/Source/umfpack_load_numeric.c
  UMFPACK/Source/umfpack_load_symbolic.c
  UMFPACK/Source/umfpack_numeric.c
  UMFPACK/Source/umfpack_qsymbolic.c
  UMFPACK/Source/umfpack_report_control.c
  UMFPACK/Source/umfpack_report_info.c
  UMFPACK/Source/umfpack_report_matrix.c
  UMFPACK/Source/umfpack_report_numeric.c
  UMFPACK/Source/umfpack_report_perm.c
  UMFPACK/Source/umfpack_report_status.c
  UMFPACK/Source/umfpack_report_symbolic.c
  UMFPACK/Source/umfpack_report_triplet.c
  UMFPACK/Source/umfpack_report_vector.c
  UMFPACK/Source/umfpack_save_numeric.c
  UMFPACK/Source/umfpack_save_symbolic.c
  UMFPACK/Source/umfpack_scale.c
  UMFPACK/Source/umfpack_serialize_numeric.c
  UMFPACK/Source/umfpack_serialize_symbolic.c
  UMFPACK/Source/umfpack_solve.c
  UMFPACK/Source/umfpack_symbolic.c
  UMFPACK/Source/umfpack_transpose.c
  UMFPACK/Source/umfpack_triplet_to_col.c
)

set (UMFPACK_ZINT_SRCS
  #UMFPACK/Source/umf_apply_order.c
  UMFPACK/Source/umf_assemble.c
  UMFPACK/Source/umf_blas3_update.c
  UMFPACK/Source/umf_build_tuples.c
  #UMFPACK/Source/umf_cholmod.c
  #UMFPACK/Source/umf_colamd.c
  UMFPACK/Source/umf_create_element.c
  UMFPACK/Source/umf_dump.c
  UMFPACK/Source/umf_extend_front.c
  #UMFPACK/Source/umf_free.c
  #UMFPACK/Source/umf_fsize.c
  UMFPACK/Source/umf_garbage_collection.c
  UMFPACK/Source/umf_get_memory.c
  UMFPACK/Source/umf_grow_front.c
  UMFPACK/Source/umf_init_front.c
  #UMFPACK/Source/umf_is_permutation.c
  UMFPACK/Source/umf_kernel.c
  UMFPACK/Source/umf_kernel_init.c
  UMFPACK/Source/umf_kernel_wrapup.c
  UMFPACK/Source/umf_local_search.c
  UMFPACK/Source/umf_lsolve.c
  UMFPACK/Source/umf_ltsolve.c
  #UMFPACK/Source/umf_malloc.c
  UMFPACK/Source/umf_mem_alloc_element.c
  UMFPACK/Source/umf_mem_alloc_head_block.c
  UMFPACK/Source/umf_mem_alloc_tail_block.c
  UMFPACK/Source/umf_mem_free_tail_block.c
  UMFPACK/Source/umf_mem_init_memoryspace.c
  #UMFPACK/Source/umf_realloc.c
  #UMFPACK/Source/umf_report_perm.c
  UMFPACK/Source/umf_report_vector.c
  UMFPACK/Source/umf_row_search.c
  UMFPACK/Source/umf_scale.c
  UMFPACK/Source/umf_scale_column.c
  UMFPACK/Source/umf_set_stats.c
  #UMFPACK/Source/umf_singletons.c
  UMFPACK/Source/umf_solve.c
  UMFPACK/Source/umf_start_front.c
  UMFPACK/Source/umf_store_lu.c
  UMFPACK/Source/umf_symbolic_usage.c
  UMFPACK/Source/umf_transpose.c
  UMFPACK/Source/umf_triplet.c
  UMFPACK/Source/umf_tuple_lengths.c
  UMFPACK/Source/umf_usolve.c
  UMFPACK/Source/umf_utsolve.c
  UMFPACK/Source/umf_valid_numeric.c
  UMFPACK/Source/umf_valid_symbolic.c
  UMFPACK/Source/umfpack_col_to_triplet.c
  UMFPACK/Source/umfpack_copy_numeric.c
  UMFPACK/Source/umfpack_copy_symbolic.c
  UMFPACK/Source/umfpack_defaults.c
  UMFPACK/Source/umfpack_deserialize_numeric.c
  UMFPACK/Source/umfpack_deserialize_symbolic.c
  UMFPACK/Source/umfpack_free_numeric.c
  UMFPACK/Source/umfpack_free_symbolic.c
  UMFPACK/Source/umfpack_get_determinant.c
  UMFPACK/Source/umfpack_get_lunz.c
  UMFPACK/Source/umfpack_get_numeric.c
  UMFPACK/Source/umfpack_get_symbolic.c
  UMFPACK/Source/umfpack_load_numeric.c
  UMFPACK/Source/umfpack_load_symbolic.c
  UMFPACK/Source/umfpack_numeric.c
  UMFPACK/Source/umfpack_qsymbolic.c
  UMFPACK/Source/umfpack_report_control.c
  UMFPACK/Source/umfpack_report_info.c
  UMFPACK/Source/umfpack_report_matrix.c
  UMFPACK/Source/umfpack_report_numeric.c
  UMFPACK/Source/umfpack_report_perm.c
  UMFPACK/Source/umfpack_report_status.c
  UMFPACK/Source/umfpack_report_symbolic.c
  UMFPACK/Source/umfpack_report_triplet.c
  UMFPACK/Source/umfpack_report_vector.c
  UMFPACK/Source/umfpack_save_numeric.c
  UMFPACK/Source/umfpack_save_symbolic.c
  UMFPACK/Source/umfpack_scale.c
  UMFPACK/Source/umfpack_serialize_numeric.c
  UMFPACK/Source/umfpack_serialize_symbolic.c
  UMFPACK/Source/umfpack_solve.c
  UMFPACK/Source/umfpack_symbolic.c
  UMFPACK/Source/umfpack_transpose.c
  UMFPACK/Source/umfpack_triplet_to_col.c
)

set (KLU_HDRS
  KLU/Include/klu.h
  KLU/Include/klu_internal.h
  KLU/Include/klu_version.h
)

set (KLU_CHOLMOD_HDRS
  KLU/User/klu_cholmod.h
)

set (KLU_SRCS
  KLU/Source/klu.c
  KLU/Source/klu_diagnostics.c
  KLU/Source/klu_dump.c
  KLU/Source/klu_extract.c
  KLU/Source/klu_factor.c
  KLU/Source/klu_free_numeric.c
  KLU/Source/klu_kernel.c
  KLU/Source/klu_refactor.c
  KLU/Source/klu_scale.c
  KLU/Source/klu_solve.c
  KLU/Source/klu_sort.c
  KLU/Source/klu_tsolve.c
)

set (KLU
  KLU/Include/klu.h
  KLU/Include/klu_version.h
)

set (LDL_HDRS
  LDL/Include/ldl.h
)

set (LDL_SRCS
  LDL/Source/ldl.c
)

set (GPURUNTIME_HDRS
  SPQR/GPURuntime/Include/SuiteSparse_GPURuntime.hpp
  SPQR/GPURuntime/Include/SuiteSparseGPU_Workspace.hpp
  SPQR/GPURuntime/Include/SuiteSparseGPU_debug.hpp
  SPQR/GPURuntime/Include/SuiteSparseGPU_internal.hpp
  SPQR/GPURuntime/Include/SuiteSparseGPU_macros.hpp
  SPQR/GPURuntime/Include/SuiteSparseGPU_workspace_macros.hpp
)

set (GPUQRENGINE_HDRS
  SPQR/GPUQREngine/Include/GPUQREngine.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_BucketList.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Common.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Front.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_FrontState.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_GraphVizHelper.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Internal.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_LLBundle.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Scheduler.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_SEntry.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_SparseMeta.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Stats.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_TaskDescriptor.hpp
  SPQR/GPUQREngine/Include/GPUQREngine_Timing.hpp
)

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

# output directory for binaries
set (OUTPUT_DIR ${CMAKE_BINARY_DIR})

# put binaries in a single directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR})

add_library (amd_int OBJECT ${AMD_SRCS})
target_include_directories (amd_int PRIVATE $<TARGET_PROPERTY:amd,INCLUDE_DIRECTORIES>)
target_compile_definitions (amd_int PRIVATE DINT)
target_compile_definitions (amd_int PRIVATE amd_EXPORTS)

add_library (amd_long OBJECT ${AMD_SRCS})
target_include_directories (amd_long PRIVATE $<TARGET_PROPERTY:amd,INCLUDE_DIRECTORIES>)
target_compile_definitions (amd_long PRIVATE DLONG)
target_compile_definitions (amd_long PRIVATE amd_EXPORTS)

add_library (amd
  $<TARGET_OBJECTS:amd_int>
  $<TARGET_OBJECTS:amd_long>
  AMD/Source/amd_version.c
  ${AMD_HDRS}
)
add_library (SuiteSparse::AMD ALIAS amd)

add_library (amd_umfpack_int OBJECT ${AMD_UMFPACK_SRCS})
target_include_directories (amd_umfpack_int PRIVATE
  $<TARGET_PROPERTY:amd,INCLUDE_DIRECTORIES>)
target_compile_definitions (amd_umfpack_int PRIVATE DINT)

add_library (amd_umfpack_long OBJECT ${AMD_UMFPACK_SRCS})
target_include_directories (amd_umfpack_long PRIVATE
  $<TARGET_PROPERTY:amd,INCLUDE_DIRECTORIES>)
target_compile_definitions (amd_umfpack_long PRIVATE DLONG)

if (WITH_LICENSE MATCHES "GPL")
  add_library (btf_int OBJECT ${BTF_SRCS})
  target_include_directories (btf_int PRIVATE $<TARGET_PROPERTY:btf,INCLUDE_DIRECTORIES>)
  target_compile_definitions (btf_int PRIVATE DINT)
  target_compile_definitions (btf_int PRIVATE btf_EXPORTS)

  add_library (btf_long OBJECT ${BTF_SRCS})
  target_include_directories (btf_long PRIVATE $<TARGET_PROPERTY:btf,INCLUDE_DIRECTORIES>)
  target_compile_definitions (btf_long PRIVATE DLONG)
  target_compile_definitions (btf_long PRIVATE btf_EXPORTS)

  add_library (btf
    $<TARGET_OBJECTS:btf_int>
    $<TARGET_OBJECTS:btf_long>
    BTF/Source/btf_version.c
    ${BTF_HDRS}
  )
  add_library (SuiteSparse::BTF ALIAS btf)

  target_include_directories (btf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
  )
endif (WITH_LICENSE MATCHES "GPL")

add_library (camd_int OBJECT ${CAMD_SRCS})
target_include_directories (camd_int PRIVATE $<TARGET_PROPERTY:camd,INCLUDE_DIRECTORIES>)
target_compile_definitions (camd_int PRIVATE DINT)
target_compile_definitions (camd_int PRIVATE camd_EXPORTS)

add_library (camd_long OBJECT ${CAMD_SRCS})
target_include_directories (camd_long PRIVATE $<TARGET_PROPERTY:camd,INCLUDE_DIRECTORIES>)
target_compile_definitions (camd_long PRIVATE DLONG)
target_compile_definitions (camd_long PRIVATE camd_EXPORTS)

add_library (camd
  $<TARGET_OBJECTS:camd_int>
  $<TARGET_OBJECTS:camd_long>
  CAMD/Source/camd_version.c
  ${CAMD_HDRS}
)
add_library (SuiteSparse::CAMD ALIAS camd)

if (WITH_CHOLMOD)
  add_library (cholmod_int OBJECT ${CHOLMOD_SRCS})
  target_include_directories (cholmod_int PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_int PRIVATE CHOLMOD_INT32)
  target_compile_definitions (cholmod_int PRIVATE cholmod_EXPORTS)

  set (CHOLMOD_BOUND_SRC CHOLMOD/Utility/t_cholmod_bound.c)
  foreach (CHOLMOD_BOUND_VARIANT
      dbound_int sbound_int dbound_long sbound_long)
    add_library (cholmod_${CHOLMOD_BOUND_VARIANT} OBJECT
      ${CHOLMOD_BOUND_SRC})
    target_include_directories (cholmod_${CHOLMOD_BOUND_VARIANT} PRIVATE
      $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_${CHOLMOD_BOUND_VARIANT} PRIVATE
      cholmod_EXPORTS)
  endforeach ()
  target_compile_definitions (cholmod_dbound_int PRIVATE CHOLMOD_INT32)
  target_compile_definitions (cholmod_sbound_int PRIVATE CHOLMOD_INT32)
  target_compile_definitions (cholmod_dbound_long PRIVATE CHOLMOD_INT64)
  target_compile_definitions (cholmod_sbound_long PRIVATE CHOLMOD_INT64)
  target_compile_options (cholmod_dbound_int PRIVATE
    "-DCHOLMOD_BOUND_FUNCTION=cholmod_dbound"
    "-DCOMMON_BOUND=(Common->dbound)"
    "-DCOMMON_BOUNDS_HIT=(Common->ndbounds_hit)" "-DReal=double")
  target_compile_options (cholmod_sbound_int PRIVATE
    "-DCHOLMOD_BOUND_FUNCTION=cholmod_sbound"
    "-DCOMMON_BOUND=(Common->sbound)"
    "-DCOMMON_BOUNDS_HIT=(Common->nsbounds_hit)" "-DReal=float")
  target_compile_options (cholmod_dbound_long PRIVATE
    "-DCHOLMOD_BOUND_FUNCTION=cholmod_l_dbound"
    "-DCOMMON_BOUND=(Common->dbound)"
    "-DCOMMON_BOUNDS_HIT=(Common->ndbounds_hit)" "-DReal=double")
  target_compile_options (cholmod_sbound_long PRIVATE
    "-DCHOLMOD_BOUND_FUNCTION=cholmod_l_sbound"
    "-DCOMMON_BOUND=(Common->sbound)"
    "-DCOMMON_BOUNDS_HIT=(Common->nsbounds_hit)" "-DReal=float")

  set (CHOLMOD_MEMORY_SRC CHOLMOD/Utility/t_cholmod_malloc.c)
  foreach (CHOLMOD_MEMORY_VARIANT
      malloc_int calloc_int malloc_long calloc_long)
    add_library (cholmod_${CHOLMOD_MEMORY_VARIANT} OBJECT
      ${CHOLMOD_MEMORY_SRC})
    target_include_directories (cholmod_${CHOLMOD_MEMORY_VARIANT} PRIVATE
      $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_${CHOLMOD_MEMORY_VARIANT} PRIVATE
      cholmod_EXPORTS)
  endforeach ()
  target_compile_definitions (cholmod_malloc_int PRIVATE CHOLMOD_INT32)
  target_compile_definitions (cholmod_calloc_int PRIVATE CHOLMOD_INT32)
  target_compile_definitions (cholmod_malloc_long PRIVATE CHOLMOD_INT64)
  target_compile_definitions (cholmod_calloc_long PRIVATE CHOLMOD_INT64)
  target_compile_options (cholmod_malloc_int PRIVATE
    "-DCHOLMOD_ALLOC_FUNCTION=cholmod_malloc"
    "-DSUITESPARSE_ALLOC_FUNCTION=SuiteSparse_malloc")
  target_compile_options (cholmod_calloc_int PRIVATE
    "-DCHOLMOD_ALLOC_FUNCTION=cholmod_calloc"
    "-DSUITESPARSE_ALLOC_FUNCTION=SuiteSparse_calloc")
  target_compile_options (cholmod_malloc_long PRIVATE
    "-DCHOLMOD_ALLOC_FUNCTION=cholmod_l_malloc"
    "-DSUITESPARSE_ALLOC_FUNCTION=SuiteSparse_malloc")
  target_compile_options (cholmod_calloc_long PRIVATE
    "-DCHOLMOD_ALLOC_FUNCTION=cholmod_l_calloc"
    "-DSUITESPARSE_ALLOC_FUNCTION=SuiteSparse_calloc")

  add_library (cholmod_allocate_work_long OBJECT
    CHOLMOD/Utility/cholmod_allocate_work.c)
  target_include_directories (cholmod_allocate_work_long PRIVATE
    $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_allocate_work_long PRIVATE
    CHOLMOD_INT64
    cholmod_EXPORTS)
  target_compile_options (cholmod_allocate_work_long PRIVATE
    "-Dcholmod_alloc_work=cholmod_l_alloc_work"
    "-Dcholmod_allocate_work=cholmod_l_allocate_work")

  add_library (cholmod_allocate_factor_long OBJECT
    CHOLMOD/Utility/cholmod_allocate_factor.c)
  target_include_directories (cholmod_allocate_factor_long PRIVATE
    $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_allocate_factor_long PRIVATE
    CHOLMOD_INT64
    cholmod_EXPORTS)

  add_library (cholmod_long OBJECT ${CHOLMOD_LONG_SRCS})
  target_include_directories (cholmod_long PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_long PRIVATE CHOLMOD_INT64)
  target_compile_definitions (cholmod_long PRIVATE cholmod_EXPORTS)

  add_library (cholmod
    $<TARGET_OBJECTS:cholmod_int>
    $<TARGET_OBJECTS:cholmod_long>
    $<TARGET_OBJECTS:cholmod_dbound_int>
    $<TARGET_OBJECTS:cholmod_sbound_int>
    $<TARGET_OBJECTS:cholmod_dbound_long>
    $<TARGET_OBJECTS:cholmod_sbound_long>
    $<TARGET_OBJECTS:cholmod_malloc_int>
    $<TARGET_OBJECTS:cholmod_calloc_int>
    $<TARGET_OBJECTS:cholmod_malloc_long>
    $<TARGET_OBJECTS:cholmod_calloc_long>
    $<TARGET_OBJECTS:cholmod_allocate_work_long>
    $<TARGET_OBJECTS:cholmod_allocate_factor_long>
    ${CHOLMOD_HDRS}
  )

  if (METIS_FOUND)
    set (_METIS_IN ${SuiteSparse_SOURCE_DIR}/cmake/FindMETIS.cmake)
    set (_METIS_OUT ${SuiteSparse_BINARY_DIR}/FindMETIS.cmake)

    add_custom_command (
      OUTPUT ${_METIS_OUT}
      COMMAND ${CMAKE_COMMAND} -E create_symlink ${_METIS_IN} ${_METIS_OUT}
      DEPENDS ${_METIS_IN}
      COMMENT "Creating FindMETIS.cmake symlink in build directory"
    )

    target_sources (cholmod PRIVATE ${_METIS_OUT})
  endif (METIS_FOUND)

  add_library (SuiteSparse::CHOLMOD ALIAS cholmod)

  target_include_directories (cholmod PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include>
  )

  if (OpenMP_FOUND)
    target_compile_definitions (cholmod PRIVATE HAVE_OPENMP)
    target_link_libraries (cholmod PRIVATE
      OpenMP::OpenMP_C
    )

    set (OpenMP_DEPENDENCY "find_dependency (OpenMP ${OpenMP_C_VERSION} COMPONENTS C)")
  endif (OpenMP_FOUND)

  if (HAVE_FEENABLEEXCEPT)
    target_compile_definitions (cholmod PRIVATE _GNU_SOURCE)
  endif (HAVE_FEENABLEEXCEPT)

  if (METIS_FOUND)
    target_link_libraries (cholmod PRIVATE METIS::METIS)

    add_library (cholmod_partition INTERFACE)
    add_library (SuiteSparse::Partition ALIAS cholmod_partition)
    target_link_libraries (cholmod_partition INTERFACE cholmod)

    set (METIS_DEPENDENCY "find_dependency (METIS ${METIS_VERSION})")
  endif (METIS_FOUND)

  target_link_libraries (cholmod
    PRIVATE
      BLAS::BLAS
      LAPACK::LAPACK
      SuiteSparse::AMD
      SuiteSparse::CAMD
      SuiteSparse::CCOLAMD
      SuiteSparse::COLAMD
    PRIVATE
      SuiteSparse::Config
  )
endif (WITH_CHOLMOD)

add_library (ccolamd_int OBJECT ${CCOLAMD_SRCS})
target_include_directories (ccolamd_int PRIVATE $<TARGET_PROPERTY:ccolamd,INCLUDE_DIRECTORIES>)
target_compile_definitions (ccolamd_int PRIVATE DINT)
target_compile_definitions (ccolamd_int PRIVATE ccolamd_EXPORTS)

add_library (ccolamd_long OBJECT ${CCOLAMD_SRCS})
target_include_directories (ccolamd_long PRIVATE $<TARGET_PROPERTY:ccolamd,INCLUDE_DIRECTORIES>)
target_compile_definitions (ccolamd_long PRIVATE DLONG)
target_compile_definitions (ccolamd_long PRIVATE ccolamd_EXPORTS)

add_library (ccolamd
  $<TARGET_OBJECTS:ccolamd_int>
  $<TARGET_OBJECTS:ccolamd_long>
  CCOLAMD/Source/ccolamd_version.c
  ${CCOLAMD_HDRS}
)
add_library (SuiteSparse::CCOLAMD ALIAS ccolamd)

add_library (colamd_int OBJECT ${COLAMD_SRCS})
target_include_directories (colamd_int PRIVATE $<TARGET_PROPERTY:colamd,INCLUDE_DIRECTORIES>)
target_compile_definitions (colamd_int PRIVATE DINT)
target_compile_definitions (colamd_int PRIVATE colamd_EXPORTS)

add_library (colamd_long OBJECT ${COLAMD_SRCS})
target_include_directories (colamd_long PRIVATE $<TARGET_PROPERTY:colamd,INCLUDE_DIRECTORIES>)
target_compile_definitions (colamd_long PRIVATE DLONG)
target_compile_definitions (colamd_long PRIVATE colamd_EXPORTS)

add_library (colamd
  $<TARGET_OBJECTS:colamd_int>
  $<TARGET_OBJECTS:colamd_long>
  COLAMD/Source/colamd_version.c
  ${COLAMD_HDRS}
)
add_library (SuiteSparse::COLAMD ALIAS colamd)

if (WITH_SUPERNODAL)
  set (SPQR_HDRS
    SPQR/Include/SuiteSparseQR.hpp
    SPQR/Include/SuiteSparseQR_C.h
    SPQR/Include/SuiteSparseQR_definitions.h
  )

  add_library (spqr
    ${SPQR_HDRS}
    SPQR/Source/spqr_1colamd.cpp
    SPQR/Source/spqr_1factor.cpp
    SPQR/Source/spqr_1fixed.cpp
    SPQR/Source/spqr_analyze.cpp
    SPQR/Source/spqr_append.cpp
    SPQR/Source/spqr_assemble.cpp
    SPQR/Source/spqr_cpack.cpp
    SPQR/Source/spqr_csize.cpp
    SPQR/Source/spqr_cumsum.cpp
    SPQR/Source/spqr_cholmod_wrappers.cpp
    SPQR/Source/spqr_debug.cpp
    SPQR/Source/spqr_factorize.cpp
    SPQR/Source/spqr_fcsize.cpp
    SPQR/Source/spqr_freefac.cpp
    SPQR/Source/spqr_freenum.cpp
    SPQR/Source/spqr_freesym.cpp
    SPQR/Source/spqr_front.cpp
    SPQR/Source/spqr_fsize.cpp
    SPQR/Source/spqr_happly.cpp
    SPQR/Source/spqr_happly_work.cpp
    SPQR/Source/spqr_hpinv.cpp
    SPQR/Source/spqr_kernel.cpp
    SPQR/Source/spqr_larftb.cpp
    SPQR/Source/spqr_maxcolnorm.cpp
    SPQR/Source/spqr_panel.cpp
    SPQR/Source/spqr_parallel.cpp
    SPQR/Source/spqr_rconvert.cpp
    SPQR/Source/spqr_rcount.cpp
    SPQR/Source/spqr_rhpack.cpp
    SPQR/Source/spqr_rmap.cpp
    SPQR/Source/spqr_rsolve.cpp
    SPQR/Source/spqr_shift.cpp
    SPQR/Source/spqr_stranspose1.cpp
    SPQR/Source/spqr_stranspose2.cpp
    SPQR/Source/spqr_tol.cpp
    SPQR/Source/spqr_trapezoidal.cpp
    SPQR/Source/spqr_type.cpp
    SPQR/Source/SuiteSparseQR.cpp
    SPQR/Source/SuiteSparseQR_C.cpp
    SPQR/Source/SuiteSparseQR_expert.cpp
    SPQR/Source/SuiteSparseQR_qmult.cpp
  )
  add_library (SuiteSparse::SPQR ALIAS spqr)

  if (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)
    target_sources (spqr PRIVATE
      SPQR/SPQRGPU/spqrgpu_buildAssemblyMaps.cpp
      SPQR/SPQRGPU/spqrgpu_computeFrontStaging.cpp
      SPQR/SPQRGPU/spqrgpu_kernel.cpp
    )

    target_link_libraries (spqr PRIVATE SuiteSparse::GPURuntime)
  endif (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)

  if (TBB_FOUND)
    target_link_libraries (spqr PRIVATE TBB::tbb)
    target_compile_definitions (spqr PRIVATE HAVE_TBB)
    target_compile_features (spqr PRIVATE cxx_std_11)

    set (TBB_DEPENDENCY "find_dependency (TBB ${TBB_VERSION})")
  endif (TBB_FOUND)
endif (WITH_SUPERNODAL)

add_library (suitesparseconfig
  $<TARGET_OBJECTS:suitesparsebase>
)
add_library (SuiteSparse::Config ALIAS suitesparseconfig)

if (WITH_LICENSE STREQUAL "GPL")
  add_library (umfpack_int OBJECT ${UMFPACK_INT_SRCS})
  target_include_directories (umfpack_int PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int PRIVATE DINT)
  target_compile_definitions (umfpack_int PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long OBJECT ${UMFPACK_INT_SRCS})
  target_include_directories (umfpack_long PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long PRIVATE DLONG)
  target_compile_definitions (umfpack_long PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint OBJECT ${UMFPACK_ZINT_SRCS})
  target_include_directories (umfpack_zint PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint PRIVATE ZINT)
  target_compile_definitions (umfpack_zint PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong OBJECT ${UMFPACK_ZINT_SRCS})
  target_include_directories (umfpack_zlong PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong PRIVATE ZLONG)
  target_compile_definitions (umfpack_zlong PRIVATE umfpack_EXPORTS)

  set (UMFPACK_SOLVE_SRCS
    UMFPACK/Source/umf_ltsolve.c
    UMFPACK/Source/umf_utsolve.c
  )

  add_library (umfpack_int_solve OBJECT ${UMFPACK_SOLVE_SRCS})
  target_include_directories (umfpack_int_solve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_solve PRIVATE DINT CONJUGATE_SOLVE)
  target_compile_definitions (umfpack_int_solve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_solve OBJECT ${UMFPACK_SOLVE_SRCS})
  target_include_directories (umfpack_long_solve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_solve PRIVATE DLONG CONJUGATE_SOLVE)
  target_compile_definitions (umfpack_long_solve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_solve OBJECT ${UMFPACK_SOLVE_SRCS})
  target_include_directories (umfpack_zint_solve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_solve PRIVATE ZINT CONJUGATE_SOLVE)
  target_compile_definitions (umfpack_zint_solve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_solve OBJECT ${UMFPACK_SOLVE_SRCS})
  target_include_directories (umfpack_zlong_solve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_solve PRIVATE ZLONG CONJUGATE_SOLVE)
  target_compile_definitions (umfpack_zlong_solve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_wsolve OBJECT UMFPACK/Source/umfpack_solve.c)
  target_include_directories (umfpack_int_wsolve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_wsolve PRIVATE DINT WSOLVE)
  target_compile_definitions (umfpack_int_wsolve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_wsolve OBJECT UMFPACK/Source/umfpack_solve.c)
  target_include_directories (umfpack_long_wsolve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_wsolve PRIVATE DLONG WSOLVE)
  target_compile_definitions (umfpack_long_wsolve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_wsolve OBJECT UMFPACK/Source/umfpack_solve.c)
  target_include_directories (umfpack_zint_wsolve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_wsolve PRIVATE ZINT WSOLVE)
  target_compile_definitions (umfpack_zint_wsolve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_wsolve OBJECT UMFPACK/Source/umfpack_solve.c)
  target_include_directories (umfpack_zlong_wsolve PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_wsolve PRIVATE ZLONG WSOLVE)
  target_compile_definitions (umfpack_zlong_wsolve PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_fixq OBJECT UMFPACK/Source/umf_assemble.c)
  target_include_directories (umfpack_int_fixq PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_fixq PRIVATE DINT FIXQ)
  target_compile_definitions (umfpack_int_fixq PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_fixq OBJECT UMFPACK/Source/umf_assemble.c)
  target_include_directories (umfpack_long_fixq PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_fixq PRIVATE DLONG FIXQ)
  target_compile_definitions (umfpack_long_fixq PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_fixq OBJECT UMFPACK/Source/umf_assemble.c)
  target_include_directories (umfpack_zint_fixq PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_fixq PRIVATE ZINT FIXQ)
  target_compile_definitions (umfpack_zint_fixq PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_fixq OBJECT UMFPACK/Source/umf_assemble.c)
  target_include_directories (umfpack_zlong_fixq PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_fixq PRIVATE ZLONG FIXQ)
  target_compile_definitions (umfpack_zlong_fixq PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_drop OBJECT UMFPACK/Source/umf_store_lu.c)
  target_include_directories (umfpack_int_drop PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_drop PRIVATE DINT DROP)
  target_compile_definitions (umfpack_int_drop PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_drop OBJECT UMFPACK/Source/umf_store_lu.c)
  target_include_directories (umfpack_long_drop PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_drop PRIVATE DLONG DROP)
  target_compile_definitions (umfpack_long_drop PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_drop OBJECT UMFPACK/Source/umf_store_lu.c)
  target_include_directories (umfpack_zint_drop PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_drop PRIVATE ZINT DROP)
  target_compile_definitions (umfpack_zint_drop PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_drop OBJECT UMFPACK/Source/umf_store_lu.c)
  target_include_directories (umfpack_zlong_drop PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_drop PRIVATE ZLONG DROP)
  target_compile_definitions (umfpack_zlong_drop PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_map_nox OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_int_map_nox PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_map_nox PRIVATE DINT DO_MAP)
  target_compile_definitions (umfpack_int_map_nox PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_map_nox OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_long_map_nox PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_map_nox PRIVATE DLONG DO_MAP)
  target_compile_definitions (umfpack_long_map_nox PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_map_nox OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zint_map_nox PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_map_nox PRIVATE ZINT DO_MAP)
  target_compile_definitions (umfpack_zint_map_nox PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_map_nox OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zlong_map_nox PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_map_nox PRIVATE ZLONG DO_MAP)
  target_compile_definitions (umfpack_zlong_map_nox PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_nomap_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_int_nomap_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_nomap_x PRIVATE DINT DO_VALUES)
  target_compile_definitions (umfpack_int_nomap_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_nomap_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_long_nomap_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_nomap_x PRIVATE DLONG DO_VALUES)
  target_compile_definitions (umfpack_long_nomap_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_nomap_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zint_nomap_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_nomap_x PRIVATE ZINT DO_VALUES)
  target_compile_definitions (umfpack_zint_nomap_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_nomap_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zlong_nomap_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_nomap_x PRIVATE ZLONG DO_VALUES)
  target_compile_definitions (umfpack_zlong_nomap_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_int_map_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_int_map_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_int_map_x PRIVATE DINT DO_MAP DO_VALUES)
  target_compile_definitions (umfpack_int_map_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_long_map_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_long_map_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_long_map_x PRIVATE DLONG DO_MAP DO_VALUES)
  target_compile_definitions (umfpack_long_map_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zint_map_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zint_map_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zint_map_x PRIVATE ZINT DO_MAP DO_VALUES)
  target_compile_definitions (umfpack_zint_map_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack_zlong_map_x OBJECT UMFPACK/Source/umf_triplet.c)
  target_include_directories (umfpack_zlong_map_x PRIVATE $<TARGET_PROPERTY:umfpack,INCLUDE_DIRECTORIES>)
  target_compile_definitions (umfpack_zlong_map_x PRIVATE ZLONG DO_MAP DO_VALUES)
  target_compile_definitions (umfpack_zlong_map_x PRIVATE umfpack_EXPORTS)

  add_library (umfpack
    $<TARGET_OBJECTS:amd_umfpack_int>
    $<TARGET_OBJECTS:amd_umfpack_long>
    $<TARGET_OBJECTS:umfpack_int>
    $<TARGET_OBJECTS:umfpack_int_drop>
    $<TARGET_OBJECTS:umfpack_int_fixq>
    $<TARGET_OBJECTS:umfpack_int_map_nox>
    $<TARGET_OBJECTS:umfpack_int_map_x>
    $<TARGET_OBJECTS:umfpack_int_nomap_x>
    $<TARGET_OBJECTS:umfpack_int_solve>
    $<TARGET_OBJECTS:umfpack_int_wsolve>
    $<TARGET_OBJECTS:umfpack_long>
    $<TARGET_OBJECTS:umfpack_long_drop>
    $<TARGET_OBJECTS:umfpack_long_fixq>
    $<TARGET_OBJECTS:umfpack_long_map_nox>
    $<TARGET_OBJECTS:umfpack_long_map_x>
    $<TARGET_OBJECTS:umfpack_long_nomap_x>
    $<TARGET_OBJECTS:umfpack_long_solve>
    $<TARGET_OBJECTS:umfpack_long_wsolve>
    $<TARGET_OBJECTS:umfpack_zint>
    $<TARGET_OBJECTS:umfpack_zint_drop>
    $<TARGET_OBJECTS:umfpack_zint_fixq>
    $<TARGET_OBJECTS:umfpack_zint_map_nox>
    $<TARGET_OBJECTS:umfpack_zint_map_x>
    $<TARGET_OBJECTS:umfpack_zint_nomap_x>
    $<TARGET_OBJECTS:umfpack_zint_solve>
    $<TARGET_OBJECTS:umfpack_zint_wsolve>
    $<TARGET_OBJECTS:umfpack_zlong>
    $<TARGET_OBJECTS:umfpack_zlong_drop>
    $<TARGET_OBJECTS:umfpack_zlong_fixq>
    $<TARGET_OBJECTS:umfpack_zlong_map_nox>
    $<TARGET_OBJECTS:umfpack_zlong_map_x>
    $<TARGET_OBJECTS:umfpack_zlong_nomap_x>
    $<TARGET_OBJECTS:umfpack_zlong_solve>
    $<TARGET_OBJECTS:umfpack_zlong_wsolve>
    UMFPACK/Source/umfpack_version.c
    ${UMFPACK_HDRS}
    UMFPACK/Source/umfpack_tictoc.c
    UMFPACK/Source/umfpack_timer.c
  )
  add_library (SuiteSparse::UMFPACK ALIAS umfpack)

  target_link_libraries (umfpack PRIVATE
    BLAS::BLAS
    LAPACK::LAPACK
  )

  if (TARGET cholmod)
    # NOTE cholmod already depends on LAPACK
    target_link_libraries (umfpack PRIVATE SuiteSparse::CHOLMOD)
  else (TARGET cholmod)
    set_property (SOURCE
      UMFPACK/Source/umf_cholmod.c
      APPEND PROPERTY COMPILE_DEFINITIONS NCHOLMOD
    )
  endif (TARGET cholmod)
endif (WITH_LICENSE STREQUAL "GPL")

set (KLU_NO_COMPLEX_SRCS
  KLU/Source/klu_defaults.c
  KLU/Source/klu_analyze.c
  KLU/Source/klu_analyze_given.c
  KLU/Source/klu_free_symbolic.c
  KLU/Source/klu_memory.c
)

if (WITH_LICENSE MATCHES "GPL")
  add_library (klu_int OBJECT ${KLU_SRCS} ${KLU_NO_COMPLEX_SRCS})
  target_include_directories (klu_int PRIVATE $<TARGET_PROPERTY:klu,INCLUDE_DIRECTORIES>)
  target_compile_definitions (klu_int PRIVATE DINT)
  target_compile_definitions (klu_int PRIVATE klu_EXPORTS)

  add_library (klu_long OBJECT ${KLU_SRCS} ${KLU_NO_COMPLEX_SRCS})
  target_include_directories (klu_long PRIVATE $<TARGET_PROPERTY:klu,INCLUDE_DIRECTORIES>)
  target_compile_definitions (klu_long PRIVATE DLONG)
  target_compile_definitions (klu_long PRIVATE klu_EXPORTS)

  add_library (klu_zint OBJECT ${KLU_SRCS})
  target_include_directories (klu_zint PRIVATE $<TARGET_PROPERTY:klu,INCLUDE_DIRECTORIES>)
  target_compile_definitions (klu_zint PRIVATE DINT COMPLEX)
  target_compile_definitions (klu_zint PRIVATE klu_EXPORTS)

  add_library (klu_zlong OBJECT ${KLU_SRCS})
  target_include_directories (klu_zlong PRIVATE $<TARGET_PROPERTY:klu,INCLUDE_DIRECTORIES>)
  target_compile_definitions (klu_zlong PRIVATE DLONG COMPLEX)
  target_compile_definitions (klu_zlong PRIVATE klu_EXPORTS)

  add_library (klu
    $<TARGET_OBJECTS:klu_int>
    $<TARGET_OBJECTS:klu_long>
    $<TARGET_OBJECTS:klu_zint>
    $<TARGET_OBJECTS:klu_zlong>
    KLU/Source/klu_version.c
    ${KLU_HDRS}
  )
  add_library (SuiteSparse::KLU ALIAS klu)

  target_link_libraries (klu PRIVATE
    SuiteSparse::AMD
    SuiteSparse::BTF
    SuiteSparse::COLAMD
  )

  if (TARGET cholmod)
    target_link_libraries (klu PRIVATE
      SuiteSparse::CHOLMOD
    )

    add_library (klu_cholmod
      KLU/User/klu_cholmod.c
      KLU/User/klu_l_cholmod.c
      ${KLU_CHOLMOD_HDRS}
    )
    add_library (SuiteSparse::KLU_CHOLMOD ALIAS klu_cholmod)
    set_property (TARGET klu_cholmod PROPERTY NAME_ALIAS klu_cholmod)
    set_property (TARGET klu_cholmod PROPERTY HEADERS_VARIABLE
      KLU_CHOLMOD_HDRS)
    set_property (TARGET klu_cholmod PROPERTY INSTALL_DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/KLU)

    target_include_directories (klu_cholmod PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/User>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_KLU_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/KLU>
    )
    target_link_libraries (klu_cholmod PRIVATE
      SuiteSparse::CHOLMOD
      SuiteSparse::KLU
    )
  endif (TARGET cholmod)
endif (WITH_LICENSE MATCHES "GPL")

if (WITH_LICENSE MATCHES "GPL")
  add_library (ldl_int OBJECT ${LDL_SRCS})
  target_include_directories (ldl_int PRIVATE $<TARGET_PROPERTY:ldl,INCLUDE_DIRECTORIES>)
  target_compile_definitions (ldl_int PRIVATE DINT)
  target_compile_definitions (ldl_int PRIVATE ldl_EXPORTS)

  add_library (ldl_long OBJECT ${LDL_SRCS})
  target_include_directories (ldl_long PRIVATE $<TARGET_PROPERTY:ldl,INCLUDE_DIRECTORIES>)
  target_compile_definitions (ldl_long PRIVATE LDL_LONG)
  target_compile_definitions (ldl_long PRIVATE ldl_EXPORTS)

  add_library (ldl
    $<TARGET_OBJECTS:ldl_int>
    $<TARGET_OBJECTS:ldl_long>
    ${LDL_HDRS}
  )
  add_library (SuiteSparse::LDL ALIAS ldl)

endif (WITH_LICENSE MATCHES "GPL")

set (SuiteSparse_AMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/AMD)
set (SuiteSparse_BTF_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/BTF)
set (SuiteSparse_CAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CAMD)
set (SuiteSparse_CCOLAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CCOLAMD)
set (SuiteSparse_CHOLMOD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CHOLMOD)
set (SuiteSparse_COLAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/COLAMD)
set (SuiteSparse_KLU_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/KLU)
set (SuiteSparse_LDL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/LDL)
set (SuiteSparse_SPQR_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/SPQR)
set (SuiteSparse_UMFPACK_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/UMFPACK)
set (SuiteSparse_GPURuntime_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/GPURuntime)
set (SuiteSparse_GPUQREngine_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/GPUQREngine)

set (_SuiteSparse_CONFIG_INCLUDE_DIRS
  "$<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>")

target_include_directories (suitesparseconfig PUBLIC
  ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
)

if (HAVE_BLAS_UNDERSCORE)
  add_compile_definitions (BLAS_UNDERSCORE)
elseif (HAVE_BLAS_NO_UNDERSCORE)
  add_compile_definitions (BLAS_NO_UNDERSCORE)
endif (HAVE_BLAS_UNDERSCORE)

if (HAVE_M)
  target_link_libraries (suitesparseconfig PRIVATE m)
endif (HAVE_M)

target_include_directories (amd PUBLIC
  ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>)
target_include_directories (camd PUBLIC
  ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:${SuiteSparse_CAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CAMD_INCLUDE_DIR}>)
target_include_directories (ccolamd PUBLIC
  ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:${SuiteSparse_CCOLAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CCOLAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CCOLAMD_INCLUDE_DIR}>)
target_include_directories (colamd PUBLIC
  ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>)

if (WITH_LICENSE MATCHES "GPL")
  target_include_directories (btf PUBLIC
    ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>)
endif (WITH_LICENSE MATCHES "GPL")

if (TARGET cholmod)
  target_include_directories (cholmod PUBLIC
    $<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_CHOLMOD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Cholesky>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/MatrixOps>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Modify>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Supernodal>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Utility>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CHOLMOD_INCLUDE_DIR}>
  )
endif (TARGET cholmod)

if (WITH_LICENSE STREQUAL "GPL")
  target_include_directories (umfpack PUBLIC
    ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_UMFPACK_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Source>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_UMFPACK_INCLUDE_DIR}>
  )
endif (WITH_LICENSE STREQUAL "GPL")

if (WITH_LICENSE MATCHES "GPL")
  target_include_directories (klu PUBLIC
    $<INSTALL_INTERFACE:${SuiteSparse_KLU_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_KLU_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>)
  target_include_directories (ldl PUBLIC
    ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${SuiteSparse_LDL_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LDL/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_LDL_INCLUDE_DIR}>)
endif (WITH_LICENSE MATCHES "GPL")

target_link_libraries (amd PRIVATE SuiteSparse::Config)
target_link_libraries (camd PRIVATE SuiteSparse::Config)
target_link_libraries (ccolamd PRIVATE SuiteSparse::Config)
target_link_libraries (colamd PRIVATE SuiteSparse::Config)

if (WITH_LICENSE MATCHES "GPL")
  target_link_libraries (klu PRIVATE SuiteSparse::Config)
endif (WITH_LICENSE MATCHES "GPL")

if (WITH_LICENSE STREQUAL "GPL")
  target_link_libraries (umfpack PRIVATE
    SuiteSparse::Config
    amd
  )
endif (WITH_LICENSE STREQUAL "GPL")

if (TARGET spqr)
  target_include_directories (spqr PUBLIC
    ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${SuiteSparse_SPQR_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SPQR_INCLUDE_DIR}>
  )

  target_link_libraries (spqr
    PUBLIC
      SuiteSparse::CHOLMOD
    PRIVATE
      BLAS::BLAS
      LAPACK::LAPACK
      SuiteSparse::Config
  )
endif (TARGET spqr)

# CUDA support
if (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)
  add_library (gpuruntime
    SPQR/GPURuntime/Source/SuiteSparseGPU_Workspace.cpp
    SPQR/GPURuntime/Source/SuiteSparseGPU_Workspace_cpuAllocators.cpp
    SPQR/GPURuntime/Source/SuiteSparseGPU_Workspace_gpuAllocators.cpp
    SPQR/GPURuntime/Source/SuiteSparseGPU_Workspace_transfer.cpp
  )
  add_library (SuiteSparse::GPURuntime ALIAS gpuruntime)

  target_include_directories (gpuruntime PUBLIC
    ${_SuiteSparse_CONFIG_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/GPURuntime/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_GPURuntime_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_GPURuntime_INCLUDE_DIR}>
  )

  target_link_libraries (gpuruntime
    PRIVATE
      SuiteSparse::Config
      CUDA::cudart
  )

  target_compile_definitions (gpuruntime PUBLIC SPQR_HAS_CUDA)

  # GPUQREngine

  add_library (gpuqrengine
    SPQR/GPUQREngine/Include/Kernel/sharedMemory.hpp
    SPQR/GPUQREngine/Source/BucketList/BucketList.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_AdvanceBundles.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_CreateBundles.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_FillWorkQueue.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_GrowBundles.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_Manage.cpp
    SPQR/GPUQREngine/Source/BucketList/BucketList_PostProcessing.cpp
    SPQR/GPUQREngine/Source/GPUQREngine_ExpertDense.cpp
    SPQR/GPUQREngine/Source/GPUQREngine_ExpertSparse.cpp
    SPQR/GPUQREngine/Source/GPUQREngine_GraphVizHelper.cpp
    SPQR/GPUQREngine/Source/GPUQREngine_Internal.cpp
    SPQR/GPUQREngine/Source/GPUQREngine_UberKernel.cu
    SPQR/GPUQREngine/Source/LLBundle/LLBundle.cpp
    SPQR/GPUQREngine/Source/LLBundle/LLBundle_AddTiles.cpp
    SPQR/GPUQREngine/Source/LLBundle/LLBundle_Advance.cpp
    SPQR/GPUQREngine/Source/LLBundle/LLBundle_GPUPack.cpp
    SPQR/GPUQREngine/Source/LLBundle/LLBundle_PipelinedRearrange.cpp
    SPQR/GPUQREngine/Source/LLBundle/LLBundle_UpdateSecondMinIndex.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_FillWorkQueue.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_Front.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_LaunchKernel.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_PostProcess.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_Render.cpp
    SPQR/GPUQREngine/Source/Scheduler/Scheduler_TransferData.cpp
    SPQR/GPUQREngine/Source/ssgpu_maxQueueSize.cpp
    SPQR/GPUQREngine/Source/TaskDescriptor/TaskDescriptor_flops.cpp
  )

  target_include_directories (gpuqrengine PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/GPUQREngine/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_GPUQREngine_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_GPUQREngine_INCLUDE_DIR}>
  )

  target_include_directories (gpuqrengine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/SPQR/GPUQREngine/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/SPQR/GPURuntime/Include
  )

  target_compile_definitions (gpuqrengine PUBLIC SPQR_HAS_CUDA)

  target_link_libraries (gpuqrengine
    PRIVATE
      CUDA::cublas
      gpuruntime
  )

  if (WITH_CHOLMOD)
    set (CHOLMOD_GPU_SRCS CHOLMOD/GPU/cholmod_gpu.c)

    add_library (cholmod_gpu_int OBJECT ${CHOLMOD_GPU_SRCS})
    target_include_directories (cholmod_gpu_int PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_gpu_int PRIVATE
      CHOLMOD_HAS_CUDA
      CHOLMOD_INT32)
    target_compile_definitions (cholmod_gpu_int PRIVATE cholmod_EXPORTS)

    add_library (cholmod_gpu_long OBJECT ${CHOLMOD_GPU_SRCS})
    target_include_directories (cholmod_gpu_long PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_gpu_long PRIVATE
      CHOLMOD_HAS_CUDA
      CHOLMOD_INT64)
    target_compile_definitions (cholmod_gpu_long PRIVATE cholmod_EXPORTS)

    target_sources (cholmod PRIVATE
      $<TARGET_OBJECTS:cholmod_gpu_int>
      $<TARGET_OBJECTS:cholmod_gpu_long>
      CHOLMOD/GPU/cholmod_gpu_kernels.cu
    )

    target_compile_definitions (cholmod PUBLIC
      CHOLMOD_HAS_CUDA
      GPU_BLAS)
    target_include_directories (cholmod PUBLIC
      $<TARGET_PROPERTY:CUDA::cublas,INTERFACE_INCLUDE_DIRECTORIES>
    )
    target_link_libraries (cholmod
      PRIVATE
        CUDA::cudart
        gpuruntime
        CUDA::cublas
    )
  endif (WITH_CHOLMOD)

  # SPQR GPU
  if (TARGET spqr)
    target_compile_definitions (spqr PUBLIC SPQR_HAS_CUDA)
    target_sources (spqr PRIVATE
      SPQR/SPQRGPU/spqrgpu_buildAssemblyMaps.cpp
      SPQR/SPQRGPU/spqrgpu_computeFrontStaging.cpp
      SPQR/SPQRGPU/spqrgpu_kernel.cpp
    )

    target_link_libraries (spqr PRIVATE gpuqrengine)
  endif (TARGET spqr)
endif (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)

set (_SuiteSparse_TARGETS amd camd ccolamd colamd suitesparseconfig)

if (WITH_LICENSE MATCHES "GPL")
  list (APPEND _SuiteSparse_TARGETS btf klu ldl)
  if (TARGET klu_cholmod)
    list (APPEND _SuiteSparse_TARGETS klu_cholmod)
  endif (TARGET klu_cholmod)
endif (WITH_LICENSE MATCHES "GPL")

if (WITH_LICENSE STREQUAL "GPL")
  list (APPEND _SuiteSparse_TARGETS umfpack)
endif (WITH_LICENSE STREQUAL "GPL")

set_property (TARGET suitesparseconfig PROPERTY NAME_ALIAS SuiteSparse)
set_property (TARGET suitesparseconfig PROPERTY HEADERS_VARIABLE
  SUITESPARSECONFIG_HDRS)
set_property (TARGET suitesparseconfig PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR})

if (TARGET cholmod)
  list (APPEND _SuiteSparse_TARGETS cholmod)

  if (TARGET cholmod_partition)
    list (APPEND _SuiteSparse_TARGETS cholmod_partition)
  endif (TARGET cholmod_partition)
endif (TARGET cholmod)

if (TARGET spqr)
  list (APPEND _SuiteSparse_TARGETS spqr)
  set_property (TARGET spqr PROPERTY NAME_ALIAS SuiteSparseQR)
  set_property (TARGET spqr PROPERTY INSTALL_DESTINATION
    ${SuiteSparse_SPQR_INCLUDE_DIR})
endif (TARGET spqr)

if (TARGET gpuruntime)
  list (APPEND _SuiteSparse_TARGETS gpuruntime)
  set_property (TARGET gpuruntime PROPERTY NAME_ALIAS GPURuntime)
endif (TARGET gpuruntime)

if (TARGET gpuqrengine)
  list (APPEND _SuiteSparse_TARGETS gpuqrengine)
  set_property (TARGET gpuqrengine PROPERTY NAME_ALIAS GPUQREngine)
endif (TARGET gpuqrengine)

set_property (TARGET amd PROPERTY SPDX_LICENSE BSD-3-Clause)
set_property (TARGET camd PROPERTY SPDX_LICENSE BSD-3-Clause)
set_property (TARGET ccolamd PROPERTY SPDX_LICENSE BSD-3-Clause)
set_property (TARGET colamd PROPERTY SPDX_LICENSE BSD-3-Clause)
set_property (TARGET suitesparseconfig PROPERTY SPDX_LICENSE NOASSERTION)

if (WITH_LICENSE MATCHES "GPL")
  set_property (TARGET btf PROPERTY SPDX_LICENSE LGPL-2.1-or-later)
  set_property (TARGET klu PROPERTY SPDX_LICENSE LGPL-2.1-or-later)
  set_property (TARGET ldl PROPERTY SPDX_LICENSE LGPL-2.1-or-later)
  if (TARGET klu_cholmod)
    set_property (TARGET klu_cholmod PROPERTY SPDX_LICENSE LGPL-2.1-or-later)
  endif (TARGET klu_cholmod)
endif (WITH_LICENSE MATCHES "GPL")

if (WITH_LICENSE STREQUAL "GPL")
  set_property (TARGET umfpack PROPERTY SPDX_LICENSE GPL-2.0-or-later)
endif (WITH_LICENSE STREQUAL "GPL")

if (TARGET cholmod)
  if (WITH_LICENSE STREQUAL "GPL")
    set (_SuiteSparse_CHOLMOD_SPDX_LICENSE
      "LGPL-2.1-or-later AND GPL-2.0-or-later")
  else (WITH_LICENSE STREQUAL "GPL")
    set (_SuiteSparse_CHOLMOD_SPDX_LICENSE LGPL-2.1-or-later)
  endif (WITH_LICENSE STREQUAL "GPL")

  set_property (TARGET cholmod PROPERTY SPDX_LICENSE
    "${_SuiteSparse_CHOLMOD_SPDX_LICENSE}")

  if (TARGET cholmod_partition)
    set_property (TARGET cholmod_partition PROPERTY SPDX_LICENSE
      "${_SuiteSparse_CHOLMOD_SPDX_LICENSE}")
  endif (TARGET cholmod_partition)
endif (TARGET cholmod)

if (TARGET spqr)
  set_property (TARGET spqr PROPERTY SPDX_LICENSE GPL-2.0-or-later)
endif (TARGET spqr)

if (TARGET gpuruntime)
  set_property (TARGET gpuruntime PROPERTY SPDX_LICENSE GPL-2.0-or-later)
endif (TARGET gpuruntime)

if (TARGET gpuqrengine)
  set_property (TARGET gpuqrengine PROPERTY SPDX_LICENSE GPL-2.0-or-later)
endif (TARGET gpuqrengine)

foreach (_target ${_SuiteSparse_TARGETS})
  if (_target STREQUAL cholmod_partition)
    # cholmod_partition is meta library that does provide its own headers.
    continue ()
  endif (_target STREQUAL cholmod_partition)

  get_property (_alias TARGET ${_target} PROPERTY NAME_ALIAS)
  string (TOUPPER "${_target}" _upper_target)

  if ("${_alias}" STREQUAL "")
    set (_alias ${_upper_target})
    set (_export_prefix ${_target})
  else ("${_alias}" STREQUAL "")
    set (_export_prefix ${_alias})
  endif ("${_alias}" STREQUAL "")

  get_property (_target_include_dir TARGET ${_target} PROPERTY
    INSTALL_DESTINATION)

  if ("${_target_include_dir}" STREQUAL "")
    set (_target_include_dir ${SuiteSparse_${_alias}_INCLUDE_DIR})
  endif ("${_target_include_dir}" STREQUAL "")

  set (_export_header
    ${CMAKE_CURRENT_BINARY_DIR}/${_target_include_dir}/${_export_prefix}_export.h)

  # Generate an export header with the name <target>_export.h for each mandatory
  # target.
  generate_export_header (${_target} EXPORT_FILE_NAME ${_export_header})

  # Add the export header to the target
  target_sources (${_target} PRIVATE ${_export_header})

  get_property (_HDRS_VAR TARGET ${_target} PROPERTY HEADERS_VARIABLE)

  if ("${_HDRS_VAR}" STREQUAL "")
    set (_HDRS_VAR ${_upper_target}_HDRS)
  endif ("${_HDRS_VAR}" STREQUAL "")

  if (NOT DEFINED ${_HDRS_VAR})
    message (FATAL_ERROR "CMake variable ${_HDRS_VAR} is not defined")
  endif (NOT DEFINED ${_HDRS_VAR})

  # Install headers
  install (FILES ${${_HDRS_VAR}} ${_export_header} DESTINATION
    ${_target_include_dir} COMPONENT Development)
endforeach (_target)

set (SuiteSparse_CMake_LIBDIR ${CMAKE_INSTALL_LIBDIR}/SuiteSparse)
# suitesparse-config.cmake relative install directory
set (SuiteSparse_CMake_INSTALLDIR ${SuiteSparse_CMake_LIBDIR}/cmake)

target_sources (suitesparseconfig PRIVATE ${_output_CMake_MODULES})

set (SuiteSparse_NAMESPACE SuiteSparse::)

configure_package_config_file (suitesparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  INSTALL_DESTINATION ${SuiteSparse_CMake_INSTALLDIR})
write_basic_package_version_file (suitesparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY AnyNewerVersion)

export (TARGETS ${_SuiteSparse_TARGETS} NAMESPACE ${SuiteSparse_NAMESPACE} FILE
  suitesparse-targets.cmake)
export (PACKAGE SuiteSparse)

install (EXPORT suitesparse-targets NAMESPACE ${SuiteSparse_NAMESPACE}
  DESTINATION ${SuiteSparse_CMake_INSTALLDIR})
install (FILES
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-version.cmake
  DESTINATION ${SuiteSparse_CMake_INSTALLDIR}
  COMPONENT Development)

if (TARGET SuiteSparse::Partition)
  install (FILES
    cmake/FindMETIS.cmake
    DESTINATION ${SuiteSparse_CMake_INSTALLDIR}
    COMPONENT Development)
endif (TARGET SuiteSparse::Partition)

set_component_version (amd 3.3.4)
set_component_version (camd 3.3.5)
set_component_version (ccolamd 3.3.5)
set_component_version (colamd 3.3.5)

if (WITH_LICENSE MATCHES "GPL")
  set_component_version (btf 2.3.3)
  set_component_version (klu 2.3.6)
  if (TARGET klu_cholmod)
    set_component_version (klu_cholmod 2.3.6)
  endif (TARGET klu_cholmod)
  set_component_version (ldl 3.3.3)
endif (WITH_LICENSE MATCHES "GPL")

set_component_version (suitesparseconfig ${SuiteSparse_VERSION})

if (WITH_LICENSE STREQUAL "GPL")
  set_component_version (umfpack 6.3.8)
endif (WITH_LICENSE STREQUAL "GPL")

if (TARGET cholmod)
  set_component_version (cholmod 5.3.5)
endif (TARGET cholmod)

if (TARGET spqr)
  set_component_version (spqr 4.3.6)
endif (TARGET spqr)

if (TARGET gpuruntime)
  set_component_version (gpuruntime 1.0.5)
endif (TARGET gpuruntime)

if (TARGET gpuqrengine)
  set_component_version (gpuqrengine 1.0.5)
endif (TARGET gpuqrengine)

if (HAVE_QSTD_C99)
  target_compile_options (suitesparseconfig PUBLIC $<$<C_COMPILER_ID:Intel>:-Qstd=c99>)
endif (HAVE_QSTD_C99)

set_target_properties (amd PROPERTIES EXPORT_NAME AMD OUTPUT_NAME amd)
set_target_properties (camd PROPERTIES EXPORT_NAME CAMD OUTPUT_NAME camd)
set_target_properties (ccolamd PROPERTIES EXPORT_NAME CCOLAMD OUTPUT_NAME ccolamd)
set_target_properties (colamd PROPERTIES EXPORT_NAME COLAMD OUTPUT_NAME colamd)
set_target_properties (suitesparseconfig PROPERTIES EXPORT_NAME Config)

if (WITH_LICENSE MATCHES "GPL")
  set_target_properties (btf PROPERTIES EXPORT_NAME BTF OUTPUT_NAME btf)
  set_target_properties (klu PROPERTIES EXPORT_NAME KLU OUTPUT_NAME klu)
  if (TARGET klu_cholmod)
    set_target_properties (klu_cholmod PROPERTIES
      EXPORT_NAME klu_cholmod OUTPUT_NAME klu_cholmod)
  endif (TARGET klu_cholmod)
  set_target_properties (ldl PROPERTIES EXPORT_NAME LDL OUTPUT_NAME ldl)
endif (WITH_LICENSE MATCHES "GPL")

if (WITH_LICENSE STREQUAL "GPL")
  set_target_properties (umfpack PROPERTIES EXPORT_NAME UMFPACK OUTPUT_NAME umfpack)
endif (WITH_LICENSE STREQUAL "GPL")

if (TARGET cholmod)
  set_target_properties (cholmod PROPERTIES EXPORT_NAME CHOLMOD OUTPUT_NAME cholmod)
endif (TARGET cholmod)

if (TARGET cholmod_partition)
  set_target_properties (cholmod_partition PROPERTIES EXPORT_NAME Partition)
endif (TARGET cholmod_partition)

if (TARGET spqr)
  set_target_properties (spqr PROPERTIES EXPORT_NAME SPQR OUTPUT_NAME spqr)
endif (TARGET spqr)

if (WITH_DEMOS)
  add_executable (amd_demo
    AMD/Demo/amd_demo.c
  )
  target_link_libraries (amd_demo PRIVATE SuiteSparse::AMD)

  add_executable (amd_demo2
    AMD/Demo/amd_demo2.c
  )
  target_link_libraries (amd_demo2 PRIVATE SuiteSparse::AMD)

  add_executable (amd_l_demo
    AMD/Demo/amd_l_demo.c
  )
  target_link_libraries (amd_l_demo PRIVATE SuiteSparse::AMD)

  add_executable (amd_simple
    AMD/Demo/amd_simple.c
  )
  target_link_libraries (amd_simple PRIVATE SuiteSparse::AMD)

  add_executable (camd_demo
    CAMD/Demo/camd_demo.c
  )
  target_link_libraries (camd_demo PRIVATE SuiteSparse::CAMD)

  add_executable (camd_demo2
    CAMD/Demo/camd_demo2.c
  )
  target_link_libraries (camd_demo2 PRIVATE SuiteSparse::CAMD)

  add_executable (camd_l_demo
    CAMD/Demo/camd_l_demo.c
  )
  target_link_libraries (camd_l_demo PRIVATE SuiteSparse::CAMD)

  add_executable (camd_simple
    CAMD/Demo/camd_simple.c
  )
  target_link_libraries (camd_simple PRIVATE SuiteSparse::CAMD)

  if (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)
    add_library (amd_f77_amd_objects OBJECT
      AMD/Source/amd.f
    )

    add_library (amd_f77_amdbar_objects OBJECT
      AMD/Source/amdbar.f
    )

    add_executable (amd_f77cross
      AMD/Demo/amd_f77wrapper.c
      AMD/Demo/amd_f77cross.f
    )
    target_link_libraries (amd_f77cross PRIVATE SuiteSparse::AMD)

    add_executable (amd_f77demo
      AMD/Demo/amd_f77demo.f
      $<TARGET_OBJECTS:amd_f77_amd_objects>
      $<TARGET_OBJECTS:amd_f77_amdbar_objects>
    )

    add_executable (amd_f77simple
      AMD/Demo/amd_f77simple.f
      $<TARGET_OBJECTS:amd_f77_amd_objects>
    )
  endif (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)

  add_executable (ccolamd_example
    CCOLAMD/Demo/ccolamd_example.c
  )
  target_link_libraries (ccolamd_example PRIVATE SuiteSparse::CCOLAMD)

  add_executable (ccolamd_l_example
    CCOLAMD/Demo/ccolamd_l_example.c
  )
  target_link_libraries (ccolamd_l_example PRIVATE SuiteSparse::CCOLAMD)

  if (WITH_CHOLESKY AND WITH_CHECK)
    add_executable (cholmod_di_demo
      CHOLMOD/Demo/cholmod_demo.h
      CHOLMOD/Demo/cholmod_di_demo.c
    )
    target_link_libraries (cholmod_di_demo PRIVATE
      SuiteSparse::AMD
      SuiteSparse::CAMD
      SuiteSparse::CCOLAMD
      SuiteSparse::CHOLMOD
      SuiteSparse::COLAMD
      SuiteSparse::Config
    )

    add_executable (cholmod_dl_demo CHOLMOD/Demo/cholmod_dl_demo.c)
    target_link_libraries (cholmod_dl_demo PRIVATE
      SuiteSparse::AMD
      SuiteSparse::CAMD
      SuiteSparse::CCOLAMD
      SuiteSparse::CHOLMOD
      SuiteSparse::COLAMD
      SuiteSparse::Config
    )

    add_executable (cholmod_si_demo CHOLMOD/Demo/cholmod_si_demo.c)
    target_link_libraries (cholmod_si_demo PRIVATE
      SuiteSparse::AMD
      SuiteSparse::CAMD
      SuiteSparse::CCOLAMD
      SuiteSparse::CHOLMOD
      SuiteSparse::COLAMD
      SuiteSparse::Config
    )

    add_executable (cholmod_sl_demo CHOLMOD/Demo/cholmod_sl_demo.c)
    target_link_libraries (cholmod_sl_demo PRIVATE
      SuiteSparse::AMD
      SuiteSparse::CAMD
      SuiteSparse::CCOLAMD
      SuiteSparse::CHOLMOD
      SuiteSparse::COLAMD
      SuiteSparse::Config
    )

    add_executable (cholmod_di_simple CHOLMOD/Demo/cholmod_di_simple.c)
    target_link_libraries (cholmod_di_simple PRIVATE
      SuiteSparse::Config
      SuiteSparse::CHOLMOD
    )

    add_executable (cholmod_dl_simple CHOLMOD/Demo/cholmod_dl_simple.c)
    target_link_libraries (cholmod_dl_simple PRIVATE
      SuiteSparse::Config
      SuiteSparse::CHOLMOD
    )

    add_executable (cholmod_si_simple CHOLMOD/Demo/cholmod_si_simple.c)
    target_link_libraries (cholmod_si_simple PRIVATE
      SuiteSparse::Config
      SuiteSparse::CHOLMOD
    )

    add_executable (cholmod_sl_simple CHOLMOD/Demo/cholmod_sl_simple.c)
    target_link_libraries (cholmod_sl_simple PRIVATE
      SuiteSparse::Config
      SuiteSparse::CHOLMOD
    )
  endif (WITH_CHOLESKY AND WITH_CHECK)

  if (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)
    add_executable (cholmod_reade CHOLMOD/Demo/reade.f)
    add_executable (cholmod_readhb CHOLMOD/Demo/readhb.f)
    add_executable (cholmod_readhb2 CHOLMOD/Demo/readhb2.f)
    set_property (TARGET cholmod_reade PROPERTY LINK_WHAT_YOU_USE OFF)
    set_property (TARGET cholmod_readhb PROPERTY LINK_WHAT_YOU_USE OFF)
    set_property (TARGET cholmod_readhb2 PROPERTY LINK_WHAT_YOU_USE OFF)
  endif (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)

  add_executable (colamd_example
    COLAMD/Demo/colamd_example.c
  )
  target_link_libraries (colamd_example PRIVATE SuiteSparse::COLAMD)

  add_executable (colamd_l_example
    COLAMD/Demo/colamd_l_example.c
  )
  target_link_libraries (colamd_l_example PRIVATE SuiteSparse::COLAMD)

  if (TARGET klu)
    add_executable (klu_simple
      KLU/Demo/klu_simple.c
    )
    target_link_libraries (klu_simple PRIVATE SuiteSparse::KLU)
  endif (TARGET klu)

  if (TARGET cholmod)
    add_executable (kludemo
      KLU/Demo/kludemo.c
    )
    target_link_libraries (kludemo PRIVATE
      SuiteSparse::CHOLMOD
      SuiteSparse::KLU
      SuiteSparse::KLU_CHOLMOD
    )

    add_executable (kluldemo
      KLU/Demo/kluldemo.c
    )
    target_link_libraries (kluldemo PRIVATE
      SuiteSparse::CHOLMOD
      SuiteSparse::KLU
      SuiteSparse::KLU_CHOLMOD
    )

    if (HAVE_M)
      target_link_libraries (kludemo PRIVATE m)
      target_link_libraries (kluldemo PRIVATE m)
    endif (HAVE_M)

    install (TARGETS
      kludemo
      kluldemo
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
    )
  endif (TARGET cholmod)

  if (TARGET ldl)
    add_executable (ldlamd
      LDL/Demo/ldlmain.c
    )
    target_compile_definitions (ldlamd PRIVATE USE_AMD)
    target_link_libraries (ldlamd PRIVATE
      SuiteSparse::AMD
      SuiteSparse::LDL
    )

    add_executable (ldlmain
      LDL/Demo/ldlmain.c
    )
    target_link_libraries (ldlmain PRIVATE SuiteSparse::LDL)

    add_executable (ldllmain
      LDL/Demo/ldlmain.c
    )
    target_compile_definitions (ldllmain PRIVATE LDL_LONG)
    target_link_libraries (ldllmain PRIVATE SuiteSparse::LDL)

    add_executable (ldllamd
      LDL/Demo/ldlmain.c
    )
    target_compile_definitions (ldllamd PRIVATE LDL_LONG USE_AMD)
    target_link_libraries (ldllamd PRIVATE
      SuiteSparse::AMD
      SuiteSparse::LDL
    )

    add_executable (ldlsimple
      LDL/Demo/ldlsimple.c
    )
    target_link_libraries (ldlsimple PRIVATE
      SuiteSparse::Config
      SuiteSparse::LDL
    )

    add_executable (ldllsimple
      LDL/Demo/ldlsimple.c
    )
    target_compile_definitions (ldllsimple PRIVATE LDL_LONG)
    target_link_libraries (ldllsimple PRIVATE
      SuiteSparse::Config
      SuiteSparse::LDL
    )
  endif (TARGET ldl)

  if (TARGET umfpack)
    add_executable (umf4 UMFPACK/Demo/umf4.c)
    target_link_libraries (umf4 PRIVATE
      SuiteSparse::AMD
      SuiteSparse::UMFPACK
    )

    add_executable (umfpack_di_demo UMFPACK/Demo/umfpack_di_demo.c)
    target_link_libraries (umfpack_di_demo PRIVATE SuiteSparse::UMFPACK)

    add_executable (umfpack_dl_demo UMFPACK/Demo/umfpack_dl_demo.c)
    target_link_libraries (umfpack_dl_demo PRIVATE SuiteSparse::UMFPACK)

    add_executable (umfpack_zi_demo UMFPACK/Demo/umfpack_zi_demo.c)
    target_link_libraries (umfpack_zi_demo PRIVATE SuiteSparse::UMFPACK)

    add_executable (umfpack_zl_demo UMFPACK/Demo/umfpack_zl_demo.c)
    target_link_libraries (umfpack_zl_demo PRIVATE SuiteSparse::UMFPACK)

    add_executable (umfpack_simple UMFPACK/Demo/umfpack_simple.c)
    target_link_libraries (umfpack_simple PRIVATE
      SuiteSparse::Config
      SuiteSparse::UMFPACK
    )
  endif (TARGET umfpack)

  set (SuiteSparse_BIN_TARGETS
    amd_demo
    amd_demo2
    amd_l_demo
    ccolamd_example
    ccolamd_l_example
    camd_demo
    camd_demo2
    camd_l_demo
    camd_simple
    colamd_example
    colamd_l_example
  )

  if (TARGET amd_f77cross)
    list (APPEND SuiteSparse_BIN_TARGETS
      amd_f77cross
      amd_f77demo
      amd_f77simple
    )
  endif (TARGET amd_f77cross)

  if (TARGET cholmod_di_simple)
    list (APPEND SuiteSparse_BIN_TARGETS
      cholmod_di_simple
      cholmod_dl_simple
      cholmod_si_simple
      cholmod_sl_simple
    )
  endif (TARGET cholmod_di_simple)

  if (TARGET cholmod_reade)
    list (APPEND SuiteSparse_BIN_TARGETS
      cholmod_reade
      cholmod_readhb
      cholmod_readhb2
    )
  endif (TARGET cholmod_reade)

  if (TARGET ldl)
    list (APPEND SuiteSparse_BIN_TARGETS
      ldlamd
      ldlmain
      ldlsimple
    )
  endif (TARGET ldl)

  if (TARGET umfpack)
    list (APPEND SuiteSparse_BIN_TARGETS
      umf4
      umfpack_di_demo
      umfpack_dl_demo
      umfpack_simple
      umfpack_zi_demo
      umfpack_zl_demo
    )
  endif (TARGET umfpack)

  if (TARGET umfpack)
    if (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)
      add_library (umf4_f77wrapper64_objects OBJECT
        UMFPACK/Demo/umf4_f77wrapper.c
      )
      target_compile_definitions (umf4_f77wrapper64_objects PRIVATE DLONG)
      target_link_libraries (umf4_f77wrapper64_objects PRIVATE
        SuiteSparse::UMFPACK
      )

      add_library (umf4_f77zwrapper64_objects OBJECT
        UMFPACK/Demo/umf4_f77zwrapper.c
      )
      target_compile_definitions (umf4_f77zwrapper64_objects PRIVATE ZLONG)
      target_link_libraries (umf4_f77zwrapper64_objects PRIVATE
        SuiteSparse::UMFPACK
      )

      add_executable (umf4hb
        UMFPACK/Demo/umf4_f77wrapper.c
        UMFPACK/Demo/umf4hb.f
      )
      target_link_libraries (umf4hb PRIVATE SuiteSparse::UMFPACK)

      add_executable (umfpack_readhb UMFPACK/Demo/readhb.f)
      add_executable (umfpack_readhb_nozeros
        UMFPACK/Demo/readhb_nozeros.f)
      add_executable (umfpack_readhb_size UMFPACK/Demo/readhb_size.f)
      set_property (TARGET umfpack_readhb PROPERTY LINK_WHAT_YOU_USE OFF)
      set_property (TARGET umfpack_readhb_nozeros
        PROPERTY LINK_WHAT_YOU_USE OFF)
      set_property (TARGET umfpack_readhb_size PROPERTY LINK_WHAT_YOU_USE OFF)

      add_executable (umf4hb64
        UMFPACK/Demo/umf4hb64.f
        $<TARGET_OBJECTS:umf4_f77wrapper64_objects>
      )
      target_link_libraries (umf4hb64 PRIVATE SuiteSparse::UMFPACK)

      add_executable (umf4zhb
        UMFPACK/Demo/umf4_f77zwrapper.c
        UMFPACK/Demo/umf4zhb.f
      )
      target_link_libraries (umf4zhb PRIVATE SuiteSparse::UMFPACK)

      add_executable (umf4zhb64
        UMFPACK/Demo/umf4zhb64.f
        $<TARGET_OBJECTS:umf4_f77zwrapper64_objects>
      )
      target_link_libraries (umf4zhb64 PRIVATE SuiteSparse::UMFPACK)

      list (APPEND SuiteSparse_BIN_TARGETS
        umf4hb
        umf4hb64
        umf4zhb
        umf4zhb64
        umfpack_readhb
        umfpack_readhb_nozeros
        umfpack_readhb_size
      )
    endif (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)
  endif (TARGET umfpack)

  if (TARGET spqr)
    if (WITH_CHECK)
      add_executable (qrdemo SPQR/Demo/qrdemo.cpp)
      add_executable (qrdemo_int32 SPQR/Demo/qrdemo_int32.cpp)
      add_executable (qrdemoc SPQR/Demo/qrdemoc.c)
      add_executable (qrdemoc_int32 SPQR/Demo/qrdemoc_int32.c)
      add_executable (qrsimple SPQR/Demo/qrsimple.cpp)
      add_executable (qrsimplec SPQR/Demo/qrsimplec.c)
      add_executable (qrsimplec_int32 SPQR/Demo/qrsimplec_int32.c)

      target_link_libraries (qrdemo PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrdemo_int32 PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrdemoc PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrdemoc_int32 PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrsimple PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrsimplec PRIVATE SuiteSparse::SPQR)
      target_link_libraries (qrsimplec_int32 PRIVATE SuiteSparse::SPQR)

      list (APPEND SuiteSparse_BIN_TARGETS
        qrdemo
        qrdemo_int32
        qrdemoc
        qrdemoc_int32
        qrsimple
        qrsimplec
        qrsimplec_int32
      )
    endif (WITH_CHECK)

    if (TARGET gpuqrengine)
      add_executable (qrdemo_gpu SPQR/Demo/qrdemo_gpu.cpp)
      add_executable (qrdemo_gpu2 SPQR/Demo/qrdemo_gpu2.cpp)
      add_executable (qrdemo_gpu3 SPQR/Demo/qrdemo_gpu3.cpp)
      add_executable (gpuqrengine_demo
        SPQR/GPUQREngine/Demo/gpuqrengine_demo.cpp)
      set_property (TARGET gpuqrengine_demo PROPERTY LINK_WHAT_YOU_USE OFF)

      target_link_libraries (qrdemo_gpu PRIVATE
        SuiteSparse::Config
        SuiteSparse::CHOLMOD
        SuiteSparse::SPQR
      )
      target_link_libraries (qrdemo_gpu2 PRIVATE
        SuiteSparse::Config
        SuiteSparse::CHOLMOD
        SuiteSparse::SPQR
      )
      target_link_libraries (qrdemo_gpu3 PRIVATE
        SuiteSparse::Config
        SuiteSparse::CHOLMOD
        SuiteSparse::SPQR
      )
      list (APPEND SuiteSparse_BIN_TARGETS
        qrdemo_gpu
        qrdemo_gpu2
        qrdemo_gpu3
        gpuqrengine_demo
      )
    endif (TARGET gpuqrengine)
  endif (TARGET spqr)

  if (TARGET cholmod_di_demo)
    foreach (_SuiteSparse_CHOLMOD_DEMO_TARGET IN ITEMS
      cholmod_di_demo
      cholmod_dl_demo
      cholmod_si_demo
      cholmod_sl_demo
      cholmod_di_simple
      cholmod_dl_simple
      cholmod_si_simple
      cholmod_sl_simple
    )
      target_include_directories (${_SuiteSparse_CHOLMOD_DEMO_TARGET} PRIVATE
        $<TARGET_PROPERTY:amd,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:camd,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:ccolamd,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:colamd,INTERFACE_INCLUDE_DIRECTORIES>
      )
    endforeach (_SuiteSparse_CHOLMOD_DEMO_TARGET)
  endif (TARGET cholmod_di_demo)
endif (WITH_DEMOS)

if (BUILD_TESTING)
  if (TARGET cholmod AND WITH_CHOLESKY AND WITH_CHECK)
    if (TARGET cholmod_di_demo)
      set (_SuiteSparse_CHOLMOD_DI_TEST_TARGET cholmod_di_demo)
      set (_SuiteSparse_CHOLMOD_DL_TEST_TARGET cholmod_dl_demo)
      set (_SuiteSparse_CHOLMOD_SI_TEST_TARGET cholmod_si_demo)
      set (_SuiteSparse_CHOLMOD_SL_TEST_TARGET cholmod_sl_demo)
    else ()
      add_executable (cholmod_di_test CHOLMOD/Demo/cholmod_di_demo.c)
      add_executable (cholmod_dl_test CHOLMOD/Demo/cholmod_dl_demo.c)
      add_executable (cholmod_si_test CHOLMOD/Demo/cholmod_si_demo.c)
      add_executable (cholmod_sl_test CHOLMOD/Demo/cholmod_sl_demo.c)

      foreach (_SuiteSparse_CHOLMOD_TEST_TARGET IN ITEMS
        cholmod_di_test
        cholmod_dl_test
        cholmod_si_test
        cholmod_sl_test
      )
        target_link_libraries (${_SuiteSparse_CHOLMOD_TEST_TARGET} PRIVATE
          SuiteSparse::AMD
          SuiteSparse::CAMD
          SuiteSparse::CCOLAMD
          SuiteSparse::CHOLMOD
          SuiteSparse::COLAMD
          SuiteSparse::Config
        )
      endforeach (_SuiteSparse_CHOLMOD_TEST_TARGET)

      set (_SuiteSparse_CHOLMOD_DI_TEST_TARGET cholmod_di_test)
      set (_SuiteSparse_CHOLMOD_DL_TEST_TARGET cholmod_dl_test)
      set (_SuiteSparse_CHOLMOD_SI_TEST_TARGET cholmod_si_test)
      set (_SuiteSparse_CHOLMOD_SL_TEST_TARGET cholmod_sl_test)
    endif ()

  add_test ( NAME CHOLMOD_int32_double_bcsstk01
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk01.tri)
  add_test ( NAME CHOLMOD_int64_double_bcsstk01
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk01.tri)
  add_test ( NAME CHOLMOD_int32_single_bcsstk01
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk01.tri)
  add_test ( NAME CHOLMOD_int64_single_bcsstk01
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk01.tri)

  add_test ( NAME CHOLMOD_int32_double_lp_afiro
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/lp_afiro.tri)
  add_test ( NAME CHOLMOD_int64_double_lp_afiro
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/lp_afiro.tri)
  add_test ( NAME CHOLMOD_int32_single_lp_afiro
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/lp_afiro.tri)
  add_test ( NAME CHOLMOD_int64_single_lp_afiro
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/lp_afiro.tri)

  add_test ( NAME CHOLMOD_int32_double_can24
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/can___24.mtx)
  add_test ( NAME CHOLMOD_int64_double_can24
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/can___24.mtx)
  add_test ( NAME CHOLMOD_int32_single_can24
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/can___24.mtx)
  add_test ( NAME CHOLMOD_int64_single_can24
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/can___24.mtx)

  add_test ( NAME CHOLMOD_int32_double_complex
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/c.tri)
  add_test ( NAME CHOLMOD_int64_double_complex
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/c.tri)
  add_test ( NAME CHOLMOD_int32_single_complex
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/c.tri)
  add_test ( NAME CHOLMOD_int64_single_complex
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/c.tri)

  add_test ( NAME CHOLMOD_int32_double_supernodal
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk02.tri)
  add_test ( NAME CHOLMOD_int64_double_supernodal
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_DL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk02.tri)
  add_test ( NAME CHOLMOD_int32_single_supernodal
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SI_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk02.tri)
  add_test ( NAME CHOLMOD_int64_single_supernodal
    COMMAND $<TARGET_FILE:${_SuiteSparse_CHOLMOD_SL_TEST_TARGET}>
      ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Demo/Matrix/bcsstk02.tri)

    if (WIN32)
      set (_SuiteSparse_CHOLMOD_TEST_PATH_MODIFICATION)
      foreach (_SuiteSparse_CHOLMOD_TEST_DEPENDENCY IN ITEMS
        amd
        camd
        ccolamd
        cholmod
        colamd
        suitesparseconfig
      )
        list (APPEND _SuiteSparse_CHOLMOD_TEST_PATH_MODIFICATION
          "PATH=path_list_prepend:$<TARGET_FILE_DIR:${_SuiteSparse_CHOLMOD_TEST_DEPENDENCY}>")
      endforeach (_SuiteSparse_CHOLMOD_TEST_DEPENDENCY)
      set_tests_properties (
      CHOLMOD_int32_double_bcsstk01
      CHOLMOD_int64_double_bcsstk01
      CHOLMOD_int32_single_bcsstk01
      CHOLMOD_int64_single_bcsstk01
      CHOLMOD_int32_double_lp_afiro
      CHOLMOD_int64_double_lp_afiro
      CHOLMOD_int32_single_lp_afiro
      CHOLMOD_int64_single_lp_afiro
      CHOLMOD_int32_double_can24
      CHOLMOD_int64_double_can24
      CHOLMOD_int32_single_can24
      CHOLMOD_int64_single_can24
      CHOLMOD_int32_double_complex
      CHOLMOD_int64_double_complex
      CHOLMOD_int32_single_complex
      CHOLMOD_int64_single_complex
      CHOLMOD_int32_double_supernodal
      CHOLMOD_int64_double_supernodal
      CHOLMOD_int32_single_supernodal
      CHOLMOD_int64_single_supernodal
        PROPERTIES ENVIRONMENT_MODIFICATION
          "${_SuiteSparse_CHOLMOD_TEST_PATH_MODIFICATION}")
    endif ()
  endif ()
endif ()

install (TARGETS ${_SuiteSparse_TARGETS} ${SuiteSparse_BIN_TARGETS}
  EXPORT suitesparse-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  LIBRARY DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development
  ARCHIVE DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development)

if (BUILD_CXSPARSE)
  add_subdirectory (CXSparse)
endif (BUILD_CXSPARSE)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.3
  AND CMAKE_EXPERIMENTAL_GENERATE_SBOM)
  set (_SuiteSparse_SBOM_EXPORTS suitesparse-targets)
  set (_SuiteSparse_SPDX_LICENSE "Apache-2.0 AND BSD-3-Clause")

  if (WITH_LICENSE MATCHES "GPL")
    string (APPEND _SuiteSparse_SPDX_LICENSE
      " AND LGPL-2.1-or-later")
  endif (WITH_LICENSE MATCHES "GPL")

  if (WITH_LICENSE STREQUAL "GPL")
    string (APPEND _SuiteSparse_SPDX_LICENSE
      " AND GPL-2.0-or-later")
  endif (WITH_LICENSE STREQUAL "GPL")

  if (BUILD_CXSPARSE AND CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
    list (APPEND _SuiteSparse_SBOM_EXPORTS cxsparse-targets)
  endif (BUILD_CXSPARSE AND CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)

  if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
    set (_SuiteSparse_SBOM_EXPORT_ARGUMENTS
      EXPORTS ${_SuiteSparse_SBOM_EXPORTS})
  else (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
    set (_SuiteSparse_SBOM_EXPORT_ARGUMENTS EXPORT suitesparse-targets)
  endif (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)

  set (_SuiteSparse_SBOM_ARGUMENTS
    PROJECT SuiteSparse
    DESTINATION ${CMAKE_INSTALL_DATADIR}/sbom
    VERSION ${SuiteSparse_VERSION}
    LICENSE ${_SuiteSparse_SPDX_LICENSE}
    DESCRIPTION "A suite of sparse matrix packages"
    HOMEPAGE_URL http://faculty.cse.tamu.edu/davis/suitesparse.html
    FORMAT spdx-3.0+json
  )

  if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
    list (APPEND _SuiteSparse_SBOM_ARGUMENTS
      PACKAGE_URL https://github.com/sergiud/SuiteSparse)
  endif (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)

  install (SBOM SuiteSparse
    ${_SuiteSparse_SBOM_EXPORT_ARGUMENTS}
    ${_SuiteSparse_SBOM_ARGUMENTS}
  )
endif (CMAKE_VERSION VERSION_GREATER_EQUAL 4.3
  AND CMAKE_EXPERIMENTAL_GENERATE_SBOM)

if (MSVC)
  foreach (bin IN LISTS SuiteSparse_BIN_TARGETS)
    install (CODE
"
include (BundleUtilities)
fixup_bundle (\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/${bin}${CMAKE_EXECUTABLE_SUFFIX} \"\" \"\")
"
      COMPONENT Runtime)
  endforeach (bin)
endif (MSVC)

add_feature_info (CUDA WITH_CUDA "GPU accelerated CHOLMOD using CUDA")
add_feature_info (CAMD WITH_CAMD "CAMD, CCOLAMD, CSYMAMD interfaces")
add_feature_info (Check WITH_CHECK "Check module")
add_feature_info (Cholesky WITH_CHOLESKY "Cholesky module")
add_feature_info (Fortran WITH_FORTRAN "Fortran compilation")
add_feature_info (MatrixOps WITH_MATRIXOPS "MatrixOps module")
add_feature_info (Partition WITH_METIS "Partition module")
add_feature_info (Print WITH_PRINT "Additional diagnostic messages")
add_feature_info (Supernodal WITH_SUPERNODAL "Supernodal module")
add_feature_info (Threading WITH_TBB "Parallelization using oneTBB")

set_package_properties (TBB PROPERTIES
  URL "https://www.threadingbuildingblocks.org"
  DESCRIPTION "C++ template library for task parallelism"
  TYPE RECOMMENDED
  PURPOSE "Enables multithreading support for SPQR"
)

set_package_properties (METIS PROPERTIES
  URL "http://glaros.dtc.umn.edu/gkhome/views/metis"
  DESCRIPTION "Serial programs for partitioning graphs, partitioning finite element meshes, and producing fill reducing orderings for sparse matrices"
)

set_package_properties (BLAS PROPERTIES
  URL "http://www.netlib.org/blas"
  DESCRIPTION "Basic Linear Algebra Subprograms"
)

set_package_properties (LAPACK PROPERTIES
  URL "http://www.netlib.org/lapack"
  DESCRIPTION "Linear Algebra PACKage"
)

set_package_properties (OpenMP PROPERTIES
  DESCRIPTION "Multi-platform shared memory multiprocessing"
  URL "http://www.openmp.org"
  PURPOSE "Enables multithreading support"
)

set_package_properties (Threads PROPERTIES
  DESCRIPTION "Multithreading support"
)

feature_summary (WHAT ALL)
