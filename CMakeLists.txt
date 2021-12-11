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

cmake_minimum_required (VERSION 3.22)

if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set (CMAKE_CUDA_ARCHITECTURES 75)
endif (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)

project (SuiteSparse LANGUAGES C CXX VERSION 5.10.1
  DESCRIPTION "A suite of sparse matrix packages"
  HOMEPAGE_URL http://faculty.cse.tamu.edu/davis/suitesparse.html
)

set (CMAKE_DEBUG_POSTFIX _debug)
set (CMAKE_C_VISIBILITY_PRESET hidden)
set (CMAKE_CXX_VISIBILITY_PRESET hidden)
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_VISIBILITY_INLINES_HIDDEN 1)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# SuiteSparse version
set (MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (VERSION_PATCH ${PROJECT_VERSION_PATCH})

set (SHORT_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

if (NOT VERSION_PATCH EQUAL 0)
  set (SHORT_VERSION "${SHORT_VERSION}.${VERSION_PATCH}")
endif (NOT VERSION_PATCH EQUAL 0)

set (VERSION ${PROJECT_VERSION})

set (CPACK_PACKAGE_NAME SuiteSparse)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SuiteSparse: an extended version of CSparse")
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_PACKAGE_VERSION ${SHORT_VERSION})

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
include (ComponentVersion)
include (CPack)
include (FeatureSummary)
include (GenerateExportHeader)
include (GNUInstallDirs)
include (SuiteSparseBase)

option (BUILD_CXSPARSE "Build CXSparse" ON)
option (WITH_CUDA "Enable CUDA support" ON)
option (WITH_DEMOS "Build demos" ON)
option (WITH_FORTRAN "Enables Fortran support" ON)
option (WITH_METIS "Enables METIS support" ON)
option (WITH_OPENMP "Enable OpenMP support" ON)
option (WITH_PRINT "Print diagnostic messages" OFF)
option (WITH_TBB "Enables Intel Threading Building Blocks support" ON)

set (WITH_LICENSE "GPL" CACHE STRING "Software license the binary distribution should adhere")
set_property (CACHE WITH_LICENSE PROPERTY STRINGS "Minimal;GPL;LGPL")

cmake_dependent_option (BUILD_SHARED_LIBS "Build shared libraires" OFF "WITH_LICENSE STREQUAL \"GPL\"" ON)

cmake_dependent_option (WITH_LGPL "Enable GNU LGPL modules" ON "WITH_LICENSE MATCHES \"GPL\"" OFF)
cmake_dependent_option (WITH_GPL "Enable GNU GPL modules" ON "WITH_LICENSE STREQUAL \"GPL\"" OFF)

cmake_dependent_option (WITH_CHOLMOD "Enable CHOLMOD" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_CHOLESKY "Enable the Cholesky module" ON "WITH_CHOLMOD" OFF)

cmake_dependent_option (WITH_CAMD "Enable interfaces to CAMD, CCOLAMD, CSYMAMD in Partition module" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_CHECK "Enable the Check module" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_MATRIXOPS "Enable the MatrixOps module" ON "WITH_GPL AND WITH_CHOLMOD" OFF)
cmake_dependent_option (WITH_MODIFY "Enable the Modify module" ON "WITH_GPL AND WITH_CHOLMOD" OFF)
cmake_dependent_option (WITH_SUPERNODAL "Enable the Supernodal module" ON "WITH_GPL" OFF)

if (WITH_CUDA)
  check_language (CUDA)

  if (CMAKE_CUDA_COMPILER)
    enable_language (CUDA)
  endif (CMAKE_CUDA_COMPILER)
else (WITH_CUDA)
  set (CMAKE_FIND_PACKAGE_DISABLE_CUDAToolkit ON)
endif (WITH_CUDA)

if (WITH_FORTRAN)
  check_language (Fortran)

  if (CMAKE_Fortran_COMPILER)
    enable_language (Fortran)
  endif (CMAKE_Fortran_COMPILER)
endif (WITH_FORTRAN)

check_c_compiler_flag (-Qstd=c99 HAVE_QSTD_C99)

if (NOT WITH_METIS)
  set (CMAKE_FIND_PACKAGE_DISABLE_METIS ON)
endif (NOT WITH_METIS)

if (NOT WITH_OPENMP)
  set (CMAKE_FIND_PACKAGE_DISABLE_OpenMP ON)
endif (NOT WITH_OPENMP)

if (NOT WITH_TBB)
  set (CMAKE_FIND_PACKAGE_DISABLE_TBB ON)
endif (NOT WITH_TBB)

find_package (BLAS REQUIRED)
find_package (CUDAToolkit)
find_package (LAPACK REQUIRED)
find_package (METIS 5.0)
find_package (OpenMP COMPONENTS C)
find_package (TBB 4.2...<2021.4 NO_MODULE)

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
  GPL
  MATRIXOPS
  MODIFY
  PARTITION
  PRINT
  SUPERNODAL
)

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
  AMD/Source/amd_global.c
  AMD/Source/amd_info.c
  AMD/Source/amd_order.c
  AMD/Source/amd_post_tree.c
  AMD/Source/amd_postorder.c
  AMD/Source/amd_preprocess.c
  AMD/Source/amd_valid.c
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
  CAMD/Source/camd_global.c
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
    CHOLMOD/Include/cholmod_blas.h
    CHOLMOD/Include/cholmod_camd.h
    CHOLMOD/Include/cholmod_check.h
    CHOLMOD/Include/cholmod_cholesky.h
    CHOLMOD/Include/cholmod_complexity.h
    CHOLMOD/Include/cholmod_core.h
    CHOLMOD/Include/cholmod_function.h
    CHOLMOD/Include/cholmod_gpu.h
    CHOLMOD/Include/cholmod_gpu_kernels.h
    CHOLMOD/Include/cholmod_internal.h
    CHOLMOD/Include/cholmod_io64.h
    CHOLMOD/Include/cholmod_matrixops.h
    CHOLMOD/Include/cholmod_modify.h
    CHOLMOD/Include/cholmod_partition.h
    CHOLMOD/Include/cholmod_supernodal.h
    CHOLMOD/Include/cholmod_template.h
  )

  set (CHOLMOD_SRCS
    CHOLMOD/Check/cholmod_check.c
    CHOLMOD/Check/cholmod_read.c
    CHOLMOD/Check/cholmod_write.c
    CHOLMOD/Cholesky/cholmod_amd.c
    CHOLMOD/Cholesky/cholmod_analyze.c
    CHOLMOD/Cholesky/cholmod_colamd.c
    CHOLMOD/Cholesky/cholmod_etree.c
    CHOLMOD/Cholesky/cholmod_factorize.c
    CHOLMOD/Cholesky/cholmod_postorder.c
    CHOLMOD/Cholesky/cholmod_rcond.c
    CHOLMOD/Cholesky/cholmod_resymbol.c
    CHOLMOD/Cholesky/cholmod_rowcolcounts.c
    CHOLMOD/Cholesky/cholmod_rowfac.c
    CHOLMOD/Cholesky/cholmod_solve.c
    CHOLMOD/Cholesky/cholmod_spsolve.c
    #CHOLMOD/Cholesky/t_cholmod_lsolve.c
    #CHOLMOD/Cholesky/t_cholmod_ltsolve.c
    #CHOLMOD/Cholesky/t_cholmod_rowfac.c
    #CHOLMOD/Cholesky/t_cholmod_solve.c
    CHOLMOD/Core/cholmod_aat.c
    CHOLMOD/Core/cholmod_add.c
    CHOLMOD/Core/cholmod_band.c
    CHOLMOD/Core/cholmod_change_factor.c
    CHOLMOD/Core/cholmod_common.c
    CHOLMOD/Core/cholmod_complex.c
    CHOLMOD/Core/cholmod_copy.c
    CHOLMOD/Core/cholmod_dense.c
    CHOLMOD/Core/cholmod_error.c
    CHOLMOD/Core/cholmod_factor.c
    CHOLMOD/Core/cholmod_memory.c
    CHOLMOD/Core/cholmod_sparse.c
    CHOLMOD/Core/cholmod_transpose.c
    CHOLMOD/Core/cholmod_triplet.c
    CHOLMOD/Core/cholmod_version.c
    #CHOLMOD/Core/t_cholmod_change_factor.c
    #CHOLMOD/Core/t_cholmod_dense.c
    #CHOLMOD/Core/t_cholmod_transpose.c
    #CHOLMOD/Core/t_cholmod_triplet.c
    #CHOLMOD/Demo/cholmod_demo.c
    #CHOLMOD/Demo/cholmod_l_demo.c
    #CHOLMOD/Demo/cholmod_simple.c
    #CHOLMOD/GPU/cholmod_gpu.c
    #CHOLMOD/GPU/cholmod_gpu_kernels.c
    CHOLMOD/MatrixOps/cholmod_drop.c
    CHOLMOD/MatrixOps/cholmod_horzcat.c
    CHOLMOD/MatrixOps/cholmod_norm.c
    CHOLMOD/MatrixOps/cholmod_scale.c
    CHOLMOD/MatrixOps/cholmod_sdmult.c
    CHOLMOD/MatrixOps/cholmod_ssmult.c
    CHOLMOD/MatrixOps/cholmod_submatrix.c
    CHOLMOD/MatrixOps/cholmod_symmetry.c
    CHOLMOD/MatrixOps/cholmod_vertcat.c
    #CHOLMOD/MatrixOps/t_cholmod_sdmult.c
    CHOLMOD/Modify/cholmod_rowadd.c
    CHOLMOD/Modify/cholmod_rowdel.c
    CHOLMOD/Modify/cholmod_updown.c
    #CHOLMOD/Modify/t_cholmod_updown.c
    #CHOLMOD/Modify/t_cholmod_updown_numkr.c
    CHOLMOD/Partition/cholmod_camd.c
    CHOLMOD/Partition/cholmod_ccolamd.c
    CHOLMOD/Partition/cholmod_csymamd.c
    CHOLMOD/Partition/cholmod_metis.c
    CHOLMOD/Partition/cholmod_nesdis.c
    CHOLMOD/Supernodal/cholmod_super_numeric.c
    CHOLMOD/Supernodal/cholmod_super_solve.c
    CHOLMOD/Supernodal/cholmod_super_symbolic.c
    #CHOLMOD/Supernodal/t_cholmod_super_numeric.c
    #CHOLMOD/Supernodal/t_cholmod_super_solve.c
  )

  if (NOT METIS_FOUND)
    set_property (SOURCE
      CHOLMOD/Check/cholmod_check.c
      CHOLMOD/Cholesky/cholmod_analyze.c
      CHOLMOD/Partition/cholmod_camd.c
      CHOLMOD/Partition/cholmod_ccolamd.c
      CHOLMOD/Partition/cholmod_csymamd.c
      CHOLMOD/Partition/cholmod_m.c
      CHOLMOD/Partition/cholmod_metis.c
      CHOLMOD/Partition/cholmod_nesdis.c
      SPQR/Source/spqr_1colamd.cpp
      SPQR/Source/spqr_analyze.cpp
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
  UMFPACK/Include/umfpack_col_to_triplet.h
  UMFPACK/Include/umfpack_defaults.h
  UMFPACK/Include/umfpack_free_numeric.h
  UMFPACK/Include/umfpack_free_symbolic.h
  UMFPACK/Include/umfpack_get_determinant.h
  UMFPACK/Include/umfpack_get_lunz.h
  UMFPACK/Include/umfpack_get_numeric.h
  UMFPACK/Include/umfpack_get_symbolic.h
  UMFPACK/Include/umfpack_global.h
  UMFPACK/Include/umfpack_load_numeric.h
  UMFPACK/Include/umfpack_load_symbolic.h
  UMFPACK/Include/umfpack_numeric.h
  UMFPACK/Include/umfpack_qsymbolic.h
  UMFPACK/Include/umfpack_report_control.h
  UMFPACK/Include/umfpack_report_info.h
  UMFPACK/Include/umfpack_report_matrix.h
  UMFPACK/Include/umfpack_report_numeric.h
  UMFPACK/Include/umfpack_report_perm.h
  UMFPACK/Include/umfpack_report_status.h
  UMFPACK/Include/umfpack_report_symbolic.h
  UMFPACK/Include/umfpack_report_triplet.h
  UMFPACK/Include/umfpack_report_vector.h
  UMFPACK/Include/umfpack_save_numeric.h
  UMFPACK/Include/umfpack_save_symbolic.h
  UMFPACK/Include/umfpack_scale.h
  UMFPACK/Include/umfpack_solve.h
  UMFPACK/Include/umfpack_symbolic.h
  UMFPACK/Include/umfpack_tictoc.h
  UMFPACK/Include/umfpack_timer.h
  UMFPACK/Include/umfpack_transpose.h
  UMFPACK/Include/umfpack_triplet_to_col.h
  UMFPACK/Include/umfpack_wsolve.h
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
  UMFPACK/Source/umfpack_defaults.c
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
  UMFPACK/Source/umfpack_defaults.c
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
  SuiteSparse_GPURuntime/Include/SuiteSparseGPU_Runtime.hpp
  SuiteSparse_GPURuntime/Include/SuiteSparseGPU_Workspace.hpp
  SuiteSparse_GPURuntime/Include/SuiteSparseGPU_debug.hpp
  SuiteSparse_GPURuntime/Include/SuiteSparseGPU_macros.hpp
  SuiteSparse_GPURuntime/Include/SuiteSparseGPU_workspace_macros.hpp
)

set (GPUQRENGINE_HDRS
  GPUQREngine/Include/GPUQREngine.hpp
  GPUQREngine/Include/GPUQREngine_BucketList.hpp
  GPUQREngine/Include/GPUQREngine_Common.hpp
  GPUQREngine/Include/GPUQREngine_Front.hpp
  GPUQREngine/Include/GPUQREngine_FrontState.hpp
  GPUQREngine/Include/GPUQREngine_GraphVizHelper.hpp
  GPUQREngine/Include/GPUQREngine_Internal.hpp
  GPUQREngine/Include/GPUQREngine_LLBundle.hpp
  GPUQREngine/Include/GPUQREngine_Scheduler.hpp
  GPUQREngine/Include/GPUQREngine_SEntry.hpp
  GPUQREngine/Include/GPUQREngine_SparseMeta.hpp
  GPUQREngine/Include/GPUQREngine_Stats.hpp
  GPUQREngine/Include/GPUQREngine_TaskDescriptor.hpp
  GPUQREngine/Include/GPUQREngine_Timing.hpp
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
  ${AMD_HDRS}
)
add_library (SuiteSparse::AMD ALIAS amd)

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
  ${BTF_HDRS}
)
add_library (SuiteSparse::BTF ALIAS btf)

target_include_directories (btf PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
)

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
  ${CAMD_HDRS}
)
add_library (SuiteSparse::CAMD ALIAS camd)

if (WITH_CHOLMOD)
  add_library (cholmod_int OBJECT ${CHOLMOD_SRCS})
  target_include_directories (cholmod_int PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_int PRIVATE DINT)
  target_compile_definitions (cholmod_int PRIVATE cholmod_EXPORTS)

  add_library (cholmod_long OBJECT ${CHOLMOD_SRCS})
  target_include_directories (cholmod_long PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
  target_compile_definitions (cholmod_long PRIVATE DLONG)
  target_compile_definitions (cholmod_long PRIVATE cholmod_EXPORTS)

  add_library (cholmod
    $<TARGET_OBJECTS:cholmod_int>
    $<TARGET_OBJECTS:cholmod_long>
    ${CHOLMOD_HDRS}
  )
  add_library (SuiteSparse::CHOLMOD ALIAS cholmod)

  target_include_directories (cholmod PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include>
  )

  if (OpenMP_FOUND)
    target_compile_definitions (cholmod PRIVATE HAVE_OPENMP)
    target_link_libraries (cholmod PRIVATE OpenMP::OpenMP_C)

    set (OpenMP_DEPENDENCY "find_dependency (OpenMP ${OpenMP_C_VERSION} COMPONENTS C)")
  endif (OpenMP_FOUND)

  if (HAVE_FEENABLEEXCEPT)
    target_compile_definitions (cholmod PRIVATE _GNU_SOURCE)
  endif (HAVE_FEENABLEEXCEPT)

  if (METIS_FOUND)
    target_include_directories (cholmod SYSTEM PRIVATE ${METIS_INCLUDE_DIRS})
    target_include_directories (cholmod SYSTEM PRIVATE
      $<BUILD_INTERFACE:${METIS_INCLUDE_DIRS}>)
    target_link_libraries (cholmod PRIVATE ${METIS_LIBRARIES})
  endif (METIS_FOUND)

  target_link_libraries (cholmod
    PRIVATE BLAS::BLAS
            LAPACK::LAPACK
            SuiteSparse::AMD
            SuiteSparse::CAMD
            SuiteSparse::CCOLAMD
            SuiteSparse::COLAMD
    PUBLIC SuiteSparse::Config
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

  if (WITH_TBB AND TBB_FOUND)
    target_link_libraries (spqr PRIVATE TBB::tbb)
    target_compile_definitions (spqr PRIVATE HAVE_TBB)

    set (TBB_DEPENDENCY "find_dependency (TBB ${TBB_VERSION})")
  endif (WITH_TBB AND TBB_FOUND)
endif (WITH_SUPERNODAL)

add_library (suitesparseconfig
  $<TARGET_OBJECTS:suitesparsebase>
)
add_library (SuiteSparse::Config ALIAS suitesparseconfig)

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
  ${UMFPACK_HDRS}
  UMFPACK/Source/umfpack_global.c
  UMFPACK/Source/umfpack_tictoc.c
  UMFPACK/Source/umfpack_timer.c
)
add_library (SuiteSparse::UMFPACK ALIAS umfpack)

target_link_libraries (umfpack PRIVATE BLAS::BLAS LAPACK::LAPACK)

if (TARGET cholmod)
  # NOTE cholmod already depends on LAPACK
  target_link_libraries (umfpack PRIVATE SuiteSparse::CHOLMOD)
else (TARGET cholmod)
  set_property (SOURCE
    UMFPACK/Source/umf_cholmod.c
    APPEND PROPERTY COMPILE_DEFINITIONS NCHOLMOD
  )
endif (TARGET cholmod)

set (KLU_NO_COMPLEX_SRCS
  KLU/Source/klu_defaults.c
  KLU/Source/klu_analyze.c
  KLU/Source/klu_analyze_given.c
  KLU/Source/klu_free_symbolic.c
  KLU/Source/klu_memory.c
)

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
  ${KLU_HDRS}
)
add_library (SuiteSparse::KLU ALIAS klu)

target_link_libraries (klu PUBLIC
  SuiteSparse::AMD
  SuiteSparse::BTF
  SuiteSparse::COLAMD
)

if (TARGET cholmod)
  target_link_libraries (klu PRIVATE SuiteSparse::CHOLMOD)
endif (TARGET cholmod)

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

target_link_libraries (ldl PRIVATE SuiteSparse::AMD)

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

target_include_directories (suitesparseconfig PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}>"
)

if (HAVE_BLAS_UNDERSCORE)
  add_compile_definitions (BLAS_UNDERSCORE)
elseif (HAVE_BLAS_NO_UNDERSCORE)
  add_compile_definitions (BLAS_NO_UNDERSCORE)
endif (HAVE_BLAS_UNDERSCORE)

if (HAVE_M)
  target_link_libraries (suitesparseconfig INTERFACE m)
endif (HAVE_M)

target_include_directories (amd PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>)
target_include_directories (btf PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>)
target_include_directories (camd PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_CAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CAMD_INCLUDE_DIR}>)
target_include_directories (ccolamd PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_CCOLAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CCOLAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CCOLAMD_INCLUDE_DIR}>)
target_include_directories (colamd PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>)

if (TARGET cholmod)
  target_include_directories (cholmod PUBLIC
    $<INSTALL_INTERFACE:${SuiteSparse_CHOLMOD_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Cholesky>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Core>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/MatrixOps>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Modify>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Supernodal>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CHOLMOD_INCLUDE_DIR}>
  )
endif (TARGET cholmod)

target_include_directories (umfpack PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_UMFPACK_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Source>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_UMFPACK_INCLUDE_DIR}>
)
target_include_directories (klu PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_KLU_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_KLU_INCLUDE_DIR}>)
target_include_directories (ldl PUBLIC
  $<INSTALL_INTERFACE:${SuiteSparse_LDL_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LDL/Include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_LDL_INCLUDE_DIR}>)

target_link_libraries (amd PUBLIC SuiteSparse::Config)
target_link_libraries (btf PUBLIC SuiteSparse::Config)
target_link_libraries (camd PUBLIC SuiteSparse::Config)
target_link_libraries (ccolamd PUBLIC SuiteSparse::Config)
target_link_libraries (colamd PUBLIC SuiteSparse::Config)
target_link_libraries (klu PUBLIC SuiteSparse::Config)
target_link_libraries (ldl PUBLIC SuiteSparse::Config)
target_link_libraries (umfpack PUBLIC amd SuiteSparse::Config)

if (TARGET spqr)
  target_include_directories (spqr PUBLIC
    $<INSTALL_INTERFACE:${SuiteSparse_SPQR_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SPQR_INCLUDE_DIR}>
  )

  target_link_libraries (spqr
    PUBLIC SuiteSparse::CHOLMOD SuiteSparse::Config
    PRIVATE BLAS::BLAS LAPACK::LAPACK
  )
endif (TARGET spqr)

# CUDA support
if (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)
  add_library (gpuruntime
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_cpuAllocators.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_gpuAllocators.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_transfer.cpp
  )
  add_library (SuiteSparse::GPURuntime ALIAS gpuruntime)

  target_include_directories (gpuruntime PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_GPURuntime/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_GPURuntime_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_GPURuntime_INCLUDE_DIR}>
  )

  target_link_libraries (gpuruntime
    PUBLIC SuiteSparse::Config
    PRIVATE CUDA::cudart
  )

  # GPUQREngine

  add_library (gpuqrengine
    #GPUQREngine/Include/Kernel/Apply/block_apply.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_1.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_1_by_1.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_2.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_2_by_1.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_3.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_3_by_1.cu
    #GPUQREngine/Include/Kernel/Apply/block_apply_chunk.cu
    #GPUQREngine/Include/Kernel/Apply/cevta_tile.cu
    #GPUQREngine/Include/Kernel/Apply/params_apply.hpp
    #GPUQREngine/Include/Kernel/Apply/pipelined_rearrange.cu
    #GPUQREngine/Include/Kernel/Assemble/packAssemble.cu
    #GPUQREngine/Include/Kernel/Assemble/sAssemble.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_3_by_1.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_1_by_1.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_1_by_1_edge.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_2_by_1.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_2_by_1_edge.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_3_by_1.cu
    #GPUQREngine/Include/Kernel/Factorize/factorize_vt_3_by_1_edge.cu
    #GPUQREngine/Include/Kernel/qrKernel.cu
    GPUQREngine/Include/Kernel/sharedMemory.hpp
    #GPUQREngine/Include/Kernel/uberKernel.cu
    GPUQREngine/Source/BucketList/BucketList.cpp
    GPUQREngine/Source/BucketList/BucketList_AdvanceBundles.cpp
    GPUQREngine/Source/BucketList/BucketList_CreateBundles.cpp
    GPUQREngine/Source/BucketList/BucketList_FillWorkQueue.cpp
    GPUQREngine/Source/BucketList/BucketList_GrowBundles.cpp
    GPUQREngine/Source/BucketList/BucketList_Manage.cpp
    GPUQREngine/Source/BucketList/BucketList_PostProcessing.cpp
    GPUQREngine/Source/GPUQREngine_ExpertDense.cpp
    GPUQREngine/Source/GPUQREngine_ExpertSparse.cpp
    GPUQREngine/Source/GPUQREngine_GraphVizHelper.cpp
    GPUQREngine/Source/GPUQREngine_Internal.cpp
    GPUQREngine/Source/GPUQREngine_UberKernel.cu
    GPUQREngine/Source/LLBundle/LLBundle.cpp
    GPUQREngine/Source/LLBundle/LLBundle_AddTiles.cpp
    GPUQREngine/Source/LLBundle/LLBundle_Advance.cpp
    GPUQREngine/Source/LLBundle/LLBundle_GPUPack.cpp
    GPUQREngine/Source/LLBundle/LLBundle_PipelinedRearrange.cpp
    GPUQREngine/Source/LLBundle/LLBundle_UpdateSecondMinIndex.cpp
    GPUQREngine/Source/Scheduler/Scheduler.cpp
    GPUQREngine/Source/Scheduler/Scheduler_FillWorkQueue.cpp
    GPUQREngine/Source/Scheduler/Scheduler_Front.cpp
    GPUQREngine/Source/Scheduler/Scheduler_LaunchKernel.cpp
    GPUQREngine/Source/Scheduler/Scheduler_PostProcess.cpp
    GPUQREngine/Source/Scheduler/Scheduler_Render.cpp
    GPUQREngine/Source/Scheduler/Scheduler_TransferData.cpp
    GPUQREngine/Source/Scheduler/ssgpu_maxQueueSize.cpp
    GPUQREngine/Source/TaskDescriptor/TaskDescriptor_flops.cpp
  )

  target_include_directories (gpuqrengine PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GPUQREngine/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_GPUQREngine_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_GPUQREngine_INCLUDE_DIR}>
  )

  target_include_directories (gpuqrengine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/GPUQREngine/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_GPURuntime/Include
  )

  target_link_libraries (gpuqrengine
    PUBLIC CUDA::cublas
    PRIVATE gpuruntime
  )

  if (WITH_CHOLMOD)
    set (CHOLMOD_GPU_SRCS CHOLMOD/GPU/cholmod_gpu.c)

    add_library (cholmod_gpu_int OBJECT ${CHOLMOD_GPU_SRCS})
    target_include_directories (cholmod_gpu_int PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_gpu_int PRIVATE DINT)
    target_compile_definitions (cholmod_gpu_int PRIVATE cholmod_EXPORTS)

    add_library (cholmod_gpu_long OBJECT ${CHOLMOD_GPU_SRCS})
    target_include_directories (cholmod_gpu_long PRIVATE $<TARGET_PROPERTY:cholmod,INCLUDE_DIRECTORIES>)
    target_compile_definitions (cholmod_gpu_long PRIVATE DLONG)
    target_compile_definitions (cholmod_gpu_long PRIVATE cholmod_EXPORTS)

    target_sources (cholmod PRIVATE
      $<TARGET_OBJECTS:cholmod_gpu_int>
      $<TARGET_OBJECTS:cholmod_gpu_long>
      CHOLMOD/GPU/cholmod_gpu_kernels.cu
    )

    target_compile_definitions (cholmod PUBLIC GPU_BLAS)
    target_link_libraries (cholmod
      PRIVATE CUDA::cudart gpuruntime
      PUBLIC CUDA::cublas
    )
  endif (WITH_CHOLMOD)

  # SPQR GPU
  if (TARGET spqr)
    target_sources (spqr PRIVATE
      SPQR/SPQRGPU/spqrgpu_buildAssemblyMaps.cpp
      SPQR/SPQRGPU/spqrgpu_computeFrontStaging.cpp
      SPQR/SPQRGPU/spqrgpu_kernel.cpp
    )

    target_link_libraries (spqr PRIVATE gpuqrengine)
  endif (TARGET spqr)
endif (WITH_CUDA AND CMAKE_CUDA_COMPILER_LOADED)

set (_SuiteSparse_TARGETS amd btf camd ccolamd colamd umfpack
  klu ldl suitesparseconfig)

set_property (TARGET suitesparseconfig PROPERTY NAME_ALIAS SuiteSparse)
set_property (TARGET suitesparseconfig PROPERTY HEADERS_VARIABLE
  SUITESPARSECONFIG_HDRS)
set_property (TARGET suitesparseconfig PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR})

if (TARGET cholmod)
  list (APPEND _SuiteSparse_TARGETS cholmod)
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

foreach (_target ${_SuiteSparse_TARGETS})
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

install (TARGETS ${_SuiteSparse_TARGETS}
  EXPORT suitesparse-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  LIBRARY DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development
  ARCHIVE DESTINATION ${SuiteSparse_CMake_LIBDIR} COMPONENT Development)

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

set_component_version (amd 2.4.6)
set_component_version (btf 1.2.6)
set_component_version (camd 2.4.6)
set_component_version (ccolamd 2.9.6)
set_component_version (colamd 2.9.6)
set_component_version (klu 1.3.9)
set_component_version (ldl 2.2.6)
set_component_version (suitesparseconfig "${VERSION}")
set_component_version (umfpack 5.7.9)

if (TARGET cholmod)
  set_component_version (cholmod 3.0.13)
endif (TARGET cholmod)

if (TARGET spqr)
  set_component_version (spqr 2.0.9)
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

set_target_properties (amd PROPERTIES EXPORT_NAME AMD)
set_target_properties (btf PROPERTIES EXPORT_NAME BTF)
set_target_properties (camd PROPERTIES EXPORT_NAME CAMD)
set_target_properties (ccolamd PROPERTIES EXPORT_NAME CCOLAMD)
set_target_properties (colamd PROPERTIES EXPORT_NAME COLAMD)
set_target_properties (klu PROPERTIES EXPORT_NAME KLU)
set_target_properties (ldl PROPERTIES EXPORT_NAME LDL)
set_target_properties (suitesparseconfig PROPERTIES EXPORT_NAME Config)
set_target_properties (umfpack PROPERTIES EXPORT_NAME UMFPACK)

if (TARGET cholmod)
  set_target_properties (cholmod PROPERTIES EXPORT_NAME CHOLMOD)
endif (TARGET cholmod)

if (TARGET spqr)
  set_target_properties (spqr PROPERTIES EXPORT_NAME SPQR)
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

  if (TARGET cholmod AND WITH_CHECK)
    add_executable (cholmod_demo
      CHOLMOD/Demo/cholmod_demo.h
      CHOLMOD/Demo/cholmod_demo.c
    )
    target_link_libraries (cholmod_demo PRIVATE SuiteSparse::CHOLMOD)

    add_executable (cholmod_l_demo CHOLMOD/Demo/cholmod_l_demo.c)
    target_link_libraries (cholmod_l_demo PRIVATE SuiteSparse::CHOLMOD)

    add_executable (cholmod_simple CHOLMOD/Demo/cholmod_simple.c)
    target_link_libraries (cholmod_simple PRIVATE SuiteSparse::CHOLMOD)
  endif (TARGET cholmod AND WITH_CHECK)

  add_executable (colamd_example
    COLAMD/Demo/colamd_example.c
  )
  target_link_libraries (colamd_example PRIVATE SuiteSparse::COLAMD)

  add_executable (colamd_l_example
    COLAMD/Demo/colamd_l_example.c
  )
  target_link_libraries (colamd_l_example PRIVATE SuiteSparse::COLAMD)

  add_executable (klu_simple
    KLU/Demo/klu_simple.c
  )
  target_link_libraries (klu_simple PRIVATE SuiteSparse::KLU)

  if (TARGET cholmod)
    add_executable (kludemo
      KLU/Demo/kludemo.c
    )
  target_link_libraries (kludemo PRIVATE SuiteSparse::KLU SuiteSparse::CHOLMOD)

    add_executable (kluldemo
      KLU/Demo/kluldemo.c
    )
    target_link_libraries (kluldemo PRIVATE SuiteSparse::KLU SuiteSparse::CHOLMOD)
  endif (TARGET cholmod)

  add_executable (ldlamd
    LDL/Demo/ldlmain.c
  )
  target_compile_definitions (ldlamd PRIVATE USE_AMD)
  target_link_libraries (ldlamd PRIVATE SuiteSparse::AMD SuiteSparse::LDL)

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
  target_link_libraries (ldllamd PRIVATE SuiteSparse::AMD SuiteSparse::LDL)

  add_executable (ldlsimple
    LDL/Demo/ldlsimple.c
  )
  target_link_libraries (ldlsimple PRIVATE SuiteSparse::LDL)

  add_executable (ldllsimple
    LDL/Demo/ldlsimple.c
  )
  target_compile_definitions (ldllsimple PRIVATE LDL_LONG)
  target_link_libraries (ldllsimple PRIVATE SuiteSparse::LDL)

  add_executable (umf4 UMFPACK/Demo/umf4.c)
  target_link_libraries (umf4 PRIVATE SuiteSparse::UMFPACK)

  add_executable (umfpack_di_demo UMFPACK/Demo/umfpack_di_demo.c)
  target_link_libraries (umfpack_di_demo PRIVATE SuiteSparse::UMFPACK)

  add_executable (umfpack_dl_demo UMFPACK/Demo/umfpack_dl_demo.c)
  target_link_libraries (umfpack_dl_demo PRIVATE SuiteSparse::UMFPACK)

  add_executable (umfpack_zi_demo UMFPACK/Demo/umfpack_zi_demo.c)
  target_link_libraries (umfpack_zi_demo PRIVATE SuiteSparse::UMFPACK)

  add_executable (umfpack_zl_demo UMFPACK/Demo/umfpack_zl_demo.c)
  target_link_libraries (umfpack_zl_demo PRIVATE SuiteSparse::UMFPACK)

  add_executable (umfpack_simple UMFPACK/Demo/umfpack_simple.c)
  target_link_libraries (umfpack_simple PRIVATE SuiteSparse::UMFPACK)

  if (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)
    add_executable (umf4hb UMFPACK/Demo/umf4_f77wrapper.c UMFPACK/Demo/umf4hb.f)
    target_link_libraries (umf4hb PRIVATE SuiteSparse::UMFPACK)
  endif (WITH_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED)

  if (TARGET spqr)
    add_executable (qrdemo SPQR/Demo/qrdemo.cpp)
    add_executable (qrdemoc SPQR/Demo/qrdemoc.c)
    add_executable (qrsimple SPQR/Demo/qrsimple.cpp)
    add_executable (qrsimplec SPQR/Demo/qrsimplec.c)

    target_link_libraries (qrdemo PRIVATE SuiteSparse::SPQR)
    target_link_libraries (qrdemoc PRIVATE SuiteSparse::SPQR)
    target_link_libraries (qrsimple PRIVATE SuiteSparse::SPQR)
    target_link_libraries (qrsimplec PRIVATE SuiteSparse::SPQR)

    if (TARGET gpuqrengine)
      add_executable (qrdemo_gpu SPQR/Demo/qrdemo_gpu.cpp)
      target_link_libraries (qrdemo_gpu PRIVATE
        SuiteSparse::CHOLMOD
        SuiteSparse::GPURuntime
        SuiteSparse::SPQR
      )
    endif (TARGET gpuqrengine)
  endif (TARGET spqr)
endif (WITH_DEMOS)

if (BUILD_CXSPARSE)
  add_subdirectory (CXSparse)
endif (BUILD_CXSPARSE)

add_feature_info (CUDA WITH_CUDA "GPU accelerated CHOLMOD using CUDA")
add_feature_info (GPL WITH_GPL "GNU GPL modules")
add_feature_info (LGPL WITH_LGPL "GNU LGPL modules")
add_feature_info (CAMD WITH_CAMD "CAMD, CCOLAMD, CSYMAMD interfaces")
add_feature_info (Check WITH_CHECK "Check module")
add_feature_info (Cholesky WITH_CHOLESKY "Cholesky module")
add_feature_info (Fortran WITH_FORTRAN "Fortran compilation")
add_feature_info (MatrixOps WITH_MATRIXOPS "MatrixOps module")
add_feature_info (Partition WITH_METIS "Partition module")
add_feature_info (Print WITH_PRINT "Additional diagnostic messages")
add_feature_info (Supernodal WITH_SUPERNODAL "Supernodal module")
add_feature_info (Threading WITH_TBB "Parallelization using Intel Threading Building Blocks")

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
