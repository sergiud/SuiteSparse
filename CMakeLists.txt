# Copyright (c) 2016, 2017 Sergiu Deitsch
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

if (POLICY CMP0007)
  cmake_policy (SET CMP0007 NEW)
endif (POLICY CMP0007)

if (POLICY CMP0022)
  cmake_policy (SET CMP0022 NEW)
endif (POLICY CMP0022)

if (POLICY CMP0023)
  cmake_policy (SET CMP0023 NEW)
endif (POLICY CMP0023)

if (POLICY CMP0042)
  cmake_policy (SET CMP0042 NEW)
endif (POLICY CMP0042)

if (POLICY CMP0054)
  cmake_policy (SET CMP0054 NEW)
endif (POLICY CMP0054)

if (POLICY CMP0063)
  cmake_policy (SET CMP0063 NEW)
endif (POLICY CMP0063)

include (CheckLanguage)

option (BUILD_CXSPARSE "Build CXSparse" ON)

option (WITH_FORTRAN "Enables Fortran support" ON)

if (WITH_FORTRAN)
  check_language (Fortran)

  if (CMAKE_Fortran_COMPILER)
    enable_language (Fortran)
  endif (CMAKE_Fortran_COMPILER)
endif (WITH_FORTRAN)

project (SuiteSparse C CXX)

set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX _debug)

include (CheckCCompilerFlag)
include (CheckFunctionExists)
include (CMakeDependentOption)
include (FeatureSummary)
include (GenerateExportHeader)
include (GNUInstallDirs)

set (CMAKE_C_VISIBILITY_PRESET hidden)
set (CMAKE_CXX_VISIBILITY_PRESET hidden)
set (CMAKE_VISIBILITY_INLINES_HIDDEN 1)

check_c_compiler_flag (-fPIC HAVE_FPIC)
check_c_compiler_flag (-Qstd=c99 HAVE_QSTD_C99)

option (WITH_CUDA "Enable CUDA support" ON)
option (WITH_OPENMP "Enable OpenMP support" ON)

set (CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

if (WITH_CUDA)
  find_package (CUDA 7.0)
endif (WITH_CUDA)

if (WITH_OPENMP)
  find_package (OpenMP)
endif (WITH_OPENMP)

find_package (METIS 5.0)

option (WITH_DEMOS "Build demos" ON)
option (WITH_PRINT "Print diagnostic messages" OFF)
option (WITH_TBB "Enables Intel Threading Building Blocks support" ON)
option (WITH_LGPL "Enable GNU LGPL modules" ON)
cmake_dependent_option (WITH_GPL "Enable GNU GPL modules" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_CAMD "Enable interfaces to CAMD, CCOLAMD, CSYMAMD
in Partition module" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_PARTITION "Enable the Partition module" ON
  "WITH_LGPL AND METIS_FOUND" OFF)
cmake_dependent_option (WITH_CHOLESKY "Enable the Cholesky module" ON
  "WITH_LGPL" OFF)
cmake_dependent_option (WITH_CHECK "Enable the Check module" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_MODIFY "Enable the Modify module" ON
  "WITH_GPL" OFF)
cmake_dependent_option (WITH_SUPERNODAL "Enable the Supernodal module" ON
  "WITH_GPL" OFF)
cmake_dependent_option (WITH_MATRIXOPS "Enable the MatrixOps module" ON
  "WITH_GPL" OFF)

function (check_options)
  foreach (_OPT ${ARGV})
    if (NOT DEFINED WITH_${_OPT})
      message (FATAL_ERROR "Option WITH_${_OPT} is not  defined")
    endif (NOT DEFINED WITH_${_OPT})

    if (NOT WITH_${_OPT})
      set (N${_OPT} 1 PARENT_SCOPE)
    else (NOT WITH_${_OPT})
      set (N${_OPT} 0 PARENT_SCOPE)
    endif (NOT WITH_${_OPT})
  endforeach (_OPT)
endfunction (check_options)

check_options (CAMD CHECK CHOLESKY GPL MATRIXOPS MODIFY PARTITION PRINT
  SUPERNODAL)

# SuiteSparse version
set (MAJOR_VERSION 4)
set (MINOR_VERSION 5)
set (VERSION_PATCH 5)

set (SHORT_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

if (NOT VERSION_PATCH EQUAL 0)
  set (SHORT_VERSION "${SHORT_VERSION}.${VERSION_PATCH}")
endif (NOT VERSION_PATCH EQUAL 0)

set (VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${VERSION_PATCH}")

set (CPACK_PACKAGE_NAME SuiteSparse)
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SuiteSparse: an extended version of CSparse")
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set (CPACK_PACKAGE_VERSION ${SHORT_VERSION})

include (CMakePackageConfigHelpers)
include (CMakePushCheckState)
include (CPack)
include (CheckIncludeFile)
include (CheckLibraryExists)
include (CheckSymbolExists)
include (CheckTypeSize)
include (ComponentVersion)
include (GenerateSourceFiles)

find_package (LAPACK QUIET NO_MODULE)

if (NOT LAPACK_FOUND)
  find_package (LAPACK MODULE REQUIRED)
endif (NOT LAPACK_FOUND)

cmake_push_check_state (RESET)

set (CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

check_function_exists (sgemm HAVE_BLAS_NO_UNDERSCORE)
check_function_exists (sgemm_ HAVE_BLAS_UNDERSCORE)

cmake_pop_check_state ()

if (WITH_TBB)
  find_package (TBB 3.0)
endif (WITH_TBB)

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

source_group ("Generated Files" REGULAR_EXPRESSION ".+_(di|dl|zi|zl)(_.+)?\\.c$")

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS
    _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

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

set (AMD_HDRS
  AMD/Include/amd.h
)

set (BTF_SRCS
  BTF/Source/btf_maxtrans.c
  BTF/Source/btf_order.c
  BTF/Source/btf_strongcomp.c
)

set (BTF_HDRS
  BTF/Include/btf.h
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

set (CAMD_HDRS
  CAMD/Include/camd.h
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

if (NOT WITH_PARTITION OR NOT METIS_FOUND)
  set_property (SOURCE
    CHOLMOD/Partition/cholmod_camd.c
    CHOLMOD/Partition/cholmod_ccolamd.c
    CHOLMOD/Partition/cholmod_csymamd.c
    CHOLMOD/Partition/cholmod_m.c
    CHOLMOD/Partition/cholmod_metis.c
    CHOLMOD/Partition/cholmod_nesdis.c
    APPEND PROPERTY COMPILE_DEFINITIONS NPARTITION)
endif (NOT WITH_PARTITION OR NOT METIS_FOUND)

configure_file (CHOLMOD/Include/cholmod_config.h.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include/cholmod_config.h)

set (CHOLMOD_HDRS
  CHOLMOD/Include/cholmod.h
  CHOLMOD/Include/cholmod_blas.h
  CHOLMOD/Include/cholmod_camd.h
  CHOLMOD/Include/cholmod_check.h
  CHOLMOD/Include/cholmod_cholesky.h
  CHOLMOD/Include/cholmod_complexity.h
  ${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include/cholmod_config.h
  CHOLMOD/Include/cholmod_core.h
  CHOLMOD/Include/cholmod_function.h
  CHOLMOD/Include/cholmod_gpu.h
  CHOLMOD/Include/cholmod_gpu_kernels.h
  CHOLMOD/Include/cholmod_io64.h
  CHOLMOD/Include/cholmod_matrixops.h
  CHOLMOD/Include/cholmod_modify.h
  CHOLMOD/Include/cholmod_partition.h
  CHOLMOD/Include/cholmod_supernodal.h
  CHOLMOD/Include/cholmod_template.h
)

set (CCOLAMD_SRCS
  CCOLAMD/Source/ccolamd.c
)

set (CCOLAMD_HDRS
  CCOLAMD/Include/ccolamd.h
)

set (COLAMD_SRCS
  COLAMD/Source/colamd.c
)

set (COLAMD_HDRS
  COLAMD/Include/colamd.h
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

set (KLU_SRCS
  KLU/Source/klu.c
  KLU/Source/klu_analyze.c
  KLU/Source/klu_analyze_given.c
  KLU/Source/klu_defaults.c
  KLU/Source/klu_diagnostics.c
  KLU/Source/klu_dump.c
  KLU/Source/klu_extract.c
  KLU/Source/klu_factor.c
  KLU/Source/klu_free_numeric.c
  KLU/Source/klu_free_symbolic.c
  KLU/Source/klu_kernel.c
  KLU/Source/klu_memory.c
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

set (LDL_SRCS
  LDL/Source/ldl.c
)

generate_source_files (GEN_AMD_SRCS_DINT FILES "${AMD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_AMD_SRCS_DLONG FILES "${AMD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_BTF_SRCS_DINT FILES "${BTF_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_BTF_SRCS_DLONG FILES "${BTF_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_CAMD_SRCS_DINT FILES "${CAMD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_CAMD_SRCS_DLONG FILES "${CAMD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_CCOLAMD_SRCS_DINT FILES "${CCOLAMD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_CCOLAMD_SRCS_DLONG FILES "${CCOLAMD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_COLAMD_SRCS_DINT FILES "${COLAMD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_COLAMD_SRCS_DLONG FILES "${COLAMD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_CHOLMOD_SRCS_DINT FILES "${CHOLMOD_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_CHOLMOD_SRCS_DLONG FILES "${CHOLMOD_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES "${UMFPACK_INT_SRCS}" SUFFIX
  _di COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES "${UMFPACK_INT_SRCS}" SUFFIX
  _dl COMPILE_DEFINITIONS DLONG)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES "${UMFPACK_ZINT_SRCS}"
  SUFFIX _zi COMPILE_DEFINITIONS ZINT)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES "${UMFPACK_ZINT_SRCS}"
  SUFFIX _zl COMPILE_DEFINITIONS ZLONG)

generate_source_files (GEN_KLU_SRCS_DINT FILES "${KLU_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_KLU_SRCS_DLONG FILES "${KLU_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS DLONG)

generate_source_files (GEN_LDL_SRCS_DINT FILES "${LDL_SRCS}" SUFFIX _di
  COMPILE_DEFINITIONS DINT)
generate_source_files (GEN_LDL_SRCS_DLONG FILES "${LDL_SRCS}" SUFFIX _dl
  COMPILE_DEFINITIONS LDL_LONG)

if (HAVE_COMPLEX_H)
  #list (APPEND PROJECT_SRCS Source/cs_convert.c)
endif (HAVE_COMPLEX_H)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_ltsolve.c UMFPACK/Source/umf_utsolve.c SUFFIX _di_c
  COMPILE_DEFINITIONS DINT CONJUGATE_SOLVE)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_ltsolve.c UMFPACK/Source/umf_utsolve.c SUFFIX _dl_c
  COMPILE_DEFINITIONS DLONG CONJUGATE_SOLVE)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_ltsolve.c UMFPACK/Source/umf_utsolve.c SUFFIX _zi_c
  COMPILE_DEFINITIONS ZINT CONJUGATE_SOLVE)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_ltsolve.c UMFPACK/Source/umf_utsolve.c SUFFIX _zl_c
  COMPILE_DEFINITIONS ZLONG CONJUGATE_SOLVE)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umfpack_solve.c SUFFIX _di_w COMPILE_DEFINITIONS DINT WSOLVE)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umfpack_solve.c SUFFIX _dl_w COMPILE_DEFINITIONS DLONG WSOLVE)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umfpack_solve.c SUFFIX _zi_w COMPILE_DEFINITIONS ZINT WSOLVE)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umfpack_solve.c SUFFIX _zl_w COMPILE_DEFINITIONS ZLONG WSOLVE)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_assemble.c SUFFIX _di_fixq COMPILE_DEFINITIONS DINT FIXQ)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_assemble.c SUFFIX _dl_fixq COMPILE_DEFINITIONS DLONG FIXQ)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_assemble.c SUFFIX _zi_fixq COMPILE_DEFINITIONS ZINT FIXQ)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_assemble.c SUFFIX _zl_fixq COMPILE_DEFINITIONS ZLONG FIXQ)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_store_lu.c SUFFIX _di_drop COMPILE_DEFINITIONS DINT DROP)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_store_lu.c SUFFIX _dl_drop COMPILE_DEFINITIONS DLONG DROP)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_store_lu.c SUFFIX _zi_drop COMPILE_DEFINITIONS ZINT DROP)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_store_lu.c SUFFIX _zl_drop COMPILE_DEFINITIONS ZLONG DROP)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_assemble.c SUFFIX _di_fixq COMPILE_DEFINITIONS DINT FIXQ)

generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _di_map_nox COMPILE_DEFINITIONS DINT DO_MAP)
generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _di_nomap_x COMPILE_DEFINITIONS DINT DO_VALUES)
generate_source_files (GEN_UMFPACK_INT_SRCS_DINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _di_map_x COMPILE_DEFINITIONS DINT DO_MAP DO_VALUES)

generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _dl_map_nox COMPILE_DEFINITIONS DLONG DO_MAP)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _dl_nomap_x COMPILE_DEFINITIONS DLONG DO_VALUES)
generate_source_files (GEN_UMFPACK_INT_SRCS_DLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _dl_map_x COMPILE_DEFINITIONS DLONG
  DO_MAP DO_VALUES)

generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zi_map_nox COMPILE_DEFINITIONS ZINT DO_MAP)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zi_nomap_x COMPILE_DEFINITIONS ZINT DO_VALUES)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZINT FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zi_map_x COMPILE_DEFINITIONS ZINT DO_MAP
  DO_VALUES)

generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zl_map_nox COMPILE_DEFINITIONS ZLONG DO_MAP)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zl_nomap_x COMPILE_DEFINITIONS ZLONG DO_VALUES)
generate_source_files (GEN_UMFPACK_ZINT_SRCS_DZLONG FILES
  UMFPACK/Source/umf_triplet.c SUFFIX _zl_map_x COMPILE_DEFINITIONS ZLONG DO_MAP
  DO_VALUES)

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_WARNINGS)
endif (MSVC)

# output directory for binaries
set (OUTPUT_DIR ${CMAKE_BINARY_DIR})

# put binaries in a single directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR})

add_library (amd
  AMD/Include/amd_internal.h
  ${AMD_HDRS}
  ${GEN_AMD_SRCS_DINT}
  ${GEN_AMD_SRCS_DLONG}
)

add_library (btf
  BTF/Include/btf_internal.h
  ${BTF_HDRS}
  ${GEN_BTF_SRCS_DINT}
  ${GEN_BTF_SRCS_DLONG}
)

add_library (camd
  CAMD/Include/camd_internal.h
  ${CAMD_HDRS}
  ${GEN_CAMD_SRCS_DINT}
  ${GEN_CAMD_SRCS_DLONG}
)

add_library (cholmod
  CHOLMOD/Include/cholmod_internal.h
  ${CHOLMOD_HDRS}
  ${GEN_CHOLMOD_SRCS_DINT}
  ${GEN_CHOLMOD_SRCS_DLONG}
)

target_include_directories (cholmod PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include>
)

if (WITH_OPENMP AND OPENMP_FOUND)
  target_compile_definitions (cholmod PRIVATE HAVE_OPENMP)
  target_compile_options (cholmod PRIVATE ${OpenMP_C_FLAGS})

  if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set_property (TARGET cholmod APPEND PROPERTY LINK_LIBRARIES
      ${OpenMP_C_FLAGS})
    set_property (TARGET cholmod APPEND PROPERTY INTERFACE_LINK_LIBRARIES
      ${OpenMP_C_FLAGS})
  endif (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
endif (WITH_OPENMP AND OPENMP_FOUND)

if (HAVE_FEENABLEEXCEPT)
  target_compile_definitions (cholmod PRIVATE _GNU_SOURCE)
endif (HAVE_FEENABLEEXCEPT)

if (TARGET blas AND TARGET lapack)
  set (LAPACK_DEPENDENCY "find_dependency (LAPACK ${LAPACK_VERSION})")
  target_link_libraries (cholmod PUBLIC blas lapack)
else (TARGET blas AND TARGET lapack)
  target_link_libraries (cholmod PUBLIC ${LAPACK_LIBRARIES} colamd)
endif (TARGET blas AND TARGET lapack)

if (WITH_PARTITION AND METIS_FOUND)
  target_include_directories (cholmod SYSTEM PRIVATE ${METIS_INCLUDE_DIRS})
  target_include_directories (cholmod SYSTEM PUBLIC
    $<BUILD_INTERFACE:${METIS_INCLUDE_DIRS}>)
  target_link_libraries (cholmod PUBLIC ${METIS_LIBRARIES})
endif (WITH_PARTITION AND METIS_FOUND)

add_library (ccolamd
  ${CCOLAMD_HDRS}
  ${GEN_CCOLAMD_SRCS_DINT}
  ${GEN_CCOLAMD_SRCS_DLONG}
)

add_library (colamd
  ${COLAMD_HDRS}
  ${GEN_COLAMD_SRCS_DINT}
  ${GEN_COLAMD_SRCS_DLONG}
)

set (SPQR_SRCS
  SPQR/Source/SuiteSparseQR.cpp
  SPQR/Source/SuiteSparseQR_C.cpp
  SPQR/Source/SuiteSparseQR_expert.cpp
  SPQR/Source/SuiteSparseQR_qmult.cpp
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
)

set (SPQR_HDRS
  SPQR/Include/SuiteSparseQR.hpp
  SPQR/Include/SuiteSparseQR_C.h
  SPQR/Include/SuiteSparseQR_definitions.h
)

add_library (spqr
  ${SPQR_SRCS}
  ${SPQR_HDRS}
)

if (WITH_CUDA AND CUDA_FOUND)
  target_sources (spqr PRIVATE
    SPQR/SPQRGPU/spqrgpu_buildAssemblyMaps.cpp
    SPQR/SPQRGPU/spqrgpu_computeFrontStaging.cpp
    SPQR/SPQRGPU/spqrgpu_kernel.cpp
  )

  target_include_directories (spqr PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries (spqr PUBLIC ${CUDA_LIBRARIES})
  target_link_libraries (spqr PUBLIC suitesparsecuda)
endif (WITH_CUDA AND CUDA_FOUND)

if (WITH_TBB AND TBB_FOUND)
  target_include_directories (spqr SYSTEM PRIVATE ${TBB_INCLUDE_DIRS})
  target_include_directories (spqr SYSTEM PUBLIC
    $<BUILD_INTERFACE:${TBB_INCLUDE_DIRS}>)
  target_link_libraries (spqr PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions (spqr PRIVATE HAVE_TBB)
endif (WITH_TBB AND TBB_FOUND)

set (SUITESPARSECONFIG_HDRS
  SuiteSparse_config/SuiteSparse_config.h
)

add_library (suitesparseconfig
  SuiteSparse_config/SuiteSparse_config.c
  ${SUITESPARSECONFIG_HDRS}
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

add_library (umfpack
  ${GEN_UMFPACK_INT_SRCS_DINT}
  ${GEN_UMFPACK_INT_SRCS_DLONG}
  ${GEN_UMFPACK_ZINT_SRCS_DZINT}
  ${GEN_UMFPACK_ZINT_SRCS_DZLONG}
  UMFPACK/Source/umfpack_global.c
  UMFPACK/Source/umfpack_tictoc.c
  UMFPACK/Source/umfpack_timer.c
  ${UMFPACK_HDRS}
)

# NOTE cholmod already dependens on LAPACK
target_link_libraries (umfpack PUBLIC cholmod)

set (KLU_HDRS
  KLU/Include/klu.h
  KLU/Include/klu_version.h
)

add_library (klu
  KLU/Include/klu_internal.h
  ${GEN_KLU_SRCS_DINT}
  ${GEN_KLU_SRCS_DLONG}
  ${KLU_HDRS}
)

target_link_libraries (klu PUBLIC amd btf colamd)

set (LDL_HDRS
  LDL/Include/ldl.h
)

add_library (ldl
  ${GEN_LDL_SRCS_DINT}
  ${GEN_LDL_SRCS_DLONG}
  ${LDL_HDRS}
)

set (SuiteSparse_AMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/AMD)
set (SuiteSparse_BTF_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/BTF)
set (SuiteSparse_CAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CAMD)
set (SuiteSparse_CCOLAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CCOLAMD)
set (SuiteSparse_CHOLMOD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/CHOLMOD)
set (SuiteSparse_COLAMD_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse/COLAMD)
set (SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse)
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
  target_compile_definitions (cholmod PRIVATE BLAS_UNDERSCORE)
  target_compile_definitions (umfpack PRIVATE BLAS_UNDERSCORE)
  target_compile_definitions (spqr PRIVATE BLAS_UNDERSCORE)
elseif (HAVE_BLAS_NO_UNDERSCORE)
  target_compile_definitions (cholmod PRIVATE BLAS_NO_UNDERSCORE)
  target_compile_definitions (umfpack PRIVATE BLAS_NO_UNDERSCORE)
  target_compile_definitions (spqr PRIVATE BLAS_NO_UNDERSCORE)
endif (HAVE_BLAS_UNDERSCORE)

if (HAVE_M)
  target_link_libraries (suitesparseconfig PUBLIC m)
endif (HAVE_M)

target_include_directories (amd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_AMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_AMD_INCLUDE_DIR}>")
target_include_directories (btf PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_BTF_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_BTF_INCLUDE_DIR}>")
target_include_directories (camd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_CAMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CAMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CAMD_INCLUDE_DIR}>")
target_include_directories (ccolamd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_CCOLAMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CCOLAMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CCOLAMD_INCLUDE_DIR}>")
target_include_directories (colamd PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_COLAMD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_COLAMD_INCLUDE_DIR}>")
target_include_directories (cholmod PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_CHOLMOD_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Cholesky>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Core>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/MatrixOps>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Modify>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Supernodal>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_CHOLMOD_INCLUDE_DIR}>"
)
target_include_directories (spqr PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_SPQR_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_SPQR_INCLUDE_DIR}>")
target_include_directories (umfpack PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_UMFPACK_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Source>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_UMFPACK_INCLUDE_DIR}>"
)
target_include_directories (klu PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_KLU_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_KLU_INCLUDE_DIR}>")
target_include_directories (ldl PUBLIC
  "$<INSTALL_INTERFACE:${SuiteSparse_LDL_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LDL/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_LDL_INCLUDE_DIR}>")

target_link_libraries (amd PUBLIC suitesparseconfig)
target_link_libraries (btf PUBLIC suitesparseconfig)
target_link_libraries (camd PUBLIC suitesparseconfig)
target_link_libraries (ccolamd PUBLIC suitesparseconfig)
target_link_libraries (colamd PUBLIC suitesparseconfig)
target_link_libraries (cholmod PUBLIC amd camd ccolamd colamd
  suitesparseconfig)
target_link_libraries (spqr PUBLIC cholmod suitesparseconfig)
target_link_libraries (klu PUBLIC suitesparseconfig)
target_link_libraries (ldl PUBLIC suitesparseconfig)
target_link_libraries (umfpack PUBLIC amd suitesparseconfig)

#add_library (SuiteSparse ${PROJECT_SRCS})

#target_include_directories (SuiteSparse PUBLIC
#  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>")

if (NOT HAVE_COMPLEX_H)
  #target_compile_definitions (SuiteSparse PUBLIC NCOMPLEX)
endif (NOT HAVE_COMPLEX_H)

if (HAVE_C_DOUBLE_COMPLEX)
  #target_compile_definitions (SuiteSparse PUBLIC HAVE_C_DOUBLE_COMPLEX)
endif (HAVE_C_DOUBLE_COMPLEX)

if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  #target_compile_definitions (SuiteSparse PRIVATE _CRT_SECURE_NO_WARNINGS)
endif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")

# CUDA support
if (WITH_CUDA AND CUDA_FOUND)
  list (APPEND CUDA_NVCC_FLAGS -std=c++11)

  # Workaround for GCC >= 6.2.x (and CUDA 8.0)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
      VERSION_GREATER "6.1")
    list (APPEND CUDA_NVCC_FLAGS -Xcompiler -D__CORRECT_ISO_CPP11_MATH_H_PROTO)
  endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
      VERSION_GREATER "6.1")

  if (HAVE_FPIC)
    list (APPEND CUDA_NVCC_FLAGS -Xcompiler -fPIC)
  endif (HAVE_FPIC)

  cuda_add_library (suitesparsecuda
    SuiteSparse_GPURuntime/Include/SuiteSparseGPU_Runtime.hpp
    SuiteSparse_GPURuntime/Include/SuiteSparseGPU_Workspace.hpp
    SuiteSparse_GPURuntime/Include/SuiteSparseGPU_debug.hpp
    SuiteSparse_GPURuntime/Include/SuiteSparseGPU_macros.hpp
    SuiteSparse_GPURuntime/Include/SuiteSparseGPU_workspace_macros.hpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_cpuAllocators.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_gpuAllocators.cpp
    SuiteSparse_GPURuntime/Source/SuiteSparseGPU_Workspace_transfer.cpp
  )

  target_include_directories (suitesparsecuda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_GPURuntime/Include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${SuiteSparse_GPURuntime_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${SuiteSparse_GPURuntime_INCLUDE_DIR}>
  )

  target_link_libraries (suitesparsecuda PUBLIC suitesparseconfig)

  # GPUQREngine

  cuda_include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}/GPUQREngine/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_GPURuntime/Include
  )

  cuda_add_library (gpuqrengine
    #GPUQREngine/Demo/gpuqrengine_demo.cpp
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

  target_link_libraries (gpuqrengine PUBLIC suitesparsecuda)

  set (CHOLMOD_GPU_SRCS CHOLMOD/GPU/cholmod_gpu.c)

  generate_source_files (GEN_CHOLMOD_GPU_SRCS_DLONG FILES "${CHOLMOD_GPU_SRCS}" SUFFIX _dl
    COMPILE_DEFINITIONS DLONG)

  cuda_wrap_srcs (cholmod OBJ GEN_CHOLMOD_GPU_KERNELS
    CHOLMOD/GPU/cholmod_gpu_kernels.cu)

  target_sources (cholmod PRIVATE
    ${CHOLMOD_GPU_SRCS}
    CHOLMOD/GPU/cholmod_gpu_kernels.cu
    ${GEN_CHOLMOD_GPU_KERNELS}
    ${GEN_CHOLMOD_GPU_SRCS_DLONG}
  )

  target_compile_definitions (cholmod PUBLIC GPU_BLAS)
  target_include_directories (cholmod PUBLIC ${CUDA_INCLUDE_DIRS})
  target_link_libraries (cholmod PUBLIC ${CUDA_LIBRARIES}
    ${CUDA_CUBLAS_LIBRARIES} suitesparsecuda)

  # SPQR GPU
  target_sources (spqr PRIVATE
    SPQR/SPQRGPU/spqrgpu_buildAssemblyMaps.cpp
    SPQR/SPQRGPU/spqrgpu_computeFrontStaging.cpp
    SPQR/SPQRGPU/spqrgpu_kernel.cpp
  )

  target_link_libraries (spqr PUBLIC gpuqrengine)
endif (WITH_CUDA AND CUDA_FOUND)

set (_SuiteSparse_TARGETS amd btf camd ccolamd colamd cholmod spqr umfpack
  klu ldl suitesparseconfig)

set_property (TARGET suitesparseconfig PROPERTY NAME_ALIAS SuiteSparse)
set_property (TARGET suitesparseconfig PROPERTY HEADERS_VARIABLE
  SUITESPARSECONFIG_HDRS)
set_property (TARGET spqr PROPERTY NAME_ALIAS SuiteSparseQR)
set_property (TARGET spqr PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SPQR_INCLUDE_DIR})
set_property (TARGET suitesparseconfig PROPERTY INSTALL_DESTINATION
  ${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR})

if (TARGET suitesparsecuda)
  list (APPEND _SuiteSparse_TARGETS suitesparsecuda)
  set_property (TARGET suitesparsecuda PROPERTY NAME_ALIAS GPURuntime)
endif (TARGET suitesparsecuda)

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

  # Install headers
  install (FILES ${${_HDRS_VAR}} ${_export_header} DESTINATION
    ${_target_include_dir} COMPONENT Development)
endforeach (_target)


set (_SuiteSparse_CMake_LIBDIR ${CMAKE_INSTALL_LIBDIR}/SuiteSparse)
set (_SuiteSparse_CMake_INSTALLDIR ${_SuiteSparse_CMake_LIBDIR}/cmake)

install (TARGETS ${_SuiteSparse_TARGETS}
  EXPORT suitesparse-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  LIBRARY DESTINATION ${_SuiteSparse_CMake_LIBDIR} COMPONENT Development
  ARCHIVE DESTINATION ${_SuiteSparse_CMake_LIBDIR} COMPONENT Development)

set (SuiteSparse_NAMESPACE SuiteSparse::)

configure_package_config_file (suitesparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  INSTALL_DESTINATION ${_SuiteSparse_CMake_INSTALLDIR})
write_basic_package_version_file (suitesparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY SameMajorVersion)

export (TARGETS ${_SuiteSparse_TARGETS} NAMESPACE ${SuiteSparse_NAMESPACE} FILE
  suitesparse-targets.cmake)
export (PACKAGE SuiteSparse)

install (EXPORT suitesparse-targets NAMESPACE ${SuiteSparse_NAMESPACE}
  DESTINATION ${_SuiteSparse_CMake_INSTALLDIR})
install (FILES
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-version.cmake
  DESTINATION ${_SuiteSparse_CMake_INSTALLDIR}
  COMPONENT Development)

set_component_version (amd 2.4.6)
set_component_version (btf 1.2.6)
set_component_version (camd 2.4.6)
set_component_version (ccolamd 2.9.6)
set_component_version (cholmod 3.0.11)
set_component_version (colamd 2.9.6)
set_component_version (klu 1.3.8)
set_component_version (ldl 2.2.6)
set_component_version (spqr 2.0.8)
set_component_version (suitesparseconfig "${VERSION}")
set_component_version (umfpack 5.7.6)

if (TARGET suitesparsecuda)
  set_component_version (suitesparsecuda 1.0.5)
endif (TARGET suitesparsecuda)

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
set_target_properties (cholmod PROPERTIES EXPORT_NAME CHOLMOD)
set_target_properties (colamd PROPERTIES EXPORT_NAME COLAMD)
set_target_properties (klu PROPERTIES EXPORT_NAME KLU)
set_target_properties (ldl PROPERTIES EXPORT_NAME LDL)
set_target_properties (spqr PROPERTIES EXPORT_NAME SPQR)
set_target_properties (suitesparseconfig PROPERTIES EXPORT_NAME Config)
set_target_properties (umfpack PROPERTIES EXPORT_NAME UMFPACK)

set_target_properties (${_SuiteSparse_TARGETS} PROPERTIES
  POSITION_INDEPENDENT_CODE ON)

if (WITH_DEMOS)
  add_executable (cholmod_demo
    CHOLMOD/Demo/cholmod_demo.h
    CHOLMOD/Demo/cholmod_demo.c
  )
  target_link_libraries (cholmod_demo PRIVATE cholmod)

  add_executable (cholmod_l_demo CHOLMOD/Demo/cholmod_l_demo.c)
  target_link_libraries (cholmod_l_demo PRIVATE cholmod)

  add_executable (cholmod_simple CHOLMOD/Demo/cholmod_simple.c)
  target_link_libraries (cholmod_simple PRIVATE cholmod)

  add_executable (umf4 UMFPACK/Demo/umf4.c)
  target_link_libraries (umf4 PRIVATE umfpack)

  add_executable (umfpack_di_demo UMFPACK/Demo/umfpack_di_demo.c)
  target_link_libraries (umfpack_di_demo PRIVATE umfpack)

  add_executable (umfpack_dl_demo UMFPACK/Demo/umfpack_dl_demo.c)
  target_link_libraries (umfpack_dl_demo PRIVATE umfpack)

  add_executable (umfpack_zi_demo UMFPACK/Demo/umfpack_zi_demo.c)
  target_link_libraries (umfpack_zi_demo PRIVATE umfpack)

  add_executable (umfpack_zl_demo UMFPACK/Demo/umfpack_zl_demo.c)
  target_link_libraries (umfpack_zl_demo PRIVATE umfpack)

  add_executable (umfpack_simple UMFPACK/Demo/umfpack_simple.c)
  target_link_libraries (umfpack_simple PRIVATE umfpack)

  if (WITH_FORTRAN AND CMAKE_Fortran_COMPILER)
    add_executable (umf4hb UMFPACK/Demo/umf4_f77wrapper.c UMFPACK/Demo/umf4hb.f)
    target_link_libraries (umf4hb PRIVATE umfpack)
  endif (WITH_FORTRAN AND CMAKE_Fortran_COMPILER)

  add_executable (qrdemo SPQR/Demo/qrdemo.cpp)
  add_executable (qrdemoc SPQR/Demo/qrdemoc.c)
  add_executable (qrsimple SPQR/Demo/qrsimple.cpp)
  add_executable (qrsimplec SPQR/Demo/qrsimplec.c)

  target_link_libraries (qrdemo PRIVATE spqr)
  target_link_libraries (qrdemoc PRIVATE spqr)
  target_link_libraries (qrsimple PRIVATE spqr)
  target_link_libraries (qrsimplec PRIVATE spqr)

  if (TARGET gpuqrengine)
    add_executable (qrdemo_gpu SPQR/Demo/qrdemo_gpu.cpp)
    target_link_libraries (qrdemo_gpu PRIVATE spqr cholmod)
  endif (TARGET gpuqrengine)
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
add_feature_info (Partition WITH_PARTITION "Partition module")
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
