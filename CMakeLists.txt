cmake_minimum_required (VERSION 2.8.11)

if (POLICY CMP0022)
  cmake_policy (SET CMP0022 NEW)
endif (POLICY CMP0022)

if (POLICY CMP0054)
  cmake_policy (SET CMP0054 NEW)
endif (POLICY CMP0054)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include (CMakeDependentOption)
include (FeatureSummary)

find_package (METIS 5.0)

option (WITH_FORTRAN "Enables Fortran support" OFF)
option (WITH_PRINT "Print diagnostic messages" OFF)
option (WITH_TBB "Enables Intel Threading Building Blocks support" ON)
option (WITH_LGPL "Enable GNU LGPL modules" ON)
cmake_dependent_option (WITH_GPL "Enable GNU GPL modules" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_CAMD "Enable interfaces to CAMD, CCOLAMD, CSYMAMD
in Partition module" ON "WITH_LGPL" OFF)
cmake_dependent_option (WITH_PARTITION "Enable the Partition module" OFF
  "WITH_LGPL AND METIS_FOUND" ON)
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

if (WITH_FORTRAN)
  enable_language (Fortran)
endif (WITH_FORTRAN)

project (SuiteSparse C CXX)

set (MAJOR_VERSION 4)
set (MINOR_VERSION 4)
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
include (CPack)
include (CheckIncludeFile)
include (CheckLibraryExists)
include (CheckTypeSize)
include (GenerateSourceFiles)

find_package (LAPACK REQUIRED)

if (WITH_TBB)
  find_package (TBB 3.0)
endif (WITH_TBB)

check_include_file (complex.h HAVE_COMPLEX_H)
check_library_exists (m ceil "" HAVE_M)

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
  AMD/Include/amd_internal.h
)

set (BTF_SRCS
  BTF/Source/btf_maxtrans.c
  BTF/Source/btf_order.c
  BTF/Source/btf_strongcomp.c
)

set (BTF_HDRS
  BTF/Include/btf.h
  BTF/Include/btf_internal.h
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
  CAMD/Include/camd_internal.h
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

if (NOT WITH_PARTITION)
  set_property (SOURCE CHOLMOD/Partition/cholmod_metis.c
    APPEND PROPERTY COMPILE_DEFINITIONS NPARTITION)
endif (NOT WITH_PARTITION)

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
  CHOLMOD/Include/cholmod_internal.h
  CHOLMOD/Include/cholmod_io64.h
  CHOLMOD/Include/cholmod_matrixops.h
  CHOLMOD/Include/cholmod_modify.h
  CHOLMOD/Include/cholmod_partition.h
  CHOLMOD/Include/cholmod_supernodal.h
  CHOLMOD/Include/cholmod_template.h
)

set (CCOLAMD_SRCS
  CCOLAMD/Source/ccolamd.c
  CCOLAMD/Source/ccolamd_global.c
)

set (CCOLAMD_HDRS
  CCOLAMD/Include/ccolamd.h
)

set (COLAMD_SRCS
  COLAMD/Source/colamd.c
  COLAMD/Source/colamd_global.c
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
  KLU/Include/klu_internal.h
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

add_library (amd
  ${GEN_AMD_SRCS_DINT}
  ${GEN_AMD_SRCS_DLONG}
)

add_library (btf
  ${GEN_BTF_SRCS_DINT}
  ${GEN_BTF_SRCS_DLONG}
)

add_library (camd
  ${GEN_CAMD_SRCS_DINT}
  ${GEN_CAMD_SRCS_DLONG}
)

add_library (cholmod
  ${GEN_CHOLMOD_SRCS_DINT}
  ${GEN_CHOLMOD_SRCS_DLONG}
)

target_include_directories (cholmod PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CHOLMOD/Include>
)

target_compile_definitions (cholmod PUBLIC BLAS_UNDERSCORE)

target_link_libraries (cholmod PUBLIC ${LAPACK_LIBRARIES}
  colamd)

if (WITH_PARTITION AND METIS_FOUND)
  target_include_directories (cholmod SYSTEM PUBLIC ${METIS_INCLUDE_DIR})
  target_link_libraries (cholmod PUBLIC ${METIS_LIBRARIES})
endif (WITH_PARTITION AND METIS_FOUND)

add_library (ccolamd
  ${GEN_CCOLAMD_SRCS_DINT}
  ${GEN_CCOLAMD_SRCS_DLONG}
)

add_library (colamd
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
  SPQR/Include/SuiteSparseQR_C.h
  SPQR/Include/SuiteSparseQR_definitions.h
)

add_library (spqr
  ${SPQR_SRCS}
  ${SPQR_HDRS}
)

set (CONFIG_HDRS
  SuiteSparse_config/SuiteSparse_config.h
)

add_library (suitesparseconfig
  SuiteSparse_config/SuiteSparse_config.c
  ${COLAMD_HDRS}
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

target_include_directories (umfpack PUBLIC ${LAPACK_INCLUDE_DIRS})
target_link_libraries (umfpack PUBLIC ${LAPACK_LIBRARIES} cholmod)

set (KLU_HDRS
  KLU/Include/klu.h
  KLU/Include/klu_internal.h
  KLU/Include/klu_version.h
)

add_library (klu
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

target_include_directories (suitesparseconfig PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config>")

if (HAVE_M)
  target_link_libraries (suitesparseconfig PUBLIC m)
endif (HAVE_M)

if (WITH_TBB AND TBB_FOUND)
  target_include_directories (spqr PUBLIC ${TBB_INCLUDE_DIRS})
  target_link_libraries (spqr PUBLIC ${TBB_LIBRARIES})
  target_compile_definitions (spqr PUBLIC WITH_TBB)
endif (WITH_TBB AND TBB_FOUND)

target_include_directories (amd PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/AMD>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include>")
target_include_directories (btf PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/BTF>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include>")
target_include_directories (camd PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/CAMD>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CAMD/Include>")
target_include_directories (ccolamd PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/CCOLAMD>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CCOLAMD/Include>")
target_include_directories (colamd PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/COLAMD>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include>")
target_include_directories (cholmod PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/CHOLMOD>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Cholesky>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Core>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/MatrixOps>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Modify>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Supernodal>"
)
target_include_directories (spqr PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/SPQR>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SPQR/Include>")
target_include_directories (umfpack PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/UMFPACK>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Source>"
)
target_include_directories (klu PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/KLU>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include>")
target_include_directories (ldl PUBLIC
  "$<INSTALL_INTERFACE:include/SuiteSparse/LDL>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LDL/Include>")

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

install (TARGETS amd btf camd ccolamd colamd cholmod spqr umfpack klu ldl
  suitesparseconfig
  EXPORT suitesparse-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

set (SuiteSparse_AMD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/AMD/Include)
set (SuiteSparse_BTF_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BTF/Include)
set (SuiteSparse_CAMD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CAMD/Include)
set (SuiteSparse_CCOLAMD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CCOLAMD/Include)
set (SuiteSparse_CHOLMOD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CHOLMOD/Include)
set (SuiteSparse_COLAMD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/COLAMD/Include)
set (SuiteSparse_CONFIG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse_config)
set (SuiteSparse_KLU_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/KLU/Include)
set (SuiteSparse_SPQR_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SPQR/Include)
set (SuiteSparse_UMFPACK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/UMFPACK/Include)

configure_package_config_file (suitesparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config.cmake
  INSTALL_DESTINATION lib/cmake/SuiteSparse
  PATH_VARS SuiteSparse_AMD_INCLUDE_DIR SuiteSparse_BTF_INCLUDE_DIR
  SuiteSparse_CAMD_INCLUDE_DIR SuiteSparse_CCOLAMD_INCLUDE_DIR
  SuiteSparse_CHOLMOD_INCLUDE_DIR SuiteSparse_COLAMD_INCLUDE_DIR
  SuiteSparse_CONFIG_INCLUDE_DIR SuiteSparse_KLU_INCLUDE_DIR
  SuiteSparse_SPQR_INCLUDE_DIR SuiteSparse_UMFPACK_INCLUDE_DIR)
write_basic_package_version_file (suitesparse-config-version.cmake
  VERSION ${VERSION} COMPATIBILITY SameMajorVersion)

set (SuiteSparse_AMD_INCLUDE_DIR include/SuiteSparse/AMD)
set (SuiteSparse_BTF_INCLUDE_DIR include/SuiteSparse/BTF)
set (SuiteSparse_CAMD_INCLUDE_DIR include/SuiteSparse/CAMD)
set (SuiteSparse_CCOLAMD_INCLUDE_DIR include/SuiteSparse/CCOLAMD)
set (SuiteSparse_CHOLMOD_INCLUDE_DIR include/SuiteSparse/CHOLMOD)
set (SuiteSparse_COLAMD_INCLUDE_DIR include/SuiteSparse/COLAMD)
set (SuiteSparse_CONFIG_INCLUDE_DIR include/SuiteSparse)
set (SuiteSparse_KLU_INCLUDE_DIR include/SuiteSparse/KLU)
set (SuiteSparse_SPQR_INCLUDE_DIR include/SuiteSparse/SPQR)
set (SuiteSparse_UMFPACK_INCLUDE_DIR include/SuiteSparse/UMFPACK)

# install version of suitesparse-config.cmake
configure_package_config_file (suitesparse-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-install.cmake
  INSTALL_DESTINATION lib/cmake/SuiteSparse
  PATH_VARS SuiteSparse_AMD_INCLUDE_DIR SuiteSparse_BTF_INCLUDE_DIR
  SuiteSparse_CAMD_INCLUDE_DIR SuiteSparse_CCOLAMD_INCLUDE_DIR
  SuiteSparse_CHOLMOD_INCLUDE_DIR SuiteSparse_COLAMD_INCLUDE_DIR
  SuiteSparse_CONFIG_INCLUDE_DIR SuiteSparse_KLU_INCLUDE_DIR
  SuiteSparse_SPQR_INCLUDE_DIR SuiteSparse_UMFPACK_INCLUDE_DIR)

export (TARGETS amd btf camd ccolamd colamd cholmod spqr umfpack klu ldl
  suitesparseconfig FILE suitesparse-targets.cmake)
export (PACKAGE SuiteSparse)

install (FILES ${AMD_HDRS} DESTINATION ${SuiteSparse_AMD_INCLUDE_DIR})
install (FILES ${BTF_HDRS} DESTINATION ${SuiteSparse_BTF_INCLUDE_DIR})
install (FILES ${CAMD_HDRS} DESTINATION ${SuiteSparse_CAMD_INCLUDE_DIR})
install (FILES ${CCOLAMD_HDRS} DESTINATION ${SuiteSparse_CCOLAMD_INCLUDE_DIR})
install (FILES ${CHOLMOD_HDRS} DESTINATION ${SuiteSparse_CHOLMOD_INCLUDE_DIR})
install (FILES ${COLAMD_HDRS} DESTINATION ${SuiteSparse_COLAMD_INCLUDE_DIR})
install (FILES ${CONFIG_HDRS} DESTINATION ${SuiteSparse_CONFIG_INCLUDE_DIR})
install (FILES ${KLU_HDRS} DESTINATION ${SuiteSparse_KLU_INCLUDE_DIR})
install (FILES ${SPQR_HDRS} DESTINATION ${SuiteSparse_SPQR_INCLUDE_DIR})
install (FILES ${UMFPACK_HDRS} DESTINATION ${SuiteSparse_UMFPACK_INCLUDE_DIR})

install (EXPORT suitesparse-targets DESTINATION lib/cmake/SuiteSparse)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-install.cmake
  DESTINATION lib/cmake/SuiteSparse RENAME suitesparse-config.cmake)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/suitesparse-config-version.cmake
  DESTINATION lib/cmake/SuiteSparse
  COMPONENT Development)

macro (set_component_version _COMPONENT _VERSION)
  set_target_properties (${_COMPONENT} PROPERTIES VERSION ${_VERSION})
  string (REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" _MAJOR_VERSION
    "${_VERSION}")
  set_target_properties (${_COMPONENT} PROPERTIES SOVERSION ${_MAJOR_VERSION})

  if (WIN32)
    set_target_properties (${_COMPONENT} PROPERTIES RUNTIME_OUTPUT_NAME
      ${_COMPONENT}-${_VERSION})
    set_target_properties (${_COMPONENT} PROPERTIES RUNTIME_OUTPUT_NAME_DEBUG
      ${_COMPONENT}-debug-${_VERSION})
  endif (WIN32)
endmacro (set_component_version)

set_component_version (amd 2.4.1)
set_component_version (btf 1.2.1)
set_component_version (camd 2.4.1)
set_component_version (ccolamd 2.9.1)
set_component_version (cholmod 3.0.6)
set_component_version (klu 1.3.3)
set_component_version (ldl 2.2.1)
set_component_version (spqr 2.0.1)
set_component_version (suitesparseconfig 4.4.5)
set_component_version (umfpack 5.7.1)

set_target_properties (amd btf camd ccolamd colamd cholmod spqr umfpack klu ldl
  suitesparseconfig PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable (umf4
  UMFPACK/Demo/umf4.c
)

target_link_libraries (umf4 PRIVATE umfpack)

add_executable (umfpack_di_demo
  UMFPACK/Demo/umfpack_di_demo.c
)

target_link_libraries (umfpack_di_demo PRIVATE umfpack)

add_executable (umfpack_dl_demo
  UMFPACK/Demo/umfpack_dl_demo.c
)

target_link_libraries (umfpack_dl_demo PRIVATE umfpack)

add_executable (umfpack_zi_demo
  UMFPACK/Demo/umfpack_zi_demo.c
)

target_link_libraries (umfpack_zi_demo PRIVATE umfpack)

add_executable (umfpack_zl_demo
  UMFPACK/Demo/umfpack_zl_demo.c
)

target_link_libraries (umfpack_zl_demo PRIVATE umfpack)

add_executable (umfpack_simple
  UMFPACK/Demo/umfpack_simple.c
)

target_link_libraries (umfpack_simple PRIVATE umfpack)


if (WITH_FORTRAN)
  add_executable (umf4hb
    UMFPACK/Demo/umf4_f77wrapper.c
    UMFPACK/Demo/umf4hb.f
  )

  target_link_libraries (umf4hb PRIVATE umfpack)
endif (WITH_FORTRAN)

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
add_feature_info (Threading "WITH_TBB AND TBB_FOUND" "Intel Threading Building Blocks support")

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

set_package_properties (GPL PROPERTIES
  PURPOSE "Enables modules with the GNU General Public License")
set_package_properties (LGPL PROPERTIES
  PURPOSE "Enables modules with the GNU Lesser General Public License")
set_package_properties (Partition PROPERTIES
  PURPOSE "Enables the partitioning module requiring METIS")

set_feature_info (BLAS "Basic Linear Algebra Subprograms"
  "http://www.netlib.org/blas")
set_feature_info (LAPACK "Linear Algebra PACKage"
  "http://www.netlib.org/lapack")

feature_summary (WHAT ALL)
