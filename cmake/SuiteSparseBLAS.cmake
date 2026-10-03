# SPDX-FileCopyrightText: 2026 Sergiu Deitsch
# SPDX-License-Identifier: Apache-2.0

# Detect the integer ABI by executing routines from the linked libraries.
function (suitesparse_detect_blas_integer_size)
  if (CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    return ()
  endif (CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)

  foreach (_blas_symbol dcopy_ dcopy)
    foreach (_lapack_symbol dgetrf_ dgetrf)
      unset (_SuiteSparse_BLAS_RUN_RESULT CACHE)
      unset (_SuiteSparse_BLAS_COMPILE_RESULT CACHE)
      try_run (_SuiteSparse_BLAS_RUN_RESULT _SuiteSparse_BLAS_COMPILE_RESULT
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/SuiteSparseBLAS
        ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/check_blas_integer.c
        COMPILE_DEFINITIONS
          -DSUITESPARSE_DCOPY=${_blas_symbol}
          -DSUITESPARSE_DGETRF=${_lapack_symbol}
        LINK_LIBRARIES LAPACK::LAPACK BLAS::BLAS
        C_STANDARD 99
        C_STANDARD_REQUIRED TRUE
        COMPILE_OUTPUT_VARIABLE _SuiteSparse_BLAS_COMPILE_OUTPUT
        RUN_OUTPUT_VARIABLE _SuiteSparse_BLAS_RUN_OUTPUT)
      if (_SuiteSparse_BLAS_COMPILE_RESULT)
        if (NOT _SuiteSparse_BLAS_RUN_RESULT MATCHES "^0$" OR
            NOT _SuiteSparse_BLAS_RUN_OUTPUT MATCHES "^[48]$")
          message (FATAL_ERROR "BLAS/LAPACK integer ABI probe failed: "
            "actual exit ${_SuiteSparse_BLAS_RUN_RESULT}, output "
            "${_SuiteSparse_BLAS_RUN_OUTPUT}, expected exit 0 and integer size 4 or 8")
        endif (NOT _SuiteSparse_BLAS_RUN_RESULT MATCHES "^0$" OR
            NOT _SuiteSparse_BLAS_RUN_OUTPUT MATCHES "^[48]$")
        if (DEFINED BLA_SIZEOF_INTEGER AND
            NOT BLA_SIZEOF_INTEGER EQUAL _SuiteSparse_BLAS_RUN_OUTPUT)
          message (FATAL_ERROR "BLAS integer size mismatch: "
            "actual ${_SuiteSparse_BLAS_RUN_OUTPUT}, expected ${BLA_SIZEOF_INTEGER}")
        endif (DEFINED BLA_SIZEOF_INTEGER AND
            NOT BLA_SIZEOF_INTEGER EQUAL _SuiteSparse_BLAS_RUN_OUTPUT)
        set (BLA_SIZEOF_INTEGER ${_SuiteSparse_BLAS_RUN_OUTPUT} PARENT_SCOPE)
        return ()
      endif (_SuiteSparse_BLAS_COMPILE_RESULT)
    endforeach (_lapack_symbol)
  endforeach (_blas_symbol)

  message (FATAL_ERROR "BLAS/LAPACK integer ABI probe did not compile: "
    "actual compiler output:\n${_SuiteSparse_BLAS_COMPILE_OUTPUT}\n"
    "expected callable dcopy and dgetrf routines")
endfunction (suitesparse_detect_blas_integer_size)
