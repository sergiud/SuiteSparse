if (NOT DEFINED PROJECT_SOURCE_DIR)
  message (FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif (NOT DEFINED PROJECT_SOURCE_DIR)

if (NOT DEFINED PROJECT_BINARY_DIR)
  message (FATAL_ERROR "PROJECT_BINARY_DIR is required")
endif (NOT DEFINED PROJECT_BINARY_DIR)

set (_common_arguments
  -DBUILD_TESTING=OFF
  -DWITH_CUDA=OFF
  -DWITH_DEMOS=OFF
  -DWITH_FORTRAN=OFF
  -DWITH_METIS=OFF
  -DWITH_OPENMP=OFF
  -DWITH_TBB=OFF
)

foreach (_license IN ITEMS Minimal LGPL GPL)
  set (_build_dir ${PROJECT_BINARY_DIR}/license-options-${_license})
  file (REMOVE_RECURSE ${_build_dir})

  set (_license_arguments -DWITH_LICENSE=${_license})

  execute_process (
    COMMAND ${CMAKE_COMMAND} -E env CCACHE_DISABLE=1
      ${CMAKE_COMMAND} -S ${PROJECT_SOURCE_DIR} -B ${_build_dir}
      ${_license_arguments} ${_common_arguments}
    RESULT_VARIABLE _configure_result
    OUTPUT_VARIABLE _configure_output
    ERROR_VARIABLE _configure_error
  )

  if (_configure_result)
    message (FATAL_ERROR
      "${_license} configuration failed with ${_configure_result}: "
      "${_configure_error}")
  endif (_configure_result)

  file (READ ${_build_dir}/suitesparse-targets.cmake _exports)
  file (READ ${_build_dir}/CMakeCache.txt _cache)

  foreach (_redundant_option IN ITEMS WITH_GPL WITH_LGPL)
    string (FIND "${_cache}" "${_redundant_option}:" _option_position)

    if (NOT _option_position EQUAL -1)
      message (FATAL_ERROR
        "${_license} configuration exposes redundant option "
        "${_redundant_option}")
    endif (NOT _option_position EQUAL -1)
  endforeach (_redundant_option)

  if (_license STREQUAL Minimal)
    if (IS_DIRECTORY ${_build_dir}/CXSparse)
      message (FATAL_ERROR
        "${_license} configuration enabled the LGPL CXSparse target")
    endif (IS_DIRECTORY ${_build_dir}/CXSparse)
  else (_license STREQUAL Minimal)
    if (NOT IS_DIRECTORY ${_build_dir}/CXSparse)
      message (FATAL_ERROR
        "${_license} configuration disabled the permitted CXSparse target")
    endif (NOT IS_DIRECTORY ${_build_dir}/CXSparse)
  endif (_license STREQUAL Minimal)

  if (_license STREQUAL Minimal)
    set (_expected_targets AMD CAMD CCOLAMD COLAMD Config)
    set (_forbidden_targets BTF CHOLMOD KLU LDL UMFPACK SPQR)
  elseif (_license STREQUAL LGPL)
    set (_expected_targets
      AMD BTF CAMD CCOLAMD CHOLMOD COLAMD KLU LDL Config)
    set (_forbidden_targets UMFPACK SPQR)
  else (_license STREQUAL Minimal)
    set (_expected_targets
      AMD BTF CAMD CCOLAMD CHOLMOD COLAMD KLU LDL Config UMFPACK SPQR)
    set (_forbidden_targets)
  endif (_license STREQUAL Minimal)

  foreach (_target IN LISTS _expected_targets)
    string (FIND "${_exports}" "SuiteSparse::${_target}" _target_position)

    if (_target_position EQUAL -1)
      message (FATAL_ERROR
        "${_license} configuration is missing target ${_target}")
    endif (_target_position EQUAL -1)
  endforeach (_target)

  foreach (_target IN LISTS _forbidden_targets)
    string (FIND "${_exports}" "SuiteSparse::${_target}" _target_position)

    if (NOT _target_position EQUAL -1)
      message (FATAL_ERROR
        "${_license} configuration contains forbidden target ${_target}")
    endif (NOT _target_position EQUAL -1)
  endforeach (_target)
endforeach (_license)

set (_demo_build_dir ${PROJECT_BINARY_DIR}/license-options-demos)
file (REMOVE_RECURSE ${_demo_build_dir})

execute_process (
  COMMAND ${CMAKE_COMMAND} -E env CCACHE_DISABLE=1
    ${CMAKE_COMMAND} -S ${PROJECT_SOURCE_DIR} -B ${_demo_build_dir}
    -DWITH_LICENSE=LGPL ${_common_arguments}
    -DWITH_DEMOS=ON
  RESULT_VARIABLE _configure_result
  OUTPUT_VARIABLE _configure_output
  ERROR_VARIABLE _configure_error
)

if (_configure_result)
  message (FATAL_ERROR
    "LGPL demo configuration failed with ${_configure_result}: "
    "${_configure_error}")
endif (_configure_result)

file (READ ${_demo_build_dir}/build.ninja _demo_build)
string (FIND "${_demo_build}"
  "CMakeFiles/ldllsimple.dir/LDL/Demo/ldlsimple.c.o:"
  _demo_source_position)

if (_demo_source_position EQUAL -1)
  message (FATAL_ERROR
    "The ldllsimple target does not use LDL/Demo/ldlsimple.c")
endif (_demo_source_position EQUAL -1)
