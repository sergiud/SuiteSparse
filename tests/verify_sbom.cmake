if (NOT DEFINED PROJECT_BINARY_DIR)
  message (FATAL_ERROR "PROJECT_BINARY_DIR is required")
endif (NOT DEFINED PROJECT_BINARY_DIR)

if (NOT DEFINED LICENSE_FILE)
  message (FATAL_ERROR "LICENSE_FILE is required")
endif (NOT DEFINED LICENSE_FILE)

if (NOT DEFINED BUILD_CXSPARSE)
  message (FATAL_ERROR "BUILD_CXSPARSE is required")
endif (NOT DEFINED BUILD_CXSPARSE)

if (NOT DEFINED WITH_LICENSE)
  message (FATAL_ERROR "WITH_LICENSE is required")
endif (NOT DEFINED WITH_LICENSE)

file (READ ${LICENSE_FILE} _target_licenses)

set (_expected_target_licenses
  amd=BSD-3-Clause
  camd=BSD-3-Clause
  ccolamd=BSD-3-Clause
  colamd=BSD-3-Clause
  suitesparseconfig=NOASSERTION
)

if (WITH_LICENSE STREQUAL GPL)
  list (APPEND _expected_target_licenses
    umfpack=GPL-2.0-or-later
    "cholmod=LGPL-2.1-or-later AND GPL-2.0-or-later"
    spqr=GPL-2.0-or-later
  )
elseif (WITH_LICENSE STREQUAL LGPL)
  list (APPEND _expected_target_licenses
    btf=LGPL-2.1-or-later
    klu=LGPL-2.1-or-later
    ldl=LGPL-2.1-or-later
    cholmod=LGPL-2.1-or-later
  )
endif (WITH_LICENSE STREQUAL GPL)

if (BUILD_CXSPARSE)
  list (APPEND _expected_target_licenses cxsparse=LGPL-2.1-or-later)
endif (BUILD_CXSPARSE)

foreach (_expected IN LISTS _expected_target_licenses)
  string (FIND "${_target_licenses}" "${_expected}\n" _license_position)

  if (_license_position EQUAL -1)
    message (FATAL_ERROR
      "Target license metadata is missing or incorrect: ${_expected}")
  endif (_license_position EQUAL -1)
endforeach (_expected)

if (BUILD_CXSPARSE)
  string (FIND "${_target_licenses}" "cxsparse=LGPL-2.1-or-later\n"
    _cxsparse_license_position)

  if (_cxsparse_license_position EQUAL -1)
    message (FATAL_ERROR "CXSparse license metadata is missing or incorrect")
  endif (_cxsparse_license_position EQUAL -1)
endif (BUILD_CXSPARSE)

set (_install_dir ${PROJECT_BINARY_DIR}/sbom-test)
file (REMOVE_RECURSE ${_install_dir})

if (BUILD_CXSPARSE)
  execute_process (
    COMMAND ${CMAKE_COMMAND} -E env CCACHE_DISABLE=1
      ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target cxsparse
    RESULT_VARIABLE _build_result
    OUTPUT_VARIABLE _build_output
    ERROR_VARIABLE _build_error
  )

  if (_build_result)
    message (FATAL_ERROR
      "SBOM test build failed with ${_build_result}: ${_build_error}")
  endif (_build_result)
endif (BUILD_CXSPARSE)

execute_process (
  COMMAND ${CMAKE_COMMAND} --install ${PROJECT_BINARY_DIR}
    --component Unspecified
    --prefix ${_install_dir}
  RESULT_VARIABLE _install_result
  OUTPUT_VARIABLE _install_output
  ERROR_VARIABLE _install_error
)

if (_install_result)
  message (FATAL_ERROR
    "SBOM installation failed with ${_install_result}: ${_install_error}")
endif (_install_result)

file (GLOB_RECURSE _sbom_files ${_install_dir}/*.json)
list (LENGTH _sbom_files _sbom_count)

if (NOT _sbom_count EQUAL 1)
  message (FATAL_ERROR
    "Expected one installed SBOM, found ${_sbom_count}: ${_sbom_files}")
endif (NOT _sbom_count EQUAL 1)

list (GET _sbom_files 0 _sbom_file)
file (READ ${_sbom_file} _sbom)
string (JSON _json_type ERROR_VARIABLE _json_error TYPE "${_sbom}")

if (_json_error OR NOT _json_type STREQUAL OBJECT)
  message (FATAL_ERROR "Installed SBOM is not a JSON object: ${_json_error}")
endif (_json_error OR NOT _json_type STREQUAL OBJECT)

foreach (_value IN ITEMS
  SuiteSparse
  5.13.0
  rdf/3.0.1
  Apache-2.0
  BSD-3-Clause
  http://faculty.cse.tamu.edu/davis/suitesparse.html
)
  string (FIND "${_sbom}" "${_value}" _value_position)

  if (_value_position EQUAL -1)
    message (FATAL_ERROR "SBOM does not contain expected value: ${_value}")
  endif (_value_position EQUAL -1)
endforeach (_value)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
  string (FIND "${_sbom}" "https://github.com/sergiud/SuiteSparse"
    _package_url_position)

  if (_package_url_position EQUAL -1)
    message (FATAL_ERROR
      "SBOM does not contain expected package URL")
  endif (_package_url_position EQUAL -1)
endif (CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)

if (WITH_LICENSE STREQUAL GPL)
  set (_license_values LGPL-2.1-or-later GPL-2.0-or-later)
elseif (WITH_LICENSE STREQUAL LGPL)
  set (_license_values LGPL-2.1-or-later)
else (WITH_LICENSE STREQUAL GPL)
  set (_license_values)
endif (WITH_LICENSE STREQUAL GPL)

foreach (_value IN LISTS _license_values)
  string (FIND "${_sbom}" "${_value}" _value_position)

  if (_value_position EQUAL -1)
    message (FATAL_ERROR "SBOM does not contain expected license: ${_value}")
  endif (_value_position EQUAL -1)
endforeach (_value)
