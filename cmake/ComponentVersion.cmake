macro (set_component_version _COMPONENT _VERSION)
  set_target_properties (${_COMPONENT} PROPERTIES VERSION ${_VERSION})
  string (REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" _MAJOR_VERSION
    "${_VERSION}")
  set_target_properties (${_COMPONENT} PROPERTIES SOVERSION ${_MAJOR_VERSION})

  if (WIN32)
    set_target_properties (${_COMPONENT} PROPERTIES RUNTIME_OUTPUT_NAME
      ${_COMPONENT}-${_VERSION})
  endif (WIN32)
endmacro (set_component_version)
