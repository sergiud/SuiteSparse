include (CMakeParseArguments)

function (generate_source_files PROJECT_SRCS ARGS)
  set (GSF_SINGLE SUFFIX)
  set (GSF_MULTI FILES COMPILE_DEFINITIONS)

  cmake_parse_arguments (GSF_ARGS "" "${GSF_SINGLE}" "${GSF_MULTI}" ${ARGS}
      ${ARGN})

  set (_GEN_FILES ${${PROJECT_SRCS}})

  foreach (_SRC ${GSF_ARGS_FILES})
    get_filename_component (_ABS_PATH ${_SRC} ABSOLUTE)
    file (RELATIVE_PATH _RELPATH ${CMAKE_CURRENT_SOURCE_DIR} ${_ABS_PATH})
    get_filename_component (_PATH ${_RELPATH} DIRECTORY)

    get_filename_component (_NAME ${_SRC} NAME_WE)
    get_filename_component (_EXT ${_SRC} EXT)

    set (_NEW_FILE ${CMAKE_CURRENT_BINARY_DIR}/${_PATH}/${_NAME}${GSF_ARGS_SUFFIX}${_EXT})

    add_custom_command (OUTPUT ${_NEW_FILE} COMMAND ${CMAKE_COMMAND} -E copy
      ${_ABS_PATH} ${_NEW_FILE} DEPENDS ${_SRC})

    get_property (_COMPILE_DEFINITIONS SOURCE ${_SRC} PROPERTY
      COMPILE_DEFINITIONS)

    set_property (SOURCE ${_NEW_FILE} APPEND PROPERTY COMPILE_DEFINITIONS
      "${GSF_ARGS_COMPILE_DEFINITIONS};${_COMPILE_DEFINITIONS}")

    list (APPEND _GEN_FILES ${_NEW_FILE})
    # Add original source file to the generate list of files. However, exclude
    # it from compilation.
    list (APPEND _GEN_FILES ${_SRC})

    set_source_files_properties (${_SRC} PROPERTIES HEADER_FILE_ONLY ON)
  endforeach (_SRC)

  set (${PROJECT_SRCS} ${_GEN_FILES} PARENT_SCOPE)
endfunction (generate_source_files)
