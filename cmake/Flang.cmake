macro(windows_compiler_clang_abi lang)
  if(CMAKE_NO_GNUtoMS)
    set(CMAKE_GNUtoMS 0)
  else()
    option(CMAKE_GNUtoMS "Convert GNU import libraries to MS format (requires Visual Studio)" OFF)
  endif()

  message (STATUS "FOO: ${CMAKE_SYSTEM_PROCESSOR}")

  unset (CMAKE_GNUtoMS_VCVARS CACHE)
  unset (CMAKE_GNUtoMS_LIB)
  set (CMAKE_GNUtoMS 1)

  if(CMAKE_GNUtoMS AND NOT CMAKE_GNUtoMS_LIB)
    # Find MS development environment setup script for this architecture.
    # We need to use the MS Librarian tool (lib.exe).
    # Find the most recent version available.

    # Query the VS Installer tool for locations of VS 2017 and above.
    set(_vs_installer_paths "")
    foreach(vs RANGE 17 15 -1) # change the first number to the largest supported version
      cmake_host_system_information(RESULT _vs_dir QUERY VS_${vs}_DIR)
      if(_vs_dir)
        list(APPEND _vs_installer_paths "${_vs_dir}/VC/Auxiliary/Build")

        file (GLOB _vcvars "${_vs_dir}/VC/Auxiliary/Build/*.bat")

        foreach (vc IN LISTS _vcvars)
          message (STATUS "VC: ${vc}")
        endforeach (vc)
      endif()
    endforeach()

    message (STATUS "${CMAKE_SIZEOF_VOID_P}")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL 4)
      message (STATUS "x86")
      find_program(CMAKE_GNUtoMS_VCVARS NAMES vcvars32.bat
        DOC "Visual Studio vcvars32.bat"
        PATHS
        ${_vs_installer_paths}
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\14.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\12.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\11.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\7.1\\Setup\\VC;ProductDir]/bin"
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\6.0\\Setup\\Microsoft Visual C++;ProductDir]/bin"
        )
      set(CMAKE_GNUtoMS_ARCH x86)
    elseif("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
      if(CMAKE_SYSTEM_PROCESSOR STREQUAL AMD64)
        message (STATUS "AMD64")
        find_program(CMAKE_GNUtoMS_VCVARS NAMES vcvars64.bat vcvarsamd64.bat
          DOC "Visual Studio vcvarsamd64.bat"
          PATHS
          ${_vs_installer_paths}
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\14.0\\Setup\\VC;ProductDir]/bin/amd64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\12.0\\Setup\\VC;ProductDir]/bin/amd64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\11.0\\Setup\\VC;ProductDir]/bin/amd64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir]/bin/amd64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC;ProductDir]/bin/amd64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC;ProductDir]/bin/amd64"
          )
        set(CMAKE_GNUtoMS_ARCH amd64)
      elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL ARM64)
        message (STATUS "ARM64")
        find_program(CMAKE_GNUtoMS_VCVARS NAMES vcvarsarm64.bat
          DOC "Visual Studio vcvarsarm64.bat"
          PATHS
          ${_vs_installer_paths}
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\14.0\\Setup\\VC;ProductDir]/bin/arm64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\12.0\\Setup\\VC;ProductDir]/bin/arm64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\11.0\\Setup\\VC;ProductDir]/bin/arm64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir]/bin/arm64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0\\Setup\\VC;ProductDir]/bin/arm64"
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC;ProductDir]/bin/arm64"
          )
        set(CMAKE_GNUtoMS_ARCH arm64)
      endif()
    endif()
    unset(_vs_installer_paths)
    set_property(CACHE CMAKE_GNUtoMS_VCVARS PROPERTY ADVANCED 1)
    if(CMAKE_GNUtoMS_VCVARS)
      # Create helper script to run lib.exe from MS environment.
      string(REPLACE "/" "\\" CMAKE_GNUtoMS_BAT "${CMAKE_GNUtoMS_VCVARS}")
      set(CMAKE_GNUtoMS_LIB ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeGNUtoMS_lib.bat)
      configure_file(${CMAKE_ROOT}/Modules/Platform/GNUtoMS_lib.bat.in ${CMAKE_GNUtoMS_LIB})
    else()
      message(WARNING "Disabling CMAKE_GNUtoMS option because CMAKE_GNUtoMS_VCVARS is not set.")
      set(CMAKE_GNUtoMS 0)
    endif()
  endif()

  message (STATUS "VCVARS: ${CMAKE_GNUtoMS_VCVARS} ${CMAKE_GNUtoMS_ARCH}")
  message (STATUS "RULE: ${lang} ${CMAKE_${lang}_GNUtoMS_RULE}")

  if(CMAKE_GNUtoMS)
    # Teach CMake how to create a MS import library at link time.
    set(CMAKE_${lang}_GNUtoMS_RULE " -Wl,--output-def,<TARGET_NAME>.def"
      "<CMAKE_COMMAND> -Dlib=\"${CMAKE_GNUtoMS_LIB}\" -Ddef=<TARGET_NAME>.def -Ddll=<TARGET> -Dimp=<TARGET_IMPLIB> -P \"${CMAKE_ROOT}/Modules/Platform/GNUtoMS_lib.cmake\""
      )
  endif()
  execute_process (COMMAND ${CMAKE_GNUtoMS_VCVARS})

endmacro()

if (WIN32)
  if (CMAKE_C_COMPILER_LOADED)
    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
      windows_compiler_clang_abi (C)
    endif (CMAKE_C_COMPILER_ID MATCHES "Clang")
  endif (CMAKE_C_COMPILER_LOADED)

  if (CMAKE_CXX_COMPILER_LOADED)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      windows_compiler_clang_abi (CXX)
    endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  endif (CMAKE_CXX_COMPILER_LOADED)

  if (CMAKE_Fortran_COMPILER_LOADED)
    message (STATUS "Compiler: ${CMAKE_Fortran_COMPILER_ID}")

    if (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
      windows_compiler_clang_abi (Fortran)
    endif (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  endif (CMAKE_Fortran_COMPILER_LOADED)
endif (WIN32)

