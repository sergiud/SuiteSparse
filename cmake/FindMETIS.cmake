# Module for locating METIS.
#
# Read-only variables:
#   METIS_FOUND
#     Indicates whether the library has been found.
#
#   METIS_INCLUDE_DIRS
#      Specifies METIS' include directory.
#
#   METIS_LIBRARIES
#     Specifies METIS libraries that should be passed to target_link_libararies.
#
# Copyright (c) 2016 Sergiu Deitsch
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

cmake_policy (VERSION 2.8)
cmake_policy (PUSH)

include (FindPackageHandleStandardArgs)

find_path (METIS_INCLUDE_DIR NAMES metis.h
  PATH_SUFFIXES include
  DOC "METIS include directory")
find_library (METIS_LIBRARY_DEBUG NAMES metis
  PATH_SUFFIXES Debug
  DOC "METIS debug library")
find_library (METIS_LIBRARY_RELEASE NAMES metis
  PATH_SUFFIXES Release
  DOC "METIS release library")

if (METIS_LIBRARY_RELEASE)
  if (METIS_LIBRARY_DEBUG)
    set (METIS_LIBRARY debug ${METIS_LIBRARY_RELEASE} optimized
      ${METIS_LIBRARY_DEBUG} CACHE STRING "METIS library")
  else (METIS_LIBRARY_DEBUG)
    set (METIS_LIBRARY ${METIS_LIBRARY_RELEASE} CACHE FILEPATH "METIS library")
  endif (METIS_LIBRARY_DEBUG)
elseif (METIS_LIBRARY_DEBUG)
  set (METIS_LIBRARY ${METIS_LIBRARY_DEBUG} CACHE FILEPATH "METIS library")
endif (METIS_LIBRARY_RELEASE)

set (_METIS_VERSION_HEADER ${METIS_INCLUDE_DIR}/metis.h)

if (EXISTS ${_METIS_VERSION_HEADER})
  file (READ ${_METIS_VERSION_HEADER} _METIS_VERSION_CONTENTS)

  string (REGEX REPLACE ".*#define METIS_VER_MAJOR[ \t]+([0-9]+).*" "\\1"
    METIS_VERSION_MAJOR "${_METIS_VERSION_CONTENTS}")
  string (REGEX REPLACE ".*#define METIS_VER_MINOR[ \t]+([0-9]+).*" "\\1"
    METIS_VERSION_MINOR "${_METIS_VERSION_CONTENTS}")
  string (REGEX REPLACE ".*#define METIS_VER_SUBMINOR[ \t]+([0-9]+).*" "\\1"
    METIS_VERSION_PATCH "${_METIS_VERSION_CONTENTS}")

  set (METIS_VERSION
    ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
  set (METIS_VERSION_COMPONENTS 3)
endif (EXISTS ${_METIS_VERSION_HEADER})

mark_as_advanced (METIS_INCLUDE_DIR METIS_LIBRARY_DEBUG METIS_LIBRARY_RELEASE
  METIS_LIBRARY)

set (METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
set (METIS_LIBRARIES ${METIS_LIBRARY})

find_package_handle_standard_args (METIS REQUIRED_VARS
  METIS_INCLUDE_DIR METIS_LIBRARY VERSION_VAR METIS_VERSION)

cmake_policy (POP)
