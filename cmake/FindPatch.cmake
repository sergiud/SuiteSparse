# Module for locating the patch executable.
#
# Read-only variables:
#   PATCH_FOUND
#     Indicates whether the program has been found.
#
#   PATCH_EXECUTABLE
#     The location of the patch binary.
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

find_program (PATCH_EXECUTABLE NAMES patch DOC "patch executable")
mark_as_advanced (PATCH_EXECUTABLE)

if (PATCH_EXECUTABLE)
  execute_process(COMMAND ${PATCH_EXECUTABLE} -v
    RESULT_VARIABLE _NUMPY_SEARCH_SUCCESS
    OUTPUT_VARIABLE _OUTPUT
    ERROR_VARIABLE _NUMPY_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  string (REGEX MATCH "patch ([0-9]+)\\.([0-9]+)\\.([0-9]+)" _VERSION
    "${_OUTPUT}")

  string (REGEX REPLACE ".*([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1"
    PATCH_VERSION_MAJOR "${_VERSION}")
  string (REGEX REPLACE ".*([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2"
    PATCH_VERSION_MINOR "${_VERSION}")
  string (REGEX REPLACE ".*([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3"
    PATCH_VERSION_PATCH "${_VERSION}")

  set (PATCH_VERSION
    ${PATCH_VERSION_MAJOR}.${PATCH_VERSION_MINOR}.${PATCH_VERSION_PATCH})
endif (PATCH_EXECUTABLE)

find_package_handle_standard_args (Patch REQUIRED_VARS
  PATCH_EXECUTABLE VERSION_VAR PATCH_VERSION)

cmake_policy (POP)
