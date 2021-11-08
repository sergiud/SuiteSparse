#
# CMake support layer for SuiteSparse
#
# Copyright 2021 Sergiu Deitsch <sergiu.deitsch@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if (TARGET suitesparsebase)
  return ()
endif (TARGET suitesparsebase)

set (SUITESPARSECONFIG_HDRS
  ${SuiteSparse_SOURCE_DIR}/SuiteSparse_config/SuiteSparse_config.h
)

add_library (suitesparsebase OBJECT
  ${SuiteSparse_SOURCE_DIR}/SuiteSparse_config/SuiteSparse_config.c
  ${SUITESPARSECONFIG_HDRS}
)

set (SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/SuiteSparse)

target_include_directories (suitesparsebase PRIVATE
  ${SuiteSparse_SOURCE_DIR}/SuiteSparse_config
  ${SuiteSparse_BINARY_DIR}/${SuiteSparse_SUITESPARSECONFIG_INCLUDE_DIR}
)
