cmake_minimum_required (VERSION 3.1)
project (SuiteSparse)

include (ExternalProject)

ExternalProject_Add (METIS
  URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
  URL_MD5 5465e67079419a69e0116de24fce58fe
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  CMAKE_CACHE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DGKLIB_PATH:PATH=${CMAKE_CURRENT_BINARY_DIR}/METIS-prefix/src/METIS/GKlib
  INSTALL_COMMAND ""
)

ExternalProject_Add (TBB
  URL https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb2017_20160916oss_src.tgz
  URL_MD5 8d15135e0bc71a82fc60446cfd78ae70
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make tbb_build_prefix=external
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
)

set (_PREFIX_PATH
  ${CMAKE_CURRENT_BINARY_DIR}/METIS-prefix/src/METIS-build/libmetis
  ${CMAKE_CURRENT_BINARY_DIR}/TBB-prefix/src/TBB/build/external_debug
  ${CMAKE_CURRENT_BINARY_DIR}/TBB-prefix/src/TBB/build/external_release
)

ExternalProject_Add (SuiteSparse
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  CMAKE_CACHE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DCMAKE_LIBRARY_PATH:STRING=${_PREFIX_PATH}
    -DMETIS_INCLUDE_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/METIS-prefix/src/METIS/include
    -DTBB_INCLUDE_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/TBB-prefix/src/TBB/include
    -DTBB_ROOT_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/TBB-prefix/src/TBB
  INSTALL_COMMAND ""
  DEPENDS METIS TBB
)
