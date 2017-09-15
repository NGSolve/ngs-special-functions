include (ExternalProject)
project(ngs_special_functions_super)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


################################
# Build f2c utility and library
ExternalProject_Add(
  project_f2c
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/f2c
  GIT_REPOSITORY git://github.com/juanjosegarciaripoll/f2c.git
  GIT_TAG d01fe07
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ./autogen.sh && ./configure --prefix=${NGSOLVE_INSTALL_DIR}
  BUILD_COMMAND make -j4
  INSTALL_COMMAND make install
  LOG_DOWNLOAD 0
  LOG_BUILD 0
  LOG_INSTALL 0
  )

################################
# Build actual project
ExternalProject_Add(
  project_ngsolve_specialcfs
  DEPENDS project_f2c
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/specialcf
  CMAKE_ARGS -DNGSolve_DIR=${NGSOLVE_DIR} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  BUILD_ALWAYS 1
  )
