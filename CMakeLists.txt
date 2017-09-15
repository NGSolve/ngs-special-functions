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
  GIT_REPOSITORY git://github.com/mhochsteger/f2c.git
  GIT_TAG bbb7423945ff224f1d69ee3d267dcad6859972f0
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${NGSOLVE_INSTALL_DIR}
)

################################
# Build actual project
ExternalProject_Add(
  project_ngsolve_specialcfs
  DEPENDS project_f2c
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/specialcf
  CMAKE_ARGS -DNGSolve_DIR=${NGSolve_DIR} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  BUILD_ALWAYS 1
  )
