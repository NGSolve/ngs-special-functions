include (ExternalProject)

# get command line cmake arguments
# MUST be done before call to 'project'
get_cmake_property(vars CACHE_VARIABLES)
foreach(var ${vars})
  get_property(currentHelpString CACHE "${var}" PROPERTY HELPSTRING)
    if("${currentHelpString}" MATCHES "No help, variable specified on the command line." OR "${currentHelpString}" STREQUAL "")
        # message("${var} = [${${var}}]  --  ${currentHelpString}") # uncomment to see the variables being processed
        list(APPEND CL_ARGS "-D${var}=${${var}}")
    endif()
endforeach()

project(ngs_special_functions_super)

cmake_minimum_required(VERSION 3.8)

find_package(NGSolve CONFIG REQUIRED
  HINTS /usr $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
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
  GIT_REPOSITORY https://github.com/mhochsteger/f2c.git
  GIT_TAG 747683f
  CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=${NGSOLVE_INSTALL_DIR}
  -DF2C_INSTALL_DIR_BIN=${NETGEN_INSTALL_DIR_BIN}
  -DF2C_INSTALL_DIR_LIB=${NETGEN_INSTALL_DIR_LIB}
  -DF2C_INSTALL_DIR_INCLUDE=${NETGEN_INSTALL_DIR_INCLUDE}
)

ExternalProject_Add(
  project_fetch_slatec
      DEPENDS project_f2c
      DOWNLOAD_COMMAND ""
      CONFIGURE_COMMAND ${CMAKE_COMMAND}
        -DBIN_DIR=${CMAKE_CURRENT_BINARY_DIR}/specialcf
        -DNGSOLVE_BINARY_DIR=${NGSOLVE_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/src/fetch_and_convert_slatec.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/specialcf/slatec/src
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
)

################################
# Build actual project
ExternalProject_Add(
  project_ngsolve_specialcfs
  DEPENDS project_f2c project_fetch_slatec
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/specialcf
  CMAKE_ARGS -DNGSolve_DIR=${NGSolve_DIR} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${CL_ARGS}
  BUILD_ALWAYS 1
  )
