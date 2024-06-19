cmake_minimum_required(VERSION 3.28)

include(FetchContent)

project(ngsolve_special_functions)
include(ngsolve_addon.cmake)

# Fetch f2c utility and library
FetchContent_Declare(
  f2c
  GIT_REPOSITORY https://github.com/mhochsteger/f2c.git
  GIT_TAG 96c2fbf
  EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(f2c)

if(UNIX)
  target_compile_options(f2c_exe PRIVATE -Wno-implicit-function-declaration)
endif()

get_target_property(F2C_DIR f2c_exe BINARY_DIR)
if(WIN32)
  set(F2C_COMMAND ${F2C_DIR}/Release/f2c.exe)
else()
  set(F2C_COMMAND ${F2C_DIR}/f2c)
endif()

set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/specialcf)

macro( fetch_slatec_function slatec_function )
  file(DOWNLOAD http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=slatec%2Fsrc%2F${slatec_function}.f ${BIN_DIR}/${slatec_function}.tgz)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ${BIN_DIR}/${slatec_function}.tgz WORKING_DIRECTORY ${BIN_DIR})
endmacro( fetch_slatec_function )

fetch_slatec_function(zbesi)
fetch_slatec_function(zbesj)
fetch_slatec_function(zbesk)
fetch_slatec_function(gamln)
fetch_slatec_function(zbesh)

file(DOWNLOAD http://www.netlib.org/blas/i1mach.f ${BIN_DIR}/slatec/src/i1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/r1mach.f ${BIN_DIR}/slatec/src/r1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/d1mach.f ${BIN_DIR}/slatec/src/d1mach.f)

file(DOWNLOAD http://ab-initio.mit.edu/Faddeeva.cc ${BIN_DIR}/Faddeeva/Faddeeva.cc)
file(DOWNLOAD http://ab-initio.mit.edu/Faddeeva.hh ${BIN_DIR}/Faddeeva/Faddeeva.hh)

file(GLOB SLATEC_SOURCES_FORTRAN ${BIN_DIR}/slatec/src/*.f )

set(SLATEC_SOURCES_C)
foreach(FORTRAN_FILE ${SLATEC_SOURCES_FORTRAN})
    string(REGEX REPLACE "\\.f$" ".c" C_FILE ${FORTRAN_FILE})
    list(APPEND SLATEC_SOURCES_C ${C_FILE})
endforeach()
message(STATUS "C files: ${SLATEC_SOURCES_C}")

add_custom_target(convert_sources ALL
  DEPENDS ${SLATEC_SOURCES_FORTRAN} f2c_exe
  COMMAND ${CMAKE_COMMAND} -E make_directory ${BIN_DIR}/slatec/src
  COMMAND ${F2C_COMMAND} -a ${SLATEC_SOURCES_FORTRAN} WORKING_DIRECTORY ${BIN_DIR}/slatec/src
  BYPRODUCTS ${SLATEC_SOURCES_C}
  USES_TERMINAL
)

add_ngsolve_addon(ngsolve_special_functions src/specialcf.cpp ${SLATEC_SOURCES_C}  ${CMAKE_CURRENT_BINARY_DIR}/specialcf/Faddeeva/Faddeeva.cc)
add_custom_target(ngsolve_special_functions_convert ALL DEPENDS convert_sources)
target_sources(ngsolve_special_functions PRIVATE ${SLATEC_SOURCES_C})
target_include_directories(ngsolve_special_functions PRIVATE ${BIN_DIR}/slatec/src ${BIN_DIR}/Faddeeva inc ${f2c_SOURCE_DIR}/include ${f2c_BINARY_DIR})
target_compile_definitions(ngsolve_special_functions PRIVATE -DSPECIALCF_LIBRARY_NAME="$<TARGET_LINKER_FILE_NAME:ngsolve_special_functions>")
if(UNIX)
  target_compile_options(ngsolve_special_functions PRIVATE -Wno-implicit-function-declaration)
endif()
target_link_libraries(ngsolve_special_functions PRIVATE f2c)
install(TARGETS ngsolve_special_functions DESTINATION ${ADDON_INSTALL_DIR_PYTHON})
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(ngsolve_special_functions PRIVATE -Wno-return-type-c-linkage)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

