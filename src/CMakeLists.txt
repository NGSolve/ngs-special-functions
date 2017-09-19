project(ngs_special_functions)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

find_package(NGSolve CONFIG REQUIRED)

find_library(F2C_LIBRARY f2c REQUIRED HINTS ${NGSOLVE_LIBRARY_DIR})
find_program(F2C_COMMAND f2c REQUIRED HINTS ${NGSOLVE_BINARY_DIR})

macro( fetch_and_convert_slatec_sources slatec_function )
    file(DOWNLOAD http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=slatec%2Fsrc%2F${slatec_function}.f ${CMAKE_CURRENT_BINARY_DIR}/${slatec_function}.tgz)
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ${CMAKE_CURRENT_BINARY_DIR}/${slatec_function}.tgz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro( fetch_and_convert_slatec_sources )

fetch_and_convert_slatec_sources(zbesi)
fetch_and_convert_slatec_sources(zbesj)
fetch_and_convert_slatec_sources(zbesk)
fetch_and_convert_slatec_sources(gamln)
file(DOWNLOAD http://www.netlib.org/blas/i1mach.f ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/i1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/r1mach.f ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/r1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/d1mach.f ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/d1mach.f)

file(GLOB SLATEC_SOURCES_FORTRAN ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/*.f )
execute_process(COMMAND ${F2C_COMMAND} -a ${SLATEC_SOURCES_FORTRAN} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/slatec/src)
file(GLOB SLATEC_SOURCES_C ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/*.c )

add_ngsolve_python_module(special_functions
  specialcf.cpp ${SLATEC_SOURCES_C}
)

target_include_directories(special_functions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
target_link_libraries(special_functions ${F2C_LIBRARY})

install(TARGETS special_functions DESTINATION ngsolve)
