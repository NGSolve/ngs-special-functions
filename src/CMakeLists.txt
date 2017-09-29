project(ngs_special_functions)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

find_package(NGSolve CONFIG REQUIRED)

find_library(F2C_LIBRARY f2c REQUIRED HINTS ${NGSOLVE_LIBRARY_DIR})

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/i1mach.c)
    execute_process(COMMAND ${CMAKE_COMMAND}
      -DBIN_DIR=${CMAKE_CURRENT_BINARY_DIR}
      -DNGSOLVE_BINARY_DIR=${NGSOLVE_BINARY_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/fetch_and_convert_slatec.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/slatec/src)
endif()

file(GLOB SLATEC_SOURCES_C ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/*.c )

add_ngsolve_python_module(special_functions
  specialcf.cpp ${SLATEC_SOURCES_C}
)

target_include_directories(special_functions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
target_link_libraries(special_functions ${F2C_LIBRARY})

install(TARGETS special_functions DESTINATION ngsolve)
