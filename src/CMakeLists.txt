project(ngs_special_functions)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

find_package(NGSolve CONFIG REQUIRED)

find_library(F2C_LIBRARY f2c REQUIRED HINTS ${NGSOLVE_LIBRARY_DIR})

file(GLOB SLATEC_SOURCES_C ${CMAKE_CURRENT_BINARY_DIR}/slatec/src/*.c )

add_ngsolve_python_module(special_functions
  specialcf.cpp ${SLATEC_SOURCES_C}
)

target_include_directories(special_functions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
target_link_libraries(special_functions PRIVATE ${F2C_LIBRARY})
target_compile_definitions(special_functions PRIVATE -DSPECIALCF_LIBRARY_NAME="${NETGEN_PYTHON_DIR}/ngsolve/$<TARGET_LINKER_FILE_NAME:special_functions>")

install(TARGETS special_functions DESTINATION ngsolve)
