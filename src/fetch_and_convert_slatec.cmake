find_program(F2C_COMMAND f2c REQUIRED HINTS ${NGSOLVE_BINARY_DIR})

macro( fetch_and_convert_slatec_sources slatec_function )
  file(DOWNLOAD http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=slatec%2Fsrc%2F${slatec_function}.f ${BIN_DIR}/${slatec_function}.tgz)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ${BIN_DIR}/${slatec_function}.tgz WORKING_DIRECTORY ${BIN_DIR})
endmacro( fetch_and_convert_slatec_sources )

fetch_and_convert_slatec_sources(zbesi)
fetch_and_convert_slatec_sources(zbesj)
fetch_and_convert_slatec_sources(zbesk)
fetch_and_convert_slatec_sources(gamln)
fetch_and_convert_slatec_sources(zbesh)
file(DOWNLOAD http://www.netlib.org/blas/i1mach.f ${BIN_DIR}/slatec/src/i1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/r1mach.f ${BIN_DIR}/slatec/src/r1mach.f)
file(DOWNLOAD http://www.netlib.org/blas/d1mach.f ${BIN_DIR}/slatec/src/d1mach.f)

file(DOWNLOAD http://ab-initio.mit.edu/Faddeeva.cc ${BIN_DIR}/Faddeeva/Faddeeva.cc)
file(DOWNLOAD http://ab-initio.mit.edu/Faddeeva.hh ${BIN_DIR}/Faddeeva/Faddeeva.hh)

file(GLOB SLATEC_SOURCES_FORTRAN ${BIN_DIR}/slatec/src/*.f )
execute_process(COMMAND ${F2C_COMMAND} -a ${SLATEC_SOURCES_FORTRAN} WORKING_DIRECTORY ${BIN_DIR}/slatec/src)
