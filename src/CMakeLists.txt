cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(bdry3d)
add_subdirectory(bie3d)
add_subdirectory(common)
add_subdirectory(deps)
add_subdirectory(fmm3d)


set(Hedgehog_dirs
    ${Hedgehog_quadrature}
    ${Hedgehog_geometry}
    ${Hedgehog_common}
    ${Hedgehog_dep}
    ${Hedgehog_fmm}
    )

add_library(Hedgehog STATIC ${Hedgehog_dirs})

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
link_directories($ENV{VTK_DIR}/lib)
target_link_libraries(Hedgehog PUBLIC ${VTK_LIBRARIES} ${FACEMAP_LIB} ${BLENDSURF_LIB} 
    ${P4EST_LIB} ${SC_LIB} ${PVFMM_INTERFACE_LIB} ${PVFMM_LIB} ${GEOGRAM_LIB} 
    ${FFTW_LIB}  ${PETSC_LIB} ${MPI_LIB}  lapack blas expat)
target_compile_options(Hedgehog PUBLIC "-DFFTW3")
