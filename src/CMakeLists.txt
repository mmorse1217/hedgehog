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
target_link_libraries(Hedgehog PUBLIC ${BLENDSURF_LIB} ${P4EST_LIB} ${SC_LIB}
    ${PVFMM_INTERFACE_LIB} ${PVFMM_DIR} ${GEOGRAM_LIB} ${FACEMAP_LIB} ${MPI_LIB}  lapack blas )
