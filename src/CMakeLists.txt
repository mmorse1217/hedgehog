
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

set(CMAKE_CXX_STANDARD 11)
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
target_link_libraries(Hedgehog PUBLIC ${VTK_LIBRARIES}) 
target_link_libraries(Hedgehog PUBLIC ${Patchwork_LIBRARIES})
target_link_libraries(Hedgehog PUBLIC Blendsurf::Blendsurf)
target_link_libraries(Hedgehog PUBLIC Geogram::geogram)
target_link_libraries(Hedgehog PUBLIC p4est::p4est p4est::sc)
#target_link_libraries(Hedgehog PUBLIC PETSc::PETSc)
target_link_libraries(Hedgehog PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(Hedgehog PUBLIC PVFMM::PVFMM)
if(WITH_STKFMM)
target_link_libraries(Hedgehog PUBLIC STKFMM::STKFMM)
target_compile_definitions(Hedgehog PRIVATE HAS_STKFMM=1)
endif()
target_link_libraries(Hedgehog PUBLIC Nanospline::Nanospline)
target_link_libraries(Hedgehog PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries(Hedgehog PUBLIC ${MPI_LIBRARIES})
#target_link_libraries(Hedgehog PUBLIC ${FFTW_LIBRARIES})
target_link_libraries(Hedgehog PUBLIC FFTW::FFTW)
target_link_libraries(Hedgehog PUBLIC ${VTK_LIBRARIES})
target_link_libraries(Hedgehog PUBLIC Eigen3::Eigen)

#target_include_directories(Hedgehog PUBLIC ${EIGEN3_INCLUDE_DIR})
target_include_directories(Hedgehog PUBLIC ${VTK_INCLUDE}) 
target_include_directories(Hedgehog PUBLIC ${Patchwork_INCLUDE_DIRS})
target_include_directories(Hedgehog PUBLIC ${PETSC_INCLUDES})
target_include_directories(Hedgehog PUBLIC ${LAPACK_INCLUDE} )
target_include_directories(Hedgehog PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(Hedgehog PUBLIC ${VTK_INCLUDE_DIRS})
#target_include_directories(Hedgehog PUBLIC ${FFTW_INCLUDE_DIRS})

target_include_directories(Hedgehog PUBLIC ${Hedgehog_SOURCE_DIR}/src)
target_compile_options(Hedgehog PUBLIC "-DFFTW3")
