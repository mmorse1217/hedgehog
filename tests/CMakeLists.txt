
add_subdirectory(bdry3d)
add_subdirectory(bie3d)
add_subdirectory(common)
add_subdirectory(examples)
add_subdirectory(results)
add_subdirectory(fmm3d)
add_subdirectory(utils)

set(Hedgehog_test_dirs
    ${Hedgehog_test_quadrature}
    ${Hedgehog_test_geometry}
    ${Hedgehog_test_common}
    ${Hedgehog_test_dep}
    ${Hedgehog_test_fmm}
    ${Hedgehog_test_utils}
    ${Hedgehog_test_results}
    )
add_executable(test_hedgehog test_bis.cpp ${Hedgehog_test_dirs})
target_link_libraries(test_hedgehog Hedgehog)
