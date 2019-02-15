
set(Hedgehog_test_utils_src
    evaluation_utils.cpp
    evaluation_utils.hpp
    regression_utils.cpp
    regression_utils.hpp
    test_benchmark.cpp
    view_qbkix_points.cpp
    write_patch_data_to_vtk.cpp
)
PREPEND(Hedgehog_test_utils_src)

set(Hedgehog_test_utils ${Hedgehog_test_utils_src} PARENT_SCOPE)
