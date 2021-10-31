set(Hedgehog_test_common_src
     test_boundary_data.cpp
     test_stats.cpp
     test_fft_refinement.cpp
     test_utils.cpp
     test_interpolate.cpp
)
PREPEND(Hedgehog_test_common_src)

set(Hedgehog_test_common ${Hedgehog_test_common_src} PARENT_SCOPE)
