set(Hedgehog_test_fmm3d_src
    test_fmm_interface.cpp
)
PREPEND(Hedgehog_test_fmm_src)

set(Hedgehog_test_fmm ${Hedgehog_test_fmm_src} PARENT_SCOPE)
