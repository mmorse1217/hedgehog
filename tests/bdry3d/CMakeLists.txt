set(Hedgehog_test_geometry_src
    test_face_map.cpp
    test_face_map_vs_blendsurf.cpp
    test_p4est_interface.cpp
    test_p4est_patch_refinement.cpp
    test_patch_samples.cpp
    test_patch_surf_blended.cpp
)
set(Hedgehog_test_geometry_header
)
PREPEND(Hedgehog_test_geometry_src)
PREPEND(Hedgehog_test_geometry_header)

set(Hedgehog_test_geometry ${Hedgehog_test_geometry_src} 
    ${Hedgehog_test_geometry_header} PARENT_SCOPE)
