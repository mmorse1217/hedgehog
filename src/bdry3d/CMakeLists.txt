set(Hedgehog_geometry_src
    analytic_surf.cpp
    patch_samples.cpp
    face_map_subpatch.cpp
    patch_surf.cpp
    geogram_interface.cpp
    patch_surf_analytic.cpp
    p4est_interface.cpp
    patch_surf_blended.cpp
    p4est_refinement.cpp
    patch_surf_face_map.cpp
)
set(Hedgehog_geometry_header
    analytic_surf.hpp
    patch_samples.hpp
    face_map_subpatch.hpp
    patch_surf.hpp
    geogram_interface.hpp
    patch_surf_analytic.hpp
    on_surface_point.hpp
    patch_surf_blended.hpp
    p4est_interface.hpp
    patch_surf_face_map.hpp
    p4est_refinement.hpp
)

PREPEND(Hedgehog_geometry_src)
PREPEND(Hedgehog_geometry_header)

set(Hedgehog_geometry ${Hedgehog_geometry_src} 
    ${Hedgehog_geometry_header} PARENT_SCOPE)
