set(Hedgehog_quadrature_src
    collocation_patch_samples.cpp
    error_estimate.cpp
    evaluator.cpp
    evaluator_far.cpp
    evaluator_near.cpp
    evaluator_near_interpolate.cpp
    evaluator_on_surface.cpp
    evaluator_qbkix.cpp
    evaluator_qbkix_average.cpp
    markgrid.cpp
    solver_gmres_double_layer.cpp
    solver_interface.cpp
    solver_utils.cpp
    spatial_grid.cpp
)
set(Hedgehog_quadrature_header
    collocation_patch_samples.hpp
    error_estimate.hpp
    evaluator.hpp
    evaluator_far.hpp
    evaluator_near.hpp
    evaluator_near_interpolate.hpp
    evaluator_on_surface.hpp
    evaluator_qbkix.hpp
    evaluator_qbkix_average.hpp
    markgrid.hpp
    solver_gmres_double_layer.hpp
    solver_interface.hpp
    solver_utils.hpp
    spatial_grid.hpp
)


PREPEND(Hedgehog_quadrature_src)
PREPEND(Hedgehog_quadrature_header)

set(Hedgehog_quadrature ${Hedgehog_quadrature_src} 
    ${Hedgehog_quadrature_header} PARENT_SCOPE)
