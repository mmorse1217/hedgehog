set(Hedgehog_test_quadrature_src
    test_error_estimate.cpp
    test_evaluation.cpp
    test_evaluator_near.cpp
    test_evaluator_near_extrapolate.cpp
    test_evaluator_near_interpolate.cpp
    test_lagrange_interpolation.cpp
    test_markgrid.cpp
    test_patch_refinement.cpp
    test_point_marking.cpp
    test_solver.cpp
    test_solver_gmres_double_layer.cpp
)
set(Hedgehog_test_quadrature_header
    test_lagrange_interpolation.hpp
)
PREPEND(Hedgehog_test_quadrature_src)
PREPEND(Hedgehog_test_quadrature_header)

set(Hedgehog_test_quadrature ${Hedgehog_test_quadrature_src} 
    ${Hedgehog_test_quadrature_header} PARENT_SCOPE)
