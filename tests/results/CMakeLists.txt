
set(Hedgehog_test_results_src
    generate_paraview_files.cpp
    test_adaptive.cpp
    test_expansion_error.cpp
    test_gmres_convergence_qbkix.cpp
    test_greens_identity.cpp
    test_high_curvature.cpp
    test_qbkix_eval_convergence.cpp
    test_qbkix_vs_singular_quad.cpp
    test_quadrature_error.cpp
    test_reproduce_ybz_2006.cpp
    test_singular_eval_convergence.cpp
)
PREPEND(Hedgehog_test_results_src)

set(Hedgehog_test_results ${Hedgehog_test_results_src} PARENT_SCOPE)
