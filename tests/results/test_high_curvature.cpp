#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "bdry3d/patch_surf_face_map.hpp"

using namespace hedgehog;
void curvature_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".25");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".25");
    Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    Options::set_value_petsc_opts("-bis3d_spacing", ".2");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".2");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "0");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".7");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "5e-3");
}


TEST_CASE("Test prop qbkix laplace ", "[qbkix-curve][results][laplace]"){
    int patch_order = 6;  // 6th degree patches
    int patch_refinement_factor = 1; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "squished_cube.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_high_curvature/";
    int num_iters = 4;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_far_field_eval,
            &curvature_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}
