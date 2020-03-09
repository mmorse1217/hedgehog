#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"

using namespace Ebi;



TEST_CASE("Test solver convergence with qbkix: torus ", "[results][qbkix-solver-conv][laplace][newtorus]"){
    int patch_order = 3;  // cubic patches
    int patch_refinement_factor = 0; // no coarse grid refinement; geometry is resolved
    int kernel_enum = 111; // laplace problem
    string domain = "newtorus.wrl"; // domain mesh
    string polynomial_patch_filename = "explicit_torus_patches.poly";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::POLYNOMIAL,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test laplace solver convergence with qbkix: cube", "[results][qbkix-solver-conv][laplace][cube]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 2; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "cube.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("Test laplace solver convergence with qbkix: pipe", "[results][qbkix-solver-conv][laplace][pipe]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "pipe.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test laplace solver convergence with qbkix: propeller", "[results][qbkix-solver-conv][laplace][ttorus]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "ttorus2.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test navier solver convergence with qbkix: cube", "[results][qbkix-solver-conv][navier][cube]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 2; // slightly refine geometry
    int kernel_enum = 511; // navier problem
    string domain = "cube.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("Test navier solver convergence with qbkix: pipe", "[results][qbkix-solver-conv][navier][pipe]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 511; // navier problem
    string domain = "pipe.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test navier solver convergence with qbkix: propeller", "[results][qbkix-solver-conv][navier][ttorus]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 511; // navier problem
    string domain = "ttorus2.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test stokes solver convergence with qbkix: cube", "[results][qbkix-solver-conv][stokes][cube]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 2; // slightly refine geometry
    int kernel_enum = 311; // stokes problem
    string domain = "cube.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("Test stokes solver convergence with qbkix: pipe", "[results][qbkix-solver-conv][stokes][pipe]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 311; // stokes problem
    string domain = "pipe.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test stokes solver convergence with qbkix: propeller", "[results][qbkix-solver-conv][stokes][ttorus]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 311; // stokes problem
    string domain = "ttorus2.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}



// ----------------------------------------------------------------------------------
// old tests
// ----------------------------------------------------------------------------------
TEST_CASE("Test evaluation convergence with qbkix", "[results][qbkix-eval]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "3");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
   
    string output_folder = "output/test_solver_convergence_qbkix/";
    
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =1; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
            surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
            surface->_coarse = true;
            surface->setFromOptions();
            surface->setup();
            //surface->refine();
            //surface->refine_test();
            surface->refine_uniform(i);
            stats.print_results(); 


            int num_coarse_patches = surface->num_patches();
            stats.add_result("num coarse patches", num_coarse_patches);

            TestConfig test;
            // data generated from singularities on a sphere outside the domain
            test.bc_type = BoundaryDataType::HARMONIC;
            //test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
            test.singularity_type = SingularityType::SPHERE;
            test.sphere_radius_bc= 2;

            // target on a sphere in the far field
            test.target_type   = TargetType::RING;
            test.sphere_radius_targets = .5;
            //test.target_type   = TargetType::COLLOCATION_POINTS;

            // Solve with GMRES using two-sided qbkix 
            test.solver_matvec_type = EvaluationType::EXTRAPOLATION_AVERAGE;
            //test.solution_scheme = SolutionScheme::GMRES_SOLVE;
            test.solution_scheme = SolutionScheme::EXPLICIT_DENSITY;
            
            // evaluate at targets with smooth quad
            test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
            //test.evaluation_scheme = EvaluationScheme::ON_QBKIX_AVERAGE;
            //test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
            
            //string filename = to_string(num_coarse_patches)+string("_patches_abs_error_singularity_radius_")+
                //to_string(test.sphere_radius_bc)+string(".vtp");
            string filename = string("ref_lvl_")+to_string(i);
            string full_path= output_folder + filename + string(".vtp");

            /*stats._file_prefix = output_folder +to_string(num_coarse_patches)+
                string("_patches_singularity_radius_") +
                to_string(test.sphere_radius_bc);*/
            
            stats._file_prefix = output_folder + filename;
            test.error_filename = full_path;

            test.dump_values = true;
            run_test(surface.get(),test);
            //cout << string("output/"+test_name+filename << endl;
            stats.dump_key_values("", stats._file_prefix);
                    //to_string(num_coarse_patches)+string("_patches"));
        }
    }
}

TEST_CASE("Test stokes solver convergence with qbkix: hourglass", "[results][qbkix-solver-conv][stokes][hourglass]"){
    int patch_order = 12;  // 6th degree patches
    int patch_refinement_factor = 1; // slightly refine geometry
    int kernel_enum = 311; // stokes problem
    string domain = "hourglass.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval,
            &gmres_test_base_options,
            polynomial_patch_filename,
            num_iters);
}


TEST_CASE("Test face-map + qbkix on propeller...",
        "[qbkix][face-map][prop][laplace]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/ppp.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/ppp.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");

    /*
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    */
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "3");
    
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive","1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","2");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","0");
    Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-2");
    //Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", ".1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-bdtype", "2"); // face-map surface
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");


    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "20");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));

    double face_map_spacing= .1;
    string output_folder = "output/test_qbkix_eval_convergence/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =4; 
        //int num_levels_refinement =1; 
        for(int i = 1; i < num_levels_refinement; i++){

            double inflation_factor =  .5 + double(i)/double(2*num_levels_refinement-1);

            Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", 
                    std::to_string(inflation_factor));

            
            unique_ptr<PatchSurfFaceMap> face_map_surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
            face_map_surface->_surface_type = PatchSurfFaceMap::BLENDED;
            face_map_surface->_coarse = true;
            face_map_surface->setFromOptions();
            face_map_surface->setup();
            face_map_surface->refine_test();

            TestConfig test;
            // data generated from Y^1_1 i.e. (1,1)-spherical harmonic
            test.bc_type = BoundaryDataType::HARMONIC;
            test.singularity_type= SingularityType::SPHERE;
            test.sphere_radius_bc = 2;

            // Evaluate density at the collocation points
            //test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
            test.target_type   = TargetType::COLLOCATION_POINTS;

            // Solve with GMRES using solver_matvec_type
            test.solver_matvec_type = EvaluationType::EXTRAPOLATION_AVERAGE;
            test.solution_scheme = SolutionScheme::GMRES_SOLVE;
            test.target_type   = TargetType::COLLOCATION_POINTS;
                
            test.evaluation_scheme = EvaluationScheme::ON_QBKIX;


            string filename = string("qbx__ref_lvl_")+to_string(i);
            string full_path= output_folder + filename + string(".vtp");


            stats._file_prefix = output_folder + filename;
            test.error_filename = full_path;

            test.dump_values = true;
            run_test(face_map_surface.get(),test);

            stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
            stats.add_result("bbox inflation factor", Options::get_double_from_petsc_opts("-adaptive_upsampling_bbox_inflation_factor"));
            stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
            stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
            stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
            stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
            stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
            stats.print_results();
            stats.dump_key_values("", stats._file_prefix);
            stats.clear();

        }
    }
}



void test_gmres_solve_near_eval_mpi(
        PatchSurfFaceMap* surface,
        string output_folder,
        int i){
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    //test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 2.;
    string d = Test::get_domain();
    /*if( d== "new_ppp" || d == "comb2" || d == "nearly_touching"){
        test.sphere_radius_bc = 2.;
    } else {
        test.sphere_radius_bc = 1.;
    }*/

    test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
    //test.target_type   = TargetType::COLLOCATION_POINTS;
   // ... via Green's Identity
   //test.evaluation_scheme = EvaluationScheme::ON_QBKIX_AVERAGE;
   test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
    // no solve step need
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    //test.solution_scheme   = SolutionScheme::EXPLICIT_DENSITY;
    //test.solver_matvec_type = INTERIOR_EXTRAPOLATION;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;
    
    string filename;
    int is_adaptive_test = Options::get_int_from_petsc_opts("-adaptive");
    if(is_adaptive_test){
        filename = string("qbx_adaptive_gmres_ref_lvl_")+to_string(i);
    } else {
        filename = string("qbx_gmres_ref_lvl_")+to_string(i);
    }
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface, test);

    int num_coarse_patches = surface->patches().size();
    stats.add_result("num coarse patches", num_coarse_patches);
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results(); 
    stats.clear();

}

void solver_test_base_options_mpi(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");/*
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".06");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
    Options::set_value_petsc_opts("-bis3d_spacing", ".090909");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909");
    */

    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".135");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".135");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".11");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".022");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

}



TEST_CASE("Test laplace solver convergence with qbkix MPI: cube", "[mpi-solve]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "pipe.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder = "output/test_qbkix_eval_convergence/";
    int num_iters = 1;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_gmres_solve_near_eval_mpi,
            &solver_test_base_options_mpi,
            polynomial_patch_filename,
            num_iters);
}
