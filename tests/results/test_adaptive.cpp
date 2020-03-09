#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "common/vtk_writer.hpp"
#include "../utils/evaluation_utils.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include <sampling.hpp>
#include "bie3d/markgrid.hpp"
using namespace Ebi;
using Sampling::sample_2d;
using Sampling::equispaced;

void adaptive_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".11");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".12");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".03");
    Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".076923");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".076923");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0667");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0667");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    if(Options::get_string_from_petsc_opts("-bd3d_filename") == "wrl_files/nearly_touching.wrl" ||
       Options::get_string_from_petsc_opts("-bd3d_filename") == "wrl_files/interlocking_torii_flip.wrl"){
        Options::set_value_petsc_opts("-dom", "1");
    } else {
        Options::set_value_petsc_opts("-dom", "0");
    }
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    /*
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "3");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "5e-3");
    */
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "3");
    //Options::set_value_petsc_opts("-target_accuracy", "1e-8");
    Options::set_value_petsc_opts("-target_accuracy", "1e-5");
}

void uniform_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".004");

    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".14285");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", "1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "20");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
}


void adaptive_convergence_test(string domain, FunctionHandle f){
    string output_folder= "output/test_adaptive/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    // adaptively resolve boundary data and upsampled quadrature to \eps.
    adaptive_test_base_options(); 
    
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type =  PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine();

    int range_dof = Test::get_kernel() == "laplace" ? 1 : 3;

    int it =0;
    for(int i = 2; i <= 10; i+=2){
        double target_accuracy = pow(10.,-i);
        Options::set_value_petsc_opts("-target_accuracy", to_string(target_accuracy));
        face_map->resolve_rhs(f, range_dof, target_accuracy);// TODO write something to refine from PatchSample values? or pass the FunctionHandle explicitly for each geometry. also an option
        
        double distance = Options::get_double_from_petsc_opts("-boundary_distance_ratio");
        double spacing = Options::get_double_from_petsc_opts("-interpolation_spacing_ratio");
        Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(distance/2.));
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(spacing/2.));

        //test_gmres_solve_near_eval(face_map.get(), output_folder, it++);
        stats.print_results(); 
        stats.clear();
    }
    
}

void uniform_convergence_test(string domain, FunctionHandle f){
    string output_folder= "output/test_adaptive/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";

    // compare against k levels of uniform refinment of coarse + upsampled grids
    uniform_test_base_options();

    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type =  PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();
    for (int i = 0; i < 5; i++) {
        face_map->refine_uniform(1);
        test_gmres_solve_near_eval(face_map.get(), output_folder, i);
        stats.print_results(); 
        stats.clear();
        
    }
}
TEST_CASE("Test adaptive cube", "[adaptive][cube][results]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    adaptive_convergence_test("", laplace_singluarity_propeller);
}
TEST_CASE("Test uniform cube", "[uniform][cube][results]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    uniform_convergence_test("", laplace_singluarity_propeller);
}

TEST_CASE("Test adaptive vs. uniform for qbx location placement", "[adaptive-vs-uniform][results][half_donut]"){
    int patch_order = 20;  // 6th degree patches
    int patch_refinement_factor = 1; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "half_donut.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    //int num_iters = 1;
    int num_iters = 4;

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/half_donut.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/half_donut.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-kt", "111");

    stats._file_prefix = "data/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";

    unique_ptr<PatchSurfFaceMap> base(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    base->_surface_type =  PatchSurfFaceMap::BLENDED;
    base->_coarse = true;
    base->setFromOptions();
    base->setup();

    // make surface admissible
    uniform_test_base_options();
    base->refine();

    Options::set_value_petsc_opts("-adaptive", "0");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    uniform_test_base_options();
    unique_ptr<PatchSurfFaceMap>face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = base->_surface_type;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup_from_existing_face_map(base.get());
    face_map->initialize_with_existing_p4est(base->_p4est);
    face_map->refine_test();

    for(int i = 0; i < num_iters; i++){
        // i levels of uniform upsampling 
        Options::set_value_petsc_opts("-uniform_upsampling_num_levels", to_string(i));

        test_gmres_solve_near_eval(face_map.get(), output_folder, i);
        stats.print_results(); 
        stats.clear();
    }
    
    Options::set_value_petsc_opts("-adaptive", "1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    adaptive_test_base_options();
    int it=0;
    for(int i = 3; i < 10; i+=2){
        Options::set_value_petsc_opts("-target_accuracy", to_string(pow(10.,-i)));
        test_gmres_solve_near_eval(face_map.get(), output_folder, it);
        it++;
        stats.print_results(); 
        stats.clear();
    }

}


TEST_CASE("Interlocking torii example", "[adaptive][results][torii]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "interlocking_torii_flip.wrl"; // domain mesh

    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    int num_iters = 1;
    Options::set_value_petsc_opts("-adaptive", "1");
    Options::set_value_petsc_opts("-dom", "1");
    Options::set_value_petsc_opts("-bis3d_ksp_rtol", "1e-8");
    Options::set_value_petsc_opts("-bis3d_ksp_max_it", "500");
    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            //&test_gmres_solve_far_field_eval,
            &test_gmres_solve_near_eval,
            &adaptive_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}


TEST_CASE("Test adaptive propeller", "[adaptive][results][propeller]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "new_ppp.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    int num_iters = 1;
    Options::set_value_petsc_opts("-adaptive", "1");
    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            //&test_gmres_solve_far_field_eval,
            &test_gmres_solve_near_eval,
            &adaptive_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}


TEST_CASE("Test adaptive comb", "[adaptive][results][comb]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "comb2.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    int num_iters = 1;
    Options::set_value_petsc_opts("-adaptive", "1");
    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            //&test_gmres_solve_far_field_eval,
            &test_gmres_solve_near_eval,
            &adaptive_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}

TEST_CASE("Test adaptive pinched ", "[adaptive][results][nearly_touching]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "nearly_touching.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    int num_iters = 1;
    Options::set_value_petsc_opts("-adaptive", "1");
    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            //&test_gmres_solve_far_field_eval,
            &test_gmres_solve_near_eval,
            &adaptive_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}

TEST_CASE("Test adaptive octopus", "[adaptive][results][octopus]"){
    int patch_order = 16;  // 6th degree patches
    int patch_refinement_factor = 0; // slightly refine geometry
    int kernel_enum = 111; // laplace problem
    string domain = "closed_octopus.wrl"; // domain mesh
    string polynomial_patch_filename = "";
    string output_folder= "output/test_adaptive/";
    int num_iters = 1;
    Options::set_value_petsc_opts("-adaptive", "1");
    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            //&test_gmres_solve_far_field_eval,
            &test_gmres_solve_near_eval,
            &adaptive_test_base_options,
            polynomial_patch_filename,
            num_iters,
            "adaptive");
}

TEST_CASE("Cube eye candy", "[eye-candy][results][cube]"){
/*
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
*/
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/larger_vessel_section2.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/larger_vessel_section2.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "10");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-5");
    adaptive_test_base_options();
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".22");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".044");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0166666");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".25");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".05882350");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05882350");

    Options::set_value_petsc_opts("-adaptive", "1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");

    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/eye_candy/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";


    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine();
    surface->resolve_rhs(&laplace_singluarity_propeller, 1, 1e-7 );
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
            vector<int>(surface->num_patches(), 0) ,
            surface.get(), 0, "output/");
    
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SINGLE;
    //cube
    //test.single_singularity_location = Point3(0., 0., 1.);
    // blood vessel
    test.single_singularity_location = Point3(20., 0., 0.);
    //test.single_singularity_location = Point3(-.05, .85, .45);
    //test.single_singularity_location = Point3(0., 0., 1.);

    test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
    
    //test.target_type   = TargetType::GRID;
    //test.target_plane_point= Point3(-.05,0., .45);
    //test.target_plane_vec1 = Point3(1.,0.,0.);
    //test.target_plane_vec2 = Point3(0.,1.,0.);
    //test.num_targets = 10;
    
    //test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
    test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;

    string filename = string("qbx_adaptive_rhs_solution");
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface.get(), test);

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
TEST_CASE("Propeller eye candy", "[eye-candy][results][propeller]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/new_ppp.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");

/*
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
*/
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-6");
    adaptive_test_base_options();
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0166666");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");

    Options::set_value_petsc_opts("-adaptive", "1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");

    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/eye_candy/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";


    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine();
    surface->resolve_rhs(&laplace_singluarity_propeller, 1, 1e-8 );
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
            vector<int>(surface->num_patches(), 0) ,
            surface.get(), 0, "output/");
    
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SINGLE;
    test.single_singularity_location = Point3(0., 0., 1.);
    //test.single_singularity_location = Point3(-.05, .85, .45);
    //test.single_singularity_location = Point3(0., 0., 1.);

    test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
    
    //test.target_type   = TargetType::GRID;
    //test.target_plane_point= Point3(-.05,0., .45);
    //test.target_plane_vec1 = Point3(1.,0.,0.);
    //test.target_plane_vec2 = Point3(0.,1.,0.);
    //test.num_targets = 10;
    
    //test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
    test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;

    string filename = string("qbx_adaptive_rhs_solution");
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface.get(), test);

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


TEST_CASE("Octopus eye candy", "[eye-candy][results][octopus]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/closed_octopus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/closed_octopus.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");

/*
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
*/
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-4");
    adaptive_test_base_options();
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".013333");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".12");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".02");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "50000");

    Options::set_value_petsc_opts("-adaptive", "1");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");

    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/eye_candy/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";


    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine();
    surface->resolve_rhs(&laplace_singluarity_propeller, 1, 1e-6 );
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
            vector<int>(surface->num_patches(), 0) ,
            surface.get(), 0, "output/");
    
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    test.singularity_type= SingularityType::SINGLE;
    test.single_singularity_location = Point3(0., 0., .4);
    //test.single_singularity_location = Point3(-.05, .85, .45);
    //test.single_singularity_location = Point3(0., 0., 1.);

    test.target_type   = TargetType::ON_SURFACE_NON_COLLOCATION;
    
    //test.target_type   = TargetType::GRID;
    //test.target_plane_point= Point3(-.05,0., .45);
    //test.target_plane_vec1 = Point3(1.,0.,0.);
    //test.target_plane_vec2 = Point3(0.,1.,0.);
    //test.num_targets = 10;
    
    //test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
    test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;

    string filename = string("qbx_adaptive_rhs_solution");
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface.get(), test);

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

TEST_CASE("Torus constant density", "[const-den][results][torii]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/interlocking_torii_flip.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/interlocking_torii_flip.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/ttorusflip.wrl");
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/ttorusflip.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-4");
    adaptive_test_base_options();
    Options::set_value_petsc_opts("-dom", "1");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".013333");

    //Options::set_value_petsc_opts("-adaptive", "1");
    //Options::set_value_petsc_opts("-upsampling_type", "adaptive");

    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05");
    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine();
/*
Options::set_value_petsc_opts("-bdtype", "1");
    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
*/
    //surface->resolve_rhs(&laplace_singluarity_propeller, 1, 1e-6 );
    //write_face_map_patches_to_vtk(DblNumMat(0,0), 
    //        vector<int>(surface->num_patches(), 0) ,
   //         surface.get(), 0, "output/");
    
    // Set up test case
    TestConfig test;
    
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
    //test.bc_type = BoundaryDataType::SINGLE;
    //test.bc_type = BoundaryDataType::HARMONIC;
    //test.singularity_type= SingularityType::SINGLE;
    //test.single_singularity_location = Point3(0., 0., 0.);
    //test.bc_type = BoundaryDataType::HARMONIC;
    //    test.singularity_type= SingularityType::SPHERE;
    //    test.sphere_radius_bc = .05;

    //test.target_type   = TargetType::COLLOCATION_POINTS;
    test.target_type   = TargetType::SPHERE_FAR;
    test.sphere_radius_targets   = 4.;
    
    //test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
    test.solver_matvec_type= EvaluationType::EXTRAPOLATION_AVERAGE;
    test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.dump_values = true;

    string filename = string("qbx_adaptive_rhs_solution");
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;
    run_test(surface.get(), test);

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

void singularity_torii_exterior(Vec samples, int dof,Vec& potential){
    DblNumMat samples_local(DIM, samples);
    DblNumMat normals(DIM, samples_local.n());
    Kernel3d laplace(111, vector<double>(2,1.));
    DblNumMat potential_local(laplace.get_tdof(), potential);
#pragma omp parallel for
    for (int i = 0; i < potential_local.n(); i++) {
        //Point3 x(-.05, .85, .45);
        Point3 x(0., 0.03, .875);
        Point3 y(samples_local.clmdata(i));
        potential_local(0,i) = 4.*M_PI/(x-y).l2();
        
    }
}


TEST_CASE("Dump surface mesh template: single proc torii", "[render][torii]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/interlocking_torii_flip.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/interlocking_torii_flip.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-7");
    adaptive_test_base_options();
    Options::set_value_petsc_opts("-bis3d_spacing", ".05");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    
    // Looks ok but needs more upsampling
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01333");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".025");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-dom", "1");
    
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "3");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");
    Options::set_value_petsc_opts("-dump_qbkix_points","0");

    //Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "5e-3");

    
    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    
    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine();
    surface->resolve_rhs(&singularity_torii_exterior, 1, 1e-9 );
    
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer(surface.get()));
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->_compute_refined_surface = true;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    
    solver->setFromOptions();
    solver->setup(); 

    //int num_samples = 10;
    double spacing = .05;
    int dof = 1; // 1 for Laplace, 3 for stokes
    int num_patches = surface->num_patches();
    int num_samples = int(round(1./spacing))+1;
    const int verts_per_patch= num_samples*num_samples;

    Vec target_points;
    Vec target_values;
    Vec computed_density;
    Vec true_target_values;
    Vec error;
    Petsc::create_mpi_vec(
            solver->mpiComm(),
            num_samples*num_samples*num_patches*DIM,
            target_points);
    Petsc::create_mpi_vec(
            solver->mpiComm(),
            num_samples*num_samples*num_patches*dof, 
            target_values);
    VecDuplicate(target_values, &true_target_values);
    VecCopy(target_values, true_target_values);
  
    int local_density_total_dofs, local_pole_total_dofs, local_total_dofs;
    int global_density_total_dofs, global_pole_total_dofs, global_total_dofs;
    solver->localSize(
            local_density_total_dofs,
            local_pole_total_dofs,
            local_total_dofs);
    solver->globalSize(
            global_density_total_dofs,
            global_pole_total_dofs,
            global_total_dofs);

    Petsc::create_mpi_vec(
            solver->mpiComm(),
            local_total_dofs, 
            computed_density);
    Vec boundary_condition;
    VecDuplicate(computed_density, &boundary_condition);
    VecCopy(computed_density, boundary_condition);

    
    VecDuplicate(target_values, &error);
    VecCopy(target_values, error);

    VecSet(target_values, 0.);
    
    // Setup point singularity infofor boundary condition
    Vec singularity_positions;
    Vec singularity_strengths;
    int num_singularities = 1;

    Petsc::create_mpi_vec(solver->mpiComm(),
            num_singularities*DIM, 
            singularity_positions);

    Petsc::create_mpi_vec(solver->mpiComm(),
            num_singularities*dof, 
            singularity_strengths);
    { // actually set values
        DblNumMat singularity_positions_local(DIM, singularity_positions);
        DblNumMat singularity_strengths_local(dof, singularity_strengths);
        setvalue(singularity_positions_local, 0.);
        setvalue(singularity_strengths_local, 4*M_PI);
        singularity_positions_local(2,0) = .875;
    }


    // Evaluate boundary condition
    Vec potential= Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            solver->problem_kernel(),
            singularity_positions,
            singularity_strengths,
            solver->patch_samples()->sample_point_3d_position(),
            solver->patch_samples()->sample_point_3d_position());
    // Fill a larger array of the appropriate size with the boundary condition values.
    vector<int64_t> range(solver->patch_samples()->local_num_sample_points(),0);
    std::iota(range.begin(), range.end(),0);
    double* potential_ptr;
    VecGetArray(potential, &potential_ptr);
    VecSetValues(boundary_condition, range.size(), range.data(), potential_ptr, INSERT_VALUES);

    // Solve with GMRES
    //VecSet(computed_density, 1.);
    int64_t size;
    VecGetSize(potential, &size);
    cout << "boundary_condition: " << size << endl;
    VecGetSize(computed_density, &size);
    cout << "computed_density: " << size << endl;
    solver->solve(boundary_condition, computed_density);
    VecDestroy(&boundary_condition);


    DblNumMat parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
            sample_2d<equispaced>(num_samples, Sampling::base_domain).data());

    const int triangles_per_patch= 2*(num_samples-1)*(num_samples-1);
    // Sample patch at equispaced points, restore to petsc vector
    int num_total_vertices = num_patches*num_samples*num_samples;
    int num_total_triangles = num_patches*triangles_per_patch; // 

    // global surface mesh
    //DblNumMat surface_vertices(DIM, num_total_vertices);
    DblNumMat surface_vertices(DIM, target_points);
    IntNumMat surface_triangles(DIM, num_total_triangles);
    
    vector<vector<uint> > patch_to_triangle_id_map(num_patches, vector<uint>());
    #pragma omp parallel for
    for (int pi = 0; pi < surface->num_patches(); pi++) {
        const auto& patch = surface->patch(pi);
        
        DblNumMat vertices;
        IntNumMat triangles;
        
        patch->mesh_patch(
                spacing, Rectangle(Interval(0.,1.), Interval(0.,1.)),
                vertices, triangles);
        int num_vertices_per_patch_1d = patch->mesh_num_vertices(spacing);
        int num_vertices_per_patch = patch->mesh_num_vertices(spacing);
        int num_triangles_per_patch = patch->mesh_num_triangles(spacing);

        assert(vertices.m() == DIM);
        assert(vertices.n() == num_vertices_per_patch);
        assert(triangles.m() == 3);
        assert(triangles.n() == num_triangles_per_patch);

        // copy vertices and triangle into global matrices
        for (int i = 0; i < num_vertices_per_patch; i++) {
            for (int d = 0; d < DIM; d++) {
                surface_vertices(d, pi*num_vertices_per_patch + i) = 
                    vertices(d,i);
            }
        }
        
        vector<uint> triangles_on_patch(num_triangles_per_patch);
        for (int i = 0; i < num_triangles_per_patch; i++) {
            for (int d = 0; d < 3; d++) {
                surface_triangles(d, pi*num_triangles_per_patch+ i) = 
                    pi*num_vertices_per_patch+triangles(d,i);
            }
            triangles_on_patch[i] = pi*num_triangles_per_patch +i;
        }
        patch_to_triangle_id_map[pi] = triangles_on_patch;
        setvalue(vertices, 0.);
        setvalue(triangles, 0);
    }

    // Compute the potential at the mesh points due to solved density
    // TODO: use the existing parallel marking code. should be plug and play as
    // is
    NumVec<OnSurfacePoint> closest_points(num_total_vertices);
    Vec eval_density;

    // get a pointer to the local array of the Petsc Vec
    double* density_ptr;
    VecGetArray(computed_density, &density_ptr);

    // Allocate the actual sub-Vec's for the density and the pole_coefficients
    VecCreateMPIWithArray(solver->mpiComm(), 1, local_density_total_dofs, global_density_total_dofs, density_ptr, &eval_density);
    solver->evaluate(target_points, 
            eval_density, target_values, closest_points, VAR_U);
    VecDestroy(&eval_density);
    VecRestoreArray(computed_density, &density_ptr);

    // Computed true potential: TODO: replace "0" with however the potential is computed.
    // probably will be 1./(singularity_location -sample)^2 or something
    /*{
        DblNumMat true_target_values_local(dof, true_target_values);
        DblNumMat target_values_local(dof, target_values);
        for (int pi = 0; pi < surface->num_patches(); pi++) {
            const auto& patch = surface->patch(pi);
            int num_vertices_per_patch = patch->mesh_num_vertices(spacing);
            for (int i = 0; i < num_vertices_per_patch; i++) {
                Point3 sample_point(surface_vertices.clmdata(pi*num_vertices_per_patch + i));
                for (int d = 0; d < dof; d++) {
                    true_target_values_local(d, pi*num_vertices_per_patch + i) = 0.;
                }
            }
        }

    }*/
    VecDestroy(&true_target_values);
    true_target_values= Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            solver->problem_kernel(),
            singularity_positions,
            singularity_strengths,
            target_points,
            target_points);
    VecCopy(true_target_values, error);
    VecAXPY(error, -1., target_values);
    VecAbs(error);
    double error_value;
    VecNorm(error, NORM_INFINITY, &error_value);
    cout << "Error, inf-norm: " << error_value << endl;
    VecNorm(error, NORM_2, &error_value);
    cout << "Error, 2-norm: " << error_value << endl;
    // compute inverted indexing of triangle ids -> patch id
    vector<uint> triangle_ids_to_patch_ids(num_total_triangles);
    for(int pi =0; pi < num_patches; pi++){
        const auto& triangles_on_patch = patch_to_triangle_id_map[pi];
        for(const int triangle_id : triangles_on_patch){
            triangle_ids_to_patch_ids[triangle_id] = pi;
        }
    }
    // dump triangle mesh for visualization
    vector<int> pids;
    pids.assign(triangle_ids_to_patch_ids.begin(), triangle_ids_to_patch_ids.end());

    DblNumMat surface_values(dof, error);
    write_triangle_mesh_to_vtk(surface_vertices, surface_triangles, 0, 
            "mpirank_"+to_string(surface->mpiRank()), 
            //"/scratch/mjm1030/mpirank2_"+to_string(surface->mpiRank()), 
            pids, surface_values);
    stats.add_result("num coarse patches", num_patches);
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results(); 
    VecDestroy(&singularity_positions);
    VecDestroy(&singularity_strengths);

}

TEST_CASE("Dump surface mesh template: single proc vessel", "[render][vessel]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/largest_vessel_section.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/largest_vessel_section.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-6");
    adaptive_test_base_options();
    Options::set_value_petsc_opts("-bis3d_spacing", ".8");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".8");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "0");

    
    stats._file_prefix = "data/";
    string output_folder= "output/test_adaptive/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    
    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_test();
    /*
    surface->refine();
    surface->resolve_rhs(&singularity_torii_exterior, 1, 1e-6 );
    */
    
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer(surface.get()));
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->_compute_refined_surface = true;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    
    solver->setFromOptions();
    solver->setup(); 

    //int num_samples = 10;
    double spacing = .05;
    int dof = 1; // 1 for Laplace, 3 for stokes
    int num_patches = surface->num_patches();
    int num_samples = int(round(1./spacing))+1;
    const int verts_per_patch= num_samples*num_samples;

    Vec target_points;
    Vec target_values;
    Vec computed_density;
    Vec true_target_values;
    Vec error;
    Petsc::create_mpi_vec(
            solver->mpiComm(),
            num_samples*num_samples*num_patches*DIM,
            target_points);
    Petsc::create_mpi_vec(
            solver->mpiComm(),
            num_samples*num_samples*num_patches*dof, 
            target_values);
    VecDuplicate(target_values, &true_target_values);
    VecCopy(target_values, true_target_values);
  
    int local_density_total_dofs, local_pole_total_dofs, local_total_dofs;
    int global_density_total_dofs, global_pole_total_dofs, global_total_dofs;
    solver->localSize(
            local_density_total_dofs,
            local_pole_total_dofs,
            local_total_dofs);
    solver->globalSize(
            global_density_total_dofs,
            global_pole_total_dofs,
            global_total_dofs);

    Petsc::create_mpi_vec(
            solver->mpiComm(),
            local_total_dofs, 
            computed_density);
    Vec boundary_condition;
    VecDuplicate(computed_density, &boundary_condition);
    VecCopy(computed_density, boundary_condition);

    
    VecDuplicate(target_values, &error);
    VecCopy(target_values, error);

    VecSet(target_values, 0.);
    
    // Setup point singularity infofor boundary condition
    Vec singularity_positions;
    Vec singularity_strengths;
    int num_singularities = 1;

    Petsc::create_mpi_vec(solver->mpiComm(),
            num_singularities*DIM, 
            singularity_positions);

    Petsc::create_mpi_vec(solver->mpiComm(),
            num_singularities*dof, 
            singularity_strengths);
    { // actually set values
        DblNumMat singularity_positions_local(DIM, singularity_positions);
        DblNumMat singularity_strengths_local(dof, singularity_strengths);
        setvalue(singularity_positions_local, 0.);
        setvalue(singularity_strengths_local, 4*M_PI);
        singularity_positions_local(2,0) = .875;
    }


    // Evaluate boundary condition
    Vec potential= Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            solver->problem_kernel(),
            singularity_positions,
            singularity_strengths,
            solver->patch_samples()->sample_point_3d_position(),
            solver->patch_samples()->sample_point_3d_position());
    // Fill a larger array of the appropriate size with the boundary condition values.
    vector<int64_t> range(solver->patch_samples()->local_num_sample_points(),0);
    std::iota(range.begin(), range.end(),0);
    double* potential_ptr;
    VecGetArray(potential, &potential_ptr);
    VecSetValues(boundary_condition, range.size(), range.data(), potential_ptr, INSERT_VALUES);

    // Solve with GMRES
    //VecSet(computed_density, 1.);
    int64_t size;
    VecGetSize(potential, &size);
    cout << "boundary_condition: " << size << endl;
    VecGetSize(computed_density, &size);
    cout << "computed_density: " << size << endl;
    //solver->solve(boundary_condition, computed_density);
    VecSet(computed_density, 1.);
    VecDestroy(&boundary_condition);


    DblNumMat parametric_samples = DblNumMat(2, num_samples*num_samples, true, 
            sample_2d<equispaced>(num_samples, Sampling::base_domain).data());

    const int triangles_per_patch= 2*(num_samples-1)*(num_samples-1);
    // Sample patch at equispaced points, restore to petsc vector
    int num_total_vertices = num_patches*num_samples*num_samples;
    int num_total_triangles = num_patches*triangles_per_patch; // 

    // global surface mesh
    //DblNumMat surface_vertices(DIM, num_total_vertices);
    DblNumMat surface_vertices(DIM, target_points);
    IntNumMat surface_triangles(DIM, num_total_triangles);
    
    vector<vector<uint> > patch_to_triangle_id_map(num_patches, vector<uint>());
    #pragma omp parallel for
    for (int pi = 0; pi < surface->num_patches(); pi++) {
        const auto& patch = surface->patch(pi);
        
        DblNumMat vertices;
        IntNumMat triangles;
        
        patch->mesh_patch(
                spacing, Rectangle(Interval(0.,1.), Interval(0.,1.)),
                vertices, triangles);
        int num_vertices_per_patch_1d = patch->mesh_num_vertices(spacing);
        int num_vertices_per_patch = patch->mesh_num_vertices(spacing);
        int num_triangles_per_patch = patch->mesh_num_triangles(spacing);

        assert(vertices.m() == DIM);
        assert(vertices.n() == num_vertices_per_patch);
        assert(triangles.m() == 3);
        assert(triangles.n() == num_triangles_per_patch);

        // copy vertices and triangle into global matrices
        for (int i = 0; i < num_vertices_per_patch; i++) {
            for (int d = 0; d < DIM; d++) {
                surface_vertices(d, pi*num_vertices_per_patch + i) = 
                    vertices(d,i);
            }
        }
        
        vector<uint> triangles_on_patch(num_triangles_per_patch);
        for (int i = 0; i < num_triangles_per_patch; i++) {
            for (int d = 0; d < 3; d++) {
                surface_triangles(d, pi*num_triangles_per_patch+ i) = 
                    pi*num_vertices_per_patch+triangles(d,i);
            }
            triangles_on_patch[i] = pi*num_triangles_per_patch +i;
        }
        patch_to_triangle_id_map[pi] = triangles_on_patch;
        setvalue(vertices, 0.);
        setvalue(triangles, 0);
    }

    // Compute the potential at the mesh points due to solved density
    // TODO: use the existing parallel marking code. should be plug and play as
    // is
    NumVec<OnSurfacePoint> closest_points(num_total_vertices);
    Vec eval_density;

    // get a pointer to the local array of the Petsc Vec
    double* density_ptr;
    VecGetArray(computed_density, &density_ptr);

    // Allocate the actual sub-Vec's for the density and the pole_coefficients
    VecCreateMPIWithArray(solver->mpiComm(), 1, local_density_total_dofs, global_density_total_dofs, density_ptr, &eval_density);
    //solver->evaluate(target_points, 
    //        eval_density, target_values, closest_points, VAR_U);
    VecDestroy(&eval_density);
    VecRestoreArray(computed_density, &density_ptr);

    // Computed true potential: TODO: replace "0" with however the potential is computed.
    // probably will be 1./(singularity_location -sample)^2 or something
    VecDestroy(&true_target_values);
    true_target_values= Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            solver->problem_kernel(),
            singularity_positions,
            singularity_strengths,
            target_points,
            target_points);
    VecCopy(true_target_values, error);
    VecAXPY(error, -1., target_values);
    VecAbs(error);
    double error_value;
    VecNorm(error, NORM_INFINITY, &error_value);
    cout << "Error, inf-norm: " << error_value << endl;
    VecNorm(error, NORM_2, &error_value);
    cout << "Error, 2-norm: " << error_value << endl;
    // compute inverted indexing of triangle ids -> patch id
    vector<uint> triangle_ids_to_patch_ids(num_total_triangles);
    for(int pi =0; pi < num_patches; pi++){
        const auto& triangles_on_patch = patch_to_triangle_id_map[pi];
        for(const int triangle_id : triangles_on_patch){
            triangle_ids_to_patch_ids[triangle_id] = pi;
        }
    }
    // dump triangle mesh for visualization
    vector<int> pids;
    pids.assign(triangle_ids_to_patch_ids.begin(), triangle_ids_to_patch_ids.end());

    DblNumMat surface_values(dof, error);
    write_triangle_mesh_to_vtk(surface_vertices, surface_triangles, 0, 
            "/scratch/mjm1030/mpirank_"+to_string(surface->mpiRank()), 
            pids, surface_values);
    stats.add_result("num coarse patches", num_patches);
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results(); 
    VecDestroy(&singularity_positions);
    VecDestroy(&singularity_strengths);

}
