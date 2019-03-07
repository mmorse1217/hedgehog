#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "common/vtk_writer.hpp"
#include "../utils/evaluation_utils.hpp"
#include "bdry3d/patch_surf_face_map.hpp"

using namespace Ebi;

void adaptive_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".11");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".2");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".076923");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".076923");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0625");
    Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0625");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    if(Options::get_string_from_petsc_opts("-bd3d_filename") == "wrl_files/nearly_touching.wrl"){
        Options::set_value_petsc_opts("-dom", "1");
    } else {
        Options::set_value_petsc_opts("-dom", "0");
    }
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "20");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "3");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "5e-3");
    //Options::set_value_petsc_opts("-target_accuracy", "1e-8");
    Options::set_value_petsc_opts("-target_accuracy", "1e-5");
}

void uniform_test_base_options(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");

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

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-target_accuracy", "1e-7");
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
    surface->resolve_rhs(&laplace_singluarity_propeller, 1, 1e-9 );
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

