#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"

using namespace Ebi;

void test_greens_identity(PatchSurfFaceMap* surface, string output_folder, int i){
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    /*test.singularity_type= SingularityType::SINGLE;
        test.single_singularity_location = Point3(2.,2.,2.);*/
    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 1.;

    // Evaluate solution at the collocation points
    test.target_type   = TargetType::COLLOCATION_POINTS;
    //test.target_type   = TargetType::RING;
    //test.sphere_radius_targets= .5;
    // ... via Green's Identity
    test.evaluation_scheme = EvaluationScheme::GREENS_IDENTITY;
    // no solve step need
    test.solution_scheme   = SolutionScheme::NONE;
    test.solver_matvec_type = NONE;

    test.dump_values = true;

    string filename = string("qbx_greens_id_ref_lvl_")+to_string(i);
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

void greens_identity_base_options(){
    /*
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "8");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "1");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".115");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".115");
    // one level of upsampling, good  convergence rates on cube and torus
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".090909090909");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909090909");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".058823");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".058823");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_np", "12");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    // DEBUG ONLY
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".016");
    
    // decent 
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".085");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".017");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".015");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".11");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".022");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".02");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".055");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".008125");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0125");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".065");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".02");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0125");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    //Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "2");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".047619");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".047619");
    //Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "2");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".037037");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01");
    //Options::set_value_petsc_opts("-near_interpolation_num_samples", "10");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".015");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0039");
    
    // 8 digits!
    //Options::set_value_petsc_opts("-bis3d_spacing", ".033333");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0075");
    
    // last options
    //Options::set_value_petsc_opts("-bis3d_spacing", ".02941");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".02");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    // 1e-1, 1e-2 4e-4
    //Options::set_value_petsc_opts("-bis3d_spacing", ".04");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    Options::set_value_petsc_opts("-bis3d_spacing", ".04");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "8000");
    Options::set_value_petsc_opts("-bis3d_np", "16");

    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    //Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".043748");
    Options::set_value_petsc_opts("-bis3d_spacing", ".0333");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".06666");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".02");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".00333");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    
    
    // Holy magic 
    Options::set_value_petsc_opts("-bis3d_ptsmax", "8000");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_spacing", "0.0454545");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.09090909");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "2");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.0588");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".066");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0111");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".02");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".00333");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".015");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0025");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
*/    
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".045");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".045");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".006666");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".090909");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "20");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "10000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    Options::set_value_petsc_opts("-bis3d_spacing", "0.05");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".0475");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".00225");
    // best for newtorus
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".0475");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".002");
    // best for cube
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".004");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0013");
    // Current best
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".03");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".005");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".0025");
    /*
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "3");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".06666");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".006667");
    */
}

void setup_and_run_greens_identity_test(
        int patch_order,
        int patch_refinement_factor,
        int kernel_enum,
        PatchSurfFaceMap::SurfaceType surface_type,
        string domain,
        string polynomial_patch_filename = string()){
    Options::set_value_petsc_opts("-bd3d_filename", string("wrl_files/") + string(domain));
    Options::set_value_petsc_opts("-bd3d_meshfile", string("wrl_files/") + string(domain));
    Options::set_value_petsc_opts("-poly_coeffs_file", 
            string("wrl_files/poly/") + string(polynomial_patch_filename));
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", to_string(patch_order));
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", to_string(patch_refinement_factor));
    Options::set_value_petsc_opts("-kt", to_string(kernel_enum));
    greens_identity_base_options();


    int num_levels_refinement = 6;
    string output_folder = "output/test_greens_identity/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    for(int i = 0; i < num_levels_refinement; i++){

        Options::set_value_petsc_opts("-debug_test_iter", to_string(i));
        auto surface = setup_face_map(surface_type, i);
        test_greens_identity(surface.get(), output_folder, i);
        stats.print_results(); 
    }
}

TEST_CASE("Test qbkix Green's Identity Laplace: torus", "[results][gid][laplace][newtorus]"){
    int patch_order = 3;  // cubic patches
    int patch_refinement_factor = 0; // no coarse grid refinement; geometry is resolved
    int kernel_enum = 111; // laplace problem
    string domain = "newtorus.wrl"; // domain mesh
    string polynomial_patch_filename = "explicit_torus_patches.poly";
    string output_folder = "output/test_greens_identity/";
    int num_iters = 6;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::POLYNOMIAL,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            polynomial_patch_filename,
            num_iters);
}

TEST_CASE("Test qbkix Green's Identity Laplace: cube", "[results][gid][laplace][cube]"){
    int patch_order = 20;  // higher order patches to more accurately resolve 
    //kinks in normals at patch corners
    int patch_refinement_factor = 2; // a bit of uniform fitting to better resolve patch corners
    int kernel_enum = 111; // laplace problem
    string domain = "cube.wrl"; // domain mesh
    string output_folder = "output/test_greens_identity/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            string(),
            num_iters);

}

TEST_CASE("Test qbkix Green's Identity Navier: torus", "[results][gid][navier][newtorus]"){
    int patch_order = 3;  // cubic patches
    int patch_refinement_factor = 0; // no coarse grid refinement; geometry is resolved
    int kernel_enum = 511; // elasticity problem
    string domain = "newtorus.wrl"; // domain mesh
    string polynomial_patch_filename = "explicit_torus_patches.poly";

    string output_folder = "output/test_greens_identity/";
    int num_iters = 6;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::POLYNOMIAL,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("Test qbkix Green's Identity Navier: cube", "[results][gid][navier][cube]"){
    int patch_order = 20;  // higher order patches to more accurately resolve 
    //kinks in normals at patch corners
    int patch_refinement_factor = 2; // a bit of uniform fitting to better resolve patch corners
    int kernel_enum = 511; // elasticity problem
    string domain = "cube.wrl"; // domain mesh
    string output_folder = "output/test_greens_identity/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            string(),
            num_iters);

}

TEST_CASE("Test qbkix Green's Identity Stokes: torus", "[results][gid][stokes][newtorus]"){
    int patch_order = 3;  // cubic patches
    int patch_refinement_factor = 0; // no coarse grid refinement; geometry is resolved
    int kernel_enum = 311; // stokes problem
    string domain = "newtorus.wrl"; // domain mesh
    string polynomial_patch_filename = "explicit_torus_patches.poly";

    string output_folder = "output/test_greens_identity/";
    int num_iters = 6;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::POLYNOMIAL,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("Test qbkix Green's Identity Stokes: cube", "[results][gid][stokes][cube]"){
    int patch_order = 20;  // higher order patches to more accurately resolve 
    //kinks in normals at patch corners
    int patch_refinement_factor = 2; // a bit of uniform fitting to better resolve patch corners
    int kernel_enum = 311; // stokes problem
    string domain = "cube.wrl"; // domain mesh
    string output_folder = "output/test_greens_identity/";
    int num_iters = 5;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_greens_identity,
            &greens_identity_base_options,
            string(),
            num_iters);
}

void test_greens_identity_mpi(PatchSurfFaceMap* surface, string output_folder, int i){
    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::HARMONIC;
    /*test.singularity_type= SingularityType::SINGLE;
        test.single_singularity_location = Point3(2.,2.,2.);*/
    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 1.;

    // Evaluate solution at the collocation points
    test.target_type   = TargetType::COLLOCATION_POINTS;
    //test.target_type   = TargetType::RING;
    //test.sphere_radius_targets= .5;
    // ... via Green's Identity
    test.evaluation_scheme = EvaluationScheme::GREENS_IDENTITY;
    // no solve step need
    test.solution_scheme   = SolutionScheme::NONE;
    test.solver_matvec_type = NONE;

    test.dump_values = true;

    string filename = string("qbx_greens_id_ref_lvl_")+to_string(i);
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

void greens_identity_base_options_mpi(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "8");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".04");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".06");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "1");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".115");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".115");
    // one level of upsampling, good  convergence rates on cube and torus
    /*Options::set_value_petsc_opts("-boundary_distance_ratio", ".135");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".135");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    */
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".090909090909");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".090909090909");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".058823");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".058823");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "8000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

}

TEST_CASE("Test qbkix Green's Identity Laplace mpi: cube", "[mpi-green][cube]"){
    int patch_order = 12;  // higher order patches to more accurately resolve 
    //kinks in normals at patch corners
    int patch_refinement_factor = 0; // a bit of uniform fitting to better resolve patch corners
    int kernel_enum = 111; // laplace problem
    string domain = "cube.wrl"; // domain mesh
    string output_folder = "output/test_greens_identity/";
    int num_iters = 1;

    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::BLENDED,
            domain,
            output_folder,
            &test_greens_identity_mpi,
            &greens_identity_base_options_mpi,
            string(),
            num_iters);

}
//TEST_CASE
