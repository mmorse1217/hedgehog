#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
#include "common/stats.hpp"
#include "evaluation_utils.hpp"
#include <numeric>
using namespace hedgehog;

TEST_CASE("dump qbkix points for paraview", "[qbkix-points][debug]"){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "3");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "1000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");

    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type = PatchSurfFaceMap::SurfaceType::POLYNOMIAL;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_uniform(1);

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(surface->num_patches(), 0);
    samples->bdry() = surface.get();
    samples->setup();

    vector<int> qbkix_ids(Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    vector<int> patches_to_sample(surface->num_patches());
    std::iota(qbkix_ids.begin(), qbkix_ids.end(), 0);
    std::iota(patches_to_sample.begin(), patches_to_sample.end(), 0);
    
    /*Vec qbkix_points =  samples->generate_qbkix_points_from_sample_points(
            qbkix_ids, patches_to_sample);*/
    Vec normals;
    VecDuplicate(samples->sample_point_normal(), &normals);
    VecCopy(samples->sample_point_normal(), normals);
    VecScale(normals,-1);
   Vec qbkix_points= 
        generate_interior_qbkix_points(
                samples->sample_point_3d_position(),
                samples->sample_point_3d_position(),
                qbkix_ids, 
                normals, 
                samples->sample_point_far_field(),
                samples->sample_point_interpolant_spacing()
                );

    Vec values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, Petsc::get_vec_size(qbkix_points)/3, values);
    write_general_points_to_vtk(qbkix_points, 1, 
            "qbkix_points.vtp", values, "output/");

}

TEST_CASE("Test marking a plane of targets", "[debug][sheet]"){
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".75");

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-kt", "111");
    //Options::set_value_petsc_opts("-boundary_distance_ratio", ".04");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".11");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".022");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");

    Options::set_value_petsc_opts("-adaptive", "0");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "2");

    stats._file_prefix = "data/";
    string output_folder= "output/test_plane_targets/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";


    unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    surface->_surface_type =  PatchSurfFaceMap::BLENDED;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_test();
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
            vector<int>(surface->num_patches(), 0) ,
            surface.get(), 0, output_folder);

    // Set up test case
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    test.bc_type = BoundaryDataType::CONSTANT_DENSITY;

    
    test.target_type   = TargetType::PLANE;
    //test.target_plane_point= Point3(-.05,0., .45);
    test.target_plane_point= Point3(0., 0., 0.);
    test.target_plane_vec1 = Point3(1., 0., 0.);
    test.target_plane_vec2 = Point3(0., 1., 0.);
    test.num_targets = 50;
    
    test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
    //test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
    test.solution_scheme   = SolutionScheme::EXPLICIT_DENSITY;
    test.solver_matvec_type = EXTRAPOLATION_AVERAGE;
    test.dump_values = true;

    string filename = string("qbx_target_plane_solution");
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
