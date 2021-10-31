#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bie3d/markgrid.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bdry3d/p4est_refinement.hpp"
#include "common/nummat.hpp"
#include "common/vtk_writer.hpp"
#include "../utils/evaluation_utils.hpp"
#include <algorithm>
extern "C" {
#include <p4est_vtk.h>
}
using namespace hedgehog;
using namespace Markgrid;


TEST_CASE("Patch refinement", "[geom][patch-refinement][cube]"){


    SECTION("Test refine patches until all points are inside domain"){
        //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/two_small_flat_patch.poly");
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "6");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    //Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","1.");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","0.");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");

        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        //face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        face_map->setFromOptions();
        face_map->_coarse = true;
        face_map->setup();

        // Refine until points are inside

        double start = omp_get_wtime();
        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_qbkix_refinement");
        refine_patches_for_qbkix_point_location(face_map->_p4est, face_map);
        stats.add_result("qbx ref time", omp_get_wtime() - start);



        vector<Patch*> subpatches = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        face_map->patches() = subpatches;


        // check that they're actually inside...
        vector<int> partition(face_map->patches().size(), 0);
        PatchSamples* test_samples = new PatchSamples("", "");
        test_samples->bdry() = face_map;
        test_samples->patch_partition() = partition;
        test_samples->setup();

        // generate qbkix targets from samples
        vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
        vector<int> patch_ids;
        int num_patches = face_map->patches().size();
        for(int i =0; i < num_patches; i++)
            patch_ids.push_back(i);
        Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

        NumVec<OnSurfacePoint> on_surface_points;
        /*NumVec<OnSurfacePoint> on_surface_points = 
            Markgrid::mark_target_points(qbkix_points_local, face_map);

        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).inside_domain == INSIDE);
        }*/
        stats.print_results();


        /*int num_qbkix_points_per_patch = num_qbkix_samples/num_patches;
          for(int pi =0; pi < num_patches; pi++){
          for(int qi =0; qi < num_qbkix_points_per_patch; qi++){
          int index = pi*num_qbkix_points_per_patch + qi;
          CHECK(on_surface_points(index).parent_patch == pi);
          }
          }*/
        /*
        PatchSurfFaceMap* refined_face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        //refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        refined_face_map->setFromOptions();
        refined_face_map->setup_from_existing_face_map(face_map);
        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_after_medial_axis_refinement");
        p4est_vtk_write_file(refined_face_map->_p4est, NULL, "squished_cube_coarse_face_map");
        refine_patches_for_fixed_qbkix_points(refined_face_map->_p4est, refined_face_map);
        p4est_vtk_write_file(refined_face_map->_p4est, NULL, "squished_cube_refined_face_map");


        on_surface_points = 
            Markgrid::mark_target_points(qbkix_points_local, refined_face_map);
        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).region == FAR);
        }
        stats.print_results();

*/
    }
}

TEST_CASE("Patch upsampling", "[geom][patch-upsampling][cube]"){


        //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/two_small_flat_patch.poly");
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "6");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    //Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","1.");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".5");
    Options::set_value_petsc_opts("-adaptive_upsampling","bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter","2");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
        Options::set_value_petsc_opts("-bis3d_spacing", ".333");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");

        stats._file_prefix = "data/";

        PatchSurfFaceMap* refined_face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        refined_face_map->setFromOptions();
        refined_face_map->_coarse = true;
        refined_face_map->setup();
        //refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        //refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        //refined_face_map->setFromOptions();
        //refined_face_map->setup_from_existing_face_map(face_map);
        //p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_after_medial_axis_refinement");
        //p4est_vtk_write_file(refined_face_map->_p4est, NULL, "squished_cube_coarse_face_map");
        double start = omp_get_wtime();

        refine_patches_for_fixed_qbkix_points(refined_face_map->_p4est, refined_face_map);
        stats.add_result("qbx ref time", omp_get_wtime() - start);
        
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(refined_face_map->_p4est, refined_face_map);
        refined_face_map->patches() = subpatches;

        int num_patches = refined_face_map->patches().size();
        
        // generate qbkix targets from samples
        vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
        vector<int> patch_ids;
        for(int i =0; i < num_patches; i++)
            patch_ids.push_back(i);

        stats.print_results();

        // check that they're actually inside...
        vector<int> partition(refined_face_map->patches().size(), 0);
        PatchSamples* test_samples = new PatchSamples("", "");
        test_samples->bdry() = refined_face_map;
        test_samples->patch_partition() = partition;
        test_samples->setup();

        Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        DblNumMat qbkix_points_local(DIM, qbkix_points);

        NumVec<OnSurfacePoint> on_surface_points = 
            Markgrid::mark_target_points(qbkix_points_local, refined_face_map);

        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).region == FAR);
        }
        stats.print_results();


}

TEST_CASE("Fully adaptive", "[geom][fully-adaptive][cube]"){


        //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/two_small_flat_patch.wrl");
        //Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/two_small_flat_patch.poly");
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/new_ppp.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/new_ppp.wrl");
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "6");
        //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    //Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","1.");
    //Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".25");
    Options::set_value_petsc_opts("-adaptive_upsampling","bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter","2");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
        Options::set_value_petsc_opts("-bis3d_spacing", ".25");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");

        stats._file_prefix = "data/";

        PatchSurfFaceMap* refined_face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        refined_face_map->setFromOptions();
        refined_face_map->_coarse = true;
        refined_face_map->setup();
        //refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        //refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        //refined_face_map->setFromOptions();
        //refined_face_map->setup_from_existing_face_map(face_map);
        //p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_after_medial_axis_refinement");
        //p4est_vtk_write_file(refined_face_map->_p4est, NULL, "squished_cube_coarse_face_map");
        double start = omp_get_wtime();

        refine_patches_for_qbkix_point_location(refined_face_map->_p4est, refined_face_map);
        refine_patches_for_fixed_qbkix_points(refined_face_map->_p4est, refined_face_map);
        stats.add_result("qbx ref time", omp_get_wtime() - start);
        
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(refined_face_map->_p4est, refined_face_map);
        refined_face_map->patches() = subpatches;
/*
        int num_patches = refined_face_map->patches().size();
        
        // generate qbkix targets from samples
        vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
        vector<int> patch_ids;
        for(int i =0; i < num_patches; i++)
            patch_ids.push_back(i);

        stats.print_results();

        // check that they're actually inside...
        vector<int> partition(refined_face_map->patches().size(), 0);
        PatchSamples* test_samples = new PatchSamples("", "");
        test_samples->bdry() = refined_face_map;
        test_samples->patch_partition() = partition;
        test_samples->setup();

        Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        DblNumMat qbkix_points_local(DIM, qbkix_points);

        NumVec<OnSurfacePoint> on_surface_points = 
            Markgrid::mark_target_points(qbkix_points_local, refined_face_map);

        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).region == FAR);
        }
        stats.print_results();
*/

}



TEST_CASE("Test refinement on two flat patches very close together", "[geom][patch-refinement][flat-patch]"){
    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/two_near_flat_patches.wrl");
    PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/two_near_flat_patches.poly");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/two_near_flat_patches.wrl");
    PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "2");

    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();

    // Refine until points are inside


        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_qbkix_refinement");
    refine_patches_for_qbkix_point_location(face_map->_p4est, face_map);

        vector<Patch*> subpatches = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        face_map->patches() = subpatches;
    //_patches = (vector<Patch*>)subpatches;
    face_map->patches() = subpatches;


    // check that they're actually inside...
    vector<int> partition(face_map->patches().size(), 0);
    PatchSamples* test_samples = new PatchSamples("", "");
    test_samples->bdry() = face_map;
    test_samples->patch_partition() = partition;
    test_samples->setup();

    // generate qbkix targets from samples
    vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
    vector<int> patch_ids;

    for(int i =0; i < face_map->patches().size(); i++)
        patch_ids.push_back(i);

    Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

    int64_t num_qbkix_samples;
    VecGetSize(qbkix_points, &num_qbkix_samples);
    num_qbkix_samples /= DIM;
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

    NumVec<OnSurfacePoint> on_surface_points = 
        Markgrid::mark_target_points(qbkix_points_local, face_map);

    for(int i = 0; i < on_surface_points.m(); i++){
        CHECK(on_surface_points(i).inside_domain == INSIDE);
    }

    // should only need three levels of refinement to land qbkix points in
    // the proper place for this case. check that it holds.
    for(int pi = 0; pi < face_map->patches().size(); pi++){
        FaceMapSubPatch* patch = (FaceMapSubPatch*)face_map->patches()[pi];
        CHECK(patch->_level == 3);
    }
}

TEST_CASE("Test adaptive upsampling on patches from previous test", "[geom][patch-upsampling][flat-patch]"){
    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/two_near_flat_patches.wrl");
    PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/two_near_flat_patches.poly");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/two_near_flat_patches.wrl");
    PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "4");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","1.");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter ","3");
    Options::set_value_petsc_opts("-target_accuracy","1e-3");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".25");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_uniform(3); // result of previous test

    // Refine until points are inside


    stats._file_prefix = "data/";
        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_qbkix_refinement");
        refine_patches_for_fixed_qbkix_points(face_map->_p4est, face_map);
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        face_map->patches() = subpatches;



    // check that they're actually inside...
    // sample qbkix points from refined surface 
    vector<int> partition(face_map->patches().size(), 0);
    PatchSamples* test_samples = new PatchSamples("", "");
    test_samples->bdry() = face_map;
    test_samples->patch_partition() = partition;
    test_samples->setup();

    // generate qbkix targets from samples
    vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
    vector<int> patch_ids;

    for(int i =0; i < face_map->patches().size(); i++)
        patch_ids.push_back(i);

    Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

    int64_t num_qbkix_samples;
    VecGetSize(qbkix_points, &num_qbkix_samples);
    num_qbkix_samples /= DIM;
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

    // check that they are far from the refined geometry
    NumVec<OnSurfacePoint> on_surface_points = 
        Markgrid::mark_target_points(qbkix_points_local, face_map);

    for(int i = 0; i < on_surface_points.m(); i++){
        CHECK(on_surface_points(i).region == FAR);
    }

    // should only need three levels of refinement to land qbkix points in
    // the proper place for this case. check that it holds.
    /*for(int pi = 0; pi < face_map->patches().size(); pi++){
        FaceMapSubPatch* patch = (FaceMapSubPatch*)face_map->patches()[pi];
        CHECK(patch->_level == 3);
    }*/
}


TEST_CASE("Test adaptive upsampling on curved patches ", "[geom][patch-upsampling][curved-patch]"){

    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/two_parallel_curved_patches.wrl");
    PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/two_parallel_curved_patches.poly");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/two_parallel_curved_patches.wrl");
    PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "2");
    //Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","0.");


    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","1.");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter ","3");
    Options::set_value_petsc_opts("-target_accuracy","1e-3");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();
    refine_patches_for_qbkix_point_location(face_map->_p4est, face_map);

    // Refine until points are inside


    stats._file_prefix = "data/";
        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_qbkix_refinement");
        refine_patches_for_fixed_qbkix_points(face_map->_p4est, face_map);
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        face_map->patches() = subpatches;



    // check that they're actually inside...
    // sample qbkix points from refined surface 
    vector<int> partition(face_map->patches().size(), 0);
    PatchSamples* test_samples = new PatchSamples("", "");
    test_samples->bdry() = face_map;
    test_samples->patch_partition() = partition;
    test_samples->setup();

    // generate qbkix targets from samples
    vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
    vector<int> patch_ids;

    for(int i =0; i < face_map->patches().size(); i++)
        patch_ids.push_back(i);

    Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

    int64_t num_qbkix_samples;
    VecGetSize(qbkix_points, &num_qbkix_samples);
    num_qbkix_samples /= DIM;
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

    // check that they are far from the refined geometry
    NumVec<OnSurfacePoint> on_surface_points = 
        Markgrid::mark_target_points(qbkix_points_local, face_map);

    for(int i = 0; i < on_surface_points.m(); i++){
        CHECK(on_surface_points(i).region == FAR);
    }

    // should only need three levels of refinement to land qbkix points in
    // the proper place for this case. check that it holds.
    /*for(int pi = 0; pi < face_map->patches().size(); pi++){
        FaceMapSubPatch* patch = (FaceMapSubPatch*)face_map->patches()[pi];
        CHECK(patch->_level == 3);
    }*/
}



TEST_CASE("Test refinement on two curved patches very close together", "[geom][patch-refinement][curved-patch]"){
    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/two_parallel_curved_patches.wrl");
    PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/two_parallel_curved_patches.poly");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/two_parallel_curved_patches.wrl");
    PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor","0.");

    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();

    // Refine until points are inside


        p4est_vtk_write_file(face_map->_p4est, NULL, "squished_cube_qbkix_refinement");
    refine_patches_for_qbkix_point_location(face_map->_p4est, face_map);

        vector<Patch*> subpatches = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        face_map->patches() = subpatches;
    //_patches = (vector<Patch*>)subpatches;
    face_map->patches() = subpatches;

    {
        // check that they're actually inside...
        vector<int> partition(face_map->patches().size(), 0);
        PatchSamples* test_samples = new PatchSamples("", "");
        test_samples->bdry() = face_map;
        test_samples->patch_partition() = partition;
        test_samples->setup();

        // generate qbkix targets from samples
        vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
        vector<int> patch_ids;

        for(int i =0; i < face_map->patches().size(); i++)
            patch_ids.push_back(i);

        Vec qbkix_points = test_samples->generate_qbkix_points_from_sample_points(indices, patch_ids);

        int64_t num_qbkix_samples;
        VecGetSize(qbkix_points, &num_qbkix_samples);
        num_qbkix_samples /= DIM;
        DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_samples, qbkix_points);

        NumVec<OnSurfacePoint> on_surface_points = 
            Markgrid::mark_target_points(qbkix_points_local, face_map);

        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).inside_domain == INSIDE);
        }
    }

}

TEST_CASE("right hand side refinement", "[geom][refine-rhs][cube]"){


        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "6");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "10");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");



        PatchSurfFaceMap* refined_face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        refined_face_map->_surface_type = PatchSurfFaceMap::BLENDED;
        refined_face_map->setFromOptions();
        refined_face_map->_coarse = true;
        refined_face_map->setup();
        refined_face_map->refine_test();
        double start = omp_get_wtime();
        
        


        //refine_patches_for_fixed_qbkix_points(refined_face_map->_p4est, refined_face_map);
        resolve_function(refined_face_map->_p4est, refined_face_map, &laplace_singluarity_cube, 1, 1e-10);

        stats.add_result("qbx ref time", omp_get_wtime() - start);
        
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(refined_face_map->_p4est, refined_face_map);
        PatchSamples* test_samples = new PatchSamples("", "");
        test_samples->bdry() = refined_face_map;
        vector<int> partition(refined_face_map->patches().size(), 0);
        test_samples->patch_partition() = partition;
        test_samples->setup();
                
        
        
        Vec func_values;
        
        
        Petsc::create_mpi_vec(refined_face_map->mpiComm(), test_samples->local_num_sample_points(),func_values); // laplace => tdof = 1
        laplace_singluarity_cube(test_samples->sample_point_3d_position(), 1, func_values);
write_general_points_to_vtk(test_samples->sample_point_3d_position(), 1, "target_function.vtp",
        func_values,
        "data/");
        refined_face_map->patches() = subpatches;

        int num_patches = refined_face_map->patches().size();
        /*
        // generate qbkix targets from samples
        vector<int> indices(1, Options::get_int_from_petsc_opts("-near_interpolation_num_samples") - 1);
        vector<int> patch_ids;
        for(int i =0; i < num_patches; i++)
            patch_ids.push_back(i);
        stats.print_results();

        // check that they're actually inside...

        Vec func_values;
        Petsc::create_mpi_vec(face_map->mpiComm(), test_samples->local_num_sample_points(),func_values); // laplace => tdof = 1
        Vec density;
        Petsc::create_mpi_vec(face_map->mpiComm(), test_samples->local_num_sample_points(),density); // laplace => tdof = 1
        VecSet(density,1.);
        Kernel3d laplace(111, vector<int>(2,1.);
                Vec source;
        Petsc::create_mpi_vec(face_map->mpiComm(), 3,source); 
        VecSet(source, 0.);
        VecSetValue(source, 1, .7, INSERT_VALUES);
        fmm = unique_ptr<FMM>(new PvFMM());
        fmm->initialize_fmm(source, source,test_samples->sample_point_3d_position(),
                 laplace);
        fmm->evaluate(density, func_values);


        for(int i = 0; i < on_surface_points.m(); i++){
            CHECK(on_surface_points(i).region == FAR);
        }
*/
        stats.print_results();


}

void test_greens_identity2(PatchSurfFaceMap* surface, string output_folder, int i){
    // Set up test case
    surface->resolve_rhs(laplace_singluarity_cube, 1, 1e-4);
    TestConfig test;
    // Singularity configuration defining boundary condition; sphere radius 2 of
    // point charges
    //test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
    test.bc_type = BoundaryDataType::HARMONIC;
    /*
    test.singularity_type= SingularityType::SINGLE;
        test.single_singularity_location = Point3(0.,2.,0.);
    */
    test.singularity_type= SingularityType::SPHERE;
    test.sphere_radius_bc = 2.;
    // Evaluate solution at the collocation points
    test.target_type   = TargetType::COLLOCATION_POINTS;
    //test.target_type   = TargetType::RING;
    //test.sphere_radius_targets= .5;
    // ... via Green's Identity
    test.evaluation_scheme = EvaluationScheme::GREENS_IDENTITY;
    // no solve step need
    test.solution_scheme   = SolutionScheme::NONE;
    //test.solution_scheme = SolutionScheme::EXPLICIT_DENSITY;
    //test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
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
    stats.clear();

}

void greens_identity_base_options2(){
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    /*Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts( "-uniform_upsampling_num_levels", "3");
    */
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox_closest_point");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor",".4");

    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

}





TEST_CASE("Debug adaptive rhs Green's Identity laplace : torus", "[results][gid][laplace][debug]"){
    int patch_order = 3;  // higher order patches to more accurately resolve 
    //kinks in normals at patch corners
    int patch_refinement_factor = 0; // a bit of uniform fitting to better resolve patch corners
    int kernel_enum = 111; // stokes problem
    string domain = "newtorus.wrl"; // domain mesh
    string output_folder = "output/test_greens_identity/";
    int num_iters = 1;
    string polynomial_patch_filename = "explicit_torus_patches.poly";


    setup_and_run_face_map_convergence_test(
            patch_order, 
            patch_refinement_factor, 
            kernel_enum,
            PatchSurfFaceMap::POLYNOMIAL,
            domain,
            output_folder,
            &test_greens_identity2,
            &greens_identity_base_options2,
            polynomial_patch_filename,
            num_iters);
}
TEST_CASE("mpi facemap setup", "[geom][patch-refinement][mpi]"){
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    //Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    auto surface = setup_face_map(PatchSurfFaceMap::BLENDED, 0, "uniform");

    unique_ptr<PatchSamples> samples(new PatchSamples("",""));
    samples->patch_partition() = vector<int>(surface->num_patches(), surface->mpiRank());
    samples->bdry() = surface.get();
    samples->setup();
    Vec x;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,samples->local_num_sample_points(), x);
    VecSet(x, 1.);
    string s = "samples_mpi_rank_" + to_string(surface->mpiRank()) + ".vtp";
    write_general_points_to_vtk(samples->sample_point_3d_position(), 1, 
            s, x, "data/");

    Vec targets;
    VecDuplicate(samples->sample_point_3d_position(), &targets);
    VecCopy(samples->sample_point_3d_position(), targets);
    VecAXPY(targets, -.01, samples->sample_point_normal());
    DblNumMat t(3, targets);
    DblNumMat sample_point_3d_position_local(3, samples->sample_point_3d_position());
    DblNumMat x_local(1, x);
    auto on_surface_points = Markgrid::mark_target_points( t,surface.get(),false );
    for(int i =0; i < on_surface_points.m(); i++){
        auto osp = on_surface_points(i);
        Point3 position;
        auto patch = surface->subpatch(osp.parent_patch);
        patch->xy_to_patch_coords(osp.parametric_coordinates, PatchSamples::EVAL_VL, position.array());
        Point3 sample_point(sample_point_3d_position_local.clmdata(i));
        CHECK(fabs((sample_point - position).l2()) <= 1e-8);
        x_local(0,i) = (sample_point - position).l2();

    }
    string ss = "targets_mpi_rank_" + to_string(surface->mpiRank()) + ".vtp";
    write_general_points_to_vtk(targets, 1, 
            ss, x, "data/");

}
