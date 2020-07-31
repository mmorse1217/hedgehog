#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include <p4est_interface.hpp>
#include "bie3d/solver_utils.hpp"
#include "common/ebi_uti.hpp"
#include <bitset>
using namespace hedgehog;

unique_ptr<PatchSurfFaceMap> build_surface(PatchSurfFaceMap::SurfaceType t){
        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = t;
        surface->_coarse = true;
        surface->setFromOptions();
        surface->setup();
        return surface;
}

void set_cube_config(){

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", ".5");
}

void set_torus_config(){

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/newtorus.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", ".5");
}

void set_pipe_config(){

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-upsampling_type", "adaptive");
    Options::set_value_petsc_opts("-adaptive_upsampling", "bbox");
    Options::set_value_petsc_opts("-adaptive_upsampling_switch_iter", "2");
    Options::set_value_petsc_opts("-adaptive_upsampling_bbox_inflation_factor", ".5");
}

TEST_CASE("Test p4est refinement", "[geom][p4est]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    /*
       PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
       face_map->_surface_type = PatchSurfFaceMap::BLENDED;

    PatchSamples* patch_samples = new PatchSamples("","");
    PatchSamples* refined_patch_samples = new PatchSamples("","");

    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();
    face_map->save("cube_base");

    vector<int> patch_partition(face_map->patches().size(), 0 );

    patch_samples->bdry() = face_map;
    patch_samples->patch_partition() = patch_partition;
    patch_samples->setup();
    SECTION("check subpatches cover face-map"){
        PatchSurfFaceMap* loaded_face_map = PatchSurfFaceMap::load("cube_base");
        PatchSamples* loaded_patch_samples = new PatchSamples("","");

        CHECK(loaded_face_map->patches().size() == face_map->patches().size());

        vector<int> patch_partition(loaded_face_map->patches().size(), 0 );

        loaded_patch_samples->bdry() = face_map;
        loaded_patch_samples->patch_partition() = patch_partition;
        loaded_patch_samples->setup();
        
        int num_samples = patch_samples->local_num_sample_points();
        int loaded_num_samples = loaded_patch_samples->local_num_sample_points();
        CHECK(num_samples == loaded_num_samples); 
        
        DblNumMat samples = 
            get_local_vector(DIM, num_samples, patch_samples->sample_point_3d_position());
        DblNumMat normals = 
            get_local_vector(DIM, num_samples, patch_samples->sample_point_normal());
        
        DblNumMat loaded_samples = 
            get_local_vector(DIM, loaded_num_samples, loaded_patch_samples->sample_point_3d_position());
        DblNumMat loaded_normals = 
            get_local_vector(DIM, loaded_num_samples, loaded_patch_samples->sample_point_normal());

        for(int i =0; i < num_samples; i++){
            Point3 sample(samples.clmdata(i));
            Point3 loaded_sample(loaded_samples.clmdata(i));
            Point3 normal(normals.clmdata(i));
            Point3 loaded_normal(loaded_normals.clmdata(i));
            CHECK((normal - loaded_normal).length() <= 1e-15);
            CHECK((sample - loaded_sample).length() <= 1e-15);

        }
        samples.restore_local_vector();
        loaded_samples.restore_local_vector();

        delete loaded_patch_samples;
        delete loaded_face_map;
        // Test that the domains of the p4est patches cover the domain of the
        // root level face-map patch [0,1]^2

        // before the end of the test case, form a PatchSample from the refined
        // FaceMap and dump to a file, to view in python to ensure that sampling
        // pattern seems correct.
    }
    SECTION("check 2-to-1 balancing"){
        // Test that 2-to-1 balancing is enforced by checking neighbor 2d domain 
        // sizes
        // MJM TODO expose p4est neighbor information
        // uses p4est_quadrant_face_neighbor_extra to do neighbor across a face
        // cross quadtrees if needed and
        // p4est_quadrant_all_face_neighbors to find all possible sized
        // neighbors across a face within a single octree. Need to combine the
        // two in a non-trivial way
    }
    SECTION(""){
    }
    SECTION("Check coordinate mapping from face-map to p4est quads"){
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
        face_map->_surface_type = PatchSurfFaceMap::BLENDED;


        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        
        // Refine until points are inside
        //p4est_connectivity_t* connectivity = build_connectivity_from_face_map(face_map);
        //p4est_t* p4est = p4est_new(MPI_COMM_WORLD, connectivity, sizeof(RefinementData<FaceMapSubPatch>), NULL, NULL);
        p4est_t* p4est = face_map->_p4est;
        PointInQuad piq;
        p4est->user_pointer = &piq;

        refine_patches_uniform(4, p4est, face_map);
        for(unsigned ti = 0; ti < p4est->trees->elem_count; ti++){

            p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, ti);

            for(unsigned qi = 0; qi < tree->quadrants.elem_count; qi++){
                p4est_quadrant_t* quad = p4est_quadrant_array_index(&tree->quadrants, qi);
                int8_t level = quad->level;
                double x= p4est_quad_to_face_map_coords(quad->x, level);
                double y= p4est_quad_to_face_map_coords(quad->y, level);
                int32_t qcoord = face_map_to_p4est_quad_coords(x+double(rand()/RAND_MAX), level);
                //CHECK(quad->x == qcoord);

                // need the offset; without it, the point hashes to the next
                // box.
                double quad_width = 1./pow(2, level);
                qcoord = face_map_to_p4est_quad_coords(x + quad_width - 1e-15, level);
                //CHECK(quad->x == qcoord);
                
                vector<OnSurfacePoint> points_to_process;
                for(int i =0; i < 2; i++){
                    for(int j =0; j < 2; j++){
                        OnSurfacePoint o;
                        double shift = quad_width - 1e-5;
                        o.parametric_coordinates = Point2(x + i*shift,y + j*shift);
                        o.parent_patch = ti;
                        points_to_process.push_back(o);
                    }
                }
                find_leaf_quad_containing_points(p4est, points_to_process);
                for(int i =0; i < points_to_process.size(); i++){
                    p4est_quadrant_t* quad_search = ((PointInQuad*)p4est->user_pointer)->quadrant_containing_points[i];
                    CHECK(quad->x  == quad_search->x);
                    CHECK(quad->y  == quad_search->y);
                    CHECK(quad->level  == quad_search->level);
                }
            }
        }
    }*/
}
/*
TEST_CASE("P4est refinement regression test initialization", "[geom][p4est][regression-init]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");

    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-4");
    auto s = build_surface(PatchSurfFaceMap::BLENDED);
}

TEST_CASE("p4est refinement regression test: base forest", "[geom][p4est][base][regression]"){
}
TEST_CASE("p4est refinement regression test: qbkix refinement", "[geom][p4est][qbkix-refine][regression]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
}
TEST_CASE("p4est refinement regression test: adaptive upsampling", "[geom][p4est][upsample][regression]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "");
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
}
*/
