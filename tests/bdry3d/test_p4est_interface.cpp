
#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bdry3d/p4est_refinement.hpp"
#include <p4est_interface.hpp>
#include "bie3d/solver_utils.hpp"
#include "common/ebi_uti.hpp"
#include <bitset>
using namespace hedgehog;

bool are_intervals_equal(Interval i1, Interval i2){
    return fabs(i1.first - i2.first) <= 1e-15 &&
            fabs(i1.second - i2.second) <= 1e-15;

}
bool is_interval_smaller(Interval small, Interval big){

    return small.first - big.first > 0. &&
           small.second- big.second< 0.;

}

bool are_patches_equal(FaceMapSubPatch* p1, FaceMapSubPatch* p2){
    return are_intervals_equal(p1->_x_interval, p2->_x_interval) && 
           are_intervals_equal(p1->_y_interval, p2->_y_interval) && 
           p1->_face_map_patch == p2->_face_map_patch &&
           p1->_id == p2->_id &&
           p1->_level == p2->_level &&
           p1->_parent_id == p2->_parent_id &&
           p1->_coarse_parent_patch == p2->_coarse_parent_patch &&
           p1->_quad_id_within_p4est_tree == p2->_quad_id_within_p4est_tree &&
           p1->_quadrant == p2->_quadrant &&
           p1 == p2;

}

bool is_patch_child_of_parent(FaceMapSubPatch* child, FaceMapSubPatch* parent){
    return is_interval_smaller(child->_x_interval, parent->_x_interval) && 
           is_interval_smaller(child->_y_interval, parent->_y_interval) && 
           child->_face_map_patch == parent->_face_map_patch &&
           child->_id                    != parent->_id &&
           child->_level                 > parent->_level &&
           child->_parent_id             != parent->_parent_id &&
           child->_coarse_parent_patch   == parent->_coarse_parent_patch &&
           child->_quad_id_within_p4est_tree != parent->_quad_id_within_p4est_tree &&
           child->_quadrant              != parent->_quadrant;

}


TEST_CASE("Test p4est-related functions", "[geom][p4est]"){

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    // dump the points and patches 
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");

    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->setFromOptions();
    face_map->setup();
    SECTION("Test p4est_to_face_map_subpatches on an unrefined face-map"){

        initialize_p4est_leaves_with_subpatches(face_map->_p4est, face_map);
        vector<Patch*> computed_subpatches= 
            p4est_to_face_map_subpatches(face_map->_p4est, face_map);
        vector<Patch*> true_patches = face_map->patches();
        //face_map->refine_test();
        

        // make sure the subpatch domains are all [0,1]^2
        for(auto p : computed_subpatches){
            FaceMapSubPatch* subpatch = dynamic_cast<FaceMapSubPatch*>(p);
            
            assert(subpatch != NULL);
            CHECK(are_intervals_equal(subpatch->_x_interval, Interval(0.,1.)));
            CHECK(are_intervals_equal(subpatch->_y_interval, Interval(0.,1.)));
        }
        // make sure that we made the same number of patches
        int num_true_patches = true_patches.size();
        int num_computed_patches = computed_subpatches.size();
        REQUIRE(num_true_patches == num_computed_patches);

        // now overwrite FaceMapPatch's with equal FaceMapSubPatch's.
        face_map->refine_test();
        true_patches = face_map->patches();

        // make sure they're equivalent patches.
        for(int pi=0; pi < num_true_patches; pi++){
            auto p = dynamic_cast<FaceMapSubPatch*>(true_patches[pi]);
            auto new_p = dynamic_cast<FaceMapSubPatch*>(computed_subpatches[pi]);
            CHECK(are_patches_equal(p, new_p));
        }
    } SECTION("Test update_face_map, no refinement"){

        PatchSurfFaceMap* intermediate_face_map = nullptr;
        PatchSamples* intermediate_patch_samples = nullptr;
        face_map->refine_test();

        // use update_face_map to clone patches to intermediate_face_map, and
        // sample intermediate_patch_samples from intermediate_face_map.
        update_face_map(face_map->_p4est,
                face_map,
                intermediate_face_map,
                intermediate_patch_samples);
       
        int num_updated_patches = intermediate_face_map->num_patches();
        int num_init_patches = face_map->num_patches();
        
        // They should now have the same number of patches...
        REQUIRE(num_updated_patches == num_init_patches);

        // but different p4est trees. due to deep copy
        REQUIRE(intermediate_face_map->_p4est != face_map->_p4est);
        
        // make sure all patches are equivalent
        for(int pi=0; pi < num_init_patches; pi++){
                auto p = face_map->subpatch(pi);
                auto updated_p = intermediate_face_map->subpatch(pi);
           CHECK(are_patches_equal(p, updated_p));
        }


    } SECTION("Test update_face_map, uniform refinement of base face-map"){

        PatchSurfFaceMap* intermediate_face_map = nullptr;
        PatchSamples* intermediate_patch_samples = nullptr;
        int num_init_patches = face_map->num_patches();

        // repeatedly update intermediate_face_map/patch_samples after another
        // level of refinement of the original face-map.
        for(int ref_level = 1; ref_level < 4; ref_level++){
            // refine face-map by one level 
            face_map->refine_uniform(1);
            // update intermediate_face_map with the patches given in face-map
            // after refinement
            update_face_map(face_map->_p4est,
                    face_map,
                    intermediate_face_map,
                    intermediate_patch_samples);

            int num_updated_patches = intermediate_face_map->num_patches();

            // make sure they are the same 
            int num_patches_per_ref_level = pow(4,ref_level);
            REQUIRE(num_updated_patches == num_patches_per_ref_level*num_init_patches);
            
            for(int pi=0; pi < num_updated_patches; pi++){

                auto p = face_map->subpatch(pi);
                auto updated_p = intermediate_face_map->subpatch(pi);
                CHECK(are_patches_equal(p, updated_p));
            }
        }
    } SECTION("Test update_face_map, make sure updated patches match p4est"){

        face_map->refine_test();
        PatchSurfFaceMap* intermediate_face_map = nullptr;
        PatchSamples* intermediate_patch_samples = nullptr;

        int num_init_patches = face_map->num_patches();

        p4est_t* p4est = face_map->_p4est;
        update_face_map(p4est,
                face_map,
                intermediate_face_map,
                intermediate_patch_samples);

        // repeatedly update intermediate_face_map/patch_samples after another
        // level of refinement of the previous iteration's intermediate_face_map
        for(int ref_level = 1; ref_level < 4; ref_level++){
            cout << "refinement level: " << ref_level << endl;
            // refine the p4est by one level 
            refine_patches_uniform(1, p4est, face_map);
            
            // update intermediate_face_map with the patches contained in p4est
            update_face_map(p4est,
                    face_map,
                    intermediate_face_map,
                    intermediate_patch_samples);

            // make we have the right number of patches first
            int num_updated_patches = intermediate_face_map->num_patches();
            int num_patches_per_ref_level = pow(4,ref_level);
            REQUIRE(num_updated_patches == num_patches_per_ref_level*num_init_patches);

            
            // get the actual patches from p4est
            auto p4est_leaf_patches = collect_face_map_subpatches(p4est);
            REQUIRE(p4est_leaf_patches.size() == num_updated_patches);

            // make everything matches.
            for(int pi=0; pi < num_updated_patches; pi++){

                auto p4est_patch = p4est_leaf_patches[pi];
                auto updated_face_map_patch= intermediate_face_map->subpatch(pi);
                REQUIRE(p4est_patch);
                REQUIRE(updated_face_map_patch);
                CHECK(are_patches_equal(p4est_patch, updated_face_map_patch));
            }
        }

    }

}

TEST_CASE("Test p4est-related functions, non-uniform", "[geom][p4est][fff]"){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/small_flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/small_flat_patch.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/small_flat_patch.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".05");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    // dump the points and patches 
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");

    //auto face_map= unique_ptr<PatchSurfFaceMap>(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    auto face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();
    SECTION("Test update_face_map, refine a single patch"){

        face_map->refine_test();
        PatchSurfFaceMap* intermediate_face_map = nullptr;
        PatchSamples* intermediate_patch_samples = nullptr;

        int num_init_patches = face_map->num_patches();

        p4est_t* p4est = face_map->_p4est;
        update_face_map(p4est,
                face_map,
                intermediate_face_map,
                intermediate_patch_samples);
        
        int num_new_patches = 0;
        for (int i = 0; i < 29; i++) {
            // refine one quad 
            mark_all_patches_for_refinement<FaceMapSubPatch>(p4est, false);
            auto patch = intermediate_face_map->subpatch(0);
            auto ref_data = get_refinement_data(patch->_quadrant);
            ref_data->refine = true;
            refine_p4est_quads(p4est);
            
            // update intermediate_face_map
            update_face_map(p4est,
                    face_map,
                    intermediate_face_map,
                    intermediate_patch_samples);

            // make we have the right number of patches first
            int num_updated_patches = intermediate_face_map->num_patches();
            num_new_patches += 4 -1;
            REQUIRE(num_updated_patches == num_init_patches + num_new_patches);
            check_patches(p4est);
             
        }

    }
}


