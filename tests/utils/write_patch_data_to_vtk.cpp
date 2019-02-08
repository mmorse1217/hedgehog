#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bie3d/markgrid.hpp"
#include "common/utils.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bie3d/markgrid.hpp"
#include "bie3d/spatial_grid.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "common/vtk_writer.hpp"
extern "C" {
#include <p4est_vtk.h>
}
using namespace Ebi;
using namespace Markgrid;

TEST_CASE("Visualize QBKIX point generation by dumping to vtk", "[vtk][visualize][utils]"){

        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "6");
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

        // Refine until points are inside
        int n = face_map->patches().size();
        vector<int> patches_refined_relative_ids(n,0);
        for(int i =0; i < n; i++){
            patches_refined_relative_ids[i] = i;
        }
        int it = 0;
        // geometry refinement with p4est works
            write_face_map_patches_to_vtk(DblNumMat(0,0), 
                    patches_refined_relative_ids,
                    face_map, it++);
        // renderer works with FaceMapSubPatch in p4est leaves
            initialize_p4est_leaves_with_subpatches(face_map->_p4est, face_map);
            write_face_map_patches_to_vtk(DblNumMat(0,0), 
                    patches_refined_relative_ids,
                    face_map, it++);
        // renderer works after collecting all FaceMapSubPatch's from the forest  and dumping
        // them into face_map
            face_map->patches() = p4est_to_face_map_subpatches(face_map->_p4est, face_map);
            initialize_p4est_leaves_with_subpatches(face_map->_p4est, face_map);
            write_face_map_patches_to_vtk(DblNumMat(0,0), 
                    patches_refined_relative_ids,
                    face_map, it++);
        /*
        p4est_connectivity_t* connectivity = build_connectivity_from_face_map(face_map);
        p4est_t* p4est = p4est_new(MPI_COMM_WORLD, connectivity, sizeof(RefinementData), NULL, NULL);


        store_p4est("cube_base", p4est);
        refine_patches_for_qbkix_point_location(p4est, face_map);
        store_p4est("cube_coarse_patches", p4est);
        refine_patches_for_fixed_qbkix_points(p4est, face_map);
        store_p4est("cube_refined_patches", p4est);
        p4est_t* p4est_coarse;
        p4est_t* p4est_refined;
        p4est_connectivity_t* connectivity_coarse;
        p4est_connectivity_t* connectivity_refined;

        load_p4est("cube_coarse_patches", face_map->mpiComm(),p4est_coarse, &connectivity_coarse);
        load_p4est("cube_refined_patches", face_map->mpiComm(),p4est_refined, &connectivity_refined);
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(p4est, face_map);
        //_patches = (vector<Patch*>)subpatches;
        face_map->patches() = subpatches;
        */
}

