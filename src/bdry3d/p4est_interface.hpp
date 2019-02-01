#ifndef __GPMESH_TO_P4EST_HPP__
#define __GPMESH_TO_P4EST_HPP__

#include "bdry3d/patch_surf_face_map.hpp"
#include <gpmesh.hpp>
#include "bdry3d/face_map_subpatch.hpp"
extern "C" {
#include <p4est_connectivity.h>
#include <p4est.h>
#include <p4est_iterate.h>
#include <p4est_extended.h>
#include <p4est_bits.h>
#include <p4est_vtk.h>
}
BEGIN_EBI_NAMESPACE

FaceMapSubPatch* extract_patch_from_quad(p4est_quadrant_t* quad);

void initialize_p4est_leaves_with_subpatches(p4est_t* p4est, PatchSurfFaceMap* face_map);
vector<Patch*> p4est_to_face_map_subpatches(p4est_t* p4est, PatchSurfFaceMap* face_map);


FaceMapSubPatch* extract_patch_from_quad(p4est_quadrant_t* quad);

RefinementData<FaceMapSubPatch>* get_refinement_data_from_forest(p4est_t* p4est,
        int tree_id, int quad_id_within_tree);
RefinementData<FaceMapSubPatch>* get_refinement_data_from_forest(p4est_t* p4est,
        pair<int,int> tree_and_quad_id);


RefinementData<FaceMapSubPatch>* get_refinement_data(p4est_quadrant_t* quad);
/*
void dump_vtk_data_for_paraview(DblNumMat qbkix_points,
        NumVec<OnSurfacePoint> closest_on_surface_points, int it,
        vector<int> global_patch_ids, PatchSurfFaceMap* face_map );
*/
Rectangle parameter_domain_of_child_quad(int parent_patch_level, 
        p4est_quadrant_t* quad, Rectangle parent_domain);
// refinement
void update_and_resample_face_map(p4est_t* p4est, 
        PatchSurfFaceMap* face_map,
        PatchSurfFaceMap*& intermediate_face_map,
        PatchSamples*& intermediate_patch_samples,
        Vec& qbkix_points,
        vector<int> qbkix_indices);

vector<pair<int, int> > patches_that_need_refinement(p4est_t* p4est);
// refinement
// refinement
vector<FaceMapSubPatch*> face_map_patches_to_refine(p4est_t* p4est);
vector<int> face_map_patch_ids_to_refine(p4est_t* p4est);

void store_p4est(string filename, p4est_t* p4est);

void load_p4est(string filename, 
        MPI_Comm comm,
        //size_t data_size
        p4est_t*& p4est,
        p4est_connectivity_t** connectivity);
void find_leaf_quad_containing_points(p4est_t* p4est, vector<OnSurfacePoint> points_to_process);

struct PointInQuad{
    vector<p4est_quadrant_t*> quadrant_containing_points;
    PointInQuad() {}

    int operator()(p4est_t * p4est,
            p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant,
            p4est_locidx_t local_num,
            void *point){
        pair<int,OnSurfacePoint>* ordered_osp = ((pair<int,OnSurfacePoint>*) point);
        if(which_tree != ordered_osp->second.parent_patch)
            return 0;
        else {
            Point2 p = ordered_osp->second.parametric_coordinates;
            //Interval x_int, y_int;
            Rectangle domain;
            quad_to_face_map_uv_bounds(quadrant, domain);

            Interval x_int = domain.first;
            Interval y_int = domain.second;

            int point_in_quad = 
                (p.x() >=  x_int.first)   &&
                (p.x() <=  x_int.second)  &&
                (p.y() >=  y_int.first)   &&
                (p.y() <=  y_int.second) ;

            // TODO check with Abtin about this piece
            /*
            int point_on_quad_edge = 
                (p.x() - x_int.first) >= EPS    &&
                (p.x() - x_int.second)<= EPS    &&
                (p.y() - y_int.first) >= EPS    &&
                (p.y() - y_int.second)<= EPS   ;*/
            if(point_in_quad /*|| point_on_quad_edge*/){
                quadrant_containing_points[ordered_osp->first] = quadrant;
            }
            return point_in_quad /*|| point_on_quad_edge*/;
        }

    }
};
vector<FaceMapSubPatch*> collect_face_map_subpatches(p4est_t* p4est);

// refinement
void update_quad_indicies(p4est_t* p4est);
// refinement
void set_coarse_patch_ids(p4est_t* p4est);
p4est_t* object_safe_p4est_copy(p4est_t* input);

void refine_p4est_quads(p4est_t* p4est);
/**
 * Given a p4est defining the current patch set and a face-map whose patches 
 * define the result of geometry refinement, return a face_map containing the
 * patches given by p4est and its associated PatchSamples object 
 */
// refinement
void update_face_map(p4est_t* p4est, 
        PatchSurfFaceMap* face_map,
        PatchSurfFaceMap*& intermediate_face_map,
        PatchSamples*& intermediate_patch_samples);

// refinement
void resample_qbkix_points(p4est_t* p4est,
        PatchSamples* patch_samples,
        Vec& qbkix_points, 
        vector<int> qbkix_indices);
// refinement
void dump_vtk_data_for_paraview(DblNumMat qbkix_points,
        NumVec<OnSurfacePoint> closest_on_surface_points, int it,
        vector<int> global_patch_ids, PatchSurfFaceMap* face_map,
        string file_prefix=string() );

void check_patches(p4est_t* p4est);
END_EBI_NAMESPACE
#endif
