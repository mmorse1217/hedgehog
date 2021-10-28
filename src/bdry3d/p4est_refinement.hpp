#ifndef __P4EST_REFINEMENT_HPP__
#define __P4EST_REFINEMENT_HPP__

#include "bdry3d/patch_surf_face_map.hpp"
#include <gpmesh.hpp>
#include "bdry3d/face_map_subpatch.hpp"
extern "C" {
#include <p4est_connectivity.h>
#include <p4est.h>
#include <p4est_extended.h>
#include <p4est_bits.h>
#include <p4est_vtk.h>
}
BEGIN_EBI_NAMESPACE
void refine_patches_uniform(int max_level, p4est_t* p4est, PatchSurfFaceMap*& face_map);

void refine_patches_for_qbkix_point_location(p4est_t* p4est, PatchSurfFaceMap*& face_map);

void refine_patches_point_location(p4est_t* p4est, PatchSurfFaceMap*& face_map, 
        vector<int> qbkix_indices);

void refine_patches_midpoint_near_medial_axis(p4est_t*& p4est, PatchSurfFaceMap*& face_map,
        vector<int> qbkix_indices);


void refine_patches_for_fixed_qbkix_points(p4est_t*& p4est, PatchSurfFaceMap*& face_map);
void resolve_function(p4est_t* p4est, PatchSurfFaceMap*& face_map, FunctionHandle f, int range_dim, double eps_abs /*, double eps_rel*/);

void balance(p4est_t* p4est);
END_EBI_NAMESPACE
#endif
