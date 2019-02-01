#ifndef __VTK_WRITER_HPP__
#define __VTK_WRITER_HPP__

#include "nummat.hpp"
#include "numvec.hpp"
#include "vec2t.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
extern "C" {
#include <p4est.h>
}
BEGIN_EBI_NAMESPACE
void write_general_points_to_vtk(Vec points, int degrees_of_freedom, string filename,
        Vec values,
        string file_prefix=string());

void write_qbkix_points_to_vtk(DblNumMat qbkix_points_local, 
        NumVec<OnSurfacePoint>final_closest_on_surface_points, int interation=0,
        string file_prefix=string());

void write_face_map_patches_to_vtk(DblNumMat qbkix_points, 
        vector<int> patch_refined_relative_ids,
        PatchSurfFaceMap* face_map, int iteration,
        string file_prefix=string());
void write_face_map_mesh_to_vtk(
        PatchSurfFaceMap* face_map, int iteration,
        string file_prefix=string(),
        int num_samples_1d=20);

void write_face_map_patch_bounding_boxes_to_vtk(
        vector<int> patches_refined_relative_ids,
        PatchSurfFaceMap* face_map, int iteration,
        string file_prefix=string(),
        bool inflate = false);

void write_lines_from_qbkix_to_closest_point(DblNumMat qbkix_points, 
        NumVec<OnSurfacePoint> final_closest_on_surface_points, 
        PatchSurfFaceMap* face_map, int iteration,
        string file_prefix=string());
void write_triangle_mesh_to_vtk(
        DblNumMat vertices, IntNumMat faces, int iteration, string file_prefix,
        vector<int> corresponding_patches=vector<int>());

END_EBI_NAMESPACE

#endif
