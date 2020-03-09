#ifndef _FACE_MAP_SUBPATCH_HPP_
#define _FACE_MAP_SUBPATCH_HPP_

#include "patch_surf.hpp"
#include "patch_surf_face_map.hpp"
#include "patch_samples.hpp"
#include "common/interpolate.hpp"
#include <face_map.hpp>
#include "vec3t.hpp"

BEGIN_EBI_NAMESPACE
typedef pair<double, double> Interval;
typedef pair<Interval, Interval> Rectangle;


class FaceMapSubPatch: public Patch{
    private:

        //double _characteristic_length;
    public:
        int _id;
        FaceMapPatch* _face_map_patch;
        DblNumVec* _quadrature_weights; 
        Interval _x_interval; 
        Interval _y_interval;
        int _level;
        int _parent_id;
        int _coarse_parent_patch;
        int _quad_id_within_p4est_tree;
        p4est_quadrant_t* _quadrant;
        double _near_zone_distance;

        FaceMapSubPatch(FaceMapPatch* face_map_patch, 
                Interval x_interval, 
                Interval y_interval,
                int level,
                int id,
                int parent_id,
                DblNumVec* quad_weights=NULL);


        // Pass virtual function calls to the FaceMapPatch containing the
        // FaceMapSubPatch
        int V(){ return _id;}
        double bnd();
        int group_id();
        int is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid);
        int face_point_to_xy(FacePointOverlapping* face_point, double* xy);
        void rescale(double* xy, double* xy_scaled);
        void unscale(double* xy, double* xy_scaled);
        int is_xy_valid( double* xy, bool& is_valid);
        int is_xy_dominant( double* xy, bool& dominant);
        int xy_to_face_point(double* xy, FacePointOverlapping* face_point);
        int xy_to_patch_coords(double* xy, int flag, double* ret);
        void eval_unsafe(double* xy, int flag, double* ret);
        int xy_to_patch_value(double* xy, int flag, double* ret);
        int estimate_jacobian(double* ret);
        void bounding_box(Point3& bounding_box_min, Point3& bounding_box_max);
        void inflated_bounding_box(Point3& bounding_box_min, Point3& bounding_box_max);
        DblNumMat control_points();
bool is_patch_valid(Vec function_values_at_children, 
        Vec function_values_at_parent, 
        Vec uv_coordinates_single_patch, int range_dim, double eps);

        static FaceMapSubPatch* as_subpatch(Patch* p){
            return dynamic_cast<FaceMapSubPatch*>(p);
        }
        void mesh_bounding_box(DblNumMat& vertices, IntNumMat& triangles);
        void single_bounding_box_triangle(DblNumMat& vertices, IntNumMat& triangles);
        double gaussian_curvature(Point2 xy);
        double mean_curvature(Point2 xy);
        void principal_curvatures(Point2 xy, double& k1, double& k2);
        Point3 normal(double* xy);
        void compute_near_zone_distance();
};

END_EBI_NAMESPACE

#endif
