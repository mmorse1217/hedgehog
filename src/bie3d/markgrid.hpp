#ifndef __MARKGRID_HPP__
#define __MARKGRID_HPP__

#include "common/numtns.hpp"
#include "solver_interface.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"
#include <omp.h>
#include "common/kernel3d.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "bdry3d/face_map_subpatch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "common/numtns.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "common/ebiobject.hpp"
#include <map>
#include "common/stats.hpp"
#include "spatial_grid.hpp"

BEGIN_EBI_NAMESPACE
class AABBTree;
namespace Markgrid{
    //---------------------------------------------------------------------------
    // P := some patch defining the surface
    // Note about vocabulary:
    // Sample: For each patch compute a fixed sampling, equispaced
    //      (s_ij=P(ih,jh), chebyshev (s_ij=(cos(i*pi/N), cos(j*pi/N)), etc. One
    //      s_ij  from the set of samplings of all patches is a "sample."
    // On-surface Point: A point on a patch is any such p = P(x,y) for (x,y) 
    //      in the patch's domain.
    //
    // It's trivial to compute the closest on-surface sample point from some
    // fixed sampling, but it's less trivial to compute the closest on-surface
    // point. closest_on_surface_point_on_extended_patch() and
    // closest_point_along_patch_edges
    //
    // Also input_points are the set of points to mark as Near/Far, and 
    // target_point refers to a single point to be marked. These could be QBKIX
    // or general points in R^3
    //---------------------------------------------------------------------------
    enum DescentType {
        PLANE = 0,
        X_DIRECTION = 1,
        Y_DIRECTION = 2 
    };
    // An hash map of id's to a set of OnSurfacePoint's that are in the near
    // zone of the item corresponding to a given id. This is used to produce a
    // map from patches to near points (OnSurfacePoint contains the index of the
    // target point that generated it) and from points to near patches
    // (OnSurfacePoint knows which patch it lives in)
    typedef map<int, vector<OnSurfacePoint> > NearFieldMap;
    
    Point2 newton_direction(Point3 X_u, Point3 X_v, 
            Point3 X_uu, Point3 X_uv, Point3 X_vv, 
            Point3 n, Point3 p,  DescentType descent_type);
    
    Point2 gradient_direction(Point3 X_u, Point3 X_v, Point3 n,
        Point3 p,  DescentType descent_type);
    bool stopping_criteria(Point3 u, Point3 v, Point3 p, DescentType descent_type);

    // Given a surface sampling patch_samples, return a vector  of size
    // input_points.size() indicating which sample point is closest. jth slot of
    // the returned vector
    // contains the index i=0,... (# sample points) of the sample that minimizes
    // \|s_i - p_j\|_2
    vector<int> closest_sample_points(PatchSamples* patch_samples, DblNumMat input_points);


    NearFieldMap find_points_near_patch(
            DblNumMat input_points,
            FaceMapSubPatch* patch);
    
    NearFieldMap find_patches_near_point(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            SpatialGrid* grid=NULL);
    
    NearFieldMap find_closest_patch_to_point_via_bfs(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            SpatialGrid* grid);

    NearFieldMap find_patches_closest_to_point(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            AABBTree* aabb_tree=NULL);

    NearFieldMap find_closest_patch_to_point_aabb_tree(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            AABBTree* aabb_tree);

    NearFieldMap compute_closest_points_on_patches(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            vector<uint> near_patches);

    NearFieldMap compute_closest_points_on_patches(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            vector<Point2> initial_uv_guesses,
            vector<uint> near_patches);

    void populate_closest_point_data(
            int target_index,
            Point3 target_point,
            FaceMapSubPatch* patch,
            OnSurfacePoint& closest_point);
    /**
     * Find the OnSurfacePoint in on_surface_points with the lowest value of
     * distance_from_target. Equal to 
     * min_{x \in on_surface_points}(x.distance_from_target)
     * @param   vector<OnSurfacePoint>  on_surface_points   a list of filled
     *                                                      OnSurfacePoint's
     * @ret     OnSurfacePoint                              the closest
     *                                                      OnSurfacePoint in 
     *                                                      on_surface_points
     */
    OnSurfacePoint find_closest_on_surface_point_in_list(
            vector<OnSurfacePoint> on_surface_points);

    /**
     * Find the OnSurfacePoint's in on_surface_points with distance_from_target
     * within eps of closest_point.distance_from_target. This should be either
     * 1,2, or 4 points if closest_point is closest to a QBKIX point.
     *
     * This is needed for the case of edge/corner marking of QBKIX points, where
     * technically there will 2/4 closest points on 2/4 different patches.
     * As a result, just find the closest one will always indicate that the
     * OnSurfacePoint with lowest distance_from_target and patch id is closest.
     * This is, in general, incorrect.
     * @param   vector<OnSurfacePoint>  on_surface_points   a list of filled
     *                                                      OnSurfacePoint's
     * @param   OnSurfacePoint          closest_point   the closest point to
     *                                                  some target that we
     *                                                  would like all the other
     *                                                  "nearby" points to
     * @ret     vector<OnSurfacePoint>                  the closest
     *                                                  OnSurfacePoint in 
     *                                                  on_surface_points to
     *                                                  closest_point
     */

    vector<OnSurfacePoint> collect_nearby_on_surface_points(
            vector<OnSurfacePoint> on_surface_points, 
            OnSurfacePoint closest_point, 
            double eps);
    

    NearFieldMap find_patches_closest_to_point(
            Point3 target_point,
            int target_index,
            PatchSurfFaceMap* face_map,
            SpatialGrid* grid=NULL);

    OnSurfacePoint closest_point_on_patch_to_target(Point3 target_point, 
        FaceMapSubPatch* patch);

    OnSurfacePoint closest_point_on_patch_to_target(Point3 target_point, 
            Point2 initial_uv_guess,
            FaceMapSubPatch* patch);

    // Compute the closest on surface points to input_points 
    // // @param DblNumMat     input_points        target points in R^3 that we
    //                                          need to know the closest points 
    //                                          on-surface
    // @param PatchSurfFaceMap* surface         Surface as a collection of patches
    // @ret   NumVec<OnSurfacePoint>            List of size input_points.n() of
    //                                          the closest OnSurfacePoint.
    NumVec<OnSurfacePoint> compute_closest_on_surface_points(
            DblNumMat input_points,
            PatchSurfFaceMap* surface);
    
    // Compute the closest on surface points to input_points GIVEN SOME FIXED
    // SAMPLING OF THE SURFACE
    // @param DblNumMat     input_points        target points in R^3 that we
    //                                          need to know the closest points 
    //                                          on-surface
    // @param PatchSamples* patch_samples       Surface sampling obj. per patch
    // @ret   NumVec<OnSurfacePoint>            List of size input_points.n() of
    //                                          the closest OnSurfacePoint.
    NumVec<OnSurfacePoint> compute_closest_on_surface_points(
            DblNumMat input_points,
            PatchSamples* patch_samples);

    // Used in new near marking code. Does NOT fail gracefully, instead finds
    // closest point on the entire polynoimal on [-\infty, \infty]^2 to
    // target_point
    OnSurfacePoint closest_on_surface_point_on_extended_patch(
        OnSurfacePoint closest_on_surface_sample,
        Point3 target_point,
        FaceMapSubPatch* patch, 
        DescentType descent_type);
    

    // Finds the closest point to target_point on the entire polynomial of the patch 
    // among the four 1d polynomials defined by the edges x=0, x=1, y=0, y=1
    OnSurfacePoint closest_point_along_patch_edges(
        Point3 target_point,
        FaceMapSubPatch* patch,
        Point2 initial_guess=Point2(.5,.5));
    
    // Find the closest corner of a polynomial patch to target_point
    // @param Point3        target_point        point in R^3  
    // @param FaceMapSubPatch* patch               polynomial surface patch of
    //                                          interest 
    // @return  OnSurfacePoint                  the corner of patch that is
    //                                          closest to target_point
    OnSurfacePoint closest_point_patch_corners(
        Point3 target_point,
        FaceMapSubPatch* patch);

    // Compute the step size of the optimization procedure used in
    // closest_point_along_patch_edges() and
    // closest_on_surface_point_on_extended_patch() to find closest on-surface 
    // points.
    double  backtracking_line_search(FaceMapSubPatch* patch, 
            Point2 xyc, 
            Point2 xyd, 
            Point3 target_point);
    double  backtracking_line_search(Patch* patch, 
            Point2 xyc, 
            Point2 xyd, 
            Point3 target_point);


    NumVec<OnSurfacePoint> mark_target_points(DblNumMat input_points, 
        PatchSurfFaceMap* face_map, bool compute_far_marking=true);

    void mark_target_points(
        DblNumMat input_points, 
        PatchSurfFaceMap* face_map, 
        NearFieldMap& on_surface_point_map);

    
    void mark_far_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
            NumVec<OnSurfacePoint>& on_surface_point);
    
    void mark_near_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
            NumVec<OnSurfacePoint>& on_surface_point);
    
    void mark_near_field(DblNumMat input_points, PatchSurfFaceMap* face_map, 
            NearFieldMap& on_surface_point_map);
    void mark_near_field_parallel(DblNumMat input_points, PatchSurfFaceMap* face_map, 
            NumVec<OnSurfacePoint>& on_surface_points);
    void check_all_points_are_marked(DblNumMat input_points, 
            NumVec<OnSurfacePoint>& on_surface_point);

    void check_all_points_are_marked(DblNumMat input_points, 
            NearFieldMap& on_surface_point);

    Vec evaluate_fmm_with_constant_density(DblNumMat input_points, PatchSurfFaceMap* face_map);

    OnSurfacePoint select_closest_point_biased_toward_target_patch(
            vector<OnSurfacePoint> closest_on_surface_points, 
            FaceMapSubPatch* target_patch, Point3 target_point,
            PatchSurfFaceMap* face_map);

}
END_EBI_NAMESPACE
#endif
