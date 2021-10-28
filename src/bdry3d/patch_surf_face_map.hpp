#ifndef _PATCH_SURF_FACE_MAP_HPP_
#define _PATCH_SURF_FACE_MAP_HPP_

#include "patch_surf.hpp"
#include <face_map.hpp>
#include "vec3t.hpp"
#include "on_surface_point.hpp"
#include <unordered_map>
#include <diff_geom.hpp>
extern "C" {
#include <p4est.h>
#include <p4est_connectivity.h>
}
BEGIN_EBI_NAMESPACE

typedef unordered_map<int, vector<int> > PatchChildrenMap;
class FaceMapSubPatch;
//----------------------------------------------------------
class PatchSurfFaceMap: public PatchSurf
{
    protected:
        //PARAMS
        string _filename;

        //COMPONENTS
        FaceMapSurf* _face_map;
    public:
        double _on_surface_threshold;
        p4est_t* _p4est;
        p4est_connectivity_t* _p4est_connectivity;
        bool _coarse;

        enum SurfaceType {
            NOT_SET = -1,
            BLENDED = 0,
            POLYNOMIAL = 1,
            ANALYTIC = 2,
        };

        SurfaceType _surface_type;
        PatchSurfFaceMap(const string& n, const string& p);
        ~PatchSurfFaceMap();
        void partial_teardown();

        string& filename() { return _filename; }
        FaceMapSurf& face_map() { return *_face_map; }
        int setFromOptions();
        int setup();
        Point3 ctr() { return _face_map->ctr(); }
        void bounding_box(Point3& bbmin, Point3& bbmax) { 
            _face_map->bbox(bbmin, bbmax); 
        }



        DblNumVec* _quadrature_weights;

        /*
         * For a particular FacePointOverlapping, patches_containing_face_point
         * returns the list of patches that contain it, where "contain" means that
         * the point is member of the patch's domain S. Note that for any face point
         * on the surface, there will always be exactly four patches that have a
         * non-trivial contribution at that point, so pivec.size() will always be 4.
         *
         * @param FacePointOverlapping *              face point in domain S of
         *                                            Figure 1, not necessarily in
         *                                            any particular domain S_k, per
         *                                            se.
         * @param vector<int>&            pivec       List of overlapping patches
         *                                            that contain the procided face
         *                                            point.
         */
        int patches_containing_face_point(FacePointOverlapping*, vector<int>& pivec); 
        int face_point_size_in_doubles();
        void setup_from_existing_face_map(PatchSurfFaceMap* face_map);
        void initialize_with_existing_p4est(p4est_t*& p4est);

        void refine(Vec near_points=NULL);
        void refine_test(Vec near_points=NULL);
        void refine_uniform(int level);
        void resolve_rhs(FunctionHandle function, int range_dim, double abs_err);
        void save(string filename);

        static PatchSurfFaceMap* load(string filename);
        PatchChildrenMap find_subpatches(
                PatchSurfFaceMap* refined_face_map);

        int num_patches(){
            return _patches.size();
        }
        // defined in face_map_subpatch.cpp
        FaceMapSubPatch* subpatch(int pi);


};

//----------------------------------------------------------
class FaceMapPatch: public Patch
{
    protected:

        // _V  - the global index of the face that corresponds to the patch. 
        // Global index is with respect to the entire surface mesh
        int _V;

        // _characteristic_length = square root of linear approximation to patch area
        //double _characteristic_length;
        // _bnd - unused
        double _bnd;
    public:


        // All figures referred to are in A Simple Manifold-Based Construction of
        // Surfaces of Arbitrary Smoothness - L. Ying & D. Zorin

        FaceMapPatch(PatchSurf* b, int pi, int V); 
        ~FaceMapPatch(){;}
        int& V() { 
            return _V;
        }

        double bnd() { 
            return _bnd; 
        }

        //double characteristic_length(){
            //return _characteristic_length;
        //}

        NumVec<OnSurfacePoint> sample_patch(int num_samples, SamplingPattern sampling_pattern);

        int group_id();
        int orientation(){
            FaceMapSurf& face_map = ((PatchSurfFaceMap*)bdry())->face_map();
            return face_map.get_orientation(_V);
        }

        /*
         * Identical functionality to is_xy_valid, except that is takes a
         * FacePointOverlapping as input instead of (x,y) coordinates. See
         * FacePointBlended and is_xy_valid below for more details
         *
         * @param FacePointOverlapping*   face_point    Sample point in
         *                                              coord-independent face 
         *                                              representation; lives in 
         *                                              domain S of Figure 1
         * @param bool&                   is_valid      whether face_point is
         *                                              contained in a face adjacent
         *                                              to _V
         */
        int is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid);


        int face_point_to_xy(FacePointOverlapping* face_point, double* xy);


        /*
         * Determines if a point (x,y) is in the box [0,1]^2
         *  @param double *   xy      position in (x,y) coordinates w.r.t to V
         *  @param bool&      is_valid     whether (x,y) is contained in the current face
         */
        int is_xy_valid( double* xy, bool& is_valid); // returns whether point is within chart upper boundary

        /*
         *  Unused legacy function from blendsurf
         *  @param double *   xy      position in (x,y) coordinates w.r.t to V
         *  @param bool&      dominant     whether (x,y) is dominant, i.e. contained in the 
         *                            box described above 
         */
        int is_xy_dominant( double* xy, bool& dominant);

        int xy_to_face_point(double* xy, FacePointOverlapping* face_point);

        /*
         * Terrible named function that evaluates a 2d coordinate xy and returns
         * its value on the patch
         * @param double *    xy      position in (x,y) coordinates w.r.t V
         * @param int         flag    BdSurf::EVAL_VALUE| BdSurf::EVAL_1ST_DERIV | 
         *                            BdSurf::EVAL_2ND_DERIV
         *                            indicate to Blendsurf whether to compute only
         *                            (u,v) coordinates on the surface, coordinates 
         *                            + 1st derivatives, or coordinates + 1st and 
         *                            2nd derivatives.
         * @param double *            surface positions + 1st/2nd derivatives
         *                            (if applicable).
         */
        int xy_to_patch_coords(double* xy, int flag, double*);
        void eval_unsafe(double* xy, int flag, double*);

        /*
         *  Unused legacy function from blendsurf, should return 1
         * @param double*     xy      position in (x,y) coordinates in D
         * @param int         flag    always EVAL_VALUE (TODO: remove param)        
         * @param double *            value of the patch at position c^{-1}_{_V}(x,y)
         *                            in domain S_{_V} (c_i as shown in Figure 1).
         */
        int xy_to_patch_value(double* xy, int flag, double*);

        /*
         * Estimate the jacobian of the patch at the point (.5, .5) (w.r.t to the
         * vertex-centered coordinates) using the current patch.
         *
         * @param double*             The resulting jacobian of the surface
         */
        int estimate_jacobian(double*);

        void bounding_box(Point3& min, Point3& max);

        Point3 normal(double* xy) ;
        double gaussian_curvature(Point2 xy){
            FaceMapSurf& face_map = dynamic_cast<PatchSurfFaceMap*>(bdry())->face_map();
            return Differential::gaussian_curvature(xy, &face_map, _V);
        }
        double mean_curvature(Point2 xy){
            FaceMapSurf& face_map = dynamic_cast<PatchSurfFaceMap*>(bdry())->face_map();
            return Differential::mean_curvature(xy, &face_map, _V);
        }
        void principal_curvatures(Point2 xy, double& k1, double& k2);

        //static double _UB;
};
//----------------------------------------------------------
class FacePointFaceMap: public FacePointOverlapping
/* 
* Unique representation of a sample point in the chart-independent coordinate
* system of the global face F.
*/
{
    protected:
        int _F;
        double _cd[2];
    public:
        FacePointFaceMap(int F, double* cd): FacePointOverlapping(), _F(F) {
            _cd[0]=cd[0];	 _cd[1]=cd[1];
        }
        int F() { return _F; }
        double* cd() { return _cd; }
        
        FacePointFaceMap static to_face_point(OnSurfacePoint p){
            return FacePointFaceMap(p.parent_patch, p.parametric_coordinates.array());
        }
};

END_EBI_NAMESPACE

#endif
