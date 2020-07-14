#ifndef _PATCH_SURF_FACE_MAP_HPP_
#define _PATCH_SURF_FACE_MAP_HPP_

#include "patch_surf.hpp"
#include "vec3t.hpp"
#include "on_surface_point.hpp"
extern "C" {
#include <p4est.h>
#include <p4est_connectivity.h>
}
BEGIN_EBI_NAMESPACE

//----------------------------------------------------------
class PatchSurfNanospline: public PatchSurf
{
    protected:
        string _filename;

        class NanosplineInterface;
        unique_ptr<NanosplineInterface> _surface;
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
        Point3 ctr() { return Point3(0.); };
        PatchSurfNanospline(const string& n, const string& p);
        ~PatchSurfNanospline();
    public:
        void partial_teardown();

        string& filename() { return _filename; }
        int setFromOptions(){};
        int setup();
        void bounding_box(Point3& bbmin, Point3& bbmax){return;}; 


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
        void setup_from_existing_face_map(PatchSurfNanospline* face_map);
        void initialize_with_existing_p4est(p4est_t*& p4est);

        void refine(Vec near_points=NULL);
        void refine_test(Vec near_points=NULL);
        void refine_uniform(int level);
        void resolve_rhs(FunctionHandle function, int range_dim, double abs_err);
        void save(string filename);

        static PatchSurfNanospline* load(string filename);
        /*PatchChildrenMap find_subpatches(
                PatchSurfNanospline* refined_face_map);
*/
        int num_patches(){
            return _patches.size();
        }
        // defined in face_map_subpatch.cpp
        //NanosplineSubPatch* subpatch(int pi);


};

//----------------------------------------------------------
class NanosplinePatch: public Patch
{
    protected:

        // _V  - the global index of the vertex that corresponds to the center of the
        // patch. Global index is with respect to the entire surface mesh
        int _V;

        // _characteristic_length = square root of linear approximation to patch area
        //double _characteristic_length;
        // _bnd - the half-width of the bounding box that encloses the entire patch.
        double _bnd;
    public:


        // All figures referred to are in A Simple Manifold-Based Construction of
        // Surfaces of Arbitrary Smoothness - L. Ying & D. Zorin

        /*
         * Initialize a blended patch object. Constructor really only updates member
         * variables and determined the value of _bnd, which is determined by 
         * internally by blendsurf.
         */
        NanosplinePatch(PatchSurf* b, int pi, int V); 
        ~NanosplinePatch(){;}
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
        int orientation();

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
         * Determines if a point (x,y) \in D (Figure 1) is in the region of the regular 
         * star  shaped wedge in domain S of Figure 1, where (0,0) corresponds to the
         * vertex center, that corresponds to a valid evaluation region [0,EVAL_UB]x[0,EVAL_UB],
         * where EVAL_UB is set in the options file. Might be a reminant from when things 
         * weren't quite working; appears to be rather redundant in blendsurf v3 where
         * EVAL_UB() == 1.
         *
         * a point (x,y) is considered valid (in this sense) with respect to the
         * current patch if it is contained within the evaluation domain for any of
         * the k faces that are valent to vertex _V. If this isn't the case, (x,y)
         * will evaluate to zero using the current patch and will contribute
         * negligibly. 
         *        ------------* <--- (1,1)
         *      /            /
         *     /            / <--- current face containing the point (x,y)
         *    /            /       region where is_valid == true is [0,EVAL_UB()]x[EVAL_UB]
         *   /            /
         *  *------------  
         *  ^--V = (0,0)
         *
         *  @param double *   xy      position in (x,y) coordinates w.r.t to V
         *  @param bool&      is_valid     whether (x,y) is contained in the current face
         */
        int is_xy_valid( double* xy, bool& is_valid); // returns whether point is within chart upper boundary

        /*
         * Determines if a point (x,y) \in D (Figure 1) is in a box [0,.5]x[0,.5] in the regular star 
         * shaped wedge in domain S of Figure 1, where (0,0) corresponds to the vertex
         * center.
         *         ------------* <--- (1,1)
         *       /            /
         *      /            / <--- Face containing the point (x,y)
         *     *------*     /
         *    /      /     /
         *   /      /     /
         *  *------*-----  
         *  |  ^---- region where (x,y) is dominant
         *  ^--V = (0,0)
         *
         *  @param double *   xy      position in (x,y) coordinates w.r.t to V
         *  @param bool&      dominant     whether (x,y) is dominant, i.e. contained in the 
         *                            box described above 
         */
        int is_xy_dominant( double* xy, bool& dominant);

        int xy_to_face_point(double* xy, FacePointOverlapping* face_point);

        /*
         * Interface to Blendsurf::BdSurf::eval(...). Evaluates the surface position, 
         * 1st and 2nd derivatives (depending on flags) from xy-coordinates in the 
         * coordinate system of global vertex _V, using contributions from all (4)
         * charts that share the face containing (x,y) ((x,y) \in D in Figure 1.)
         * @param double *    xy      position in (x,y) coordinates w.r.t V
         * @param int         flag    BdSurf::EVAL_VALUE| BdSurf::EVAL_1ST_DERIV | 
         *                            BdSurf::EVAL_2ND_DERIV
         *                            indicate to Blendsurf whether to compute only
         *                            (u,v) coordinates on the surface, coordinates 
         *                            + 1st derivatives, or coordinates + 1st and 
         *                            2nd derivatives.
         * @param double *            surface positions + 1st/2nd derivatives
         *                            (if applicable).
         * @param int         LL      number of legendre points (TODO: remove param)
         */
        int xy_to_patch_coords(double* xy, int flag, double*);
        void eval_unsafe(double* xy, int flag, double*);

        /*
         * Computes the value of the patch centered at vertex _V at the point
         * (x,y) ((x,y) \in D in Figure 1). First, it transforms an xy point in D to
         * cd coordinates in S (Figure 1), the regular star-shaped domain. The value of 
         * the patch centered at _V  is then evaluated.
         * @param double*     xy      position in (x,y) coordinates in D
         * @param int         flag    always EVAL_VALUE (TODO: remove param)        
         * @param double *            value of the patch at position c^{-1}_{_V}(x,y)
         *                            in domain S_{_V} (c_i as shown in Figure 1).
         */
        int xy_to_patch_value(double* xy, int flag, double*);

        /*
         * Estimate the jacobian of the surface at the point (0,0) (w.r.t to the
         * vertex-centered coordinates) using the current patch.
         *
         * @param double*             The resulting jacobian of the surface
         */
        int estimate_jacobian(double*);

        void bounding_box(Point3& min, Point3& max);

        Point3 normal(double* xy) ;
        //double gaussian_curvature(Point2 xy);
        //double mean_curvature(Point2 xy);
        //void principal_curvatures(Point2 xy, double& k1, double& k2);

        //static double _UB;
};
//----------------------------------------------------------
class FacePointNanospline: public FacePointOverlapping
/* 
* Unique representation of a sample point in the chart-independent coordinate
* system of the global face F. Note that F lives in the S-domian of Figure 1 of 
* Blendsurf paper (v1), i.e. in a single regular wedge with arbitrary origin.
*
*/
{
    protected:
        int _F;
        double _cd[2];
    public:
        FacePointNanospline(int F, double* cd): FacePointOverlapping(), _F(F) {
            _cd[0]=cd[0];	 _cd[1]=cd[1];
        }
        int F() { return _F; }
        double* cd() { return _cd; }

        FacePointNanospline static to_face_point(OnSurfacePoint p){
            return FacePointNanospline(p.parent_patch, p.parametric_coordinates.array());
        }
};

END_EBI_NAMESPACE

#endif
