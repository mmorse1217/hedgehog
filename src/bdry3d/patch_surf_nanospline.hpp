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
        int num_patches(){
            return _patches.size();
        }
        


};

//----------------------------------------------------------
class NanosplinePatch: public Patch
{
    protected:
        class NanosplineInterface;
        unique_ptr<NanosplineInterface> _surface;

        // _V  - the global index of the vertex that corresponds to the center of the
        // patch. Global index is with respect to the entire surface mesh
        int _V;

        // _characteristic_length = square root of linear approximation to patch area
        //double _characteristic_length;
        // _bnd - the half-width of the bounding box that encloses the entire patch.
        double _bnd;
    public:



        NanosplinePatch(PatchSurf* b, int pi, int V); 
        ~NanosplinePatch();
        int& V() { 
            return _V;
        }

        double bnd() { 
            return _bnd; 
        }

        NumVec<OnSurfacePoint> sample_patch(int num_samples, SamplingPattern sampling_pattern);

        int group_id();
        int orientation();

        int is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid);


        int face_point_to_xy(FacePointOverlapping* face_point, double* xy);


        int is_xy_valid( double* xy, bool& is_valid); // returns whether point is within chart upper boundary

        int is_xy_dominant( double* xy, bool& dominant);

        int xy_to_face_point(double* xy, FacePointOverlapping* face_point);

        int xy_to_patch_coords(double* xy, int flag, double*);
        void eval_unsafe(double* xy, int flag, double*);

        int xy_to_patch_value(double* xy, int flag, double*);

        /*
         * Estimate the jacobian of the patch
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
