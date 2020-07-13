#include <nanospline/NURBSPatch.h>
#include "patch_surf_nanospline.hpp"
#include "common/interpolate.hpp"
#include "common/utils.hpp"
#include <sampling.hpp>

// TODO abstract library calls out
//#include "p4est_interface.hpp"
#include "p4est_refinement.hpp"
#include "bie3d/solver_utils.hpp"
extern "C" {
#include <p4est.h>
#include <p4est_vtk.h>
}

using std::istringstream;
using Sampling::sample_1d;
using Sampling::sample_2d;
using Sampling::base_domain;
using Sampling::chebyshev2;

BEGIN_EBI_NAMESPACE

class PatchSurfNanospline::NanosplineInterface {
    public:
  NanosplineInterface(string filename) {
  }
  ~NanosplineInterface() {
  }
};

PatchSurfNanospline::PatchSurfNanospline(const string &n, const string &p)
    : PatchSurf(n, p), _surface_type(NOT_SET), _coarse(false) {
    // 1. Get the name of surface we want to load
    // For now: load a flat surface, all weights = 1, uniform knots
    // 2. Load surface data with obj reader
    // For now: explicitly generate control points, knots and weights.
    // 3. Initialize NanosplineInterface
    // 4. Perform any precomputations for patch-related data
    // surface area?
}
// ---------------------------------------------------------------------- 



// ---------------------------------------------------------------------- 
int PatchSurfNanospline::setup()
{
    return 0;
}

// ---------------------------------------------------------------------- 
int PatchSurfNanospline::face_point_size_in_doubles() {
    return (sizeof(FacePointNanospline)-1)/sizeof(double)+1;
}

// ---------------------------------------------------------------------- 
int PatchSurfNanospline::patches_containing_face_point(FacePointOverlapping* face_point, vector<int>& pivec)
{
    return 0;
}


PatchSurfNanospline* PatchSurfNanospline::load(string filename){

    return NULL;

}

NanosplinePatch::NanosplinePatch(PatchSurf* b, int pi, int V): Patch(b, pi), _V(V)
{
  
}

int NanosplinePatch::group_id()
{
    return 0;
}

NumVec<OnSurfacePoint> NanosplinePatch::sample_patch(
        int num_samples,
        SamplingPattern sampling_pattern){
    NumVec<OnSurfacePoint> ret;
    return ret;

}

// ---------------------------------------------------------------------- 
int NanosplinePatch::is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid)
{
    return 0;
}
// ---------------------------------------------------------------------- 

int NanosplinePatch::face_point_to_xy(FacePointOverlapping* face_point, double* xy)
{
    return 0;
}
// ---------------------------------------------------------------------- 
//
int NanosplinePatch::is_xy_valid(double* xy, bool& is_valid)
{
    return 0;
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::is_xy_dominant(double* xy, bool& dominant)
{
    return 0;
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point)
{
    return 0;
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::xy_to_patch_coords(double* xy, int flag, double* ret)
{
    return 0;
}
Point3 NanosplinePatch::normal(double* xy) {
    vector<Point3> position_and_derivs(3,Point3(0.));
    xy_to_patch_coords(xy, EVAL_VL|EVAL_FD, (double*)position_and_derivs.data());

    Point3 normal = cross(position_and_derivs[1],position_and_derivs[2]);
    return normal.dir();

}

void NanosplinePatch::eval_unsafe(double* xy, int flag, double* ret) {
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::estimate_jacobian(double* jac)
{
    return 0;
}
// ---------------------------------------------------------------------- 

int NanosplinePatch::xy_to_patch_value(double* xy, int flag, double* ret)
{
    return 0;
}

void NanosplinePatch::bounding_box(Point3& bounding_box_min, Point3& bounding_box_max){
}

void NanosplinePatch::principal_curvatures(Point2 xy, double& k1, double& k2){
    double H = mean_curvature(xy);
    double K = gaussian_curvature(xy);
    double discriminant = H*H-K;
    if(fabs(discriminant) <= 1e-5){
        discriminant = fabs(discriminant);
    }
    k1 = H + sqrt(discriminant);
    k2 = H - sqrt(discriminant);
}

END_EBI_NAMESPACE
