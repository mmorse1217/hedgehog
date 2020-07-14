#include <nanospline/forward_declaration.h>
#include <nanospline/BSplinePatch.h>
#include <Eigen/Core>
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
using Scalar = double;
BEGIN_EBI_NAMESPACE

Eigen::MatrixXd initialize_flat_patch_control_grid(const int num_control_points_1D, const int degree_u, 
        const int degree_v){
    const int num_control_points = num_control_points_1D*num_control_points_1D;
    Eigen::MatrixXd control_grid(num_control_points, DIM);

    // identify x,y in R^3 with (u,v) in spline parameter space
    // Place control points  equispaced on [-1,1] x [-1,1] x {0}
    for (int ci = 0; ci < num_control_points_1D; ci++) { // u
      for (int cj = 0; cj < num_control_points_1D; cj++) { // v
          const int index = ci*num_control_points_1D + cj;
          
          // x and y spacing are the same
          const double spacing = 2./(num_control_points_1D-1);
          const double x = ci*spacing - 1.;
          const double y = cj*spacing - 1.;
          const double z = 0.; 
          control_grid.row(index) << x, y, z;
      }
    }

    return control_grid;
    
}

Eigen::MatrixXd initialize_flat_patch_knots(const int num_control_points_1D,
                                            const int degree) {
  const int num_knots = num_control_points_1D + degree + 1;
  Eigen::MatrixXd knots(num_knots, 1);
  int counter = 0;
  for (int i = 0; i < degree + 1; i++) {
    knots(counter) = 0.;
    counter++;
  }

  // uniform knots for now.
  const double knot_spacing = 1. / (num_control_points_1D - degree - 1 + 1);
  for (int i = 0; i < num_control_points_1D - degree - 1; i++) {
      cout << "knot_spacing: " << knot_spacing << endl;
      cout << "i+1: " << i+1<< endl;
      cout << "knot: " << (i+1) * knot_spacing<< endl;
    knots(counter) = (i+1) * knot_spacing;
    counter++;
  }
  for (int i = 0; i < degree + 1; i++) {
    knots(counter) = 1.;
    counter++;
  }
  assert(counter == num_knots);
  return knots;
}
class PatchSurfNanospline::NanosplineInterface {
    public:
  NanosplineInterface() {
    // 1. Get the name of surface we want to load
    // For now: load a flat surface, all weights = 1, uniform knots
        const int num_control_points_1D = 5;
        const int degree = 3;
        const int degree_u = degree;
        const int degree_v = degree;
        const int num_control_points = num_control_points_1D*num_control_points_1D;
        nanospline::BSplinePatch<Scalar, DIM, -1, -1> patch;
        
    // 2. Load surface data with obj reader
    // For now: explicitly generate control points, knots and weights.
        Eigen::MatrixXd control_grid = 
            initialize_flat_patch_control_grid(num_control_points_1D, degree_u, degree_v);
        cout << "control_grid: " << endl;
        cout << control_grid << endl;
        patch.set_control_grid(control_grid);

        Eigen::MatrixXd knots_u = initialize_flat_patch_knots(num_control_points_1D, degree_u);
        Eigen::MatrixXd knots_v = initialize_flat_patch_knots(num_control_points_1D, degree_v);
        cout << "knots_u: " << endl;
        cout << knots_u << endl;
        cout << "knots_v: " << endl;
        cout << knots_v << endl;
        patch.set_knots_u(knots_u);
        patch.set_knots_v(knots_v);
        
        patch.set_degree_u(degree_u);
        patch.set_degree_v(degree_v);
        patch.initialize();
  }
  ~NanosplineInterface() {
  }
};

PatchSurfNanospline::~PatchSurfNanospline(){};
PatchSurfNanospline::PatchSurfNanospline(const string &n, const string &p)
    : PatchSurf(n, p), _surface_type(NOT_SET), _coarse(false) {
    // 1. Initialize NanosplineInterface
    // 2. Perform any precomputations for patch-related data
    // surface area?
  unique_ptr<NanosplineInterface> impl(new NanosplineInterface);
  _surface = std::move(impl);
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

/*void NanosplinePatch::principal_curvatures(Point2 xy, double& k1, double& k2){
    double H = mean_curvature(xy);
    double K = gaussian_curvature(xy);
    double discriminant = H*H-K;
    if(fabs(discriminant) <= 1e-5){
        discriminant = fabs(discriminant);
    }
    k1 = H + sqrt(discriminant);
    k2 = H - sqrt(discriminant);
}*/

END_EBI_NAMESPACE
