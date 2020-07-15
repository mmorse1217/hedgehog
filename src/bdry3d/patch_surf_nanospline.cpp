#include <nanospline/PatchBase.h>
#include <nanospline/forward_declaration.h>
#include <nanospline/BSplinePatch.h>
#include <nanospline/save_obj.h>
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
Eigen::MatrixXd get_greville_abscissae(const Eigen::MatrixXd knots, const int degree){
    const int num_control_points = knots.size() - degree - 1;
    Eigen::MatrixXd greville_abscissae(num_control_points,1);
    greville_abscissae.setZero();
    for (int k=0; k < num_control_points; k++) {
        for (int d=0; d < degree; d++) {
            greville_abscissae(k,0) += knots(k + d+1,0);
        }
    }
    greville_abscissae /= degree;
    return greville_abscissae;
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
    knots(counter) = (i+1) * knot_spacing;
    counter++;
  }
  for (int i = 0; i < degree + 1; i++) {
    knots(counter) = 1.;
    counter++;
  }
  assert(counter == num_knots);
  // uniform knots for now.
  /*const double knot_spacing = 1. / (num_knots - 1);
  for (int i = 0; i < num_knots; i++) {
    knots(i) = (i) * knot_spacing;
  }*/
  return knots;
}

Eigen::MatrixXd initialize_flat_patch_control_grid(const int num_control_points_1D, const int degree_u, 
        const int degree_v){
    const int num_control_points = num_control_points_1D*num_control_points_1D;
    Eigen::MatrixXd control_grid(num_control_points, DIM);

    Eigen::MatrixXd knots_u = initialize_flat_patch_knots(num_control_points_1D, degree_u);
    Eigen::MatrixXd knots_v = initialize_flat_patch_knots(num_control_points_1D, degree_v);
    Eigen::MatrixXd greville_u = get_greville_abscissae(knots_u, degree_u);
    Eigen::MatrixXd greville_v = get_greville_abscissae(knots_v, degree_v);

    // identify x,y in R^3 with (u,v) in spline parameter space
    // Place control points  equispaced on [-1,1] x [-1,1] x {0}
    for (int ci = 0; ci < num_control_points_1D; ci++) { // u
      for (int cj = 0; cj < num_control_points_1D; cj++) { // v
          const int index = ci*num_control_points_1D + cj;
          
          // x and y spacing are the same
          const double spacing = 2./(num_control_points_1D-1);
          const double x = 2.*greville_u(ci) - 1.;
          const double y = 2.*greville_u(cj) - 1.;
          const double z = 0.; 
          control_grid.row(index) << x, y, z;
      }
    }

    return control_grid;
    
}

class NanosplinePatch::NanosplineInterface {
    public:
        nanospline::BSplinePatch<Scalar, DIM, -1, -1>  _patch;
  NanosplineInterface() {
    // 1. Get the name of surface we want to load
    // For now: load a flat surface, all weights = 1, uniform knots
        const int num_control_points_1D = 15;
        const int degree = 3;
        const int degree_u = degree;
        const int degree_v = degree;
        const int num_control_points = num_control_points_1D*num_control_points_1D;

        // 2. Load surface data with obj reader
        // For now: explicitly generate control points, knots and weights.
        Eigen::MatrixXd control_grid = 
            initialize_flat_patch_control_grid(num_control_points_1D, degree_u, degree_v);
        _patch.set_control_grid(control_grid);

        Eigen::MatrixXd knots_u = initialize_flat_patch_knots(num_control_points_1D, degree_u);
        Eigen::MatrixXd knots_v = initialize_flat_patch_knots(num_control_points_1D, degree_v);
        _patch.set_knots_u(knots_u);
        _patch.set_knots_v(knots_v);

        _patch.set_degree_u(degree_u);
        _patch.set_degree_v(degree_v);
        _patch.initialize();
        //nanospline::save_patch_obj("flat_patch.obj", _patch);
  }
  ~NanosplineInterface() {
  }
};

PatchSurfNanospline::~PatchSurfNanospline(){};
PatchSurfNanospline::PatchSurfNanospline(const string &n, const string &p)
    : PatchSurf(n, p), _surface_type(NOT_SET), _coarse(false) {
        // TODO generalize to arbtirary patch sets
    _patches.push_back(new NanosplinePatch(this, 0, 0));


    // TODO move to patch constructor
    const int num_samples = 20;
    NumMatrix samples = sample_2d<chebyshev2>(num_samples, base_domain);

    DblNumVec quad_weight(num_samples*num_samples);


    DblNumVec quadrature_nodes(num_samples, true,
            sample_1d<chebyshev2>(num_samples,base_domain).data());
    
    DblNumVec* quadrature_weights = new DblNumVec(num_samples);
    
    for(int si = 0; si < num_samples; si++){
        (*quadrature_weights)(si) = 
            Interpolate::integrate_ith_lagrange_basis_func(
                si,0., 1., num_samples, quadrature_nodes, 50, 1.);

    }
    _quadrature_weights = quadrature_weights;

    for(int si = 0; si < num_samples; si++){
        for(int sj = 0; sj < num_samples; sj++){
            int index = si*num_samples + sj;

            quad_weight(index) = 
                (*quadrature_weights)(si)*(*quadrature_weights)(sj);

        }
    }
    for (int pi=0; pi <_patches.size(); pi++) {

        Patch* patch = _patches[pi];

        DblNumMat normals(DIM,num_samples*num_samples);
        for(int i = 0; i < num_samples*num_samples; i++){
            Point2 sample_uv(samples.clmdata(i));

            vector<Point3> values(3, Point3());
            patch->xy_to_patch_coords(sample_uv.array(), 
                    PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, 
                    (double*)values.data());

            Point3 normal = cross(values[1], values[2]);
            for(int d = 0; d < DIM; d++)
                normals(d, i) = normal(d);
        }

        patch->characteristic_length(normals, quad_weight);
    }
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

NanosplinePatch::~NanosplinePatch(){ }

NanosplinePatch::NanosplinePatch(PatchSurf* b, int pi, int V): Patch(b, pi), _V(V)
{
    // 1. Initialize NanosplineInterface
    // 2. Perform any precomputations for patch-related data
    // surface area?
  unique_ptr<NanosplineInterface> impl(new NanosplineInterface);
  _surface = std::move(impl);
  // Compute surface area of the patch
}

int NanosplinePatch::group_id()
{
    return 0;
}


// ---------------------------------------------------------------------- 
int NanosplinePatch::is_face_point_valid(FacePointOverlapping* face_point, bool& is_valid)
{
  FacePointNanospline* face_point_blended = (FacePointNanospline*) face_point;
  double* cd = face_point_blended->cd();
  is_xy_valid(cd, is_valid);
    return 0;
}
// ---------------------------------------------------------------------- 

int NanosplinePatch::face_point_to_xy(FacePointOverlapping* face_point, double* xy)
{
  FacePointNanospline* face_point_blended = (FacePointNanospline*)face_point;  
  double* cd = face_point_blended->cd();
  xy[0] = cd[0];
  xy[1] = cd[1];
    return 0;
}
// ---------------------------------------------------------------------- 
//
int NanosplinePatch::is_xy_valid(double* xy, bool& is_valid)
{
    const double x = xy[0];
    const double y = xy[1];
    const double x_lower = _surface->_patch.get_u_lower_bound();
    const double x_upper =_surface->_patch.get_u_upper_bound();
    const double y_lower =_surface->_patch.get_v_lower_bound();
    const double y_upper =_surface->_patch.get_v_upper_bound();
    is_valid = x_lower <= x && x <= x_upper &&
        y_lower <= y && y <= y_upper;
    return is_valid;
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::is_xy_dominant(double* xy, bool& dominant)
{
  // All valid points are "dominant" for face-centered maps
  is_xy_valid(xy, dominant);
    return 0;
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::xy_to_face_point(double* xy, FacePointOverlapping* face_point)
{
    return 0;
}
// ---------------------------------------------------------------------- 
void copy_eigen_to_point(Eigen::MatrixXd eigen_matrix, Point3& point){
    assert(eigen_matrix.cols() == 3);
    assert(eigen_matrix.rows() == 1);
  for (int d = 0; d < DIM; d++) {
      point(d) = eigen_matrix(0,d);
  }
}
int NanosplinePatch::xy_to_patch_coords(double* xy, int flag, double* ret)
{ 
    bool is_valid;
    is_xy_valid(xy, is_valid);
    assert(is_valid);

    Point3* point_ret = (Point3*)ret;
    Eigen::MatrixXd position = _surface->_patch.evaluate(xy[0], xy[1]);
    copy_eigen_to_point(position, point_ret[0]); 

    if( flag & EVAL_FD){

        Eigen::MatrixXd derivative_u = _surface->_patch.evaluate_derivative_u(xy[0], xy[1]);
        Eigen::MatrixXd derivative_v = _surface->_patch.evaluate_derivative_v(xy[0], xy[1]);
        copy_eigen_to_point(derivative_u, point_ret[1]);
        copy_eigen_to_point(derivative_v, point_ret[2]);
    }
    if (flag & EVAL_2ND_DERIV) {
        
        Eigen::MatrixXd derivative_uu = _surface->_patch.evaluate_2nd_derivative_uu(xy[0], xy[1]);
        Eigen::MatrixXd derivative_uv = _surface->_patch.evaluate_2nd_derivative_uv(xy[0], xy[1]);
        Eigen::MatrixXd derivative_vv = _surface->_patch.evaluate_2nd_derivative_vv(xy[0], xy[1]);
        copy_eigen_to_point(derivative_uu, point_ret[3]);
        copy_eigen_to_point(derivative_uv, point_ret[4]);
        copy_eigen_to_point(derivative_vv, point_ret[5]);
    }

    return 0;
}
Point3 NanosplinePatch::normal(double* xy) {
    vector<Point3> position_and_derivs(3,Point3(0.));
    xy_to_patch_coords(xy, EVAL_VL|EVAL_FD, (double*)position_and_derivs.data());

    Point3 normal = cross(position_and_derivs[1],position_and_derivs[2]);
    return normal.dir();

}

void NanosplinePatch::eval_unsafe(double* xy, int flag, double* ret) {
    xy_to_patch_coords(xy, flag, ret);
}
// ---------------------------------------------------------------------- 
int NanosplinePatch::estimate_jacobian(double* jac)
{
  double xy[2];  xy[0]=.5;  xy[1]=.5;
    vector<Point3> position_and_derivs(3,Point3(0.));
    xy_to_patch_coords(xy, EVAL_VL|EVAL_FD, (double*)position_and_derivs.data());

    Point3 normal = cross(position_and_derivs[1],position_and_derivs[2]);

  jac[0] = normal.length();
    return 0;
}
// ---------------------------------------------------------------------- 

int NanosplinePatch::xy_to_patch_value(double* xy, int flag, double* ret)
{
    return 0;
}

void NanosplinePatch::bounding_box(Point3 &bounding_box_min,
                                   Point3 &bounding_box_max) {
  Eigen::MatrixXd control_grid = _surface->_patch.get_control_grid();
  Eigen::MatrixXd max = control_grid.colwise().maxCoeff();
  Eigen::MatrixXd min = control_grid.colwise().minCoeff();
  copy_eigen_to_point(max, bounding_box_max);
  copy_eigen_to_point(min, bounding_box_min);
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
