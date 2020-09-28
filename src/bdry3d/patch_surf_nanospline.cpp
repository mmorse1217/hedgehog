#include <nanospline/PatchBase.h>
#include <nanospline/forward_declaration.h>
#include <nanospline/BSplinePatch.h>
#include <nanospline/save_obj.h>
#include <Eigen/Core>
#include "patch_surf_nanospline.hpp"
#include "common/interpolate.hpp"
#include "common/utils.hpp"
#include <sampling.hpp>
#include <common/vtk_writer.hpp>

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
    cout << "greville_abscissae : " << greville_abscissae << endl;
    return greville_abscissae;
}

Eigen::MatrixXd initialize_flat_patch_knots(const int num_control_points_1D,
                                            const int degree) {
  const int num_knots = num_control_points_1D + degree + 1;
  Eigen::MatrixXd knots(num_knots, 1);
  int counter = 0;
  const double knot_spacing = 1. / (num_control_points_1D - degree - 1 + 1);
  for (int i = 0; i < degree + 1; i++) {
    knots(counter) = (i-(degree) )*knot_spacing;
    counter++;
  }


  // uniform knots for now.
  for (int i = 0; i < num_control_points_1D - degree - 1; i++) {
    knots(counter) = (i+1) * knot_spacing;
    counter++;
  }

  for (int i = 0; i < degree + 1; i++) {
    knots(counter) = (num_control_points_1D - degree + i)*knot_spacing;
    counter++;
  }

  assert(counter == num_knots);
  /*int counter = 0;
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
  assert(counter == num_knots);*/
  // uniform knots for now.
  /*const double knot_spacing = 1. / (num_knots - 1);
  for (int i = 0; i < num_knots; i++) {
    knots(i) = (i) * knot_spacing;
  }*/
  cout << "knots: "<< knots << endl;
  return knots;
}

Eigen::MatrixXd initialize_flat_patch_control_grid(const int num_control_points_1D, const int degree_u, 
        const int degree_v){
    const int num_control_points = num_control_points_1D*num_control_points_1D;
    Eigen::MatrixXd control_grid(num_control_points, DIM);

    Eigen::MatrixXd knots_u = initialize_flat_patch_knots(num_control_points_1D, degree_u);
    Eigen::MatrixXd knots_v = initialize_flat_patch_knots(num_control_points_1D, degree_v);
    //Eigen::MatrixXd greville_u = get_greville_abscissae(knots_u, degree_u);
    //Eigen::MatrixXd greville_v = get_greville_abscissae(knots_v, degree_v);

    // identify x,y in R^3 with (u,v) in spline parameter space
    // Place control points  equispaced on [-1,1] x [-1,1] x {0}
    for (int ci = 0; ci < num_control_points_1D; ci++) { // u
      for (int cj = 0; cj < num_control_points_1D; cj++) { // v
          const int index = ci*num_control_points_1D + cj;
          
          // x and y spacing are the same
          const double spacing = 1./(num_control_points_1D-1);
          //const double x = 2.*greville_u(ci) - 1.;
          //const double y = 2.*greville_v(cj) - 1.;
          const double x = spacing*ci;// - 1.;
          const double y = spacing*cj;// - 1.;
          const double z = 0.; 
          control_grid.row(index) << x, y, z;
      }
    }

    return control_grid;
    
}

class NanosplinePatch::NanosplineInterface {
private:
  const int _num_control_pts_u;
  const int _num_control_pts_v;
  const int _degree_u;
  const int _degree_v;
  Eigen::MatrixXd _knots_u;
  Eigen::MatrixXd _knots_v;

public:
  nanospline::BSplinePatch<Scalar, DIM, -1, -1> _patch;
  NanosplineInterface()
      : _degree_u(2), _degree_v(2), _num_control_pts_u(6),
        _num_control_pts_v(6) {
    // 1. Get the name of surface we want to load
    // For now: load a flat surface, all weights = 1, uniform knots
    //
    assert(_num_control_pts_u == _num_control_pts_v);
    const int num_control_points = _num_control_pts_u * _num_control_pts_v;

    // 2. Load surface data with obj reader
    // For now: explicitly generate control points, knots and weights.
    Eigen::MatrixXd control_grid = initialize_flat_patch_control_grid(
        _num_control_pts_u, _degree_u, _degree_v);
    _patch.set_control_grid(control_grid);

    _knots_u = initialize_flat_patch_knots(_num_control_pts_u, _degree_u);
    _knots_v = initialize_flat_patch_knots(_num_control_pts_v, _degree_v);
    _patch.set_knots_u(_knots_u);
    _patch.set_knots_v(_knots_v);

    _patch.set_degree_u(_degree_u);
    _patch.set_degree_v(_degree_v);
    cout << "u bound " <<  _patch.get_u_lower_bound() << ", "<<_patch.get_u_upper_bound() << endl;
    cout << "v bound " <<  _patch.get_v_lower_bound() << ", "<<_patch.get_v_upper_bound() << endl;
    _patch.initialize();
    nanospline::save_patch_obj("flat_patch.obj", _patch);
  }
  

  void deform_periodic(DblNumMat _parametric_coordinates, DblNumMat _changes_in_position, 
          Eigen::Vector<double, DIM> cell_size, int periodic_dim_u, int periodic_dim_v){
    assert(_parametric_coordinates.n() == _changes_in_position.n());
    assert(_parametric_coordinates.m() == 2);
    assert(_changes_in_position.m() == DIM);
    assert(periodic_dim_u >= 0 && periodic_dim_u <= 2);
    assert(periodic_dim_v >= 0 && periodic_dim_v <= 2);
    // 0. copy into eigen matrices

    const int num_constraints = _parametric_coordinates.n();
    Eigen::MatrixXd parametric_coordinates(num_constraints, 2);
    Eigen::MatrixXd changes_in_position(num_constraints, DIM);
    
    for (int i =0; i < num_constraints; i++) {
        parametric_coordinates.row(i) << _parametric_coordinates(0,i), _parametric_coordinates(1,i);
        changes_in_position.row(i) << _changes_in_position(0,i), _changes_in_position(1,i), _changes_in_position(2,i);
    }
      
    // 1. get vanilla least squares matrix
    // TODO replace with call to nanospline when p[ull request approved
    const int num_control_pts = _num_control_pts_u * _num_control_pts_v;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basis_func_control_pts(num_control_pts, 1);
    nanospline::BSplinePatch<Scalar, 1, -1, -1> basis_function;
    basis_func_control_pts.setZero();
    basis_function.set_knots_u(_knots_u);
    basis_function.set_knots_v(_knots_v);
    basis_function.set_degree_u(_degree_u);
    basis_function.set_degree_v(_degree_v);

    // Suppose we are fitting samples p_0, ..., p_n of a function f, with
    // parameter values (u_0,v_0) ..., (u_n,v_n). If B_q^n(u,v) is the
    // jth basis function value at (u,v), then
    // least_squares_matrix(i,q) = B_q^n(u_i, v_i),
    // where q is the linearized index over basis elements:
    // q = i*num_control_pts_v + j
    Eigen::MatrixXd least_squares_matrix =
        Eigen::MatrixXd::Zero(num_constraints, num_control_pts);
    cout << "evaluate basis functions" << endl;
    for (int j = 0; j < _num_control_pts_u; j++) {
        for (int k = 0; k < _num_control_pts_v; k++) {
            int index = j * _num_control_pts_v + k; // TODO bug prone
            basis_func_control_pts(index) = 1.;
            basis_function.set_control_grid(basis_func_control_pts);
            basis_function.initialize();
            for (int i = 0; i < num_constraints; i++) {
                Scalar u = parametric_coordinates(i, 0);
                Scalar v = parametric_coordinates(i, 1);
                Scalar bezier_value = basis_function.evaluate(u, v)(0);
                least_squares_matrix(i, index) = bezier_value;
            }
            //basis_func_control_pts.row(index) << 0.;
            basis_func_control_pts.setZero();
        }
    }
    
    // 2. get selection matrices to select control points
    /* We need to extract the  control point along the edges of the patch 
     * in order to enforce periodicity. Suppose our spline patch has 
     * u-degree p_u, and v-degree p_v, and is given by 
     *      P(u,v) = \sum_i^{N_u}\sum_j^{N_v} c_ij B_ij^p(u,v), 
     *
     * Suppose that our spline looks something like this:
     *          *  *  *  *  *  *  *  *
     *
     *          *  *--*--*--*--*--*  *
     *             |              | 
     *          *  *__*__*__*__*  *  *
     *             |           |  | 
     *          *  *  *  *  *  *  *  *
     *             |           |  | 
     *          *  *  *  *  *  *  *  *
     *             |           |  | 
     *   ^      *  *  *  *  *  *  *  *
     *   |         |___________|  | 
     *   v      *  *--*--*--*--*--*  *
     *
     *          *  *  *  *  *  *  *  *
     *
     *              u ->
     * The stars (*) correspond to control points c_{ij} above. The outer square
     * (--) corresponds to the periodic cell boundary, and the inner square (__)
     * corresponds to the control points contained strictly inside the periodic
     * cell. The control points outside the cell are "dependent" control points 
     * and the control points inside the cell are "independent" control points. 
     * We want to constrain each dependent control point to uniquely match with 
     * an independent control point. In terms of the picture above:
     *
     *          *  *  *  *  *  *  *  *
     *           S6|     S7    |  S8
     *          *  *--*--*--*--*--*  *
     *             |           |  | 
     *          *__*__*__*__*__*__*__*
     *             |           |  | 
     *          *  *  *  *  *  *  *  *
     *           S4|     I     |  | S5
     *          *  *  *  *  *  *  *  *
     *             |           |  | 
     *   ^      *  *  *  *  *  *  *  *
     *   |       __|___________|__|__ 
     *   v      *  *--*--*--*--*--*  *
     *          S1 |    S2     |  S3
     *          *  *  *  *  *  *  *  *
     *
     *              u ->
     * Each control point in S1, ..., S8 needs to map to a control point in I.
     * If the periodic cell has widths d_u, d_v in the x/y directions, 
     * respectively, the equalities for each dependent control point in each section is as
     * follows:
     *  S1: c_{ij} = c_{i+p_u, j+p_v}   + (-d_u,-d_v, 0)
     *  S2: c_{ij} = c_{i, j+p_v}       + (   0,-d_v, 0)
     *  S3: c_{ij} = c_{i-p_u, j+p_v}   + ( d_u,-d_v, 0)
     *  S4: c_{ij} = c_{i+p_u, j}       + (-d_u,   0, 0)
     *  S5: c_{ij} = c_{i-p_u, j}       + ( d_u,   0, 0)
     *  S6: c_{ij} = c_{i+p_u, j-p_v}   + (-d_u, d_v, 0)
     *  S7: c_{ij} = c_{i, j-p_v}       + (   0, d_v, 0)
     *  S8: c_{ij} = c_{i-p_u, j-p_v}   + ( d_u, d_v, 0)
     *  
     * (This diagram is consistent with the image above, with xy plane in R3 in
     * line with the patch; signs can change with a change of coordinates and
     * periodic shifts correspond to the periodic portion of the cell.
     *
     * The indices for each section are as follows (note that control points
     * on the edge and corner of S1, S3, S6, S8 belong to these sections. S2,
     * S4, S5, and S7 contain the control points on the edges of the periodic
     * cell):
     *  S1: 0 <= i <=1
     *      0 <= j <=1
     *
     *  S2: 2 <= i <= num_control_pts_u - degree_u-2
     *      0 <= j <=1
     *
     *  S3: num_control_pts_u-degree_u-1 <= i <=num_control_pts_u-1
     *      0 <= j <=1
     *
     *  S4: 0 <= i <=1
     *      2 <= j <= num_control_pts_v - _degree_u - 2
     *
     *  S5: num_control_pts_u-degree_u-1 <= i <=num_control_pts_u-1
     *      2 <= j <= num_control_pts_v - _degree_u - 2
     *
     *  S6: 0 <= i <=1
     *      num_control_pts_v-degree_v -1 <= j <=num_control_pts_v -1 
     *
     *  S7: 2 <= i <= num_control_pts_u - degree_u-2
     *      num_control_pts_v-degree_v -1 <= j <=num_control_pts_v -1
     *
     *  S8: num_control_pts_u-degree_u -1 <= i <=num_control_pts_u -1
     *      num_control_pts_v-degree_v -1 <= j <=num_control_pts_v -1
     */

    cout << "construct selection matrices u " << endl;
    // Get the periodic cell size corresponding to u/v periodicity.
    Eigen::VectorXd periodic_cell_size_u = Eigen::Vector<double, DIM>::Zero();
    Eigen::VectorXd periodic_cell_size_v = Eigen::Vector<double, DIM>::Zero();
    periodic_cell_size_u(periodic_dim_u) = cell_size(periodic_dim_u);
    periodic_cell_size_v(periodic_dim_v) = cell_size(periodic_dim_v);

    //const int num_periodic_control_pts = (_degree_v+1)*_num_control_pts_u + 
    //                (_degree_u + 1)*_num_control_pts_v;
    const int num_periodic_control_pts = (_degree_u )*_num_control_pts_u + 
                    ( _degree_v )*_num_control_pts_v;
    //Eigen::MatrixXd cell_size_constraint(num_periodic_control_pts, DIM);
    
    //cell_size_constraint.setZero();
    int iter = 0;
    
    Eigen::MatrixXd cp = _patch.get_control_grid();
    const int num_internal_knots_u = _num_control_pts_u - _degree_u - 1;
    const int num_internal_knots_v = _num_control_pts_v - _degree_v - 1;
    const int num_independent_ctrl_pts = num_internal_knots_u*num_internal_knots_v;
    const int num_dependent_ctrl_pts = (num_internal_knots_u+_degree_u+1)*(num_internal_knots_v+_degree_v+1) - num_independent_ctrl_pts;
    cout << "num_independent_ctrl_pts: " << num_independent_ctrl_pts << endl;
    cout << "num_dependent_ctrl_pts: " << num_dependent_ctrl_pts << endl;

    const double d_u = cell_size(periodic_dim_u);
    const double d_v = cell_size(periodic_dim_v);
    Eigen::MatrixXd periodic_constraints(num_dependent_ctrl_pts, num_control_pts);
    Eigen::MatrixXd cell_size_constraint(num_dependent_ctrl_pts, DIM);
    cell_size_constraint.setZero();

    auto fill_constraint_matrix =
        [this, &iter, &periodic_constraints, &cell_size_constraint](
            int ui_min, int ui_max, int vj_min, int vj_max,
            std::function<int(int)> ui_index, std::function<int(int)> vj_index,
            Eigen::MatrixXd periodic_cell_size) {

            assert(int(periodic_cell_size.rows()) == 1);
            assert(int(periodic_cell_size.cols()) == 3);

          for (int ui = ui_min; ui <= ui_max; ui++) {
            for (int vj = vj_min; vj <= vj_max; vj++) {
              const int dependent_ctrl_pt_index = ui * _num_control_pts_v + vj;
              const int independent_ctrl_pt_index = ui_index(ui) * _num_control_pts_v + vj_index(vj);
              cout << "dependent: " << ui << ", "<< vj << endl;
              cout << "independent: " << ui_index(ui) << ", "<< vj_index(vj) << endl;
              periodic_constraints(iter, dependent_ctrl_pt_index) = 1.;
              periodic_constraints(iter, independent_ctrl_pt_index) = -1.;
              cell_size_constraint.row(iter) << periodic_cell_size;
              iter++;
              cout << ui << ", " << vj << ", " << iter << endl;
            }
          }
        };

    /*  
     *  S1: i = 0 
     *      j = 0 
     *  S1: c_{ij} = c_{i+p_u, j+p_v}   + (-d_u,-d_v, 0)
     */
    Eigen::MatrixXd temp(1,3);
    temp << -d_u, -d_v, 0.;
    fill_constraint_matrix(0, 0, 
            0, 0, 
            [this](int ui){return _degree_u + ui;}, 
            [this](int vj){return _degree_v + vj;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    /*
     *  S2: 1 <= i <= num_control_pts_u - degree_u-1
     *      j = 0
     *  S2: c_{ij} = c_{i, j+p_v}       + (   0,-d_v, 0)
     */
    
    temp << 0., -d_v, 0.;
    fill_constraint_matrix(1, _num_control_pts_u- _degree_u - 1, 
            0, 0, 
            [this](int ui){return ui;}, 
            [this](int vj){return _degree_v + vj;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    /*
     *  S3: num_control_pts_u-degree_u <= i <=num_control_pts_u-1
     *      j = 0
     *  S3: c_{ij} = c_{i-p_u, j+p_v}   + ( d_u,-d_v, 0)
    */
    
    temp << d_u, -d_v, 0.;
    fill_constraint_matrix(_num_control_pts_u- _degree_u,_num_control_pts_u - 1,
            0, 0, 
            [this](int ui){return ui-_degree_u;}, 
            [this](int vj){return _degree_v + vj;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    /*
     *  S4: i = 0 
     *      1 <= j <= num_control_pts_v - _degree_u - 1 
     *  S4: c_{ij} = c_{i+p_u, j}       + (-d_u,   0, 0)
     */
    
    temp << -d_u, 0., 0.;
    fill_constraint_matrix(0, 0, 
            1, _num_control_pts_v- _degree_v-1,
            [this](int ui){return ui+_degree_u;}, 
            [this](int vj){return vj;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    /*
     *  S5: num_control_pts_u-degree_u <= i <=num_control_pts_u-1
     *      1 <= j <= num_control_pts_v - _degree_u - 1
     *  S5: c_{ij} = c_{i-p_u, j}       + ( d_u,   0, 0)
     */
    
    temp << d_u, 0., 0.;
    fill_constraint_matrix(_num_control_pts_u - _degree_u, _num_control_pts_u - 1, 
            1, _num_control_pts_v - _degree_v - 1,
            [this](int ui){return ui - _degree_u;}, 
            [this](int vj){return vj;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    cout << "S6" << endl;
    /*
     *  S6:  i = 0
     *      num_control_pts_v-degree_v <= j <=num_control_pts_v -1 
     *  S6: c_{ij} = c_{i+p_u, j-p_v}   + (-d_u, d_v, 0)
     */
    
    temp << -d_u, d_v, 0.;
    fill_constraint_matrix(0, 0, 
            _num_control_pts_v - _degree_v , _num_control_pts_v - 1,
            [this](int ui){return ui + _degree_u;}, 
            [this](int vj){return vj - _degree_v;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    cout << "S7" << endl;
    /*
     *  S7: 1 <= i <= num_control_pts_u - degree_u-1
     *      num_control_pts_v-degree_v <= j <=num_control_pts_v -1
     *  S7: c_{ij} = c_{i, j-p_v}       + (   0, d_v, 0)
     */
    
    temp << 0., d_v, 0.;
    fill_constraint_matrix(1, _num_control_pts_u - _degree_u - 1, 
            _num_control_pts_v - _degree_v , _num_control_pts_v - 1,
            [this](int ui){return ui ;}, 
            [this](int vj){return vj - _degree_v;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    cout << "S8" << endl;
    /*
     *  S8: num_control_pts_u-degree_u -1 <= i <=num_control_pts_u -1
     *      num_control_pts_v-degree_v -1 <= j <=num_control_pts_v -1
     *  S8: c_{ij} = c_{i-p_u, j-p_v}   + ( d_u, d_v, 0)
     */
    
    temp << d_u, d_v, 0.;
    fill_constraint_matrix(_num_control_pts_u - _degree_u, _num_control_pts_u - 1, 
            _num_control_pts_v - _degree_v , _num_control_pts_v - 1,
            [this](int ui){return ui - _degree_u;}, 
            [this](int vj){return vj - _degree_v;}, 
            temp);
    cout << "iter: "<< iter << ", " << periodic_constraints.rows() << " x  " << periodic_constraints.cols() << endl;
    cout << "constraint matrix " << endl << periodic_constraints << endl;
   

    /*
     * The final system we want to form is 
     * [2*A^T*A  E^T ]*[c]   [2*A^Tp]
     * [E        0   ] [y] = [     d]
     * E = periodic_constraints
     *  size: num_dependent_ctrl_pts x num_control_pts
     * A = least_squares_matrix for a standard least squares solve with BSplines
     *  size: num_constraints x num_control_pts  
     *  (A^T*A = num_control_pts x *  num_control_pts)
     * p = points (or deformations) we want to fit 
     *  size: num_constraints x DIM
     * c = control points (or changes in control points) of the final patch
     *  size: num_control_pts x DIM
     * y = lagrange multiplier for the periodic constraint imposed by E
     *  size: num_dependent_ctrl_pts x DIM
     * d = vector of d_1/d_2 values from above in the appropriate place
     *  size: num_dependent_ctrl_pts x DIM
     */
    cout << "form full system" << endl;
    Eigen::MatrixXd normal_equations = least_squares_matrix.transpose()*least_squares_matrix;
    //Eigen::MatrixXd periodic_constraints = top_left_selector + bottom_right_selector;

    cout << "test selection matrix: "<<endl << periodic_constraints*_patch.get_control_grid() << endl;
    cout << "constraint: "<<endl << cell_size_constraint << endl;
    cout << "diff: "<<endl << periodic_constraints*_patch.get_control_grid() - cell_size_constraint << endl;
    
    const int constrained_system_size = num_control_pts + num_dependent_ctrl_pts;
    Eigen::MatrixXd right_hand_side(constrained_system_size, DIM);
    Eigen::MatrixXd full_least_sq_system(constrained_system_size, constrained_system_size);
    full_least_sq_system.setZero();
    right_hand_side.setZero();
    
    full_least_sq_system.topLeftCorner(num_control_pts, num_control_pts)          = 2*normal_equations;
    full_least_sq_system.topRightCorner(num_control_pts, num_dependent_ctrl_pts)  = periodic_constraints.transpose();
    full_least_sq_system.bottomLeftCorner(num_dependent_ctrl_pts,num_control_pts) = periodic_constraints;

    

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(least_squares_matrix);
    double cond = svd.singularValues()(0)
    / svd.singularValues()(svd.singularValues().size()-1);
    cout << "condition number of A: " << cond << endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd3(normal_equations);
    cond = svd3.singularValues()(0)
    / svd3.singularValues()(svd3.singularValues().size()-1);
    cout << "condition number of A^TA: " << cond << endl;
    //cout << "normal_equations: "  << endl << normal_equations<< endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(full_least_sq_system);
    cond = svd2.singularValues()(0)
    / svd2.singularValues()(svd2.singularValues().size()-1);
    cout << "condition number of full system: " << cond << endl;
    


    right_hand_side.topLeftCorner(num_control_pts, DIM) = 2*least_squares_matrix.transpose()*changes_in_position;
    right_hand_side.bottomLeftCorner(num_dependent_ctrl_pts, DIM) = cell_size_constraint;
    
    Eigen::MatrixXd solution = 
        full_least_sq_system.fullPivLu().solve(right_hand_side);
    Eigen::MatrixXd control_point_updates = solution.topLeftCorner(num_control_pts, DIM);
    cout << "deformations: " << control_point_updates << endl;
    
    Eigen::MatrixXd original_control_points = _patch.get_control_grid();
    Eigen::MatrixXd deformed_control_points = cp + control_point_updates;
    cout << "new cp's: " << deformed_control_points<< endl;

    cout << "are constraints satisfied" << endl;
    cout << periodic_constraints*deformed_control_points<< endl;
    _patch.set_control_grid(deformed_control_points);
    nanospline::save_patch_obj("deformed_flat_patch.obj", _patch);
    DblNumMat temp_control_points(DIM, num_control_pts);
    NumVec<OnSurfacePoint> temp_closest_points(num_control_pts);
    for (int i = 0; i < num_control_pts; i++) {
        temp_closest_points(i).parent_patch = i;
        for(int d = 0; d < DIM; d++){
            temp_control_points(d,i) = original_control_points(i,d);
        }
    }
    write_qbkix_points_to_vtk(temp_control_points,temp_closest_points,0,"original");
    for (int i = 0; i < num_control_pts; i++) {
        temp_closest_points(i).parent_patch = i;
        for(int d = 0; d < DIM; d++){
            temp_control_points(d,i) = deformed_control_points(i,d);
        }
    }
    write_qbkix_points_to_vtk(temp_control_points,temp_closest_points,0,"deformed");
  }
  ~NanosplineInterface() {}
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
void NanosplinePatch::deform_periodic(DblNumMat parametric_coordinates, DblNumMat changes_in_position){
    Eigen::VectorXd cell_size = Eigen::VectorXd::Constant(DIM, 1.); //TODO WHAT IS HAPPENING
    _surface->deform_periodic(parametric_coordinates, changes_in_position, cell_size, 0, 1);
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
