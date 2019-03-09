#include "interpolate.hpp"
#include "vecmatop.hpp"

BEGIN_EBI_NAMESPACE
using namespace Interpolate;


DblNumVec Interpolate::compute_barycentric_weights_chebyshev_first_kind(int n){
    DblNumVec weights(n);
    setvalue(weights, 1.0);
    for(int i = 0; i < n; i++){
        weights(i) = pow(1., i)*sin((2*i+1)*M_PI/(2*n));
    }
    return weights;
}



DblNumVec Interpolate::evaluate_lagrange_interpolant_1d(
        DblNumVec interpolation_nodes,
        DblNumVec function_values_at_nodes,
        DblNumVec evaluation_points){
    ebiAssert(interpolation_nodes.m() == function_values_at_nodes.m());

    DblNumVec interpolant_values(evaluation_points.m());
    for(int i = 0; i < evaluation_points.m(); i++){
        DblNumVec lagrange_weights(interpolation_nodes.m()); //l_j(x)
        setvalue(lagrange_weights, 1.);

        double x = evaluation_points(i);
        for(int j = 0; j < interpolation_nodes.m(); j++){
            double x_j = interpolation_nodes(j);
            double f_j = function_values_at_nodes(j);
            
            if(fabs(x - x_j) <=1e-15){
                interpolant_values(i) = f_j; 
                break;
            }

            for(int m = 0; m < interpolation_nodes.m(); m++){
                if(m != j){
                    double x_m = interpolation_nodes(m);
                    lagrange_weights(j) *= (x - x_m)/(x - x_j);
                }
            }
            interpolant_values(i) += lagrange_weights(j)*f_j;
        }
    }
    return interpolant_values;

}
DblNumMat Interpolate::evaluate_lagrange_interpolant_1d(
        DblNumVec interpolation_nodes,
        DblNumMat function_values_at_nodes,
        DblNumVec evaluation_points){
    ebiAssert(interpolation_nodes.m() == function_values_at_nodes.n());
    
    int dof = function_values_at_nodes.m();
    DblNumMat interpolant_values(dof, evaluation_points.m());
    for(int d = 0; d < dof; d++){
        for(int i = 0; i < evaluation_points.m(); i++){
            DblNumVec lagrange_weights(interpolation_nodes.m()); //l_j(x)
            setvalue(lagrange_weights, 1.);

            double x = evaluation_points(i);
            for(int j = 0; j < interpolation_nodes.m(); j++){
                double x_j = interpolation_nodes(j);
                double f_j = function_values_at_nodes(d,j);

                if(fabs(x - x_j) <=1e-15){
                    interpolant_values(d,i) = f_j; 
                    break;
                }

                for(int m = 0; m < interpolation_nodes.m(); m++){
                    if(m != j){
                        double x_m = interpolation_nodes(m);
                        lagrange_weights(j) *= (x - x_m)/(x - x_j);
                    }
                }
                interpolant_values(d,i) += lagrange_weights(j)*f_j;
            }
        }
    }
    return interpolant_values;

}



void clenshaw_curtis_compute ( int order, DblNumVec& x, DblNumVec& w )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = 1.0.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    1 <= ORDER.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  double b;
  //int i;
  //int j;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 ){
      ebiAssert(order >= 1);
  }
  else if ( order == 1 ){
    x(0) = 0.0;
    w(0) = 2.0;
  }
  else{
    for (int i = 0; i < order; i++ ){
      x(i) =  cos(double( order- i - 1) * pi/double(order-1));
    }
    x(0) = -1.0;

    if ( (order % 2) == 1 )
    {
      x((order-1)/2) = 0.0;
    }
    x(order-1) = 1.0;

    for (int i = 0; i < order; i++ ){
      theta = double(i) * pi/double( order-1 );

      w(i) = 1.0;

      for (int j = 1; j <= ( order - 1 ) / 2; j++ ){
        if ( 2*j == ( order - 1 ) ){
          b = 1.0;
        } else { 
          b = 2.0;
        }

        w(i) = w(i) - b*cos(2.0 * double(j) * theta )/ double( 4*j*j - 1);
      }
    }
    w(0) = w(0) / ( double ) ( order - 1 );

    for (int i = 1; i < order - 1; i++ ){
      w(i) = 2.0 * w(i) / ( double ) ( order - 1 );
    }

    w(order-1) = w(order-1) / ( double ) ( order - 1 );
  }

}

void rescale ( double a, double b, int n, DblNumVec& x, DblNumVec& w )

//****************************************************************************80
//
//  Purpose:
//
//    RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2009
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the new interval.
//
//    Input, int N, the order.
//
//    Input/output, double X[N], on input, the abscissas for [-1,+1].
//    On output, the abscissas for [A,B].
//
//    Input/output, double W[N], on input, the weights for [-1,+1].
//    On output, the weights for [A,B].
//
{

  for (int i = 0; i < n; i++ ){
    x(i) = ( ( a + b ) + ( b - a ) * x(i) ) / 2.0;
  }

  for (int i = 0; i < n; i++ ){
    w(i) = ( b - a ) * w(i) / 2.0;
  }
}

DblNumVec compute_chebyshev_coefficients(double a, double b, DblNumVec nodes, DblNumVec function_values){
    int n = nodes.m();
    assert(n == function_values.m());
    // Form chebyshev_vandermonde matrix
    DblNumMat chebyshev_vandermonde(n,n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double xt = 2./(b-a)*(nodes(i)-a)-1.;
            chebyshev_vandermonde(i,j) = cos(j*acos(xt));
        }
    }

    DblNumMat inverse(n,n);
    inv(chebyshev_vandermonde, inverse);
    DblNumVec chebyshev_coeffs(n);
    setvalue(chebyshev_coeffs, 0.);
    dgemv(1., inverse, function_values, 0., chebyshev_coeffs);
    return chebyshev_coeffs ;
}

double Interpolate::integrate_ith_lagrange_basis_func(int i, double a, double b,
 int basis_degree,
        DblNumVec interpolation_nodes, int quad_order,double jacobian){

    DblNumVec integration_nodes(quad_order);
    DblNumVec function_values(basis_degree);
    setvalue(function_values, 0.); 
function_values(i) = 1.*jacobian;

    /*
    DblNumVec nodes(quad_order);
    DblNumVec quad_weights(quad_order);
    clenshaw_curtis_compute(quad_order, nodes, quad_weights);
    rescale(a, b, quad_order, nodes, quad_weights);
    */

    //assert(0); // TODO BUG chebyshev quad weights are wrong!!!!
    DblNumVec cheb_coefficients = 
        compute_chebyshev_coefficients(a, b, interpolation_nodes, function_values);
    
    int n = interpolation_nodes.m();
    DblNumVec cheb_weights(n);
    setvalue(cheb_weights, 0.);
    for (int i = 0; i < n; i++) {
        if(i % 2 == 0){
            cheb_weights(i) = 2./double(1-i*i)*(b-a)/2.;
        }
    }
    
/*
    // Compute barycentric interpolation weights
    DblNumVec weights= compute_barycentric_weights_1d(interpolation_nodes);
    DblNumVec integration_node_values = evaluate_barycentric_interpolant_1d(
            interpolation_nodes, 
            weights, 
            function_values, 
            nodes);
    // Trapezoidal rule

    // Form chebyshev Vandermonde matrix to form coefficients of the basis
    // function
    
    double integral = 0;
    for(int k = 0; k < quad_order; k++){
        //cout << "(node, weight) " << k << ": " << nodes(k) << ", " << quad_weights(k) << endl;
        integral += integration_node_values(k)*quad_weights(k);
    }*/

    double exact_integral = dot(cheb_weights, cheb_coefficients);
    return exact_integral;
}

void Interpolate::evaluate_barycentric_interpolant_2d(
        int dof,
        DblNumMat& xy_coordinates, 
        DblNumMat& function_values,
        int num_samples_1d,
        DblNumMat& refined_xy_coordinates,
        DblNumMat& refined_function_values){
    ebiAssert(refined_xy_coordinates.n() == refined_function_values.n());
    ebiAssert(xy_coordinates.n() == function_values.n());
    setvalue(refined_function_values, 0.);

    /*
     * Let s_ij = (x_i, y_j) be the sample points in xy_coordinates and
     * refined_xy_coordinates. We assume the sample points are on a 
     * tensor-product grid, i.e. x_i = y_i and (s_ik)_1 = (s_ij)_1 for any j,k
     *
     */


    DblNumVec interpolation_nodes_x(num_samples_1d);
    DblNumVec interpolation_nodes_y(num_samples_1d);
    for(int i =0; i < num_samples_1d; i++){
        interpolation_nodes_x(i) = xy_coordinates(0,i);
        interpolation_nodes_y(i) = xy_coordinates(1,i*num_samples_1d); 
    }


    // TODO implement 2d barycentric interpolation

    // O(n^2) precomputation of weights; we assume that the x and y
    // interpolation nodes are the same.
    DblNumVec barycentric_weights_x = 
        compute_barycentric_weights_1d<double>(interpolation_nodes_x);
    DblNumVec barycentric_weights_y = 
        compute_barycentric_weights_1d<double>(interpolation_nodes_y);
    evaluate_barycentric_interpolant_2d(dof, xy_coordinates, function_values, num_samples_1d, interpolation_nodes_x, interpolation_nodes_y, barycentric_weights_x, barycentric_weights_y, refined_xy_coordinates, refined_function_values);    
}
void Interpolate::evaluate_barycentric_interpolant_2d(
        int dof,
        DblNumMat& xy_coordinates, 
        DblNumMat& function_values,
        int num_samples_1d,
        DblNumVec& interpolation_nodes_x,
        DblNumVec& interpolation_nodes_y,
        DblNumVec& barycentric_weights_x,
        DblNumVec& barycentric_weights_y,
        DblNumMat& refined_xy_coordinates,
        DblNumMat& refined_function_values){

    // iterate  over coefficients of function values.... TODO bring inside loop
    // iteration over desired samples of interpolant
    for(int si =0; si < refined_xy_coordinates.n(); si++){
        Point2 eval_point(refined_xy_coordinates.clmdata(si));

        DblNumVec distance_to_nodes_x(num_samples_1d);
        DblNumVec distance_to_nodes_y(num_samples_1d);
        for(int i =0; i < num_samples_1d; i++){
            double node_x= interpolation_nodes_x(i);
            double node_y= interpolation_nodes_y(i); 
            Point2 interp_point(node_x,node_y);
            Point2 difference = eval_point - interp_point;

            // specify a minimum distance to prevent swamping 
            if(fabs(difference.x()) <= 1e-15)
                difference.x() = 1e-15;
            if(fabs(difference.y()) <= 1e-15)
                difference.y() = 1e-15;

            // invert to use dot product below
            distance_to_nodes_x(i) = 1./difference.x();
            distance_to_nodes_y(i) = 1./difference.y();
        }

        double denominator = 
            dot(barycentric_weights_x, distance_to_nodes_x)*
            dot(barycentric_weights_y, distance_to_nodes_y);

        DblNumVec numerator(dof);
        for(int i =0; i < num_samples_1d; i++){
            for(int j =0; j < num_samples_1d; j++){
                    double w_i = barycentric_weights_x(i);
                    double w_j = barycentric_weights_y(j);
                    
                    double x_minus_x_i_inv = distance_to_nodes_x(j);
                    double y_minus_y_j_inv = distance_to_nodes_y(i);
                    for(int d =0; d < dof; d++){
                        double f_ij = function_values(d, i*num_samples_1d + j);
                        
                        numerator(d) += w_i * w_j * 
                            x_minus_x_i_inv * y_minus_y_j_inv * f_ij;
                    }
            }
        }
        
        for(int d =0; d < dof; d++){
            refined_function_values(d,si) = numerator(d)/denominator;
        }
    }
}

END_EBI_NAMESPACE
