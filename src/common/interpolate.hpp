/*! \file */
#ifndef _INTERPOLATE_HPP_
#define _INTERPOLATE_HPP_

//#include "ebi_namespace.hpp"
#include "nummat.hpp"
#include "vec3t.hpp"
#include "bdry3d/patch_surf.hpp"

BEGIN_EBI_NAMESPACE
namespace Interpolate{
    template<typename T>
    NumVec<T> compute_barycentric_weights_1d(NumVec<T> interpolation_nodes);
    
    template<typename T>
    NumVec<T> evaluate_barycentric_interpolant_1d(NumVec<T> interpolation_nodes,
            NumVec<T> barycentric_weights, NumVec<T> function_values_at_nodes,
            NumVec<T> evaluation_points);

    template<typename T>
    NumVec<T> evaluate_barycentric_derivative_1d(
            NumVec<T> interpolation_nodes,
            NumVec<T> barycentric_weights, 
            NumVec<T> function_values_at_nodes,
            NumVec<T> evaluation_points);
    
    // For vector valued functions
    template<typename T>
    NumMat<T> evaluate_barycentric_interpolant_1d(
            NumVec<T> interpolation_nodes,
            NumVec<T> barycentric_weights, 
            NumMat<T> function_values_at_nodes,
            NumVec<T> evaluation_points);

    DblNumVec evaluate_lagrange_interpolant_1d(
        DblNumVec interpolation_nodes,
        DblNumVec function_values_at_nodes,
        DblNumVec evaluation_points);

    DblNumMat evaluate_lagrange_interpolant_1d(
        DblNumVec interpolation_nodes,
        DblNumMat function_values_at_nodes,
        DblNumVec evaluation_points);

    /*
    //template<typename T>
    DblNumMat evaluate_barycentric_interpolant_1d(DblNumVec interpolation_nodes,
            DblNumVec barycentric_weights, DblNumMat function_values_at_nodes,
            DblNumVec evaluation_points);
    */
    double integrate_ith_lagrange_basis_func(int i, double a, double b,int basis_degree, 
            DblNumVec interpolation_nodes, int quad_order, double jacobian);

    void evaluate_barycentric_interpolant_2d(
            int dof,
            DblNumMat xy_coordindates, 
            DblNumMat function_values,
            int num_samples_1d,
            int refinement_factor,
            int axis,
            DblNumMat refined_xy_coordinates,
            DblNumMat& refined_function_values);


}

END_EBI_NAMESPACE

#include "interpolate.tpp"
#endif
