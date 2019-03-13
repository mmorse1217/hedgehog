BEGIN_EBI_NAMESPACE
template<typename T>
NumVec<T> Interpolate::compute_barycentric_weights_1d(NumVec<T> interpolation_nodes){
    int n = interpolation_nodes.m();
    //pvfmm::Profile::Add_FLOP();
  //pvfmm::Profile::Add_FLOP((long long)(n-1)*(long long)(n-1)*(26*dof));
    NumVec<T> weights(n);
    setvalue(weights, T(1.0));

    // Compute w_i = 1/prod_{j \not = i} (x_i - x_j)
    for(int i = 0; i < n; i++){

        for(int j = 0; j < n; j++){
            if(j != i){
                weights(i) *= T(1.)/(interpolation_nodes(j) - interpolation_nodes(i));
            }
        }
        /*for(int j = 0; j < i; j++){
            weights(i) *= 1./(interpolation_nodes(j) - interpolation_nodes(i));
        }
        // skip the ith product
        for(int j = i+1; j < n; j++){
            weights(i) *= 1./(interpolation_nodes(j) - interpolation_nodes(i));
        }*/
    }

    return weights;

}


template<typename T>
NumVec<T> Interpolate::evaluate_barycentric_interpolant_1d(
        NumVec<T> interpolation_nodes,
        NumVec<T> barycentric_weights, 
        NumVec<T> function_values_at_nodes,
        NumVec<T> evaluation_points){
    ebiAssert(interpolation_nodes.m() == function_values_at_nodes.m());
    ebiAssert(interpolation_nodes.m() == barycentric_weights.m());
    
    // Evaluating equation 4.2 in Barycentric Lagrange Interpolation by L.N.
    // Trefethen
    
    // stores p(x_i) for i=1,..., evaluation_points.m()
    NumVec<T> interpolant_values(evaluation_points.m());

    for(int i = 0; i < evaluation_points.m(); i++){
        T numerator = 0; // to be \sum_j=0^n (w_j/(x - x_j))*f_j
                              // j = 1..interpolation_nodes.m()
        T denominator= 0;// to be \sum_j=0^n (w_j/(x - x_j))
                              // j = 1..interpolation_nodes.m(
        
        T x = evaluation_points(i);
        for(int j = 0; j < interpolation_nodes.m(); j++){
            T x_j = interpolation_nodes(j);
            T w_j = barycentric_weights(j);
            T f_j = function_values_at_nodes(j);

            if(fabs(x - x_j) <=1e-15){
                numerator = f_j;
                denominator = 1.;
                break;

            }

            T coeff_j = w_j/(x - x_j);
            
            numerator   += coeff_j*f_j;
            denominator += coeff_j;
        }
        
        interpolant_values(i) = numerator/denominator;
    }
    return interpolant_values;
}
/*
 * compute m x n matrix D with entries D[i,j] = points[i] - nodes[j], 
 * where m = size(points) and n = size(nodes)
 */
template<typename T>
NumMat<T> compute_distance_matrix(NumVec<T> nodes, NumVec<T> points){
    int m = points.m();
    int n = nodes.m();
    NumMat<T> distance_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            distance_matrix(i,j) = points(i) - nodes(j);
        }
    }
    return distance_matrix;
}

template<typename T>
NumVec<T> Interpolate::evaluate_barycentric_derivative_1d(
        NumVec<T> interpolation_nodes,
        NumVec<T> barycentric_weights, 
        NumVec<T> function_values_at_nodes,
        NumVec<T> evaluation_points){
    int m = evaluation_points.m();
    int n = interpolation_nodes.m();
    assert(n == function_values_at_nodes.m());
    assert(n == function_values_at_nodes.m());
    
    NumMat<T> D = compute_distance_matrix<T>(interpolation_nodes, interpolation_nodes);
    for (int i = 0; i < n; i++) {
        D(i,i)  = 0;
        for (int j = 0; j < n; j++) {
            if(i != j){
                D(i,j) = barycentric_weights(j)/barycentric_weights(i)/D(i,j);
                D(i,i) -= D(i,j);
            } 
        }
    }

    NumVec<T> derivative_values(n);
    setvalue(derivative_values, 0.);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            derivative_values(i) += D(i,j)*function_values_at_nodes(j);
        }
    }

    return evaluate_barycentric_interpolant_1d<double>(
            interpolation_nodes, barycentric_weights, 
            derivative_values, evaluation_points);
}
template<typename T>
NumMat<T> Interpolate::evaluate_barycentric_interpolant_1d(
        NumVec<T> interpolation_nodes,
        NumVec<T> barycentric_weights, 
        NumMat<T> function_values_at_nodes,
        NumVec<T> evaluation_points){
    ebiAssert(interpolation_nodes.m() == function_values_at_nodes.n());
    ebiAssert(interpolation_nodes.m() == barycentric_weights.m());

    // Evaluating equation 4.2 in Barycentric Lagrange Interpolation by L.N.
    // Trefethen

    int dof = function_values_at_nodes.m();

    // stores p(x_i) for i=1,..., evaluation_points.m()
    NumMat<T> interpolant_values(dof, evaluation_points.m());
    for(int d = 0; d < dof; d++){
        for(int i = 0; i < evaluation_points.m(); i++){
            T numerator = 0; // to be \sum_j=0^n (w_j/(x - x_j))*f_j
            // j = 1..interpolation_nodes.m()
            T denominator= 0;// to be \sum_j=0^n (w_j/(x - x_j))
            // j = 1..interpolation_nodes.m(

            T x = evaluation_points(i);
            for(int j = 0; j < interpolation_nodes.m(); j++){
                T x_j = interpolation_nodes(j);
                T w_j = barycentric_weights(j);
                T f_j = function_values_at_nodes(d,j);

                if(fabs(x - x_j) <=1e-15){
                    numerator = f_j;
                    denominator = 1.;
                    break;

                }

                T coeff_j = w_j/(x - x_j);

                numerator   += coeff_j*f_j;
                denominator += coeff_j;
            }
            interpolant_values(d,i) = numerator/denominator;
        }
    }
    return interpolant_values;
}

END_EBI_NAMESPACE
