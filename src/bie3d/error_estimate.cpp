#include "error_estimate.hpp"
#include <complex>
#include <sampling.hpp>
BEGIN_EBI_NAMESPACE

ComplexNumVec ErrorEstimate::compute_pullback(
        DblNumMat target_points,
        NumVec<OnSurfacePoint> closest_points,
        CurvatureDirection curvature_direction,
        PatchSurfFaceMap* face_map){
    
    int num_targets = target_points.n();
    assert(num_targets == closest_points.m());
    ComplexNumVec pullbacks(num_targets);

#pragma omp parallel for
    for (int i = 0; i < num_targets; i++) {
        Point3 target(target_points.clmdata(i));
        OnSurfacePoint closest_point = closest_points(i);
        auto patch = face_map->subpatch(closest_point.parent_patch);
        
        pullbacks(i) = compute_pullback_under_patch(target, closest_point,
                                curvature_direction, patch);
        
    }
    return pullbacks;
}

imaginary compute_initial_guess(OnSurfacePoint on_surface_point){
    return 0.;
}

Point2 compute_principal_curvature_direction(OnSurfacePoint on_surface_point, 
        FaceMapSubPatch* patch, ErrorEstimate::CurvatureDirection curvature_direction){
    // TODO  compute eigenvectors of shape operator (analytically)
    if(curvature_direction == ErrorEstimate::CurvatureDirection::FIRST){
        return Point2(1.,0.);
        //return Point2(.5,.5);
    } else {
        return Point2(0.,1.);
        //return Point2(-.5,-.5);
    }
}

void sample_patch_along_direction(int num_samples, Point2 uv_direction, 
        Point2 point_preimage, FaceMapSubPatch* patch,
        DblNumMat& uv_values, DblNumMat& samples_along_line){
    // sample along the line in the uv_direction intersecting point_preimage
    // point_preimage = (u_0, v_0), uv_direction = (u*, v*)
    // line in parameter space is \ell(t) = (u_0, v_0) + t*(u*, v*)
    // find t values of line that intersects x=0, y=0, x=1, y=1.
    // \ell(t) = (0,z) when t = -u_0/u*; \ell(t) = (1,z) when t = (1-u_0)/u*
    // same for v
    double u_0 = point_preimage.x();
    double v_0 = point_preimage.y();
    double u_star = uv_direction.x();
    double v_star = uv_direction.y();
    
    // TODO add check for divisions by zero.
    vector<double> parameter_domain_intercepts;
    if (fabs(u_star) < 1e-12){
        // vertical vector
        parameter_domain_intercepts = {
            -v_0/v_star, (1.-v_0)/v_star
        };

    } else if (fabs(v_star) < 1e-12){
        // horizontal vector
        parameter_domain_intercepts = {
            -u_0/u_star, (1.-u_0)/u_star
        };
    } else {
        vector<double> test_values = {
            -u_0/u_star, (1.-u_0)/u_star, 
            -v_0/v_star, (1.-v_0)/v_star
        };
        for(auto t : test_values){
            Point2 uv_along_line = point_preimage + t*uv_direction;
            double u = uv_along_line.x();
            double v = uv_along_line.y();
            if(0 <=u && u <= 1. &&  0 <=v && v <= 1. ){
                parameter_domain_intercepts.push_back(t);
            }
        }
        // edge case: diagonal line, all points inside the domain
        // solution: choose any two unique points from either end
        if(parameter_domain_intercepts.size() == 4){
            vector<double> temp;
            set<double> value_set(parameter_domain_intercepts.begin(), 
                    parameter_domain_intercepts.end());
            while(!value_set.empty()){
                double fi = *value_set.begin();
                for (double fj : parameter_domain_intercepts){
                    if(fabs(fi - fj) <= 1e-7){
                        temp.push_back(fi);
                        value_set.erase(fj);
                        break;
                    }
                }
                value_set.erase(fi);

            }
            parameter_domain_intercepts = temp;
        }

        std::sort(parameter_domain_intercepts.begin(), parameter_domain_intercepts.end());
    }
    
    assert(parameter_domain_intercepts.size() == 2); 
    Point2 sampling_interval(parameter_domain_intercepts[0], 
            parameter_domain_intercepts[1]);

    // actually sample along line over specified interval
    double step_size = (sampling_interval.y() - sampling_interval.x())/double(num_samples-1);

    //DblNumMat samples_along_line(3, num_samples);
    // TODO sampling is incorrect; need to start sampling at
    // parameter_domain_intercepts[0] i think
    for (int i = 0; i < num_samples; i++) {
        double t = i*step_size;
        Point2 uv_along_line = point_preimage + (sampling_interval.x() + t)*uv_direction;
        
        Point3 position(0.);
        // need eval_unsafe in case we get a negative param value due to numerical
        // issues
        patch->eval_unsafe(uv_along_line.array(), PatchSamples::EVAL_VL, 
               position.array());
        for (int d = 0; d < 2; d++) {
            uv_values(d,i) = uv_along_line(d);
        }
        
        for (int d = 0; d < DIM; d++) {
            samples_along_line(d,i) = position(d); 
        }  
    }

}

imaginary pullback_newton(double curvature, // unneeded, assumed to be on [-1,1]
        double distance_to_target, double eps=1e-7, imaginary initial_guess=0){
    ErrorEstimate::CircleArc c(curvature);
    imaginary z_0(0, distance_to_target);
    
    // objective function to minimize
    auto function = [&c, z_0] (imaginary t)-> imaginary { return c(t) - z_0; };
    // note that derivative of objective function matches the derivative of c
    // because z_0 drops out
    
    // TODO refactor newton out of objective function set up somehow
    // issue is how to wrap function pointer to CircleArc and lambdas
    // consistently.

    imaginary ti = initial_guess;
    imaginary ti_plus_one = 1.;
    imaginary update = 1.;
    int i = 0;
    
    while ( i < 500 && std::abs(ti_plus_one - ti) >= eps && std::abs(update) > eps){
        update = -function(ti)/c.derivative(ti);
        ti_plus_one= ti;
        ti += update;
        i++;
    }
    return ti_plus_one;
    
}
imaginary ErrorEstimate::compute_pullback_under_patch(
        Point3 target,
        OnSurfacePoint closest_point,
        CurvatureDirection curvature_direction,
        FaceMapSubPatch* patch){
    double distance_from_target = closest_point.distance_from_target;
    double L = patch->characteristic_length();
    Point2 uv = closest_point.parametric_coordinates;
    
    double k1, k2;
    patch->principal_curvatures(uv, k1, k2);
    
    double k;
    if(curvature_direction == CurvatureDirection::FIRST){
        k = k1;
    } else {
        k = k2;
    }
    // normal vector is pointing outside of problem domain
    // recall that if curvature is positive, surface is bending toward the
    // normal vector. (concave); if curvature is negative, surface is bending
    // away from the normal vector.
    
    // If point is inside domain and surface curvature is negative, we want to
    // compute the estimate with a curve bending toward the target i.e. with
    // positive curvature.
    // If point is inside domain and surface curvature is positive, we want to
    // compute the estimate with a curve away from the target i.e. with
    // negative curvature.
    
    // If the point is outside the domain, we can use the same curvature sign as
    // the surface; the point is in the proper location relative to the normal.
    /*
    if(closest_point.inside_domain == INSIDE){
        k *= -1.;
    }
    */
    // If the point is outside, move the 2d projection outside the domain in
    // order to recover the correct preimage sign.
    // TODO BUG why is this required?
    /*if(closest_point.inside_domain == OUTSIDE){
        distance_from_target *= -1.;
    }*/
    return pullback_newton(k*L, distance_from_target/L, 1e-7);
}

DblNumMat ErrorEstimate::push_forward(
        ComplexNumVec pullbacks,
        NumVec<OnSurfacePoint> closest_points,
        CurvatureDirection curvature_direction,
        PatchSurfFaceMap* face_map){
    
    int num_points = pullbacks.m();
    assert(num_points == closest_points.m());
    DblNumMat points(3, num_points);

//#pragma omp parallel for
    for (int i = 0; i < num_points; i++) {
        imaginary pullback = pullbacks(i);
        OnSurfacePoint closest_point = closest_points(i);
        auto patch = face_map->subpatch(closest_point.parent_patch);

        double L = patch->characteristic_length();
        Point2 uv = closest_point.parametric_coordinates;

        double k1, k2;
        patch->principal_curvatures(uv, k1, k2);

        double k;
        if(curvature_direction == CurvatureDirection::FIRST){
            k = k1;
        } else {
            k = k2;
        }

        // push forward point in rescaled reference domain
        ErrorEstimate::CircleArc c(k*L);
        imaginary push_forward_ref = c(pullback);
        
        // add patch-relevant length scaling along the normal
        Point3 n = patch->normal(uv.array());
        double scale_factor = push_forward_ref.imag()*L;
        if(closest_point.inside_domain == INSIDE){
            scale_factor *= -1;
        }

        Point3 p;
        patch->xy_to_patch_coords(uv.array(), PatchSamples::EVAL_VL, p.array());
        Point3 push_forward = n*scale_factor + p;

        for (int d = 0; d < DIM; d++) {
            points(d,i) = push_forward(d);
        }

    }
    return points;
}


double evaluate_error_estimate_expression(double curve_derivative, 
        double density_magnitude, imaginary t, 
        int quadrature_order, int m=1){
    
    int q = max(1, quadrature_order-1);
    imaginary tstar(abs(t.real()), abs(t.imag()));
    imaginary square_root_term = sqrt(tstar*tstar - 1.);

    imaginary iq(q,0);
    double first_factor = density_magnitude/pow(curve_derivative, m);
    double second_factor = 4./pow(abs(iq*square_root_term), 3)*
        1./pow(abs(tstar + square_root_term), q)*
        pow(abs(iq/square_root_term), m);

    // recall that \Gamma(n) = (n-1)!
    //return  1./tgamma(m+1)*first_factor*second_factor;
    return  first_factor*second_factor*1./(2*fabs(t.imag()));

}

double ErrorEstimate::evaluate_near_zone_distance(double curvature,
        double density_magnitude, double target_error,
        int quadrature_order, int m, double eps){

    ErrorEstimate::CircleArc c(curvature);
    // objective function to minimize
    auto error = [&c, density_magnitude, target_error, 
                     quadrature_order, m] (double t)-> double {
            imaginary imag_t(0., t);
            double curve_deriv = fabs(c.derivative(imag_t));
            return evaluate_error_estimate_expression(curve_deriv, 
                    density_magnitude, imag_t, quadrature_order, m);
        };
    // note that derivative of objective function matches the derivative of c
    // because z_0 drops out
    
    // TODO refactor newton out of objective function set up somehow
    // issue is how to wrap function pointer to CircleArc and lambdas
    // consistently.
    double upper_bound = 1.;
    double lower_bound = 0.;
    
    for (int i = 0; i < 10; i++) {
        double midpoint = (upper_bound + lower_bound)/2.;
        double error_at_mipoint = error(midpoint);
        if(error_at_mipoint < target_error){
            upper_bound = midpoint;
        } else {
            lower_bound = midpoint;
        }

    }
    return (lower_bound+ upper_bound)/2.;
    
}




double ErrorEstimate::evaluate_error_estimate_on_patch( FaceMapSubPatch* patch,
        CurvatureDirection curvature_direction, Point3 target, 
        OnSurfacePoint closest_point,
        int quadrature_order, DblNumMat density_values){
    assert(quadrature_order == density_values.n());

    double k1, k2;
    patch->principal_curvatures(closest_point.parametric_coordinates,k1,k2);
    double k = curvature_direction == CurvatureDirection::FIRST ? k1 : k2;
    //cout << "rad of curvature: " << k << ", " << closest_point.distance_from_target<< endl;
    if(closest_point.distance_from_target >.85/fabs(k) && k < 0){
        // assume point is far if it's near the radius of curvature where the
        // pullback breaks down.
        return 1e-12;
    }

    imaginary target_pullback= ErrorEstimate::compute_pullback_under_patch(
            target, closest_point, curvature_direction, patch);
    
    //cout << "target:" << target << ", target_pullback: " << target_pullback << endl;

    int density_dof = density_values.m();
    ComplexNumMat complex_densities(density_dof, quadrature_order);
    imaginary r(1., 0.);
    for (int i = 0; i < quadrature_order; i++) {
        for (int k = 0; k < density_dof; k++) {
            complex_densities(k,i) = r*density_values(k,i);
        }
    }
    //   b. compute param values of density 
    //      i. sample as float
    DblNumVec uv_float(quadrature_order, true, 
            Sampling::sample_1d<Sampling::chebyshev1>(quadrature_order, 
                Rectangle(Interval(0.,1.), Interval(0., 1.))).data() );
            
    //      ii. convert to complex (data type are different size so no memcpy)
    ComplexNumVec uv(quadrature_order);
    for (int i = 0; i < quadrature_order; i++) {
        uv(i) = imaginary(uv_float(i), 0.);
    }
    //   c. interpolate the boundary data on the panel to the target point
    ComplexNumVec weights = 
        Interpolate::compute_barycentric_weights_1d<imaginary>(uv);

    ComplexNumMat density_at_target =
        Interpolate::evaluate_barycentric_interpolant_1d<imaginary>(
                 uv, weights, complex_densities, 
                 ComplexNumVec(1, false, &target_pullback));
    assert(density_at_target.n() == 1); 
    
    double L = patch->characteristic_length();
    CircleArc arc(k*L);
    imaginary curve_deriv_at_target = arc.derivative(target_pullback);

    /*cout << "curve_deriv_at_target: " << curve_deriv_at_target << endl;
    cout << "complex_densities: " << complex_densities << endl;
    cout << "density_at_target(i,0): " << density_at_target(0,0) << endl;
    cout << "target_pullback: " << target_pullback << endl;*/
    //cout << "curvature: " << k<< endl;
    double max_error_estimate = 0.;
    for (int i = 0; i < density_dof; i++) {
        double ith_error_estimate = evaluate_error_estimate_expression(
            abs(curve_deriv_at_target), 
            abs(density_at_target(i,0)), 
            target_pullback, quadrature_order);
        //ith_error_estimate /= closest_point.distance_from_target/L;
        max_error_estimate = max(max_error_estimate, ith_error_estimate);
    }
    return max_error_estimate;

}

double ErrorEstimate::evaluate_error_estimate( FaceMapSubPatch* patch, Point3 target, 
        OnSurfacePoint closest_point, int quadrature_order, 
        DblNumMat uv_values, DblNumMat density_values){
    // 1. evaluate the boundary data at the computed preimage
    //   a. interpolate the boundary data on the panel
    // need to prove that we can get away with interpolating |density| values 
    // instead of just density values. should be ok
    assert(uv_values.n() == quadrature_order*quadrature_order);
    assert(density_values.n() == quadrature_order*quadrature_order);
    // MJM TODO possible bug: currently assuming density values are sampled from
    // the circle over the interval [-1,1]. 
    // targets near patch boundaries will break since this is clearly not the
    // case
    // Need to properly remap [-1,1] to the domain sampled in 
    // sample_patch_along_direction(), or target param value to [-1,1].
    // Testing without this for now to see what explodes
    double max_error = 0.; // fully accurate, error will always be positive 
    
    for(auto dir_type : { CurvatureDirection::FIRST, CurvatureDirection::SECOND } ){
        //cout << closest_point.parametric_coordinates << endl;
        //cout << "compute curvature" << endl;
        Point2 curvature_dir = 
            compute_principal_curvature_direction(closest_point, patch, dir_type);

        //cout << "sample along curv. dir." << endl;
        DblNumMat samples_along_line(3, quadrature_order);
        DblNumMat uv_values_line(2, quadrature_order);
        sample_patch_along_direction(quadrature_order, curvature_dir, 
                closest_point.parametric_coordinates, patch,
                uv_values_line, samples_along_line);

        int density_dof = density_values.m();
        DblNumMat density_values_along_line(density_dof, quadrature_order);
        Interpolate::evaluate_barycentric_interpolant_2d( density_dof,
                uv_values, density_values, quadrature_order, 0, 0,
                uv_values_line, density_values_along_line);

        
        //cout << "eval principal curvatures" << endl;
        /*double k1, k2;
        patch->principal_curvatures(closest_point.parametric_coordinates,k1,k2);
        double k = dir_type == CurvatureDirection::FIRST ? k1 : k2;

        // rescale by patch size
        double L = patch->characteristic_length();
        imaginary target_imag(0., closest_point.distance_from_target/L);
        CircleArc arc(k*L);*/
        //cout << "compute error estimate on arc" << endl;
        double error_estimate = ErrorEstimate::evaluate_error_estimate_on_patch(
                patch, dir_type, target, closest_point, quadrature_order, density_values_along_line);

        max_error = max(max_error, error_estimate);
    }
    return max_error;
}

END_EBI_NAMESPACE
