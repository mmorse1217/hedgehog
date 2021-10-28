#include "error_estimate.hpp"
#include <complex>
#include <sampling.hpp>
BEGIN_EBI_NAMESPACE

imaginary characteristic_function(imaginary z, int quadrature_order){
    imaginary zstar(fabs(real(z)), fabs(imag(z)));
    auto square_root_term = sqrt(zstar*zstar -1.);
    auto square_root_term_qth_pow = pow(zstar + square_root_term, quadrature_order-1);
    return 2./square_root_term_qth_pow*(1./square_root_term_qth_pow + 
            4./pow(double(quadrature_order-1)*square_root_term,3));
}



imaginary singularity(double d, double k){
    return imaginary(0., acosh(k*k*d*d + 2*k*d + 1/(2*k*d)));
}

imaginary evaluate_error_estimate_expression_single_term(imaginary curve_derivative, 
        imaginary density_magnitude, imaginary t, 
        int quadrature_order, double curvature, double distance_from_target, double L){
    double phi = atan2(2./(curvature*curvature), 2.*distance_from_target/curvature);
    imaginary pole = curvature/2.*t + phi; // OFF BY A FACTOR OF TWO
    imaginary value = characteristic_function(t, quadrature_order)*fabs(curve_derivative)*fabs(density_magnitude)/sin(pole);
    return value;
}
double evaluate_error_estimate_expression(double L, double distance_from_target, int quadrature_order,
        double curvature,
        imaginary t, imaginary curve_derivative, imaginary density_magnitude, 
        imaginary t_conj, imaginary curve_derivative_conj, imaginary density_magnitude_conj){
    double norm_z = abs(distance_from_target)/L;
    distance_from_target /= L;
    curvature *= L;

    return 1./(4.*M_PI*M_PI*norm_z)*(
            abs( evaluate_error_estimate_expression_single_term(curve_derivative_conj, 
                density_magnitude_conj, t_conj, quadrature_order, curvature, distance_from_target, L) ) + 
            abs( evaluate_error_estimate_expression_single_term(curve_derivative, 
                density_magnitude, t, quadrature_order, curvature, distance_from_target, L)  )
            );

}


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
    } else {
        return Point2(0.,1.);
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
    bool is_vertical = fabs(u_star) < 1e-12;
    bool is_horizontal = fabs(v_star) < 1e-12;
    if (is_vertical || is_horizontal){ 
        // note only one should be zero for a given direction
        if (is_vertical){
            parameter_domain_intercepts = {
                -v_0/v_star, (1.-v_0)/v_star
            };

        } else {
            // horizontal vector
            parameter_domain_intercepts = {
                -u_0/u_star, (1.-u_0)/u_star
            };
        }
    } else {
        // check all four possible end points to see which lands in [0,1]^2
        // only keep the ones that do
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
        if(parameter_domain_intercepts.size() == 4){
            // solution: choose any two unique points from either end
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

imaginary pullback_exact(double curvature, // unneeded, assumed to be on [-1,1]
        double distance_to_target){
    return singularity(curvature, distance_to_target);
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
    //return pullback_exact(k*L, distance_from_target/L);
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

ComplexNumMat density_to_complex(DblNumMat density_values, int quadrature_order){
    int density_dof = density_values.m();
    ComplexNumMat complex_densities(density_dof, quadrature_order);
    imaginary r(1., 0.);
    for (int i = 0; i < quadrature_order; i++) {
        for (int k = 0; k < density_dof; k++) {
            complex_densities(k,i) = r*density_values(k,i);
        }
    }
    return complex_densities;
}

ComplexNumVec imaginary_chebyshev_sampling(int quadrature_order){

    DblNumVec uv_float(quadrature_order, true, 
            Sampling::sample_1d<Sampling::chebyshev1>(quadrature_order, 
                Rectangle(Interval(0.,1.), Interval(0., 1.))).data() );
            
    //      ii. convert to complex (data type are different size so no memcpy)
    ComplexNumVec uv(quadrature_order);
    for (int i = 0; i < quadrature_order; i++) {
        uv(i) = imaginary(uv_float(i), 0.);
    }
    return uv;
}



double evaluate_error_estimate_expression_legacy(double curve_derivative, 
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


double ErrorEstimate::evaluate_near_zone_distance( FaceMapSubPatch* patch){
    OnSurfacePoint closest_point(patch->V(), DBL_MAX, Point2(.5, .5), UNMARKED,  0);
    closest_point.inside_domain = INSIDE;

    double target_accuracy = Options::get_double_from_petsc_opts("-target_accuracy");
    
    double L = patch->characteristic_length();
    double spacing = Options::get_double_from_petsc_opts("-bis3d_spacing");
    int quadrature_order = int(floor(1./spacing))+1;
    
    DblNumMat density_values(1, quadrature_order*quadrature_order);
    setvalue(density_values, 1.);
    DblNumMat uv_values(2, quadrature_order*quadrature_order, true, 
            Sampling::sample_2d<Sampling::chebyshev2>(quadrature_order, Sampling::base_domain).data());
    auto error = [&patch, &closest_point, &uv_values, &density_values, 
                     quadrature_order, L] (double distance)-> double {
            Point3 sample;
            patch->xy_to_patch_coords(closest_point.parametric_coordinates.array(), 
                PatchSamples::EVAL_VL, sample.array());
            Point3 normal = patch->normal(closest_point.parametric_coordinates.array());
            Point3 target = sample - distance*normal;
            closest_point.distance_from_target = distance;
            return evaluate_error_estimate(patch, target,
                    closest_point, quadrature_order, uv_values, density_values);
        };

    // assumes target_error > 1e-11 and that 1e-11 is achieved two patch lengths
    // away. therefore error(0.) - target_error \approx 1 and 
    // error(2.) - target_error < 0. Root is between, so we apply bisection
    double upper_bound = 2.;
    double lower_bound = 0.;

    int i =0; 
    
    // random initialization to get in the while loop...
    double midpoint_prev = 10.;
    double midpoint = 3.;
    while( i < 30  && fabs(midpoint - midpoint_prev)/fabs(midpoint_prev) >= target_accuracy){
        midpoint_prev = midpoint;
        midpoint = (upper_bound + lower_bound)/2.;
        double error_at_mipoint = error(midpoint);
        if(fabs(error_at_mipoint -target_accuracy)/fabs(target_accuracy) < 1e-2){
            return midpoint;
        } else if(error_at_mipoint - target_accuracy > 0){
            lower_bound = midpoint;
        } else {
            upper_bound = midpoint;
        }
        i++;
    }
    return (lower_bound+ upper_bound)/2.;
    

}
double ErrorEstimate::evaluate_near_zone_distance_legacy(double curvature,
        double density_magnitude, double target_error,
        int quadrature_order, int m, double eps){

    ErrorEstimate::CircleArc c(curvature);
    // objective function to minimize
    auto error = [&c, density_magnitude, target_error, 
                     quadrature_order, m] (double t)-> double {
            imaginary imag_t(0., t);
            double curve_deriv = fabs(c.derivative(imag_t));
            return evaluate_error_estimate_expression_legacy(curve_deriv, 
                    density_magnitude, imag_t, quadrature_order, m);
        };
    // note that derivative of objective function matches the derivative of c
    // because z_0 drops out
    
    // TODO refactor newton out of objective function set up somehow
    // issue is how to wrap function pointer to CircleArc and lambdas
    // consistently.
    // TODO searhc in log space
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
    if(closest_point.distance_from_target >.85/fabs(k) && k < 0){
        // assume point is far if it's near the radius of curvature where the
        // pullback breaks down.
        return 1e-12;
    }
    // TODO replace with exact pullback
    imaginary target_pullback= ErrorEstimate::compute_pullback_under_patch(
            target, closest_point, curvature_direction, patch);
    imaginary target_pullback_conj = conj(target_pullback);

    ComplexNumMat complex_densities = density_to_complex(density_values, quadrature_order);
    int density_dof = complex_densities.m();
    ComplexNumVec uv = imaginary_chebyshev_sampling(quadrature_order);
    /*
    int density_dof = density_values.m();
    ComplexNumMat complex_densities(density_dof, quadrature_order);
    imaginary r(1., 0.);
    for (int i = 0; i < quadrature_order; i++) {
        for (int k = 0; k < density_dof; k++) {
            complex_densities(k,i) = r*density_values(k,i);
        }
    }*/
    //   b. compute param values of density 
    //      i. sample as float
    /*
    DblNumVec uv_float(quadrature_order, true, 
            Sampling::sample_1d<Sampling::chebyshev1>(quadrature_order, 
                Rectangle(Interval(0.,1.), Interval(0., 1.))).data() );
            
    //      ii. convert to complex (data type are different size so no memcpy)
    ComplexNumVec uv(quadrature_order);
    for (int i = 0; i < quadrature_order; i++) {
        uv(i) = imaginary(uv_float(i), 0.);
    }*/
    
    //   c. interpolate the boundary data on the panel to the target point...
    ComplexNumVec weights = 
        Interpolate::compute_barycentric_weights_1d<imaginary>(uv);

    ComplexNumMat density_at_target =
        Interpolate::evaluate_barycentric_interpolant_1d<imaginary>(
                 uv, weights, complex_densities, 
                 ComplexNumVec(1, false, &target_pullback));
    
    // ... and at its conjugate location 
    ComplexNumMat density_at_target_conj =
        Interpolate::evaluate_barycentric_interpolant_1d<imaginary>(
                 uv, weights, complex_densities, 
                 ComplexNumVec(1, false, &target_pullback_conj));
    assert(density_at_target.n() == 1); 
    
    double L = patch->characteristic_length();
    CircleArc arc(k*L);
    // Evaluate curve derivative at target point and at its conjugate
    imaginary curve_deriv_at_target = arc.derivative(target_pullback);
    imaginary curve_deriv_at_target_conj = arc.derivative(target_pullback_conj);

    double max_error_estimate = 0.;
    for (int i = 0; i < density_dof; i++) {
        
        double error = evaluate_error_estimate_expression(L, closest_point.distance_from_target, quadrature_order, k,
                target_pullback, curve_deriv_at_target, density_at_target(i,0), 
                target_pullback_conj, curve_deriv_at_target_conj, density_at_target_conj(i,0));
        max_error_estimate = max(max_error_estimate, error);
    }
    return max_error_estimate;

}

double ErrorEstimate::evaluate_error_estimate( FaceMapSubPatch* patch, Point3 target, 
        OnSurfacePoint closest_point, int quadrature_order, 
        DblNumMat uv_values, DblNumMat density_values){
    assert(uv_values.n() == quadrature_order*quadrature_order);
    assert(density_values.n() == quadrature_order*quadrature_order);
    // MJM TODO possible bug: currently assuming density values are sampled from
    // the circle over the interval [-1,1]. 
    // targets near patch boundaries will break since this is clearly not the
    // case
    // Need to properly remap [-1,1] to the domain sampled in 
    // sample_patch_along_direction(), or target param value to [-1,1].
    // Testing without this for now to see what explodes
    double total_error = 0.; // fully accurate, error will always be positive 
    
    for(auto dir_type : { CurvatureDirection::FIRST, CurvatureDirection::SECOND } ){
        Point2 curvature_dir = 
            compute_principal_curvature_direction(closest_point, patch, dir_type);

        // 1. evaluate the boundary data at the computed preimage
        //   a. interpolate the boundary data on the panel
        // need to prove that we can get away with interpolating |density| values 
        // instead of just density values. should be ok 1/20 MJM i think this
        // doesn't apply anymore
        DblNumMat samples_along_line(3, quadrature_order);
        DblNumMat uv_values_line(2, quadrature_order);
        sample_patch_along_direction(quadrature_order, curvature_dir, 
                closest_point.parametric_coordinates, patch,
                uv_values_line, samples_along_line);

        int density_dof = density_values.m();
        DblNumMat density_values_along_line(density_dof, quadrature_order);
        Interpolate::evaluate_barycentric_interpolant_2d( density_dof,
                uv_values, density_values, quadrature_order, 
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
        total_error += error_estimate;
    }
    return total_error;
}

END_EBI_NAMESPACE
