#include "../catch.hpp"
#include "common/interpolate.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/face_map_subpatch.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bdry3d/p4est_refinement.hpp"
#include <sampling.hpp>
#include "../utils/evaluation_utils.hpp"

using namespace hedgehog;
using Sampling::equispaced;
using Sampling::chebyshev1;
using Sampling::sample_1d;
using Sampling::sample_2d;

double extrapolation_eval_point_alt_2(double distance_to_closest_sample, double h){
    return -h;
}
double eval_poly_series(Point2 s, int k){
    double true_function_value = 0.;
    for(int kk = 0; kk < k; kk++)
        true_function_value += pow(double(kk),2)*pow(s.x(), max(0,kk-3))* pow(s.y(), kk);
    return true_function_value;
}


void setup_coarse_and_fine_samplings(PatchSurfFaceMap*& coarse_face_map, 
        PatchSurfFaceMap*& fine_face_map, 
        PatchSamples*& coarse_samples, 
        PatchSamples*& fine_samples,
        int coarse_refinement_level,
        int fine_refinement_level){
    assert(coarse_refinement_level >=0);
    assert(fine_refinement_level >=0);

    // Coarse grid for interpolation
    coarse_face_map =  new PatchSurfFaceMap("BD3D_", "bd3d_");
    coarse_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    coarse_face_map->_coarse = true;

    coarse_face_map->setFromOptions();
    coarse_face_map->setup();
    
    if(coarse_refinement_level > 0){
        coarse_face_map->refine_uniform(coarse_refinement_level);
    } else {
        coarse_face_map->refine_test();
    }

    cout << "coarse setup" << endl;

    // evaluation grid to interpolate to
    fine_face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    fine_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    fine_face_map->initialize_with_existing_p4est(coarse_face_map->_p4est);
    fine_face_map->_quadrature_weights = coarse_face_map->_quadrature_weights;
    
    if(fine_refinement_level > 0){
        fine_face_map->refine_uniform(fine_refinement_level);
    } else {
        fine_face_map->refine_test();
    }
    cout << "fine setup" << endl;


    int num_coarse_patches = coarse_face_map->patches().size();
    int num_fine_patches = fine_face_map->patches().size();

    // Sample coarse and fine patch sets
    coarse_samples = new PatchSamples("", "");
    vector<int> patch_partition(num_coarse_patches, 0);
    coarse_samples->bdry() = coarse_face_map;
    coarse_samples->patch_partition() = patch_partition;
    coarse_samples->setup();

    fine_samples = new PatchSamples("", "");
    vector<int> patch_partition_fine(num_fine_patches, 0);
    fine_samples->bdry() = fine_face_map;
    fine_samples->patch_partition() = patch_partition_fine;
    fine_samples->setup();

}

void setup_adaptive_samplings(PatchSurfFaceMap*& coarse_face_map, 
        PatchSurfFaceMap*& fine_face_map, 
        PatchSamples*& coarse_samples, 
        PatchSamples*& fine_samples,
        int coarse_refinement_level,
        int fine_refinement_level){
    assert(coarse_refinement_level >=0);
    assert(fine_refinement_level >=0);

    // Coarse grid for interpolation
    coarse_face_map =  new PatchSurfFaceMap("BD3D_", "bd3d_");
    coarse_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    coarse_face_map->_coarse = true;

    coarse_face_map->setFromOptions();
    coarse_face_map->setup();
    
    if(coarse_refinement_level > 0){
        coarse_face_map->refine_uniform(coarse_refinement_level);
    } else {
        coarse_face_map->refine_test();
    }

    cout << "coarse setup" << endl;

    // evaluation grid to interpolate to
    fine_face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    fine_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    fine_face_map->initialize_with_existing_p4est(coarse_face_map->_p4est);
    fine_face_map->_quadrature_weights = coarse_face_map->_quadrature_weights;
    
    resolve_function(fine_face_map->_p4est, fine_face_map, &laplace_singluarity_flat_patch, 1, 1e-5);

    cout << "fine setup" << endl;


    int num_coarse_patches = coarse_face_map->patches().size();
    int num_fine_patches = fine_face_map->patches().size();

    // Sample coarse and fine patch sets
    coarse_samples = new PatchSamples("", "");
    vector<int> patch_partition(num_coarse_patches, 0);
    coarse_samples->bdry() = coarse_face_map;
    coarse_samples->patch_partition() = patch_partition;
    coarse_samples->setup();

    fine_samples = new PatchSamples("", "");
    vector<int> patch_partition_fine(num_fine_patches, 0);
    fine_samples->bdry() = fine_face_map;
    fine_samples->patch_partition() = patch_partition_fine;
    fine_samples->setup();

}


void evaluate_function_on_patch_samples(Vec samples, int range_dim, Vec& function_values){
    // evaluate target function on coarse grid to produce values to
    // interpolate

    int num_func_values = Petsc::get_vec_size(function_values)/range_dim;
    int num_samples = Petsc::get_vec_size(samples)/3;

    assert(num_func_values ==  num_samples);
    DblNumMat func_local = get_local_vector(range_dim,
            num_samples, function_values);

    DblNumMat coarse_sample_positions = 
        get_local_vector(3, num_samples, samples);

    for(int i =0;  i < num_samples; i++){
        Point3 sample(coarse_sample_positions.clmdata(i));
        
        Point3 func_value =  grad_exp_sin_sin_sin(sample);
        for(int d = 0; d < range_dim;d++){
            func_local(d,i) = func_value(d);
        }
    }

    func_local.restore_local_vector();
    coarse_sample_positions.restore_local_vector();

}



TEST_CASE("Test common/interpolate.cpp", "[interpolate][no-fmm][common][critical]"){
    //cout.precision(16);
    SECTION("barycentric_interpolation"){
        int n = 12;
        int num_eval_points = 300;
        double a = 0; 
        double b = 1;
        Rectangle domain(Interval(a,b), Interval(a,b));
        // Interpolate at chebyshev points
        DblNumVec x = DblNumVec(n, true, sample_1d<chebyshev1>(n, domain).data());
        
        // evaluate on equispaced grid for testing
        DblNumVec eval_points = DblNumVec(num_eval_points, true,
            sample_1d<equispaced>(num_eval_points,domain).data());

        
        DblNumVec f(n);
        for(int i = 0; i < n; i++){
            f(i) = sin(x(i));
        }

        DblNumVec weights = Interpolate::compute_barycentric_weights_1d<double>(x);
        DblNumVec evaluated_interpolant = 
            Interpolate::evaluate_barycentric_interpolant_1d<double>(
                    x, weights, f, eval_points);

        for(int i =0; i < eval_points.m(); i++){
            double fi = sin(eval_points(i));
            double error = fabs(evaluated_interpolant(i) - fi);
            error /= fabs(fi);
            CHECK((error <= 2e-14 || (fi <= 1e-15)) );
        }

        // Clenshaw-curtis can integrate polys of degree 2n exactly, so nth
        // order quadrature is enough
        DblNumVec integrals_of_lagrange_polys(n);
        int quad_order = n; 
        DblNumVec jacobian(quad_order);
        setvalue(jacobian, 1.0);
        for(int i = 0; i < n; i++){
            integrals_of_lagrange_polys(i) = 
                Interpolate::integrate_ith_lagrange_basis_func(i, a, b, n, x, quad_order, jacobian(i));
        }

        
        // Approximate integral of sin(x), which is 1-cos(1)
        double computed_integral = dot(f, integrals_of_lagrange_polys); 
        double true_integral = 1. - cos(1.);

        //cout << "Approx integral: " << computed_integral << endl;
        //cout << "True integral:   " << true_integral<< endl;
        //cout <<  fabs(true_integral - computed_integral)/fabs(true_integral) << endl;
        CHECK( fabs(true_integral - computed_integral)/fabs(true_integral) <=1e-14);

        // Approximate integral of sin(3x), which is 2/3*(sin(3/2)^2
        true_integral = 2./3.*(pow(sin(3./2.), 2));
        for(int i = 0; i < n; i++){
            f(i) = sin(3*x(i));
        }

        computed_integral = dot(f, integrals_of_lagrange_polys);
        //cout << "Approx integral: " << computed_integral << endl;
        //cout << "True integral:   " << true_integral<< endl;
        //cout <<  fabs(true_integral - computed_integral)/fabs(true_integral) << endl;
        cout << fabs(true_integral - computed_integral)/fabs(true_integral) << endl;
        CHECK( fabs(true_integral - computed_integral)/fabs(true_integral) <=5e-12);


        // Approximate integral of x^7 + 9*x^5 + - 98*x^4 + 400*x^2 - 80x, which
        // is 9043/120
        true_integral = 9043./120.;
        for(int i = 0; i < n; i++){
            double xi = x(i);
            double x7 = pow(xi, 7);
            double x5 = pow(xi, 5);
            double x4 = pow(xi, 4);
            double x2 = pow(xi, 2);
            f(i) = x7 + 9*x5 - 98*x4 + 400*x2 - 80*xi;
        }

        computed_integral = dot(f, integrals_of_lagrange_polys);;
        //cout <<  fabs(true_integral - computed_integral)/fabs(true_integral) << endl;
        //cout << "Approx integral: " << computed_integral << endl;
        //cout << "True integral:   " << true_integral<< endl;
        CHECK( fabs(true_integral - computed_integral)/fabs(true_integral) <=1e-15);


    }
    SECTION("Ttest barycentric interpolation in 1d"){
        int n = 12;
        int num_eval_points = 20;
        double grid_spacing = .01;
        double a = 0; 
        double b = double(n-1)*grid_spacing;
        double eval_b = b;
        Rectangle domain(Interval(a,b), Interval(a,b));
        Rectangle domain_test(Interval(eval_b,a), Interval(eval_b,a));
        // Interpolate at chebyshev points
        DblNumVec x = DblNumVec(n, true, sample_1d<equispaced>(n, domain).data());
        
        // evaluate on equispaced grid for testing
        DblNumVec eval_points = DblNumVec(num_eval_points, true,
            sample_1d<equispaced>(num_eval_points,domain_test).data());

        DblNumVec f(n);
        for(int i = 0; i < n; i++){
            f(i) = 1./(1.+x(i));
        }

        DblNumMat target_3d(3, 1);
        DblNumMat clst_pts(3, 1);
        DblNumMat evaluated_interpolant(1,1);
        vector<double> int_nodes(n);
        DblNumMat fmat(1,n);
        for(int i =0 ; i < n; i++){
            int_nodes[i] = x(i);
            fmat(0,i) = f(i);
        }
        for(int i =0 ; i < num_eval_points; i++){
            target_3d(0,0) = eval_points(i);
            clst_pts(0,0) = target_3d(0,0);


            DblNumMat expansion_dist_from_target(1,1);
            expansion_dist_from_target(0,0) = -eval_points(i);
            DblNumMat node_spacing(1,1);
            node_spacing(0,0) = 1.;
        lagrange_extrapolation_bary(target_3d, 
                clst_pts, 
                int_nodes,
                fmat,
                node_spacing,
                expansion_dist_from_target,
                1,
                &extrapolation_eval_point_qbkix,
                //extrapolation_eval_point_alt_2,
                evaluated_interpolant);

            double fi = 1./(1.+eval_points(i));
            double error = fabs(evaluated_interpolant(0,0) - fi);
            error /= fabs(fi);
            CHECK((error <= 1e-14) );
        }

    }
    SECTION("Ttest barycentric extrapolation in 1d"){
        int n = 12;
        int num_eval_points = 500;
        double grid_spacing = .01;
        double a = 0; 
        double b = double(n-1)*grid_spacing;
        double eval_b = -grid_spacing*5;
        Rectangle domain(Interval(a,b), Interval(a,b));
        Rectangle domain_test(Interval(eval_b,a), Interval(eval_b,a));
        // Interpolate at chebyshev points
        DblNumVec x = DblNumVec(n, true, sample_1d<equispaced>(n, domain).data());
        
        // evaluate on equispaced grid for testing
        DblNumVec eval_points = DblNumVec(num_eval_points, true,
            sample_1d<equispaced>(num_eval_points,domain_test).data());

        DblNumVec f(n);
        for(int i = 0; i < n; i++){
            f(i) = 1./(1.+x(i));
        }

        DblNumMat target_3d(3, 1);
        DblNumMat clst_pts(3, 1);
        DblNumMat evaluated_interpolant(1,1);
        vector<double> int_nodes(n);
        DblNumMat fmat(1,n);
        for(int i =0 ; i < n; i++){
            int_nodes[i] = x(i);
            fmat(0,i) = f(i);
        }
        for(int i =0 ; i < num_eval_points; i++){
            target_3d(0,0) = eval_points(i);
            clst_pts(0,0) = target_3d(0,0);

            DblNumMat expansion_dist_from_target(1,1);
            expansion_dist_from_target(0,0) = fabs(eval_points(i));
            DblNumMat node_spacing(1,1);
            node_spacing(0,0) = 1.;
        lagrange_extrapolation_bary(target_3d, 
                clst_pts, 
                int_nodes,
                fmat,
                node_spacing,
                expansion_dist_from_target,
                1,
                &extrapolation_eval_point_qbkix,
                //extrapolation_eval_point_alt_2,
                evaluated_interpolant);

            double fi = 1./(1.+eval_points(i));
            double error = fabs(evaluated_interpolant(0,0) - fi);
            error /= fabs(fi);
            CHECK((error <= 1e-6) );
        }
    }
    SECTION("Test barycentric 2d interpolation on a vector-valued function (all components are equal)"){
        double a = 0.;
        double b = 1.;
        Rectangle domain(Interval(a,b), Interval(a,b));

        int num_samples = 10;
        int refinement_factor = 2;
        int range_dim = 3;
        
        DblNumMat sample_points(2, num_samples*num_samples, true,
            sample_2d<chebyshev1>(num_samples, domain).data());
        DblNumMat function_values(range_dim, num_samples*num_samples);
        
        DblNumMat refined_function_values(range_dim, pow(num_samples*refinement_factor, 2));
        DblNumMat refined_sample_points(2, pow(num_samples*refinement_factor, 2), true,
            sample_2d<chebyshev1>(num_samples*refinement_factor, domain).data());

        for(int k = 0; k < num_samples; k++){
            setvalue(function_values, 0.);
            setvalue(refined_function_values, 0.);
            for(int i = 0; i < num_samples; i++){
                for(int j = 0; j < num_samples; j++){
                    int index = num_samples*i+j;
                    Point2 s(sample_points.clmdata(index));
                    for(int d =0;  d< range_dim; d++){
                        function_values(d,index) = eval_poly_series(s,k);
                    }
                }
            }

            Interpolate::evaluate_barycentric_interpolant_2d(
                    range_dim,
                    sample_points, 
                    function_values, 
                    num_samples,
                    refined_sample_points,
                    refined_function_values);

            for(int i = 0; i < num_samples*refinement_factor; i++){
                for(int j = 0; j < num_samples*refinement_factor; j++){
                    int index = num_samples*refinement_factor*i+j;
                    Point2 s(refined_sample_points.clmdata(index));

                    double true_function_value =  eval_poly_series(s,k);
                    cout.precision(16);
                    for(int d =0;  d< range_dim; d++){
                        double diff = fabs(true_function_value - refined_function_values(d,index));
                        CHECK(diff <= 1e-11);
                    }

                }
            }


        }

    }
    SECTION("Test 2d barycentric interpolation of a weird function on a [0,1]^2 to four child samplings, to be used in refinement"){
        double a = 0.;
        double b = 1.;
        Rectangle domain(Interval(a,b), Interval(a,b));

        int num_samples = 20;
        int refinement_factor = 2;
        int range_dim = 2;

        DblNumMat sample_points(2, num_samples*num_samples, true,
                sample_2d<chebyshev1>(num_samples, domain).data());
        DblNumMat function_values(range_dim, num_samples*num_samples);

        DblNumMat refined_function_values(range_dim, 4*num_samples*num_samples);

        // create child nodes to interpolate to
        DblNumMat refined_sample_points(2, 4*num_samples*num_samples);
        for(int i=0; i < 2; i++){
            for(int j=0; j < 2; j++){
                double width = (b-a)/2.;
                double ax_child = a + width*double(i);
                double ay_child = a + width*double(j);
                
                Rectangle domain(Interval(ax_child, ax_child+width), 
                        Interval(ay_child, ay_child+width));

                DblNumMat child_samples(2, num_samples*num_samples, true,
                        sample_2d<chebyshev1>(num_samples, domain).data());

                // copy them into the array to passed to interpolate func
                // TODO kill all these terrible loops
                int stride = 2*i+j;
                for(int sj=0; sj < num_samples; sj++){
                    for(int si=0; si < num_samples; si++){
                        int index = sj*num_samples+si;
                        int block = num_samples*num_samples;
                        for(int d =0; d < 2; d++)
                            refined_sample_points(d, stride*block +index) 
                                = child_samples(d,index);
                    }
                }

            }
        }

        // evaluate original function data 
        for(int i =0; i < sample_points.n(); i++){
            Point2 s(sample_points.clmdata(i));
            Point2 f_s = grad_exp_sin_cos(s);
            for(int d =0; d < 2; d++){
                function_values(d,i) = f_s(d);
            }
        }

        // actually interpolate
        Interpolate::evaluate_barycentric_interpolant_2d(
                range_dim,
                sample_points, 
                function_values, 
                num_samples,
                refined_sample_points,
                refined_function_values);

        // compare error
        for(int i =0; i < refined_sample_points.n(); i++){
            Point2 s(refined_sample_points.clmdata(i));
            Point2 true_f_s = grad_exp_sin_cos(s);
            Point2 computed_f_s(refined_function_values.clmdata(i));
            CHECK((computed_f_s - true_f_s).length() <=1e-14);
        }

    }
    SECTION("test barycentric derivatives: polynomials"){
        int m = 150;
        DblNumVec eval_points(m, true, Sampling::sample_1d<Sampling::equispaced>(m, Sampling::base_domain).data());
        DblNumVec true_interpolant_values(m);
        for (int k = 1; k < 16; k++) {
            
            DblNumVec nodes(k, true, Sampling::sample_1d<Sampling::chebyshev1>(k+1, Sampling::base_domain).data());
            
            DblNumVec function_values(k);
            DblNumVec coeffs(k);
            setvalue(true_interpolant_values, 0.);
            setvalue(function_values, 0.);
            setvalue(coeffs, 1.);

            for (int kk = 0; kk < k; kk++) {
                for (int j = 0; j < k; j++) {
                    function_values(kk) += coeffs(j)*pow(nodes(kk),j);
                }
            }

            DblNumVec weights = Interpolate::compute_barycentric_weights_1d<double>(nodes);
            DblNumVec interpolant_values = 
                Interpolate::evaluate_barycentric_derivative_1d<double>( 
                        nodes, weights, function_values, eval_points);
            
            for (int kk = 0; kk < m; kk++) {
                for (int j = 1; j < k; j++) 
                    true_interpolant_values(kk) += coeffs(j)*j*pow(eval_points(kk),j-1);
            }

            for (int kk = 0; kk < m; kk++) {
                double computed_value = interpolant_values(kk);
                double true_value = true_interpolant_values(kk);
                CHECK(fabs(computed_value - true_value) <=fabs(true_value)*1e-11 +1e-11);
            }
        }
         
    }

}


TEST_CASE("Test patch h-refinement with varied refinement levels","[refine][interpolate][patch_samples]"){
    // Load patch + polynomial from file
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/flat_patch.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".066666");
    int num_samples_per_patch_1d = floor(1./Options::get_double_from_petsc_opts("-bis3d_spacing"))+1;
    int num_samples_per_patch = num_samples_per_patch_1d*num_samples_per_patch_1d;


    PatchSurfFaceMap* coarse_face_map;
    PatchSamples* coarse_samples;
    PatchSurfFaceMap* fine_face_map;
    PatchSamples* fine_samples;

    // initialize coarse/fine surface and patch samplings
    int coarse_ref_level = 1;
    int fine_ref_level = 3;

    setup_coarse_and_fine_samplings(coarse_face_map, fine_face_map, 
            coarse_samples, fine_samples, coarse_ref_level, fine_ref_level);

    int num_coarse_patches = coarse_face_map->patches().size();
    int num_fine_patches = fine_face_map->patches().size();
    // evaluate target function on coarse grid to produce values to
    // interpolate
    int range_dim = 3;
    Vec function;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,
            num_samples_per_patch*range_dim*num_coarse_patches,
            function);
    evaluate_function_on_patch_samples(coarse_samples->sample_point_3d_position(), range_dim, function);

    // interpolate values to refined patches
    Vec refined_function_values = refine_function(range_dim, function, coarse_samples, fine_samples);

    // evaluate true function values on refined patch grids
    Vec true_refined_function_values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,
            num_samples_per_patch*range_dim*num_fine_patches,
            true_refined_function_values);

    evaluate_function_on_patch_samples(fine_samples->sample_point_3d_position(), 
            range_dim, true_refined_function_values);

    // they better be the same size...
    int num_ref_func_values = Petsc::get_vec_size(refined_function_values)/range_dim;
    int num_true_ref_func_values = Petsc::get_vec_size(true_refined_function_values)/range_dim;
    int num_fine_samples= Petsc::get_vec_size(fine_samples->sample_point_3d_position())/3;
    assert(num_true_ref_func_values == num_ref_func_values);
    assert(num_fine_samples == num_ref_func_values);


    DblNumMat refined_function_values_local = get_local_vector(range_dim,
            num_ref_func_values, refined_function_values);

    DblNumMat true_refined_function_values_local = get_local_vector(range_dim,
            num_true_ref_func_values, true_refined_function_values);
    DblNumMat fine_sample_positions_local= get_local_vector(3,
            num_true_ref_func_values, fine_samples->sample_point_3d_position());


    for(int si = 0; si < num_ref_func_values; si++){
        for(int d = 0; d < range_dim; d++){
            double computed_value = refined_function_values_local(d,si);
            double true_value= true_refined_function_values_local(d,si);

            CHECK(fabs(computed_value - true_value) <=fabs(true_value)*1e-12 +1e-12);
        }
    }

}
TEST_CASE("Test patch h-refinement with adaptive refinement","[refine][interpolate][patch_samples][adaptive]"){
    // Load patch + polynomial from file
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/flat_patch.poly");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".066666");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    int num_samples_per_patch_1d = floor(1./Options::get_double_from_petsc_opts("-bis3d_spacing"))+1;
    int num_samples_per_patch = num_samples_per_patch_1d*num_samples_per_patch_1d;


    PatchSurfFaceMap* coarse_face_map;
    PatchSamples* coarse_samples;
    PatchSurfFaceMap* fine_face_map;
    PatchSamples* fine_samples;

    // initialize coarse/fine surface and patch samplings
    int coarse_ref_level = 1;
    int fine_ref_level = 1;

    setup_adaptive_samplings(coarse_face_map, fine_face_map, 
            coarse_samples, fine_samples, coarse_ref_level, fine_ref_level);

    int num_coarse_patches = coarse_face_map->patches().size();
    int num_fine_patches = fine_face_map->patches().size();
    // evaluate target function on coarse grid to produce values to
    // interpolate
    int range_dim = 3;
    Vec function;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,
            num_samples_per_patch*range_dim*num_coarse_patches,
            function);
    evaluate_function_on_patch_samples(coarse_samples->sample_point_3d_position(), range_dim, function);

    // interpolate values to refined patches
    Vec refined_function_values = refine_function(range_dim, function, coarse_samples, fine_samples);

    // evaluate true function values on refined patch grids
    Vec true_refined_function_values;
    Petsc::create_mpi_vec(MPI_COMM_WORLD,
            num_samples_per_patch*range_dim*num_fine_patches,
            true_refined_function_values);

    evaluate_function_on_patch_samples(fine_samples->sample_point_3d_position(), 
            range_dim, true_refined_function_values);

    // they better be the same size...
    int num_ref_func_values = Petsc::get_vec_size(refined_function_values)/range_dim;
    int num_true_ref_func_values = Petsc::get_vec_size(true_refined_function_values)/range_dim;
    int num_fine_samples= Petsc::get_vec_size(fine_samples->sample_point_3d_position())/3;
    assert(num_true_ref_func_values == num_ref_func_values);
    assert(num_fine_samples == num_ref_func_values);


    DblNumMat refined_function_values_local = get_local_vector(range_dim,
            num_ref_func_values, refined_function_values);

    DblNumMat true_refined_function_values_local = get_local_vector(range_dim,
            num_true_ref_func_values, true_refined_function_values);
    DblNumMat fine_sample_positions_local= get_local_vector(3,
            num_true_ref_func_values, fine_samples->sample_point_3d_position());


    for(int si = 0; si < num_ref_func_values; si++){
        for(int d = 0; d < range_dim; d++){
            double computed_value = refined_function_values_local(d,si);
            double true_value= true_refined_function_values_local(d,si);

            CHECK(fabs(computed_value - true_value) <=fabs(true_value)*1e-14 +1e-14);
        }
    }
    /*delete coarse_face_map;
    delete coarse_samples;
    delete fine_face_map;
    delete fine_samples;*/

}
