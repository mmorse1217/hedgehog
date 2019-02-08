#include "bie3d/evaluator_qbkix.hpp"
#include "bie3d/evaluator_near.hpp"
#include "common/nummat.hpp"
#include "../catch.hpp"

using namespace Ebi;

double degree_less_than_L_poly(double x){
    return pow(x,7) + 3*pow(x,4) + pow(x,2) + 2*x * 2;
}

TEST_CASE("Lagrange interpolation", "[lagrange]"){


    cout.precision(16);
    // Generate interpolation nodes 0, ..., L
    vector<double> interpolation_nodes = get_interpolation_nodes(INTERPOLATE_ACROSS_SURFACE);
    int L = interpolation_nodes.size();

    // spacing between interpolation nodes
    double h = .1;

    // Interpolating a 1D function f : R^1 \to R^1
    int target_dof = 1;

    // Samples of function to interpolate 
    // We're going to interpolate pow(x,7)
    
    // function values are in the interval [0, L*h], while points are
    // parameterized from [0,L]
    DblNumMat values_at_interpolation_nodes(target_dof, L);
    for(int i = 0; i < L; i++){
        // function in [0, L*h]
        double x = interpolation_nodes[i]; 
        values_at_interpolation_nodes(0,i) = pow(x,7);
    }

    SECTION("Interpolant values at nodes exactly equal original function samples"){
        // 1D points in the interval [0, L*h]
        // We're checking to make sure that p_L(x_i) = f_i, for interpolation
        // nodes x_i and function samples f_i at nodes x_i, where p_L(x) is the
        // degree Lagrange interpolant to f
        DblNumVec evaluation_points(L);
        DblNumMat interpolated_function_values(target_dof, L);
        for(int i = 0; i < L; i++){
            // evaluate interpolant at points in [0,L]
            double x = interpolation_nodes[i];
            evaluation_points(i) = x;


            DblNumVec ith_function_value(target_dof, false,
                    interpolated_function_values.clmdata(i));

            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    values_at_interpolation_nodes,
                    h, 
                    target_dof,
                    evaluation_points(i),
                    ith_function_value);

        }
        double tol = 1e-15;
        for(int i = 0; i < L; i++){
            // For each interpolation node, check that 
            // |p_L(x_i) - f_i| \approx 0
            double x = evaluation_points(i);
            CHECK( (fabs(interpolated_function_values(0,i) - values_at_interpolation_nodes(0,i)) 
                    == Approx(tol))
                 );
        }

    }
    SECTION("Interpolant reasonably approximates the underlying function"){
        // Test the interpolant for points x_i' in between the interpolation
        // nodes.
        
        int num_points_to_test = L*2;
        DblNumVec evaluation_points(num_points_to_test);

        DblNumMat interpolated_function_values(target_dof, num_points_to_test);
        for(int i = 0; i < num_points_to_test; i++){
            // points in [0,L]
            double x = double(i)/2;
            evaluation_points(i) = x;

            DblNumVec ith_function_value(target_dof, false, 
                    interpolated_function_values.clmdata(i));
            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    values_at_interpolation_nodes,
                    h, 
                    target_dof,
                    evaluation_points(i),
                    ith_function_value);
        }

        double tol = 1e-1;
        for(int i = 0; i < num_points_to_test; i++){
            // For each interpolation node, check that 
            // |p_L(x_i) - f_i| \approx 0
            double x = evaluation_points(i);
            // Remember: evaluate original function in interval [0, L*h]
            CHECK( (interpolated_function_values(0,i) == Approx(pow(x,7)).epsilon(tol)));
        }
    }
    
    SECTION("Check the the error decreases for a higher order polynomial approximation"){
        // WARNING test might be bogus. lower order is interpolating on the
        // interval [0,L_1], and the high order poly is interpolating on the
        // interval [0,L_2] with L_1 < L_2. Moreover, the first L_1 nodes of the
        // higher order interpolant coincide with the the nodes of the low order
        // one. This give misleading accuracy.

        // Interpolate
        int num_points_to_test = L*10;
        DblNumVec evaluation_points(num_points_to_test);

        for(int i = 0; i < num_points_to_test; i++){
            // points in [0,L]
            double x = double(i)/10;
            evaluation_points(i) = x;
        }
        
        DblNumMat low_order_interpolated_function_values(target_dof, num_points_to_test);
        
        for(int i = 0; i < num_points_to_test; i++){
            DblNumVec ith_function_value(target_dof, false,
                    low_order_interpolated_function_values.clmdata(i));

            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    values_at_interpolation_nodes,
                    h, 
                    target_dof,
                    evaluation_points(i),
                    ith_function_value);
        }
        

        vector<double> interpolation_nodes;
        for(int i = 0; i < 2*L; i++){
            interpolation_nodes.push_back(double(i)/2.);
        }

        int L = interpolation_nodes.size();

        DblNumMat values_at_interpolation_nodes(target_dof, L);
        for(int i = 0; i < L; i++){
            // function in [0, L*h]
            double x = interpolation_nodes[i]; 
            values_at_interpolation_nodes(0,i) = pow(x,7);
        }

        DblNumMat high_order_interpolated_function_values(target_dof, num_points_to_test);
        for(int i = 0; i < num_points_to_test; i++){
            DblNumVec ith_function_value(target_dof, false,
                    high_order_interpolated_function_values.clmdata(i));
            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    values_at_interpolation_nodes,
                    h, 
                    target_dof,
                    evaluation_points(i),
                    ith_function_value);
        }

        for(int i = 0; i < num_points_to_test; i++){
            // For each interpolation node, check that 
            // |p_L(x_i) - f_i| \approx 0
            double x = evaluation_points(i);
            // Should be roughly the same 
            CHECK( (fabs(low_order_interpolated_function_values(0,i) - pow(x,7)) ==
                    Approx(fabs(high_order_interpolated_function_values(0,i) - pow(x,7))))
                 );
        }

    }
    SECTION("Test interpolating analytic scalar function along a 1D line"){

        // Evaluate f(x,y,z) = xyz + 1 at points along a line, interpolate the
        // potential, and check error
        vector<double> interpolation_nodes = get_interpolation_nodes(INTERPOLATE_ACROSS_SURFACE);
        int L = interpolation_nodes.size();

        // Form line parametrization p = e + t*dir (arbitrarily chosen line)
        Point3 e(1.,0., 1.);
        Point3 dir(0.,10., 0.);

        // parametrized t-values
        DblNumVec interp_t_values(L);
        DblNumMat interpolation_points(DIM, L);
        DblNumMat interpolation_point_potential(target_dof, L);
        
        for(int i = 0; i < L; i++){
            // Store the t-values
            double t = i*h;
            interp_t_values(i) = t;

            // Compute the 3D position on the line
            Point3 xi = e + t*dir;

            for(int d = 0; d < DIM; d++)
                interpolation_points(d,i) = xi(d);

            // compute the potential at that 3D point
            for(int d = 0; d < target_dof; d++){
                double x = xi(0);
                double y = xi(1);
                double z = xi(2);
                interpolation_point_potential(d,i) = x*y*z + 1;
            }
        }
        int num_points_to_test = L*5;

        DblNumMat interpolated_potential(target_dof, num_points_to_test);
        DblNumVec evaluation_points(num_points_to_test);

        for(int i = 0; i < num_points_to_test; i++){
            DblNumVec ith_interpolated_potential(target_dof, false, interpolated_potential.clmdata(i));
            evaluation_points(i) = i*h/5.;

            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    interpolation_point_potential,
                    h,
                    target_dof,
                    evaluation_points(i)/h,
                    ith_interpolated_potential);
        }
        for(int i = 0; i < num_points_to_test; i++){
            Point3 xi = e + evaluation_points(i)*dir; 
            
            double x = xi(0);
            double y = xi(1);
            double z = xi(2);

            CHECK( (interpolated_potential(0,i) == Approx(x*y*z + 1.)));
        }



    }
    SECTION("Test interpolating analytic vector function along a line in 3D"){
        int target_dof = 3;
        vector<double> interpolation_nodes = get_interpolation_nodes(INTERPOLATE_ACROSS_SURFACE);
        int L = interpolation_nodes.size();

        // Form line parametrization p = e + t*dir (arbitrarily chosen line)
        Point3 e(1.,0., 1.);
        Point3 dir(0.,1., 2.);

        // parametrized t-values
        DblNumVec interp_t_values(L);
        DblNumMat interpolation_points(DIM, L);
        DblNumMat interpolation_point_potential(target_dof, L);
        for(int i = 0; i < L; i++){
            // Store the t-values
            double t = i*h;
            interp_t_values(i) = t;

            // Compute the 3D position on the line
            Point3 xi = e + t*dir;

            for(int d = 0; d < DIM; d++)
                interpolation_points(d,i) = xi(d);

            // compute the potential at that 3D point
            for(int d = 0; d < target_dof; d++){
                double x = xi(0);
                double y = xi(1);
                double z = xi(2);
                
                if(d == 0){
                    interpolation_point_potential(d,i) = pow(x,2);
                } else if (d == 1){
                    interpolation_point_potential(d,i) = pow(y,3);
                } else { 
                    interpolation_point_potential(d,i) = pow(z,4);
                }
            }
        }
        int num_points_to_test = L*5;

        DblNumMat interpolated_potential(target_dof, num_points_to_test);
        DblNumVec evaluation_points(num_points_to_test);

        for(int i = 0; i < num_points_to_test; i++){
            DblNumVec ith_interpolated_potential(target_dof, false, interpolated_potential.clmdata(i));
            evaluation_points(i) = i*h/5.;

            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    interpolation_point_potential,
                    h,
                    target_dof,
                    evaluation_points(i)/h,
                    ith_interpolated_potential);
        }
        for(int i = 0; i < num_points_to_test; i++){
            Point3 xi = e + evaluation_points(i)*dir; 
            
            double x = xi(0);
            double y = xi(1);
            double z = xi(2);
            for(int d = 0; d < target_dof; d++){
                // Check that |(x^2,y^3,z^4) - p_L(x,y,z)| <= \epsilon
                if(d == 0){
                    CHECK( (interpolated_potential(d,i) == Approx(pow(x,2))));
                } else if(d == 1){
                    CHECK( (interpolated_potential(d,i) == Approx(pow(y,3))));
                } else {
                    CHECK( (interpolated_potential(d,i) == Approx(pow(z,4))));

                }
                //
            }

            // Check that the accuracy of the interpolant at the interpolation 
            // nodes is extremely high
            if(i % 5 == 0) {
                for(int d = 0; d < target_dof; d++){
                    CHECK( fabs(interpolation_point_potential(d,i/5) - interpolated_potential(d,i)) <=1e-12);
                }
            }
        }



    }

    SECTION("Test interpolating Laplace potential along a line in 3D"){

        // Evaluate Laplace potential at points along a line induced by
        // singularities, interpolate the result, and check error
        int target_dof = 1;
        int upsample_factor = 5;
        vector<double> interpolation_nodes = get_interpolation_nodes(INTERPOLATE_ACROSS_SURFACE);
        int L = interpolation_nodes.size();

        // Form line parametrization p = e + t*dir (arbitrarily chosen line)
        Point3 e(0.,0., 1.);
        Point3 dir(0.,0.,-1.);

        // parametrized t-values
        DblNumVec interp_t_values(L);
        
        Vec interpolation_node_positions;
        VecCreateMPI(PETSC_COMM_WORLD,
                L*DIM,
                PETSC_DETERMINE,
                &interpolation_node_positions);
        
        DblNumMat interpolation_node_positions_local(DIM, L, interpolation_node_positions);

        for(int i = 0; i < L; i++){
            // Store the t-values
            double t = i*h;
            interp_t_values(i) = t;

            // Compute the 3D position on the line
            Point3 xi = e + t*dir;

            // Store the position of interpolation nodes
            for(int d = 0; d < DIM; d++)
                interpolation_node_positions_local(d,i) = xi(d);

        }
        interpolation_node_positions_local.restore_local_vector();
        
        // compute the potential at the target points 
        int num_points_to_test = L*upsample_factor;
        Kernel3d kernel(121, vector<double>(1,1));
        Vec interpolation_node_potentials = 
            evaluate_singularities_along_basis(
                                PETSC_COMM_WORLD, 
                                kernel, 
                                interpolation_node_positions);

        // Get local copy
        DblNumMat interpolation_node_potentials_local(target_dof, L, interpolation_node_potentials);

        DblNumMat interpolated_potential(target_dof, num_points_to_test);
        DblNumVec evaluation_points(num_points_to_test);
        

        // Store 3d positions of upsampled points on line to pass to FMM to
        // compute true potential 
        Vec upsampled_position;
        VecCreateMPI(PETSC_COMM_WORLD,
                L*upsample_factor*DIM,
                PETSC_DETERMINE,
                &upsampled_position);
        DblNumMat upsampled_position_local(DIM, L*upsample_factor, upsampled_position);

        for(int i = 0; i < num_points_to_test; i++){
            // Load where the solution should go 
            DblNumVec ith_interpolated_potential(target_dof, false, interpolated_potential.clmdata(i));
            
            // Compute the current t-value
            evaluation_points(i) = i*h/double(upsample_factor);
            
            //Store it
            Point3 xi = e + evaluation_points(i)*dir;
            for(int d =0; d < DIM; d++)
                upsampled_position_local(d,i) = xi(d);

            // Evaluate the interpolant
            evaluate_lagrange_interpolant(
                    interpolation_nodes,
                    interpolation_node_potentials_local,
                    h,
                    target_dof,
                    evaluation_points(i)/h,
                    ith_interpolated_potential);
        }
        upsampled_position_local.restore_local_vector();

        // Compute the true potential at the upsampled points
        Vec upsampled_potential = 
            evaluate_singularities_along_basis(
                                PETSC_COMM_WORLD, 
                                kernel, 
                                upsampled_position);

        DblNumMat upsampled_potential_local(target_dof, L*upsample_factor, upsampled_potential);

        // Compare the interpolated and true values to ensure accuracy
        for(int i = 0; i < num_points_to_test; i++){
            Point3 xi = e + evaluation_points(i)*dir; 
            
            double x = xi(0);
            double y = xi(1);
            double z = xi(2);

            // Check that the accuracy of the interpolant at the interpolation 
            // nodes is extremely high
                for(int d = 0; d < target_dof; d++){
                    CHECK( fabs(upsampled_potential_local(d,i) - interpolated_potential(d,i)) <=1e-8);
                }
        }



    }
        

}
