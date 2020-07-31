#include "bie3d/evaluator_qbkix.hpp"
#include "bie3d/evaluator_near.hpp"
#include "common/nummat.hpp"
#include "common/kernel3d.hpp"
#include "../catch.hpp"
#include "vec3t.hpp"

using namespace hedgehog;


TEST_CASE("Interpolation point generation", "[eval-near]"){
    double h;
    PetscBool err;
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &err);
    assert(err);
    cout.precision(16);
    SECTION("Test interpolation point generation for INTERPOLATE_ACROSS_SURFACE"){

        //Laplace kernel
        Kernel3d kernel(111, vector<double>(1, 1.0));

        int num_target_points = 1;
        const double one = 1.;
        const double half = .25;
        

        // Pick a target point, closest point and interpolation direction
        Vec target_point;
        VecCreateMPI(PETSC_COMM_WORLD,
                num_target_points*DIM, 
                PETSC_DETERMINE,
                &target_point);
        

        // Target point at x = (1., 1., 1.)
        VecSet(target_point, one);

        Vec closest_point;
        VecDuplicate(target_point, &closest_point);
        VecCopy(target_point, closest_point);
        
        Vec interpolation_direction;
        VecDuplicate(target_point, &interpolation_direction);
        VecCopy(target_point, interpolation_direction);

        // Interpolate in the direction (.25, .25, .25)
        VecScale(interpolation_direction, half);

        // Compute the interpolation points
        Vec interpolation_points = generate_auxiliary_interpolation_points(
                target_point,
                closest_point,
                kernel,
                INTERPOLATE_ACROSS_SURFACE,
                interpolation_direction,NULL);

        
        vector<double> interpolation_nodes;
        for(int i = 0; i < 8+1; i++){
            if(i != 8/2){
                interpolation_nodes.push_back(i);
            }
        }

        int L = interpolation_nodes.size();
        DblNumMat interpolation_points_local = get_local_vector(DIM, L, interpolation_points);
        
        double* target_point_ptr;
        double* interpolation_direction_ptr;
        
        VecGetArray(target_point, &target_point_ptr);
        VecGetArray(interpolation_direction, &interpolation_direction_ptr);

        Point3 target(target_point_ptr);
        Point3 direction(interpolation_direction_ptr);
        direction = direction/direction.length();

        cout << "Target: " << target << endl;
        cout << "Direction: " << direction << endl;
            cout << "[ ";
        for(int i = 0; i < L; i++){
                Point3 xi = (target + L/2*h*direction ) - direction*h*interpolation_nodes[i];
                //cout << "i: " << i << endl;
                cout << "[ ";
                for(int d= 0; d < DIM; d++){
                    //cout << xi(d) << ", " ;
                    cout << interpolation_points_local(d,i) << ", " ;
                    //cout <<  interpolation_points_local(d,i) << endl;
                    CHECK( (xi(d) == Approx(interpolation_points_local(d,i)) ) );
                }
                cout  << " ]" << endl;
        }
        cout << "]" << endl;
        VecRestoreArray(target_point, &target_point_ptr);
        VecRestoreArray(interpolation_direction, &interpolation_direction_ptr);
        //cout << interpolation_points_local << endl;
        //interpolation_points_local.restore_local_vector();


    }
}
