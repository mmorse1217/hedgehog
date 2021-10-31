#include "../catch.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bie3d/evaluator_near.hpp"
#include "bie3d/evaluator_near_interpolate.hpp"
#include "bie3d/evaluator_qbkix.hpp"
#include "bie3d/solver_utils.hpp"


using namespace hedgehog;

TEST_CASE("EvaluatorNearInterpolate tests", "[eval-near]"){

    SECTION("Test exact solution vs smooth quadrature at collocation points"){

        // Initiialize surface representation
        PatchSurfAnalytic* patch_surf =  new PatchSurfAnalytic("BD3D_", "bd3d_");
        patch_surf->setFromOptions();
        patch_surf->_filename = "../wrl_files/sphere.wrl"; // TODO make this less bad
        patch_surf->setup();

        // Generate refined surface samples
        PatchSamples* refined_patch_samples = new PatchSamples("", "");
        refined_patch_samples->bdry() = (PatchSurf*) patch_surf;
        vector<int> patch_partition(patch_surf->patches().size(), 0);
        refined_patch_samples->patch_partition() = patch_partition;
        refined_patch_samples->setup(true); // true => refined
        
        PatchSamples* patch_samples = new PatchSamples("", "");
        patch_samples->bdry() = (PatchSurf*) patch_surf;
        patch_samples->patch_partition() = patch_partition;
        patch_samples->setup(); // true => refined
        
        // Stokes kernel
        Kernel3d stokes_kernel(321, vector<double>(1,1));

        // Shuffle around collocation points....
        CollocationPatchSamples* refined_collocation_data = new CollocationPatchSamples();
        refined_collocation_data->distribute_collocation_points(
                refined_patch_samples->sample_point_3d_position(),
                refined_patch_samples->sample_as_face_point(),
                refined_patch_samples,
                stokes_kernel.get_sdof(),
                stokes_kernel.get_tdof());

        // Create interpolation directions:
        // interior point normals for sphere
        Vec interpolation_directions;
        VecDuplicate(patch_samples->sample_point_normal(), &interpolation_directions);
        VecCopy(patch_samples->sample_point_normal(), interpolation_directions);
        
        // Flip normals
        const double MINUS_ONE = -1.0;
        VecScale(interpolation_directions, MINUS_ONE);

        // Bump up the number of "interpolation points" to generate to produce a
        // finer sampling along each line
        PetscOptionsSetValue(NULL, "-near_interpolation_num_samples", "8");


        // Generate interpolation points
        Vec interpolation_points = 
            generate_auxiliary_interpolation_points(
                patch_samples->sample_point_3d_position(),
                patch_samples->sample_point_3d_position(),
                stokes_kernel,
                //INTERPOLATE_ACROSS_SURFACE,
                EXTRAPOLATE_ONE_SIDE,
                interpolation_directions, NULL);

        // Load density from a previous Stokes solve 
        // TODO move this file to make independent of external test
        Vec density;
        VecCreateMPI(PETSC_COMM_WORLD,
                stokes_kernel.get_sdof()*patch_samples->global_num_sample_points(),
                PETSC_DETERMINE,
                &density);
        {
            PetscViewer viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD, 
                    "stokes_solved_density_singularities.bin",
                    FILE_MODE_READ,
                    &viewer);
            VecLoad(density, viewer);
            PetscViewerDestroy(&viewer);
        }

        // Compute the refined density
        Vec refined_density = compute_refined_density(density, 
                patch_samples, 
                refined_patch_samples,
                refined_collocation_data, 
                stokes_kernel);
        
        // Evaluate smooth quadrature at interpolation points
        Vec smooth_quad_potentials = evaluate_smooth_quadrature(
                PETSC_COMM_WORLD, 
                refined_patch_samples->sample_point_3d_position(),
                refined_patch_samples->sample_point_normal(),
                interpolation_points,
                stokes_kernel,
                refined_density);
        
        // Evaluate exact solution  at interpolation points
        //Vec exact_potentials =evaluate_solution_zxy(
        Vec exact_potentials = evaluate_singularities_along_basis(
                    PETSC_COMM_WORLD,
                    stokes_kernel,
                    interpolation_points);
        /*
        Vec exact_potentials = evaluate_singularities_along_basis(
                    PETSC_COMM_WORLD,
                    stokes_kernel,
                    interpolation_points);
                    */

        // Get local copies of arrays
        int num_local_targets = patch_samples->global_num_sample_points(); 
        int L = 8;
        DblNumMat smooth_quad_potentials_local(stokes_kernel.get_tdof(), L*num_local_targets, smooth_quad_potentials);
        DblNumMat exact_potentials_local(stokes_kernel.get_tdof(), L*num_local_targets, exact_potentials);
        DblNumMat density_local(stokes_kernel.get_sdof(), num_local_targets, density);
        DblNumMat absolute_error(stokes_kernel.get_tdof(), L*num_local_targets); //PER POINT PER DIMENSION
        DblNumMat relative_error(1, L*num_local_targets); // PER POINT ONLY
        DblNumMat interpolation_points_local(DIM, L*num_local_targets, interpolation_points);
        DblNumMat targets_local(DIM, num_local_targets, patch_samples->sample_point_3d_position());

        double h = patch_samples->spacing();
        
        // Correct evaluated potentials by +/-.5*\phi(x) at each interpolation
        // point to compute the on-surface potential
        
        for(int i = 0; i < num_local_targets; i++){
            for(int d = 0; d < stokes_kernel.get_tdof(); d++){
                // Density at the on-surface target point x
                double density_correction = density_local(d,i);

                // Interior interpolation points
                // Correct potential by -.5*\phi(x)
                for(int l = 0; l  < L; l++){
                    exact_potentials_local(d, i*L + l) += -.5*density_correction;
                    smooth_quad_potentials_local(d, i*L + l) += -.5*density_correction;
                }
                // Exterior interpolation points
                // Correct potential by +.5*\phi(x)
                
            }
        }
        
        for(int i = 0; i < num_local_targets; i++){
            for(int l = 0; l < L; l++){
                double error_at_point = 0;
                double true_value_at_point = 0;
                int index = i*L+l;
                for(int d =0; d < stokes_kernel.get_tdof(); d++){
                    
                    absolute_error(d,index) = exact_potentials_local(d,index) - smooth_quad_potentials_local(d,index);
                    error_at_point += pow(absolute_error(d,index),2);
                    true_value_at_point += pow(exact_potentials_local(d,index),2);

                    
                }
                error_at_point = sqrt(error_at_point);
                true_value_at_point = sqrt(true_value_at_point);
                error_at_point /= true_value_at_point;
                relative_error(0, index) = error_at_point;
                Point3 current_target(targets_local.clmdata(i));
                    CAPTURE(i);
                    CAPTURE(l);
                    CAPTURE(current_target);
                    CHECK(error_at_point <= 1e-6);

            }
        }
        DblNumVec dummy_t(num_local_targets*L);
        write_to_file("interpolation_point_positions.txt", interpolation_points_local, dummy_t);
        write_to_file("interpolation_point_errors.txt", relative_error, dummy_t);
        DblNumVec dummy_t_targets(num_local_targets);
        write_to_file("target_point_positions.txt", targets_local, dummy_t_targets);
        smooth_quad_potentials_local.restore_local_vector();
        exact_potentials_local.restore_local_vector();
        density_local.restore_local_vector();
        targets_local.restore_local_vector();

        VecDestroy(&interpolation_directions);
        VecDestroy(&interpolation_points);
        VecDestroy(&density);
        VecDestroy(&refined_density);
        VecDestroy(&smooth_quad_potentials);
        VecDestroy(&exact_potentials);
        delete refined_collocation_data;
        delete patch_samples;
        delete patch_surf;

    }
}
