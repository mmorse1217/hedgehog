#include "../catch.hpp"
#include "common/utils.hpp"
#include "bie3d/solver_utils.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"

using namespace hedgehog;

// create list of factors to scale the on-surface samples by, in such a fashion
// that in/out is known
/*
vector<double> create_scale_factors(){
    vector<double> scale_factors;
    // scale_factor  = .1, .2, ..., .9
    int num_interior = 7;
    for(int i = 4; i < num_interior; i++){
        double scale_factor = double(i)/double(num_interior);
        scale_factors.push_back(scale_factor);
    }
    // scale_factor= .999, .9999, .99999, .999999
    int num_near_int = 5;
    for(int i = 1; i < num_near_int; i++){
        double scale_factor = 1. - pow(10, -(i) );
        scale_factors.push_back(scale_factor);
    }
    // scale_factor = 1.000001, 1.00001, 1.0001, 1.001
    int num_near_ext = 4;
    for(int i = num_near_ext; i > 0; i--){
        double scale_factor = 1. + pow(10., -(i) );
        scale_factors.push_back(scale_factor);
    }
    // scale_factor = 1.025, 1.05, 1.075, 1.1
    int num_exterior= 5;
    for(int i = 1; i < num_exterior ; i++){
        double scale_factor = 1. + double(i)/double(num_exterior);
        scale_factors.push_back(scale_factor);
    }
    return scale_factors;

}

Vec create_scaled_surface_target_points(vector<double> scale_factors, PatchSamples* patch_samples){

    int num_samples = patch_samples->local_num_sample_points();

    // returns: scale_factors.size() parallel copies of the surface, ith copy
    // scaled by scale_factors[i]
    Vec targets;
    VecCreateMPI(
            MPI_COMM_WORLD,
            num_samples*DIM*scale_factors.size(),
            PETSC_DETERMINE,
            &targets);

    DblNumMat targets_local = get_local_vector(DIM, num_samples*scale_factors.size(), targets);
    DblNumMat sample_points_local = get_local_vector(DIM, num_samples, patch_samples->sample_point_3d_position());

    for(size_t si = 0; si < scale_factors.size(); si++){
        double scaling = scale_factors[si];
        for(int i = 0; i < num_samples; i++){
            int point_index = si*num_samples + i;
            for(int d =0; d < DIM; d++){
                targets_local(d, point_index) = scaling*sample_points_local(d,i);
            }
        }
    }
    targets_local.restore_local_vector();
    sample_points_local.restore_local_vector();
    return targets;
}

Vec evaluate_potential_with_constant_density(Vec targets, PatchSamples* patch_samples){

    Options::set_value_petsc_opts("-eqn", "111"); // implies source_dof == 1
    int source_degrees_of_freedom = 1;
    int target_degrees_of_freedom = 1;
    int num_targets = Petsc::get_vec_size(targets)/DIM;
    int num_samples = patch_samples->local_num_sample_points();


    SolverGMRESDoubleLayer* solver = new SolverGMRESDoubleLayer("BIS3D_", "bis3d_");
    solver->bdry() = patch_samples->bdry();
    solver->patch_partition() = patch_samples->patch_partition();
    solver->dom() = Options::get_int_from_petsc_opts("-dom");
    solver->setFromOptions(); 
    solver->setup();
    
    Vec density;
    VecCreateMPI(
            MPI_COMM_WORLD,
            num_samples*source_degrees_of_freedom,
            PETSC_DETERMINE,
            &density);

    VecSet(density, 1.); // set constant density

    Vec potential;
    VecCreateMPI(
            MPI_COMM_WORLD,
            num_targets*target_degrees_of_freedom,
            PETSC_DETERMINE,
            &potential);
    solver->fareval(targets, VAR_U, density, potential);

    VecDestroy(&density);
    return potential;
}*/
TEST_CASE("Test point marking, in/out and near/far", "[point-marking]"){
    
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    Options::set_value_petsc_opts("-bis3d_np", "10");

    // initialize surface
    PatchSurfFaceMap* face_map = new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    // sample surface
    PatchSamples* patch_samples = new PatchSamples("", "" );
    patch_samples->bdry() = face_map;
    vector<int> patch_partition(face_map->patches().size(), 0);
    patch_samples->patch_partition() = patch_partition;
    patch_samples->setup();
/*
    SECTION("Test FMM for far point marking"){
        vector<double> scale_factors = create_scale_factors();
        Vec targets = create_scaled_surface_target_points(scale_factors, patch_samples);
        Vec potential = evaluate_potential_with_constant_density(targets,  patch_samples);
        
        int num_interior_levels = 0;
        int num_exterior_levels = 0;
        for(size_t i = 0; i < scale_factors.size(); i++){
            double scaling = scale_factors[i];
            if(scaling < 1.)
                num_interior_levels++;
            else
                num_exterior_levels++;
        }
        int num_targets = Petsc::get_vec_size(targets)/DIM;
        int num_target_levels = scale_factors.size();
        int num_targets_per_level = num_targets/num_target_levels;

        int target_dof = Petsc::get_vec_size(potential)/num_targets;
        DblNumMat potential_local = get_local_vector(target_dof, num_targets, potential);
        assert(target_dof == 1); // in case I get ambitious or something
        for(int i = 0; i < num_targets_per_level*num_interior_levels; i++){
            for(int d = 0; d < target_dof; d++){
                CHECK(fabs(potential_local(d,i) - 1.) <=1e-3);
            }
        }
        for(int i = 0; i < num_targets_per_level*num_exterior_levels; i++){
            for(int d = 0; d < target_dof; d++){
                int offset = num_targets_per_level*num_interior_levels;
                CHECK(fabs(potential_local(d,offset+i)) <=1e-3);
            }
        }

        
        VecDestroy(&targets); 
        VecDestroy(&potential); 
    }
*/
}
