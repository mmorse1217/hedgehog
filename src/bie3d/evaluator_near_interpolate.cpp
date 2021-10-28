#include "evaluator_near_interpolate.hpp"

BEGIN_EBI_NAMESPACE

int EvaluatorNearInterpolate::setup(){
    EvaluatorQBKIX::setup();
}

vector<double> EvaluatorNearInterpolate::get_interpolation_nodes(){
    PetscBool flag;
    int64_t near_interpolation_num_samples;
    PetscOptionsGetInt(NULL, "",  "-near_interpolation_num_samples", &near_interpolation_num_samples,  &flag);

    vector<double> interpolation_nodes;
    for(int i = 0; i < near_interpolation_num_samples + 1; i++){
        if(i != near_interpolation_num_samples/2)
            interpolation_nodes.push_back(i);
    }
    return interpolation_nodes;
}


double interpolation_eval_point(double distance_to_closest_sample, double h){
    vector<double> interpolation_nodes = 
        EvaluatorNearInterpolate::get_interpolation_nodes();

    int L = interpolation_nodes.size();
    return ((double) L/2) - distance_to_closest_sample/h;
}

int EvaluatorNearInterpolate::eval(Vec density, Vec potential){

    double h;
    PetscBool err = PETSC_FALSE;
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &err);
    ebiAssert(err);

    Vec refined_density = compute_refined_density(density);
    Vec interpolation_point_potential = compute_interpolation_target_potential(refined_density);


    int num_local_targets = num_local_points(_target_3d_position);

    vector<double> interpolation_nodes = 
        EvaluatorNearInterpolate::get_interpolation_nodes();
    int L = interpolation_nodes.size();


    double* target_3d_position_ptr;
    double* closest_sample_3d_position_ptr;
    double* interpolation_point_potential_ptr;
    double* final_potential_ptr;
    double* density_ptr;
    double* target_in_out_ptr;

    VecGetArray(interpolation_point_potential, &interpolation_point_potential_ptr);
    VecGetArray(_target_3d_position, &target_3d_position_ptr);
    VecGetArray(_closest_sample_3d_position, &closest_sample_3d_position_ptr);
    VecGetArray(potential, &final_potential_ptr);
    VecGetArray(density, &density_ptr);
    VecGetArray(_target_in_out, &target_in_out_ptr);

    DblNumMat target_3d_position_local(DIM, num_local_targets, false, target_3d_position_ptr);
    DblNumMat closest_sample_3d_position_local(DIM, num_local_targets, false, closest_sample_3d_position_ptr);
    DblNumMat final_potential_local(target_dof(), num_local_targets, false, final_potential_ptr);
    DblNumMat interpolation_point_potential_local(target_dof(), L*num_local_targets, false, interpolation_point_potential_ptr);
    DblNumVec target_in_out_local(num_local_targets, false, target_in_out_ptr);
    


    
    DblNumMat density_local(target_dof(), num_local_targets, false, density_ptr);

    //Correct interior interpolation point potentials by density value at
    // on surface target point
    for(int i = 0; i < num_local_targets; i++){
        // Get correction for interpolation points for the ith target
        DblNumVec density_correction(density_local.m(), false, density_local.clmdata(i));
        for(int l = 0; l < L/2; l++){
            for(int d = 0; d < target_dof(); d++){
                    
                interpolation_point_potential_local(d, i*L + l) -= .5*density_correction(d);
                interpolation_point_potential_local(d, i*L + L/2 + l) += .5*density_correction(d);
                
                    
            }
        }
    }
    

    lagrange_extrapolation(target_3d_position_local,
            closest_sample_3d_position_local,
            interpolation_nodes,
            interpolation_point_potential_local,
            h,
            &interpolation_eval_point,
            final_potential_local);
    
    // upsample along interpolation line and dump to file
    

    
    for(int i = 0; i < num_local_targets; i++){
        DblNumVec density_correction(density_local.m(), false, density_local.clmdata(i));
        for(int d = 0; d < target_dof(); d++){
            //if(target_in_out_local(i) == IN){
            if(target_in_out_local(i) == int(INSIDE)-1){
                final_potential_local(d, i) -= .5*density_correction(d);
            }
        }
    }
    
    VecRestoreArray(_target_3d_position, &target_3d_position_ptr);
    VecRestoreArray(_closest_sample_3d_position, &closest_sample_3d_position_ptr);
    VecRestoreArray(potential, &final_potential_ptr);
    VecRestoreArray(interpolation_point_potential, &interpolation_point_potential_ptr);
    VecRestoreArray(density, &density_ptr);

    VecDestroy(&refined_density);
    }
Vec evaluate_smooth_quadrature(MPI_Comm comm, Vec refined_samples, Vec refined_normals, 
        Vec targets, Kernel3d kernel, Vec refined_density){

    int num_local_targets = num_local_points(targets);
    Vec target_potential;
    VecCreateMPI(comm,
            num_local_targets*kernel.get_tdof(),
            PETSC_DETERMINE,
            &target_potential);

    PvFMM* fmm = new PvFMM(refined_samples, refined_normals, targets, kernel);
    fmm->evaluate(refined_density, target_potential);
    return target_potential;
}
END_EBI_NAMESPACE
