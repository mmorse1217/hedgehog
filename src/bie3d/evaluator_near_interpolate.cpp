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
/*
 // MJM TODO can't test this functionality until markgrid is working.
 // Need to evaluate at arbitrary targets, so we need the density correction at
 // the target point to do so.
Vec EvaluatorNearInterpolate::interpolate_density_to_targets(){
    int num_local_targets = num_local_points(_target_3d_position);

    Vec density_at_targets;
    VecCreateMPI(this->mpiComm(),
            num_local_targets*target_dof(),
            PETSC_DETERMINE,
            &density_at_targets);
}
*/

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
    
    /*
     // Dump interpolation points to a file
    Vec empty;
    VecCreateMPI(this->mpiComm(), 0, PETSC_DETERMINE, &empty);
    write_to_text_file(
            "test_interpolation_point_values.txt",
            _aux_interpolation_points,
            L*num_local_targets,
            empty,
            0,
            interpolation_point_potential,
            target_dof(),
            
            empty,
            0,
            empty,
            0,
            empty,
            0
            );*/


    
    // MJM BUG TODO fix for the case of pressure where tdof != sdof
    // MJM TODO BUG FIXME MUST INTERPOLATE DENSITY TO TARGET POINTS
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
                
                    
                //interpolation_point_potential_local(d, i*L + l) -= density_correction(d);
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
    /*
    int upsample_factor = 10;
    double interval_length = double(L); // ensure >= L
    DblNumVec subsampled_eval_points(interval_length*upsample_factor);
    DblNumMat subsampled_potential(target_dof(), interval_length*upsample_factor);


    int index = 1220; // lol
    double* interpolation_directions_ptr;
    VecGetArray(_interpolation_directions, &interpolation_directions_ptr);
    DblNumMat interp_dir_local(DIM, num_local_targets, false, interpolation_directions_ptr);

    Point3 single_target(target_3d_position_local.clmdata(index));
    Point3 single_closest_sample(closest_sample_3d_position_local.clmdata(index));
    Point3 single_interpolation_direction(interp_dir_local.clmdata(index));
    DblNumMat single_interpolation_potentials(target_dof(), L, false, interpolation_point_potential_local.clmdata(L*index));
    
    Vec upsampled_true_potential;
    VecCreateMPI(this->mpiComm(),
            target_dof()*interval_length*upsample_factor,
            PETSC_DETERMINE,
            &upsampled_true_potential);
    
    Vec upsampled_3d_position;
    VecCreateMPI(this->mpiComm(),
            DIM*interval_length*upsample_factor,
            PETSC_DETERMINE,
            &upsampled_3d_position);

    double* upsampled_3d_position_ptr;
    VecGetArray(upsampled_3d_position, &upsampled_3d_position_ptr);
    DblNumMat upsampled_3d_position_local(DIM, 
            interval_length*upsample_factor, 
            false, 
            upsampled_3d_position_ptr);

    // First interpolation point
    single_interpolation_direction = single_interpolation_direction/single_interpolation_direction.length();
    Point3 e = (single_closest_sample + (interval_length-L/2)*h*single_interpolation_direction);
    Point3 np = e; 
    Point3 step_dir =  -single_interpolation_direction*h/double(upsample_factor);
    //Point3 step_dir =  -single_interpolation_direction;
    for(int it=0; it < upsample_factor*interval_length; it++){
        DblNumVec ith_subsampled_potential(target_dof(),false, subsampled_potential.clmdata(it));
        subsampled_eval_points(it) = double(it)/upsample_factor;
        evaluate_lagrange_interpolant(
                interpolation_nodes,
                single_interpolation_potentials,
                h,
                target_dof(),
                double(it)/upsample_factor,
                ith_subsampled_potential);
                

        for(int d =0; d < DIM; d++)
            upsampled_3d_position_local(d,it) = np(d);
        //cout << "np: " << np(0) << ", " << np(1) << ", "<< np(2) << endl;
        np = np + step_dir;
        //np = e + subsampled_eval_points(it)*step_dir;
    }
    VecRestoreArray(upsampled_3d_position, &upsampled_3d_position_ptr);

    Vec exact = evaluate_singularities_along_basis(this->mpiComm(), this->knl(), upsampled_3d_position);
    Vec val = evaluate_smooth_quadrature(this->mpiComm(), 
            _refined_patch_samples->sample_point_3d_position(), 
            _refined_patch_samples->sample_point_normal(), 
            upsampled_3d_position,
            this->knl(), refined_density);
    //Vec val = evaluate_solution_x(this->mpiComm(), this->knl(), upsampled_3d_position);
    */
/*
   lagrange_extrapolation(
           interpolation_nodes,
           single_interpolation_potentials,
           h,
           subsampled_eval_points,
           subsampled_potential);
           */
    /*
    double* val_ptr;
    double*exact_ptr;
    VecGetArray(val, &val_ptr);
    VecGetArray(exact, &exact_ptr);
    DblNumVec density_correction(density_local.m(), false, density_local.clmdata(index));
    
    for(int i=0; i < interval_length*upsample_factor/2; i++){
        for(int d=0; d < target_dof(); d++){
            val_ptr[i*target_dof() + d] -= .5*density_correction(d);
            val_ptr[(L*upsample_factor/2  + i)*target_dof() + d] += .5*density_correction(d);
            exact_ptr[i*target_dof() + d] -= .5*density_correction(d);
            exact_ptr[(L*upsample_factor/2  + i)*target_dof() + d] += .5*density_correction(d);
        }
    }
    VecRestoreArray(exact, &exact_ptr);
    VecRestoreArray(val, &val_ptr);

    DblNumMat exact_local = get_local_vector(target_dof(), interval_length*upsample_factor, exact);
    DblNumMat val_local = get_local_vector(target_dof(), interval_length*upsample_factor, val);
   write_to_file("upsampled_interp_potential.txt", subsampled_potential, subsampled_eval_points);
   write_to_file("upsampled_interp_smooth_quad.txt", val_local, subsampled_eval_points);
   write_to_file("upsampled_interp_actual_potential.txt", exact_local, subsampled_eval_points);
   {
    double* upsampled_3d_position_ptr;
    VecGetArray(upsampled_3d_position, &upsampled_3d_position_ptr);
    DblNumMat upsampled_3d_position_local(DIM, interval_length*upsample_factor, false, upsampled_3d_position_ptr);

   write_to_file("upsampled_interp_point_position.txt", upsampled_3d_position_local, subsampled_eval_points);
    VecRestoreArray(upsampled_3d_position, &upsampled_3d_position_ptr);
   }
    VecRestoreArray(_interpolation_directions, &interpolation_directions_ptr);
    */
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

    PvFMM* fmm = new PvFMM();
    fmm->initialize_fmm(refined_samples, refined_normals, targets, kernel);
    fmm->evaluate(refined_density, target_potential);
    return target_potential;
}
END_EBI_NAMESPACE
