#include "evaluator_qbkix_average.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "common/utils.hpp"
#include "common/stats.hpp"
#include "common/vtk_writer.hpp"
BEGIN_EBI_NAMESPACE



int EvaluatorQBKIXAverage::setup(){
    // Use the same setup() as regular extrapolation, just use tro sided
    // extrapolation
    //_expansion_type = EXTRAPOLATE_TWO_SIDE;
    //EvaluatorQBKIX::setup(); 
    //-----------------------------------------------------------------------
    // compute interpolation points in intermediate zone
    _expansion_type = EXTRAPOLATE_TWO_SIDE;
    cout << "before generate" << endl;
    _aux_interpolation_points = construct_check_points();

    cout << "aftergenerate" << endl;
    
    // set up FMM
    //PatchSamples* refined_patch_samples = this->_refined_patch_samples;
    _fmm = unique_ptr<FMM>(
        new PvFMM(_refined_patch_samples->sample_point_3d_position(),
                  _refined_patch_samples->sample_point_normal(),
                  _aux_interpolation_points, this->knl()));

    _collocation_data->distribute_collocation_points(
        _closest_sample_3d_position, _closest_sample_as_face_point,
        this->_patch_samples, this->source_dof(), this->target_dof());
  
    _refined_collocation_data->distribute_collocation_points(
            _refined_patch_samples->sample_point_3d_position(),
            _refined_patch_samples->sample_as_face_point(),
            _refined_patch_samples,
            this->source_dof(),
            this->target_dof());
    /*stats.add_result("total density interp time", 0.);
    
    stats.add_result("total qbx time", 0.);
    stats.add_result("total fmm time", 0.);
    */

    stats.add_result("total qbx time", 0.);
    stats.add_result("total matvec time", 0.);
    stats.add_result("num. matvecs", 0.);
    stats.add_result("total density interp time", 0.);
    stats.add_result("total fmm time", 0.);
    /*
    auto fm = dynamic_cast<PatchSurfFaceMap*>(_patch_samples->bdry());
    write_face_map_patches_to_vtk(DblNumMat(0,0), 
        vector<int>(fm->num_patches(), 0),
        fm, -1, "output/surface_"+to_string(this->mpiRank())+"_");*/
    return 0;
    //-----------------------------------------------------------------------
}

Vec EvaluatorQBKIXAverage::construct_check_points(){

    vector<int> qbkix_point_index;
    int num_check_points = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");
    
    for(int i=0; i < num_check_points; i++)
        qbkix_point_index.push_back(i);
    
    // generate interior check points
    _interior_interpolation_points = 
        generate_interior_qbkix_points(
                _target_3d_position,
                _closest_sample_3d_position,
                qbkix_point_index, 
                _interpolation_directions, 
                _target_far_field,
                _target_interpolant_spacing);

    int n = Petsc::get_vec_local_size(_interior_interpolation_points)/3;
    


        /*generate_auxiliary_interpolation_points(_target_3d_position,
                _closest_sample_3d_position, 
                this->knl(),
                _expansion_type, 
                _interpolation_directions,
                _refined_patch_samples->sample_point_far_field(),
                _refined_patch_samples->sample_point_interpolant_spacing());*/
    
    // flip interpolation directions and  generate exterior check points
    VecScale(_interpolation_directions, -1.);
    _exterior_interpolation_points = 
        generate_interior_qbkix_points(
                _target_3d_position,
                _closest_sample_3d_position,
                qbkix_point_index, 
                _interpolation_directions, 
                _target_far_field,
                _target_interpolant_spacing
                );
    
    
 
    // flip back, to be safe
    VecScale(_interpolation_directions, -1.);

    // Concatenate check point vectors
    Vec check_points =  Petsc::concatenate(_interior_interpolation_points, 
            _exterior_interpolation_points);
    VecDestroy(&_interior_interpolation_points);
    VecDestroy(&_exterior_interpolation_points);
    return check_points;
    
}


int EvaluatorQBKIXAverage::eval(Vec density, Vec potential){

    // total number of check points generated
    // note that num. interior check points == num. exterior check points == 
    // num_check_points/2
    double matvec_start = omp_get_wtime();
    int num_check_points = num_local_points(_aux_interpolation_points);
    int num_check_points_one_side = num_check_points/2;
    int num_targets = Petsc::get_vec_local_size(potential)/target_dof();

    //Vec scaled_density;
    //VecDuplicate(density, &scaled_density);

    /*denscale(source_dof(), 
            this->_patch_samples->sample_point_blend_func_value(),
            density,
            scaled_density);*/
    double start = omp_get_wtime();
    Vec refined_density = EvaluatorQBKIX::compute_refined_density(density);
    //write_general_points_to_vtk(_refined_patch_samples->sample_point_3d_position(), source_dof(), 
            //"density.vtp", refined_density,"output/ref_den_"+to_string(this->mpiRank())+"_");
    stats.result_plus_equals("total density interp time", (omp_get_wtime() - start) );

    Vec check_point_potential;
    Petsc::create_mpi_vec(this->mpiComm(), 
            num_check_points*target_dof(),
            check_point_potential);

    start = omp_get_wtime();
    _fmm->evaluate(refined_density, check_point_potential);
  //((PvFMM*)fmm.get())->evaluate_direct(refined_density, check_point_potential);
    //write_general_points_to_vtk(_aux_interpolation_points, target_dof(), 
    //        "check_potential.vtp", check_point_potential,"output/check_potential_"+to_string(this->mpiRank())+"_");
    
    stats.result_plus_equals("total fmm time", (omp_get_wtime() - start) );
    start = omp_get_wtime();

    DblNumMat expansion_dist_from_boundary_local(1, _target_far_field);
    DblNumMat node_spacing_local (1, _target_interpolant_spacing);
    DblNumMat target_3d_position_local = get_local_vector(DIM, num_targets, _target_3d_position);
    DblNumMat closest_sample_3d_position_local = get_local_vector(DIM, num_targets, _closest_sample_3d_position);

    vector<double> interpolation_nodes = get_interpolation_nodes(_expansion_type);
    int L = interpolation_nodes.size();
    cout << "L:" << L << endl;

    DblNumMat interior_target_potential(target_dof(), num_targets);
    DblNumMat exterior_target_potential(target_dof(), num_targets);

    IS interior_check_index_set;
    IS exterior_check_index_set;
    int rank = this->mpiRank();
    int size= this->mpiSize();
    
    int64_t interior_idx[] = {2*rank};
    int64_t exterior_idx[] = {2*rank+1};

    // Index set for the first 1/2 of the vector
    ISCreateBlock(this->mpiComm(),
            num_check_points_one_side*target_dof(), 
            1,
            interior_idx,
            PETSC_COPY_VALUES,
            &interior_check_index_set);
    // Index set for the second 1/2 of the vector
    ISCreateBlock(this->mpiComm(),
            num_check_points_one_side*target_dof(), 
            1,
            exterior_idx,
            PETSC_COPY_VALUES,
            &exterior_check_index_set);
    
    // get the subvectors containing check point potentials

    Vec interior_check_potential;
    Vec exterior_check_potential;
    VecGetSubVector(check_point_potential,
            interior_check_index_set,
            &interior_check_potential);

    VecGetSubVector(check_point_potential,
            exterior_check_index_set,
            &exterior_check_potential);

    assert(Petsc::get_vec_local_size(interior_check_potential) 
            == num_check_points_one_side*target_dof());
    assert(Petsc::get_vec_local_size(exterior_check_potential) 
            == num_check_points_one_side*target_dof());


    DblNumMat interior_check_potential_local = get_local_vector(
            target_dof(), num_check_points_one_side, interior_check_potential);
    DblNumMat exterior_check_potential_local = get_local_vector(
            target_dof(), num_check_points_one_side, exterior_check_potential);


    // extrapolate interior values
    lagrange_extrapolation_bary(target_3d_position_local,
            closest_sample_3d_position_local,
            interpolation_nodes,
            interior_check_potential_local,
            node_spacing_local,
            expansion_dist_from_boundary_local,
            target_dof(),
            &extrapolation_eval_point_qbkix,
            interior_target_potential);

    // extrapolate exterior values
    lagrange_extrapolation_bary(target_3d_position_local,
            closest_sample_3d_position_local,
            interpolation_nodes,
            exterior_check_potential_local,
            node_spacing_local,
            expansion_dist_from_boundary_local,
            target_dof(),
            &extrapolation_eval_point_qbkix,
            exterior_target_potential);


    interior_check_potential_local.restore_local_vector(); 
    exterior_check_potential_local.restore_local_vector();

    VecRestoreSubVector(check_point_potential,
            interior_check_index_set,
            &interior_check_potential);

    VecRestoreSubVector(check_point_potential,
            exterior_check_index_set,
            &exterior_check_potential);


    DblNumMat potential_local = get_local_vector(target_dof(), num_targets, potential);
    DblNumMat density_local = get_local_vector(source_dof(), num_targets, density);
    // Average interior and exterior potentials
#pragma omp parallel for
    for (int i = 0; i < num_targets; i++) {
        for (int d = 0; d < target_dof(); d++) {
            double interior_potential = interior_target_potential(d,i);
            double exterior_potential = exterior_target_potential(d,i);
            potential_local(d,i) = .5*(interior_potential + exterior_potential);
        }
    }
    stats.result_plus_equals("total qbx time", (omp_get_wtime() - start) );
    stats.result_plus_equals("total matvec time", (omp_get_wtime() - matvec_start) );
    stats.result_plus_equals("num. matvecs", 1 );

    potential_local.restore_local_vector();
    //write_general_points_to_vtk(_patch_samples->sample_point_3d_position(), source_dof(), 
            //"target_potential.vtp", potential,"output/target_pot_"+to_string(this->mpiRank())+"_");


    target_3d_position_local.restore_local_vector();
    closest_sample_3d_position_local.restore_local_vector();
    potential_local.restore_local_vector();
    density_local.restore_local_vector();
    //VecRestoreArray(_target_in_out, &target_in_out_ptr);

    ISDestroy(&interior_check_index_set);
    ISDestroy(&exterior_check_index_set);
    Petsc::destroy_vec(check_point_potential);
    Petsc::destroy_vec(refined_density);
    //Petsc::destroy_vec(scaled_density);

}

END_EBI_NAMESPACE
