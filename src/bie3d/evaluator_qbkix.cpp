#include "evaluator_qbkix.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "common/vtk_writer.hpp"
#include "common/stats.hpp"
#include <profile.hpp>
BEGIN_EBI_NAMESPACE
int EvaluatorQBKIX::setup(){
    // compute interpolation points in intermediate zone

    vector<int> qbkix_point_index;
    int num_check_points = Options::get_int_from_petsc_opts("-near_interpolation_num_samples");
    for(int i=0; i < num_check_points; i++)
        qbkix_point_index.push_back(i);
    
    _aux_interpolation_points = 
        generate_interior_qbkix_points(
                _target_3d_position,
                _closest_sample_3d_position,
                qbkix_point_index, 
                _interpolation_directions, 
                _target_far_field,
                _target_interpolant_spacing);
    
    int n = Petsc::get_vec_local_size(_aux_interpolation_points)/3;
    // set up FMM
    _fmm = unique_ptr<FMM>(new PvFMM(
                _refined_patch_samples->sample_point_3d_position(),
                _refined_patch_samples->sample_point_normal(),
                         _aux_interpolation_points,
                         this->knl()));
    
    if(dynamic_cast<PatchSurfFaceMap*>(_patch_samples->bdry()) == NULL){
  _collocation_data->distribute_collocation_points(
                _closest_sample_3d_position,
                _closest_sample_as_face_point,
                this->_patch_samples, 
                this->source_dof(), 
                this->target_dof());
  
    _refined_collocation_data->distribute_collocation_points(
            _refined_patch_samples->sample_point_3d_position(),
            _refined_patch_samples->sample_as_face_point(),
            _refined_patch_samples,
            this->source_dof(),
            this->target_dof());
    }
    stats.add_result("total qbx time", 0.);
    stats.add_result("total matvec time", 0.);
    stats.add_result("num. matvecs", 0.);
    stats.add_result("total density interp time", 0.);
    stats.add_result("total fmm time", 0.);

    return 0;
}
Vec EvaluatorQBKIX::compute_refined_density(Vec density){

    int64_t refinement_factor = Options::get_int_from_petsc_opts("-dnref");

    // Scale density and interpolate to refined surface
    Vec scaled_density;
    VecDuplicate(density, &scaled_density);

    denscale(source_dof(), 
            this->_patch_samples->sample_point_blend_func_value(),
            density,
            scaled_density);

    vector<DblNumMat> refined_datvec;
    Vec refined_density;
    if(dynamic_cast<PatchSurfFaceMap*>(_patch_samples->bdry()) != NULL){
        cout << "face map refinement" << endl;
        refined_density = refine_function(source_dof(), 
                density,
                _patch_samples,
                _refined_patch_samples);


        Vec extra_scaled_density;
        //VecView(refined_density, PETSC_VIEWER_STDOUT_SELF);
        VecDuplicate(refined_density, &extra_scaled_density);

        denscale(source_dof(), 
                this->_refined_patch_samples->sample_point_combined_weight(),
                refined_density,
                extra_scaled_density);
        VecDestroy(&refined_density);
        refined_density = extra_scaled_density;

    } else {
        this->_patch_samples->refine_data(source_dof(),
                refinement_factor,
                scaled_density,
                refined_datvec);

        Vec refined_sample_position = _refined_patch_samples->sample_point_3d_position();
        Vec refined_sample_weight = _refined_patch_samples->sample_point_combined_weight();
        refined_density = interpolate_density_to_refined_grid( 
                this->_patch_samples, 
                refined_sample_position,
                refined_sample_weight,
                refined_datvec,
                _refined_collocation_data,
                source_dof());
    }
    int64_t temp; 
    VecGetLocalSize(refined_density, &temp);

        cout << "refined2" << endl;
    ebiAssert(temp == _refined_patch_samples->local_num_sample_points()*source_dof());
    VecDestroy(&scaled_density);
    return refined_density;
}

Vec compute_refined_density(Vec density, 
        PatchSamples* patch_samples, 
        PatchSamples* refined_patch_samples,
        CollocationPatchSamples* refined_collocation_data,
        Kernel3d kernel){

    int64_t refinement_factor = Options::get_int_from_petsc_opts("-dnref");
    int source_dof = kernel.get_sdof();
    
    // Scale density and interpolate to refined surface
    Vec scaled_density;
    VecDuplicate(density, &scaled_density);

    denscale(source_dof, 
            patch_samples->sample_point_blend_func_value(),
            density,
            scaled_density);

    vector<DblNumMat> refined_datvec;
    
    patch_samples->refine_data(source_dof,
            refinement_factor,
            scaled_density,
            refined_datvec);

    Vec refined_sample_position = refined_patch_samples->sample_point_3d_position();
    Vec refined_sample_weight = refined_patch_samples->sample_point_combined_weight();
    Vec refined_density = interpolate_density_to_refined_grid( 
            patch_samples, 
            refined_sample_position,
            refined_sample_weight,
            refined_datvec,
            refined_collocation_data,
            source_dof);
    VecDestroy(&scaled_density);
    return refined_density;
}



Vec EvaluatorQBKIX::compute_interpolation_target_potential(
        Vec refined_density){
    


    // Evaluate interpolation point potentials via FMM on refined surface 
    int num_local_targets = num_local_points(_aux_interpolation_points);
         
    
    Vec interpolation_point_potential;
    VecCreateMPI(this->mpiComm(),
            num_local_targets*target_dof(),
            PETSC_DETERMINE,
            &interpolation_point_potential);
    cout << "sizes in compute_interpolation_target_potential" << endl;
            cout <<  num_local_targets*target_dof() << ", " << Petsc::get_vec_local_size(refined_density) << endl;
    _fmm->evaluate(refined_density, interpolation_point_potential);
    

    return interpolation_point_potential;
}

double extrapolation_eval_point_blendsurf_legacy(double distance_from_target_to_closest_sample, 
        double node_spacing){
    return -(node_spacing - distance_from_target_to_closest_sample)/node_spacing;
}

double extrapolation_eval_point_blendsurf(double distance_from_target_to_closest_sample, 
        double node_spacing, double exp_distance_to_boundary){
    assert(node_spacing == exp_distance_to_boundary);
    return -(node_spacing - distance_from_target_to_closest_sample)/node_spacing;
}
double extrapolation_eval_point_qbkix(double distance_from_target_to_closest_sample, 
        double node_spacing, double exp_distance_to_boundary){
 
    return -(exp_distance_to_boundary - distance_from_target_to_closest_sample)/node_spacing;
}

double extrapolation_eval_point_blendsurf(double distance_to_closest_sample, double h){
    return -(h - distance_to_closest_sample)/h;
}
double extrapolation_eval_point_qbkix(double distance_to_target, double h){
    return -1.;
    //return -distance_to_target/h;
}
int EvaluatorQBKIX::eval(Vec density, Vec potential){
    double matvec_start = omp_get_wtime();
    //pvfmm::Profile::Tic("QBX Evaluation", &PETSC_COMM_WORLD, true);
    {

    double h;
    PetscBool err = PETSC_FALSE;
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &err);
    ebiAssert(err);
    h = Options::get_double_from_petsc_opts("-boundary_distance_ratio"); // unused


    double start = omp_get_wtime();

    //pvfmm::Profile::Tic("Density interpolation", &PETSC_COMM_WORLD, true);
    Vec refined_density = compute_refined_density(density);
    //pvfmm::Profile::Toc();

    stats.result_plus_equals("total density interp time", (omp_get_wtime() - start) );
    string s= "total fmm";
    start = omp_get_wtime();
    //pvfmm::Profile::Tic("Eval potential at checks", &PETSC_COMM_WORLD, true);
    Vec interpolation_point_potential = compute_interpolation_target_potential(refined_density);
    //pvfmm::Profile::Toc();
    stats.result_plus_equals("total fmm time", (omp_get_wtime() - start) );
    
    
    VecDestroy(&refined_density);
    


    int num_local_targets = num_local_points(_target_3d_position);
    
    //pvfmm::Profile::Tic("Extrap to target", &PETSC_COMM_WORLD, true,1);
    vector<double> interpolation_nodes; 
    int L;
    if(_expansion_type == EXTRAPOLATE_ONE_SIDE_CHEBYSHEV){
         interpolation_nodes = get_chebyshev_nodes();
         L = interpolation_nodes.size();
        for(int i = 0; i < L; i++){
            interpolation_nodes[i] /=h;

        }
    } else {
        interpolation_nodes = get_interpolation_nodes(_expansion_type);
         L = interpolation_nodes.size();
    }
    
    start = omp_get_wtime();

    DblNumMat target_3d_position_local = 
        get_local_vector(DIM, num_local_targets, _target_3d_position);
    DblNumMat closest_sample_3d_position_local = 
        get_local_vector(DIM, num_local_targets, _closest_sample_3d_position);
    DblNumMat final_potential_local = 
        get_local_vector(target_dof(), num_local_targets, potential);
    DblNumMat density_local = 
        get_local_vector(target_dof(), num_local_targets, density);
    DblNumMat interpolation_point_potential_local = 
        get_local_vector(target_dof(), num_local_targets*L, interpolation_point_potential);
    DblNumMat target_in_out_local = 
        get_local_vector(1, num_local_targets, _target_in_out);
    DblNumMat expansion_dist_from_boundary_local = 
        get_local_vector(1, num_local_targets, _target_far_field);
    DblNumMat node_spacing_local = 
        get_local_vector(1, num_local_targets, _target_interpolant_spacing);

    lagrange_extrapolation_bary(target_3d_position_local,
            closest_sample_3d_position_local,
            interpolation_nodes,
            interpolation_point_potential_local,
            node_spacing_local,
            expansion_dist_from_boundary_local,
            target_dof(),
            &extrapolation_eval_point_qbkix,
            final_potential_local);
    
#pragma omp parallel for
    for(int i = 0; i < num_local_targets; i++){
        if(target_in_out_local(0,i) == ON_SURFACE){
            for(int d = 0; d < target_dof(); d++){
                final_potential_local(d,i) -= .5*density_local(d,i);
            }
        }
    }
     
    expansion_dist_from_boundary_local.restore_local_vector(); 
    target_3d_position_local.restore_local_vector();
    closest_sample_3d_position_local.restore_local_vector();
    final_potential_local.restore_local_vector();
    density_local.restore_local_vector();
    interpolation_point_potential_local.restore_local_vector();
    target_in_out_local.restore_local_vector();
    expansion_dist_from_boundary_local.restore_local_vector();
    node_spacing_local.restore_local_vector();
    
    Petsc::destroy_vec(interpolation_point_potential);
    stats.result_plus_equals("total qbx time", (omp_get_wtime() - start) );
    stats.result_plus_equals("total matvec time", (omp_get_wtime() - matvec_start) );
    stats.result_plus_equals("num. matvecs", 1 );
    //pvfmm::Profile::Toc();
    }//pvfmm::Profile::Toc();

    return 0;
}

void EvaluatorQBKIX::lagrange_extrapolation(DblNumMat target_3d_position, 
        DblNumMat closest_sample_3d_position, 
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        double h,
        double (*eval_point)(double, double),
        DblNumMat& final_potential){ // TODO get rid of pointer access
    
    // Generate the vector of 1D t values that we need to evaluate the
    // interpolant at for each 3D target. This is given by eval_point(double, double)
    // which returns this value as a function of target distance from surface
    // and h.
    int num_local_targets = target_3d_position.n();
    DblNumVec eval_points(num_local_targets);
    int L = interpolation_nodes.size();

    for(int i = 0; i < num_local_targets; i++){
        Point3 target(target_3d_position.clmdata(i));
        Point3 point_on_boundary(closest_sample_3d_position.clmdata(i));
        Point3 ray_toward_surface(target - point_on_boundary);
        double distance_to_closest_sample = ray_toward_surface.length();
        
        DblNumVec ith_final_potential(target_dof(), false, final_potential.clmdata(i));
        DblNumMat current_target_interpolation_point_potential
            (target_dof(), L, false, interpolation_point_potential_local.clmdata(L*i));
        //eval_points(i) = eval_point(distance_to_closest_sample, h);
        evaluate_lagrange_interpolant(
                interpolation_nodes, 
                current_target_interpolation_point_potential, 
                h, 
                target_dof(),
                eval_point(distance_to_closest_sample, h),
                ith_final_potential);
    }
}

void EvaluatorQBKIX::lagrange_extrapolation(
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        double h,
        DblNumVec eval_points,
        DblNumMat& final_potential){ // TODO get rid of pointer access
    
    
    int num_local_targets = eval_points.m();
    int L = interpolation_nodes.size();

    for(int i = 0; i < num_local_targets; i++){

        DblNumVec ith_final_potential(target_dof(), false, final_potential.clmdata(i));
        evaluate_lagrange_interpolant(
                interpolation_nodes, 
                interpolation_point_potential_local,
                h, 
                target_dof(),
                eval_points(i),
                ith_final_potential);
    }
}
void evaluate_lagrange_interpolant(
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        double h,
        int target_dof,
        double eval_point,
        DblNumVec& final_potential){ // TODO get rid of pointer access
        int L = interpolation_nodes.size();
        // compute extrapolated value from interpolant 
        DblNumVec interpolation_weights(L);
        compute_interpolation_weights(L,
                &(interpolation_nodes[0]),
                eval_point,
                //-(h - distance_to_closest_sample)/h, // MJM TODO CHECK FOR BUG
                //0,
                interpolation_weights.data());
        //cout << interpolation_weights << endl;   
        for(int ll = 0; ll < L; ll++){
            for(int d = 0; d < target_dof; d++){
                final_potential(d) += interpolation_weights(ll) *
                    interpolation_point_potential_local(d, ll);
                        //current_target_interpolation_point_potential(d,ll);
            }
        }
}


void lagrange_extrapolation_bary(DblNumMat target_3d_position, 
        DblNumMat closest_sample_3d_position, 
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        DblNumMat node_spacing,
        DblNumMat expansion_distance_to_boundary,
        int target_dof,
        double (*eval_point)(double, double, double),
        DblNumMat& final_potential){ 
    
    
    // WARNING extrapolation distance is incorrect in the unequal case
        //assert(Options::get_double_from_petsc_opts("-boundary_distance_ratio") == 
        //Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    // Generate the vector of 1D t values that we need to evaluate the
    // interpolant at for each 3D target. This is given by eval_point(double, double)
    // which returns this value as a function of target distance from surface
    // and h.
    int num_local_targets = target_3d_position.n();
    DblNumVec eval_points(num_local_targets);
    int L = interpolation_nodes.size();
/*
    cout << num_local_targets << ", " << node_spacing.n() << endl;
    cout << num_local_targets << ", " << expansion_distance_to_boundary.n() << endl;
    cout << num_local_targets << ", " << closest_sample_3d_position.n() << endl;
    cout << num_local_targets*L << ", " << interpolation_point_potential_local.n() << endl;
    cout << target_dof << ", " << final_potential.m() << endl;
    cout << target_dof << ", " << interpolation_point_potential_local.m() << endl;
*/

    assert(num_local_targets == node_spacing.n());
    assert(num_local_targets == expansion_distance_to_boundary.n());
    assert(num_local_targets == closest_sample_3d_position.n());
    assert(num_local_targets*L == interpolation_point_potential_local.n());
    assert(target_dof == final_potential.m());
    assert(target_dof == interpolation_point_potential_local.m());
    
    // interpolation nodes are x =0, 1, ..., p-1, where
    // p=Options::get_double_from_petsc_opts("-near_interpolation_num_samples")
    // in 3d space they are scaled by the patch characteristic length L, i.e.
    // check point c_i = -n*L*i, where n in the outward normal at the closest
    // point on-surface.
    DblNumVec interp_nodes(L, false, interpolation_nodes.data());
#pragma omp parallel for
    for(int i = 0; i < num_local_targets; i++){
        Point3 target(target_3d_position.clmdata(i));
        Point3 point_on_boundary(closest_sample_3d_position.clmdata(i));
        Point3 ray_toward_surface(target - point_on_boundary);
        double distance_from_target_to_boundary= ray_toward_surface.length();

        // doing math to avoid passing the patches around everywhere.
        
        DblNumVec ith_final_potential(target_dof, false, final_potential.clmdata(i));
        DblNumMat current_target_interpolation_point_potential
            (target_dof, L, false, interpolation_point_potential_local.clmdata(L*i));
        
        DblNumVec weights = 
        Interpolate::compute_barycentric_weights_1d<double>(interp_nodes);
        
        DblNumVec evaluation_points(1);
        evaluation_points(0) = 
            eval_point(distance_from_target_to_boundary,
                node_spacing(0,i),
                expansion_distance_to_boundary(0,i));

        DblNumMat values = Interpolate::evaluate_barycentric_interpolant_1d(
                interp_nodes,
                weights, 
                current_target_interpolation_point_potential,
                evaluation_points);

        assert(evaluation_points.m() == values.n());
        assert(target_dof == values.m());

        for(int d =0; d < target_dof; d++){
            ith_final_potential(d) = values(d,0);
        }
    }
} 


END_EBI_NAMESPACE
