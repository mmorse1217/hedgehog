#include "evaluator_near.hpp"
#include "solver_utils.hpp"
#include "common/utils.hpp"
BEGIN_EBI_NAMESPACE

using std::cerr;

// ---------------------------------------------------------------------- 

int EvaluatorNear::setFromOptions()
{
  ebiFunctionBegin;

  PetscBool flg = PETSC_FALSE;
  PetscOptionsGetReal(NULL, this->prefix().c_str(),  "-alfcoef", &_alfcoef,  &flg);
  if(!flg){
      _alfcoef = 1.;
  }
  
  PetscOptionsGetInt(NULL, "", "-dnref", &_refinement_factor, &flg);
  if(flg){
    _refinement_factor = 16;
  }

  PetscOptionsGetInt(NULL, "", "-LL", &_surface_interpolation_num_samples, &flg);
  ebiAssert(flg);
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// The setup function
// --  separates the input array of target points for different zones
// (Omega_1, Omega_2 on p. 260 of the paper)
// --  for each point in Omega_2 zone (near) generates L points in the
// Omega_1/2 zones along the line connecting it to the closest surface point
// -- initializes FMM
// -- creates and calls setup on on-the-surface and jump evaluators
//

int EvaluatorNear::setup()
{
  ebiFunctionBegin;

  iA(_alfcoef > 0);
    
  int face_point_size_in_doubles = 
      this->_patch_samples->face_point_size_in_doubles();
  int num_local_near_targets = num_local_points(_target_3d_position);
  //1. separate points into two classes, near/on, and "far" which here
  //means Omega_1 zone in the paper (Omega_2 points are not passed to this
  //function
  // (really three, but the difference between on-the-surface and near points is minor);
  
  Vec aux_interpolation_points = generate_auxiliary_interpolation_points(
                                    _target_3d_position, 
                                    _closest_sample_3d_position, 
                                    this->knl());
  // create the large vector of intermdediate targets to pass to FMM
  // This consists of interpolation points + the intermediate targets that are
  // explicitly passed by the user for evaluation
  _target_positions_intermediate = 
      Petsc::concatenate(aux_interpolation_points, _intermediate_target_3d_position);


  // REFACTOR NOTE
  // Points are sorted externally into target_positions_intermediate (Omega_2)  
  // and target_3d_positions (Omega_1) if there are N near targets in
  // target_3d_positions, then the first N*(L-1) points of
  // target_positions_intermediate will be the auxiliary interpolation points,
  // by convention.

  //2. Setup  FMM with refined sources, initialize non-distributed arrays
  // @BUG this probably does not work in parallel
  // the "far" points passed to near are assumed to be zone Omega_1
  // (p. 260) i.e. within sqrt(h), so require finer sampling for evaluation
  PatchSamples* refined_patch_samples = this->_refined_patch_samples;
  if(_fmm_initialize){
      fmm = new PvFMM();
      fmm->initialize_fmm(refined_patch_samples->sample_point_3d_position(),
                             refined_patch_samples->sample_point_normal(),
                             _target_positions_intermediate,
                             this->knl());
  }
  
  //3. compute collocation point-related distributed vectors, using refined samples
  // Initialize regular collocation needed in jump_evaluation(...)
 
  _collocation_data->distribute_collocation_points(
                _closest_sample_3d_position,
                _closest_sample_as_face_point,
                this->_patch_samples, 
                this->source_dof(), 
                this->target_dof());
  

  // Initialize refined collocation points for near
  iC(_refined_collocation_data->distribute_collocation_points(
              refined_patch_samples->sample_point_3d_position(), 
              refined_patch_samples->sample_as_face_point(),
              _refined_patch_samples, this->source_dof(), this->target_dof()) );
  
  //4. for near points, create be3dovron and be3dovjmp; in evaluation, we
  //need the jumps for discontinuous potentials, to compute correct value for
  //interpolation on each side (see eval())
  
  _on_surface_evaluator = new EvaluatorOnSurface(this->name()+"RE_", this->prefix()+"re_");
  _on_surface_evaluator->knl() = (this->knl());
  _on_surface_evaluator->target_3d_position() = _closest_sample_3d_position;
  _on_surface_evaluator->target_as_face_point() = _closest_sample_as_face_point;
  _on_surface_evaluator->set_surface_discretization(this->_patch_samples);
 
  iC( _on_surface_evaluator->setFromOptions() );
  iC( _on_surface_evaluator->setup_no_fmm() );
  
  cout<<"nea setup re+je  "<<endl;
  
  
  ebiFunctionReturn(0);
}

// --------------------------------------------------------------------------- 
vector<double> get_interpolation_nodes(ExpansionType expansion_type){
  PetscBool flag;
  int64_t near_interpolation_num_samples;
  PetscOptionsGetInt(NULL, "",  "-near_interpolation_num_samples", &near_interpolation_num_samples,  &flag);

  vector<double> interpolation_nodes(near_interpolation_num_samples);
  
  interpolation_nodes[0] = 0.0;
  interpolation_nodes[1] = 1.0;

  for (int nn = 2; nn < near_interpolation_num_samples; nn++){
	 interpolation_nodes[nn] = interpolation_nodes[1] + (double)(nn-1)/1.0;
  }
  return interpolation_nodes;
}
// --------------------------------------------------------------------------- 
vector<double> get_chebyshev_nodes(){
  PetscBool flag;
  int64_t near_interpolation_num_samples;
  PetscOptionsGetInt(NULL, "",  "-near_interpolation_num_samples", &near_interpolation_num_samples,  &flag);
  
  double h;
  PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &flag);
  ebiAssert(flag);

  vector<double> interpolation_nodes(near_interpolation_num_samples+1);
  double L = double(near_interpolation_num_samples);
  for(int i = 0; i < near_interpolation_num_samples+1; i++){
      double t = -cos(double(i)*M_PI/L);
      t = (L*h)/2.*(t  + 1.);
      interpolation_nodes[i] = t;
  }
  return interpolation_nodes;

}
// --------------------------------------------------------------------------- 
Vec generate_auxiliary_interpolation_points(Vec near_targets, 
        Vec closest_surface_points_near, Kernel3d kernel, ExpansionType expansion_type,
        Vec interpolation_directions, 
        // TODO merge into one vector because they're generally equal
        Vec sample_point_far_field,
        Vec sample_point_interpolant_spacing){
    PetscBool err = PETSC_FALSE;
    double h;
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &err);
    ebiAssert(err);
    vector<double> interpolation_nodes;
    if(expansion_type == EXTRAPOLATE_ONE_SIDE_CHEBYSHEV){
         interpolation_nodes = get_chebyshev_nodes();
    } else {
        interpolation_nodes = get_interpolation_nodes(expansion_type);
    }
    int64_t num_aux_targets = interpolation_nodes.size();
    
    int num_local_targets = num_local_points(near_targets);
    
    // if we're extrapolating, all L points are generated here
    // if we're not, one point on-surface is generated externally, 
    //      so only generate (L-1) points here.
    int external_interp_points = 
        expansion_type != INTERPOLATE ? 0 : 1;
    Vec aux_targets;
    VecCreateMPI(PETSC_COMM_WORLD,
            num_local_targets*DIM*(num_aux_targets-external_interp_points), 
            PETSC_DETERMINE, 
            &aux_targets);

    double* near_trg_ptr;
    double* closest_surf_ptr;
    double* aux_trg_ptr;
    double* interpolation_directions_ptr;
    double* far_field_ptr;
    double* interpolant_spacing_ptr;
    
    VecGetArray(near_targets, &near_trg_ptr);
    VecGetArray(closest_surface_points_near, &closest_surf_ptr);
    VecGetArray(aux_targets, &aux_trg_ptr);

    if(sample_point_far_field != NULL){
        VecGetArray(sample_point_far_field, &far_field_ptr);
    }
    if(sample_point_far_field != NULL){
        VecGetArray(sample_point_interpolant_spacing, &interpolant_spacing_ptr);
    }
    DblNumVec far_field_local(num_local_targets, false, far_field_ptr);
    DblNumVec interpolant_spacing_local(num_local_targets, false, interpolant_spacing_ptr);

    DblNumMat interpolation_directions_local;
    if(interpolation_directions != NULL){
        VecGetArray(interpolation_directions, &interpolation_directions_ptr);
        interpolation_directions_local = DblNumMat(DIM, num_local_targets, 
                false, interpolation_directions_ptr);
    }

  // generate L-1 points in the far zone along the normal at each target point 
  // in the near
  // zone (where we interpolate); the spacing is determined by interpolation_nodes
  // constants*spacing*alpha_coeff  (which is generally 1)
  // add these points to the far list
  
  // This is important to note: this ONLY generates the L-1 interpolation
  // points in the intermediate zone that are passed to FMM. The Lth point is
  // passed explicitly as closest_sample_3d_position, evaluated on-surface via
  // singularity cancellation, and together these L points are used to
  // interpolate.
  int aux_points_counter = 0;
  int one_sided_point_counter = 0;

  for(int ti=0; ti<num_local_targets; ti++) {
        double distance_from_boundary;
        double interpolant_spacing;
        if(sample_point_far_field != NULL){
            distance_from_boundary = far_field_local(ti);
        }
        if(sample_point_interpolant_spacing!= NULL){
            interpolant_spacing = interpolant_spacing_local(ti);
        }
    //  same as above, but for these points (in Omega_2, but
    //  not on the surface ) we also generate points in the
    //  Omega_1 zone which we will use in interpolation, which we put in
    //  aux_targets
        Point3 tp(near_trg_ptr     + ti*DIM);
        Point3 bp(closest_surf_ptr + ti*DIM);
        Point3 dir;
      if(interpolation_directions == NULL){
              dir = tp-bp;
          dir = dir/dir.length(); //get direction from the closest surface point to the near target point
      } else { // we're passing the directions explicitly
        dir = interpolation_directions_local.clmdata(ti);
        dir = dir/dir.length();
      }
    if(expansion_type == EXTRAPOLATE_ONE_SIDE){ // make all L points in the intermediate zone
        for(int l=0; l<num_aux_targets; l++) {
            // interp point = (first point along ray from bp to tp that is in
            //                 the intermediate zone) + (i*h in the dir
            //                 direction) for i = 1, ..., L
            //Point3 np = (bp + h*dir) + dir * interpolation_nodes[l] * h;
            Point3 np = (bp + distance_from_boundary*dir) + 
                           dir * interpolation_nodes[l] * interpolant_spacing;
            memcpy( &aux_trg_ptr[aux_points_counter*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
            aux_points_counter++;
        }
    } else if (expansion_type == EXTRAPOLATE_TWO_SIDE){ // put L/2 points on inside and L/2 points on outside
        for(int l=0; l<num_aux_targets/2; l++) {
            // interp point = (first point along ray dir from bp to tp that is in
            //                 the intermediate zone) + (i*h in the dir
            //                 direction) for i = 1, ..., L/2
            Point3 np = (bp + distance_from_boundary*dir) + 
                        dir * interpolation_nodes[l] * interpolant_spacing;
            memcpy( &aux_trg_ptr[(aux_points_counter + l)*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
        }
        int exterior_offset = num_local_targets*num_aux_targets/2;
        for(int l=0; l<num_aux_targets/2; l++) {
            // interp point = (first point along ray dir from tp to bp that is in
            //                 the intermediate zone) + (i*h in the dir
            //                 direction) for i = 1, ..., L/2
            //                 note that points are on the opposite side of the
            //                 surface as in the previous loop (minus signs)
            //Point3 np = (bp - h*dir) - dir * interpolation_nodes[l] * h;
            Point3 np = (bp - distance_from_boundary*dir) - 
                            dir * interpolation_nodes[l] * interpolant_spacing;
            memcpy( &aux_trg_ptr[(exterior_offset + aux_points_counter + l)*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
        }
            aux_points_counter += num_aux_targets/2;

    } else if(expansion_type == EXTRAPOLATE_ONE_SIDE_CHEBYSHEV){ // make all L points in the intermediate zone
        for(int l=0; l<num_aux_targets; l++) {
            // interp point = (first point along ray from bp to tp that is in
            //                 the intermediate zone) + (i*h in the dir
            //                 direction) for i = 1, ..., L
            Point3 np = (bp + h*dir) + dir * interpolation_nodes[l];
            memcpy( &aux_trg_ptr[aux_points_counter*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
            aux_points_counter++;
        }
    } else if (expansion_type == INTERPOLATE_ACROSS_SURFACE){ // put L/2 points on inside and L/2 points on outside
        /*
         * Final interpolation point distribution looks as below:
         *                     
         * |---\Omega_1----|Omega_2|---\Omega_1----|
         * |-----L/2*h-----|--2h---|-----L/2*h-----|
         * |-h-|-h-|- ... -|-h-|-h-|-h-|-h-|-h-|-h-|
         * *---*---*- ... -*---X---*- ... -*---*---*
         *                     ^---- \Gamma (on-surface point in this case)
         *     * := interpolation node locations
         * Note that X can be anywhere in \Omega_2
         * EXTRAPOLATE_TWO_SIDE also looks like this, except that we
         * explicitly interpolate each side separately. As a result, we build
         * the interpolation point vector differently so as to allow for
         * contiguous blocks of inside/outside interpolation nodes.
         */

        //int left_most_interp_point = num_aux_targets/2;
        int left_most_interp_point = num_aux_targets/2;
        for(int l=0; l<num_aux_targets/2; l++) {
            // interp point = (furthest point along ray dir from bp to tp that is in
            //                 the intermediate zone) + (i*h in the  -dir
            //                 direction) for i = 1, ..., L/2
            //                 i.e. back towards the surface
            Point3 np = (bp + left_most_interp_point*h*dir) -
                dir * interpolation_nodes[l] * h;
            memcpy( &aux_trg_ptr[aux_points_counter*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
            aux_points_counter++;
        }
        for(int l=0; l<num_aux_targets/2; l++) {
            // interp point = (first point along ray dir from tp to bp that is in
            //                 the intermediate zone) + (i*h in the dir
            //                 direction) for i = 1, ..., L
            //                 note that points are on the opposite side of the
            //                 surface as in the previous loop (minus signs)
            Point3 np = (bp - h*dir) - dir * interpolation_nodes[l] * h;
            memcpy( &aux_trg_ptr[aux_points_counter*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
            aux_points_counter++;
        }
    } else if (expansion_type == INTERPOLATE){ //interpolate with an on-surface point, so generate (L-1) points in intermediate zone
        for(int l=1; l<num_aux_targets; l++) {
            Point3 np = bp + dir * interpolation_nodes[l] * h;
            memcpy( &aux_trg_ptr[aux_points_counter*DIM], 
                    np.array(), 
                    sizeof(double)*DIM );
            aux_points_counter++;
        }
    } else {
        assert(0);
    }
  }
  VecRestoreArray(near_targets, &near_trg_ptr);
  if(sample_point_far_field != NULL){ VecRestoreArray(sample_point_far_field, &far_field_ptr);}
  if(sample_point_interpolant_spacing != NULL){ VecRestoreArray(sample_point_interpolant_spacing, &interpolant_spacing_ptr);}
  VecRestoreArray(closest_surface_points_near, &closest_surf_ptr);
  VecRestoreArray(aux_targets, &aux_trg_ptr);

  cout<<"computed auxillary points for near interpolation"<<endl;
  return aux_targets;
}

Vec generate_interior_qbkix_points(Vec near_targets, 
        Vec closest_surface_points_near, 
        vector<int> qbkix_point_index,
        Vec interpolation_directions, 
        // TODO merge into one vector because they're generally equal
        Vec sample_point_far_field,
        Vec sample_point_interpolant_spacing
        ){

    assert(sample_point_interpolant_spacing != NULL);
    assert(sample_point_far_field != NULL);
    vector<double> interpolation_nodes = get_interpolation_nodes(EXTRAPOLATE_ONE_SIDE);

    
    // make sure the indices we asked for actually are interpolation nodes...
    for(vector<int>::iterator it = qbkix_point_index.begin(); 
            it != qbkix_point_index.end();
            it++){
        assert(find(qbkix_point_index.begin(),
                    qbkix_point_index.end(),
                    *it) != qbkix_point_index.end());
    }

    // we're only generating the qbkix points specified in qbkix_point_index

    int num_local_targets = num_local_points(near_targets);
    int64_t num_qbkix_targets = qbkix_point_index.size()*num_local_targets;
    cout << "num generated qbx points" << num_qbkix_targets << endl;

    // if we're extrapolating, all L points are generated here
    // if we're not, one point on-surface is generated externally, 
    //      so only generate (L-1) points here.
    Vec qbkix_points;
    VecCreateMPI(PETSC_COMM_WORLD,
            DIM*num_qbkix_targets, 
            PETSC_DETERMINE, 
            &qbkix_points);
    DblNumMat qbkix_points_local = get_local_vector(DIM, num_qbkix_targets, qbkix_points);

    // TODO MJM REFACTOR 
    // make get_local_vector() and restore() for DblNumVec
    double* near_trg_ptr;
    double* closest_surf_ptr;
    //double* aux_trg_ptr;
    double* interpolation_directions_ptr;
    double* far_field_ptr;
    double* interpolant_spacing_ptr;

    VecGetArray(near_targets, &near_trg_ptr);
    VecGetArray(closest_surface_points_near, &closest_surf_ptr);
    //VecGetArray(qbkix_points, &aux_trg_ptr);
    VecGetArray(sample_point_far_field, &far_field_ptr);
    VecGetArray(sample_point_interpolant_spacing, &interpolant_spacing_ptr);

    DblNumVec far_field_local(num_local_targets, false, far_field_ptr);
    DblNumVec interpolant_spacing_local(num_local_targets, false, interpolant_spacing_ptr);

    DblNumMat interpolation_directions_local;
    if(interpolation_directions != NULL){
        interpolation_directions_local = get_local_vector(DIM, num_local_targets, interpolation_directions);
    }
 



    for(int ti=0; ti<num_local_targets; ti++) {
        double distance_from_boundary = far_field_local(ti);;
        double interpolant_spacing = interpolant_spacing_local(ti);;
        //  same as above, but for these points (in Omega_2, but
        //  not on the surface ) we also generate points in the
        //  Omega_1 zone which we will use in interpolation, which we put in
        //  aux_targets
        Point3 tp(near_trg_ptr     + ti*DIM);
        Point3 bp(closest_surf_ptr + ti*DIM);
        Point3 dir;
        if(interpolation_directions == NULL){
            dir = (tp-bp);
            // if they're too close, you end up dividing by zero when normalizing, 
            // leading to an excitingly vague fmm nan potential assert.
            assert(dir.length() > 1e-13);
            dir = dir/dir.length(); //get direction from the closest surface point to the near target point
        } else { // we're passing the directions explicitly
            dir = Point3(interpolation_directions_local.clmdata(ti));
            dir = dir.dir();
        }


        for(int i = 0; i < qbkix_point_index.size(); i++){
            int index = qbkix_point_index[i];
            
            Point3 qbkix_point = 
                (bp + distance_from_boundary*dir) + 
                dir * interpolation_nodes[index] * interpolant_spacing;

            for(int d = 0; d < DIM; d++)
                qbkix_points_local(d, ti*qbkix_point_index.size() + i) = qbkix_point(d);
        }
    }
    VecRestoreArray(near_targets, &near_trg_ptr);
    VecRestoreArray(sample_point_far_field, &far_field_ptr);
    VecRestoreArray(sample_point_interpolant_spacing, &interpolant_spacing_ptr);
    VecRestoreArray(closest_surface_points_near, &closest_surf_ptr);

    qbkix_points_local.restore_local_vector();

    return qbkix_points;
}


// --------------------------------------------------------------------------- 

int EvaluatorNear::eval(Vec den, Vec val)
{
  ebiFunctionBegin;
  int source_dof = this->source_dof();
  int target_dof = this->target_dof();  
  double h = this->_patch_samples->spacing();// * alfcoef();
  
  //1. den -> detailed dstz version, and eval using fmm
  PatchSurf* bd = this->_patch_samples->bdry();
  PatchSamples* patch_samples = this->_patch_samples;
  PatchSamples* refined_patch_samples = this->_refined_patch_samples;
  
  vector<DblNumMat> _refined_datvec; 

  // make a copy of the input, scale it by blending function, and refine to upsampled grid
  Vec dat;
  VecDuplicate(den, &dat);
  denscale(source_dof, patch_samples->sample_point_blend_func_value(), den, dat);
  patch_samples->refine_data(source_dof, _refinement_factor, dat, _refined_datvec);

  Vec refined_sample_position = refined_patch_samples->sample_point_3d_position();
  Vec refined_sample_weight = refined_patch_samples->sample_point_combined_weight();
  
  // interpolate refined density to the refined surface
  Vec scaled_density = interpolate_density_to_refined_grid( 
          patch_samples, 
          refined_sample_position,
          refined_sample_weight,
          _refined_datvec,
          _refined_collocation_data,
          this->source_dof());
  
  iC( VecDestroy(&dat) );


 //RFDDN VERY IMPORTANT
  int num_local_targets = num_local_points(_target_positions_intermediate);
  
  Vec valfar;
  VecCreateMPI(this->mpiComm(),
          num_local_targets*target_dof, 
          PETSC_DETERMINE,
          &valfar);

  // Evaluate intermediate points for interpolation (targets =
  // _target_positions_intermediate in setup)
  iC( fmm->evaluate(scaled_density, valfar) );
  iC( VecDestroy(&scaled_density) );

  // copy the FMM result back 
  
  //3. den -> be3dovron, be3dovjmp
  num_local_targets = num_local_points(_closest_sample_3d_position);
  
  Vec roncls; 
  VecCreateMPI(this->mpiComm(),
              num_local_targets*target_dof, 
              PETSC_DETERMINE,
              &roncls);

  // on-the-surface evaluation at closest sample points for all points in the near zone
  fmm->initialize_fmm(this->_patch_samples->sample_point_3d_position(),
                      this->_patch_samples->sample_point_normal(),
                      _closest_sample_3d_position,
                      this->knl(),
                      false);

  fmm->set_rebuild_tree(true);

  //Scale density
  Vec unrefined_scaled_density;
  VecDuplicate(den, &unrefined_scaled_density);
  denscale(this->source_dof(),
          this->_patch_samples->sample_point_combined_weight(),
          den,
          unrefined_scaled_density);

  //FMM eval + singular_correction
  //MJM Note: replaces _on_surface_evaluator->eval() and uses same FMM as for
  //the refined evaluation
  fmm->evaluate(unrefined_scaled_density, roncls);
  _on_surface_evaluator->singular_evaluation(den, roncls);
  _on_surface_evaluator->apply_singularity_cancellation(den, roncls);

  // jump evaluation at the closest sample points
  Vec jump_potential;
  VecCreateMPI(this->mpiComm(), 
              num_local_targets*target_dof, 
              PETSC_DETERMINE, 
              &jump_potential);
          
  jump_evaluation(den, jump_potential, _refined_datvec, this->knl(), 
          this->_patch_samples, this->_collocation_data);
  
  cout<<"nea eval re+je  "<<endl;


  // MJM TODO reduce number of arguments?
  interpolation_near(_target_positions_intermediate, 
          _closest_sample_3d_position,
          _target_3d_position, 
          valfar, 
          roncls,
          jump_potential, 
          _refined_datvec, 
          _target_in_out, // MJM BUG is this arg. correct?
          val,
          target_dof,
          INTERPOLATE);
  
  iC( VecDestroy(&valfar) );
  iC( VecDestroy(&roncls) );
  iC( VecDestroy(&jump_potential) );
  ebiFunctionReturn(0);
}

void jump_evaluation(Vec den, Vec val, vector<DblNumMat>& _refined_datvec, 
        Kernel3d kernel, PatchSamples* patch_samples, 
        CollocationPatchSamples* _collocation_data) {
 
    int64_t _refinement_factor, _surface_interpolation_num_samples;
    PetscBool err = PETSC_FALSE;
    PetscOptionsGetInt(NULL, "", "-dnref", &_refinement_factor, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-LL", &_surface_interpolation_num_samples, &err);
    ebiAssert(err);
        
  double zero = 0.0;
  iC( VecSet( val, zero) );
  
  PatchSurf* bd = patch_samples->bdry();
  vector<int>& prtn = patch_samples->patch_partition();
  
  int source_dof = kernel.get_sdof();
  int target_dof = kernel.get_tdof();

  Vec _colloc_point_3d_position = _collocation_data->colloc_point_3d_position();
  int num_local_colloc_points = num_local_points(_colloc_point_3d_position);
  
  // scale density by blending function value store in dat
  Vec dat; 
  iC( VecDuplicate(den,&dat) );
  iC( denscale(source_dof,patch_samples->sample_point_blend_func_value(),den,dat) );
  
  // refine data to the upsampled grid 
  iC( patch_samples->refine_data(source_dof, _refinement_factor, dat, _refined_datvec) );
  
  //2. interpolating
  // temp vector for collocation points 
  Vec density_at_collocs;
  iC( VecCreateMPI(PETSC_COMM_WORLD, num_local_colloc_points*target_dof, 
                    PETSC_DETERMINE, &density_at_collocs) );
  //------------------------------------------------------------------
  int kt = kernel.kernelType();


  for(int pi=0; pi<bd->patches().size(); pi++) {
	 Patch* curpch = bd->patches()[pi];

	 // interface to collocation face points in tdata
	 DblNumMat colloc_point_as_face_point(
             _collocation_data->colloc_point_as_face_point(pi, patch_samples->face_point_size_in_doubles()));

	 // interface to density_at_collocs
	 DblNumMat colval(_collocation_data->coldat(pi, density_at_collocs, target_dof));

       
     int num_col_points = _collocation_data->num_collocation_points(pi);
	 if(kt == KNL_LAP_D_U || kt == KNL_MODHEL_D_U ||
             kt==KNL_STK_D_U || kt == KNL_NAV_D_U) {
#pragma omp parallel for
	   // for all collocation points on the patch
		for(int ci=0; ci<num_col_points; ci++) {
		  FacePointOverlapping* cur_face_point = (FacePointOverlapping*)(colloc_point_as_face_point.clmdata(ci));
		  double* curval = colval.clmdata(ci);

		  bool is_valid;
		  iC( curpch->is_face_point_valid(cur_face_point, is_valid) );
		  ebiAssert(is_valid==true);

		  double xy[2]; iC( curpch->face_point_to_xy(cur_face_point, xy) );
		  iC( patch_samples->interpolate_data(pi, xy, _refinement_factor,
                      _surface_interpolation_num_samples,
                      PatchSamples::EVAL_VL, _refined_datvec, curval) );
		}
	 } else if (kt == KNL_STK_D_P) {
	   // kernel is hypersingular, the jump is more complicated;
	   // p. 258, p 259, 3.3.1, and appendix B, code specific for
	   // pressure, although can be generalized
#pragma omp parallel for
		for(int ci = 0; ci < num_col_points; ci++) {
		  FacePointOverlapping* cur_face_point = (FacePointOverlapping*)
              (colloc_point_as_face_point.clmdata(ci));
		  double* curval = colval.clmdata(ci);

		  // debugging
		  bool is_valid;
		  iC( curpch->is_face_point_valid(cur_face_point, is_valid) );
		  ebiAssert(is_valid==true);

		  // get face position,  	
		  double xy[2];
          iC( curpch->face_point_to_xy(cur_face_point, xy) );

		  //-------------------------------------------
		  // Stokes pressure -specific coefficient
		  double mu = kernel.coefs()[0];

		  // value and derivatives of the surface
		  Point3 point_position_and_derivs[3];
          iC( patch_samples->interpolated_position_and_derivatives(
                  pi, /* pi-th patch */
                  xy, /*current sample point */
                  PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, /*eval flag */
                  (double*)point_position_and_derivs) /* result pos + derivatives */
              );
		  
          Point3& ga = point_position_and_derivs[1];
		  Point3& gb = point_position_and_derivs[2];

		  // normal and its length
		  Point3 normal(cross(ga,gb)); 
          Point3 unit_normal = normal/(normal.l2());

		  double nga = ga.l2();
		  // second tangent orthoginal to ga
          //
		  Point3 tm( gb - (ga*gb)/(nga*nga) * ga );
		  double ntm = tm.l2();

		  // evaluate the jump using the formula on p 258
		  // -2*mu*(a^T dphi/da + b^Tdphi/db) where a,b is an
		  // orthonormal basis in the tangent plane, here given by
		  // ga/|ga|, tm/|tm|
		  // folds everything into a hard-to-parse matvec
		  
		  double arr[9];
		  DblNumMat ARR(3,3,false,arr);
		  ARR(0,0) = 1.0;				ARR(0,1) = 0.0;	        ARR(0,2) = 0.0;
		  ARR(1,0) = 0.0;				ARR(1,1) = 1.0/nga; 	ARR(1,2) = -(ga*gb)/(nga*nga)/ntm;
		  ARR(2,0) = 0.0;				ARR(2,1) = 0.0;         ARR(2,2) = 1.0/ntm;
		  double fng[9];
           DblNumMat FNG(3,3,false,fng); //final geometry

		  for(int d=0; d<3; d++) FNG(d,0) = unit_normal(d);     //col 0
		  for(int d=0; d<3; d++) FNG(d,1) = ga(d)/nga;          //col 1
		  for(int d=0; d<3; d++) FNG(d,2) = tm(d)/ntm;          //col 2
		  
		  double old_density[9]; // MJM TODO remove me 
          DblNumMat old_density_mat(3,3,false,old_density); //old density

		  iC( patch_samples->interpolate_data(pi, xy, _refinement_factor,
                      _surface_interpolation_num_samples,
                      PatchSamples::EVAL_VL|PatchSamples::EVAL_FD,
                      _refined_datvec,
                      old_density) );
		  
          double final_density[9]; // MJM TODO remove me 
          DblNumMat final_density_mat(3,3,false,final_density); //final density

		  iC( dgemm(1.0, old_density_mat, ARR, 0.0, final_density_mat) );
		  for(int k=1; k<3; k++) //num_collocation_points
			 for(int d=0; d<3; d++)
				curval[0] += FNG(d,k)*final_density_mat(d,k);
		  curval[0] *= (-2.0*mu);
		}
	 } else {
		ebiAssert(0);
	 }
  }
  //3. vecscatter to put the summation back into val
  iC( VecScatterBegin( _collocation_data->colloc_to_target_value_scatter(), 
                        density_at_collocs,
                        val,
                        ADD_VALUES,
                        SCATTER_FORWARD) );
  iC( VecScatterEnd( _collocation_data->colloc_to_target_value_scatter(),
                      density_at_collocs,  
                      val,
                      ADD_VALUES,
                      SCATTER_FORWARD) );
  
  iC( VecDestroy(&density_at_collocs) );
  

}
// ---------------------------------------------------------------------- 

int compute_interpolation_weights(int n, double* p, double u, double* w)
{
  ebiFunctionBegin;
  for(int j=0; j<n; j++) {
	 w[j] = 1.0;
	 for(int i=0; i<n; i++) {
		if(i!=j) w[j] *= (u-p[i])/(p[j]-p[i]);
	 }
  }
  ebiFunctionReturn(0);
}

Vec interpolate_density_to_refined_grid(PatchSamples* patch_samples, 
        Vec& refined_sample_point_3d_position,
        Vec& refined_sample_point_combined_weight,
        vector<DblNumMat>& _refined_datvec,
        CollocationPatchSamples* _refined_collocation_data,
        int source_dof){

    // Load Petsc options
  int64_t refinement_factor = Options::get_int_from_petsc_opts("-dnref");
  int64_t surface_interpolation_num_samples= Options::get_int_from_petsc_opts("-LL");

  //1. den -> detailed dstz version, and eval using fmm
  PatchSurf* bd = patch_samples->bdry();
  int num_local_colloc_points = 
      num_local_points(_refined_collocation_data->colloc_point_3d_position());
  int num_local_targets = 
      num_local_points(refined_sample_point_3d_position);

  // temporary vectors for collocation and target points 
  Vec density_at_collocs;
  Vec density_at_targets;
  
  VecCreateMPI(PETSC_COMM_WORLD,
          num_local_colloc_points*source_dof, 
          PETSC_DETERMINE, 
          &density_at_collocs);

  VecCreateMPI(PETSC_COMM_WORLD, 
          num_local_targets*source_dof, 
          PETSC_DETERMINE, 
          &density_at_targets);

  
  int face_point_sze_dbl = patch_samples->face_point_size_in_doubles();
  for(int pi=0; pi<bd->patches().size(); pi++) {
      int num_colloc_points = _refined_collocation_data->num_collocation_points(pi);
		Patch* curpch = bd->patches()[pi];
		DblNumMat colloc_point_as_face_point(
            _refined_collocation_data->colloc_point_as_face_point(pi, face_point_sze_dbl));
		DblNumMat colden(_refined_collocation_data->coldat(pi, density_at_collocs, source_dof));

#pragma omp parallel for
		// for all collocation points on this patch coming from the refined grid
		for(int ci = 0; ci < num_colloc_points; ci++) {
		  FacePointOverlapping* cur_face_point =
              (FacePointOverlapping*)(colloc_point_as_face_point.clmdata(ci));

          // debugging code
          bool is_valid;
		  iC( curpch->is_face_point_valid(cur_face_point, is_valid) );
		  ebiAssert(is_valid==true); //inside

		  // interpolate density at xy, store in density_at_collocs (colden is an interface to density_at_collocs)
		  double xy[2]; iC( curpch->face_point_to_xy(cur_face_point,xy) );
		  iC( patch_samples->interpolate_data(pi, xy, refinement_factor, 
                     surface_interpolation_num_samples, 
                     PatchSamples::EVAL_VL, _refined_datvec, colden.clmdata(ci)) );
		}
  }
  
  iC( VecScatterBegin( 
              _refined_collocation_data->colloc_to_target_density_scatter(),
              density_at_collocs,  
              density_at_targets,  
              ADD_VALUES,  
              SCATTER_FORWARD) );
  iC( VecScatterEnd( 
              _refined_collocation_data->colloc_to_target_density_scatter(),   
              density_at_collocs,  
              density_at_targets,  
              ADD_VALUES,  
              SCATTER_FORWARD) );
  Vec scaled_density;
  VecDuplicate(density_at_targets, &scaled_density);
  iC( denscale(source_dof, refined_sample_point_combined_weight,
              density_at_targets, scaled_density) );
  iC( VecDestroy(&density_at_collocs) );
  iC( VecDestroy(&density_at_targets) );
  cout<<"nea eval interpolate_data  " << endl;
  return scaled_density;

}





// Vec intermediate_targets has the auxiliary interpolation points in the first
// (# of near targets)*(interpolation order - 1)*DIM elements 
void interpolation_near(Vec& intermediate_targets, Vec& closest_surface_targets,
        Vec& near_targets, Vec& intermediate_potentials, 
        Vec& closest_surface_potential, Vec& jump_potential, 
        vector<DblNumMat>& _refined_datvec, Vec& near_target_in_out,
        Vec& val, int target_dof, ExpansionType expansion_type){
// ----------------------------------------------------------------------------
    PetscBool err = PETSC_FALSE;
    double h;
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &h, &err);
    ebiAssert(err);

  vector<double> interpolation_nodes = get_interpolation_nodes(expansion_type);

  int L = interpolation_nodes.size();

  int num_local_total_intermed_targets = num_local_points(intermediate_targets); 
  int num_local_near_targets = num_local_points(near_targets);
  int num_local_intermed_targets = 
      num_local_total_intermed_targets - (L-1)*num_local_near_targets;
  
  
// ----------------------------------------------------------------------------
  double* jvarr;
  double* rvarr;
  double* fvarr;
  double* tiarr;
  double* near_targets_ptr;
  double* closest_boundary_point_ptr;
  double* vlarr;
  cout<<"nea eval re+je  "<<endl;

  VecGetArray(jump_potential, &jvarr);
  VecGetArray(closest_surface_potential, &rvarr);
  VecGetArray(intermediate_potentials, &fvarr);
  VecGetArray(near_target_in_out, &tiarr);
  VecGetArray(near_targets, &near_targets_ptr);
  VecGetArray(closest_surface_targets, &closest_boundary_point_ptr);
  VecGetArray(val, &vlarr);

  /* for each target T in \Omega_2
   *  D = distance from T to the surface \Gamma
   *  N = unit vector pointing from T to nearest point on \Gamma
   *
   *  Note: P(T) = potential at target T
   *  if D < 1e-8:
   *    Depending on whether T is inside/outside \Gamma, 
   *    P(T) = on_surface_eval(T) +/- .5* jump_potential(T) 
   *  else if 1e-8 <= D < h:
   *    P(nearest_point_on_surface) = on_surface_eval(nearest_point_on_surface)
   *                        +/- .5 * jump_potential(nearest_point_on_surface) 
   *    Using P(nearest_point_on_surface) and the intermediate potentials
   *    evaluated in FMM, interpolate to compute potential at P(T)
   *
   *  else: //D >= h
   *    already computed in FMM, copy values into output array
   */
  DblNumMat near_targets_local(DIM, num_local_near_targets, false, near_targets_ptr);
  DblNumMat closest_surface_targets_local(DIM, num_local_near_targets, false, closest_boundary_point_ptr);
  DblNumMat closest_surface_potential_local(target_dof, num_local_near_targets, false, rvarr);
  DblNumMat jump_potential_local(target_dof, num_local_near_targets, false, jvarr);
  DblNumMat final_potential_local(target_dof, num_local_near_targets + num_local_intermed_targets, false, vlarr);

  // here we assume that in the list of the intermediate zone 
  // ("far" = Omega_1 or beyond) points, the points added for 
  // near zone interpolation come first, in groups of L-1
  DblNumMat interp_point_potential_local(target_dof, num_local_near_targets*(L-1), false, fvarr);
            
  DblNumMat interpolation_point_potentials(target_dof, L);
  DblNumVec target_in_out_local(num_local_near_targets,false, tiarr);
  
    for(int i = 0; i < num_local_near_targets; i++){
        Point3 target(near_targets_local.clmdata(i));
        Point3 point_on_boundary(closest_surface_targets_local.clmdata(i));
        Point3 ray_toward_surface(target - point_on_boundary);
        double distance_to_closest_sample = ray_toward_surface.length();


        // Use on-surface evaluated potential, but correct by the jump
        // potential, since we're really taking the interior limit  to the
        // surface
		// value = that produced by the on-the-surface evaluator +/- 1/2 of 
        // the jump (depending if it is inside or out)
        DblNumVec closest_surface_potential(target_dof, false, 
                closest_surface_potential_local.clmdata(i));
        DblNumVec jump_potential(target_dof, false, 
                jump_potential_local.clmdata(i));
        double jump_correction_magnitude;
        
        if (target_in_out_local(i)==1){
            jump_correction_magnitude = 0.5;
        } else if (target_in_out_local(i) == 0){
            jump_correction_magnitude =  -0.5;
        } else {
            // Bail; only can accept 1/0 as in/out marker
            ebiAssert(false);
        }

        if(distance_to_closest_sample < 1e-8){
            // points on the surface DZ FIXME: replace 1e-8 with a parameter
            // MJM WARNING TODO this is untested; verify that this result
            // matches the pointer based loop commented out below
            
            
            // value = that produced by the on-the-surface evaluator +/- 1/2 of 
            // the jump (depending if it is inside or out)
            for(int d = 0; d < target_dof; d++){
                final_potential_local(d, i) = 
                                closest_surface_potential(d) + 
                                jump_correction_magnitude * jump_potential(d);

            }
        } else if (distance_to_closest_sample < h){
            
            // as for on-the surface points, combine on-the-surface value with the jump
            // put it as a first point into the list of points for interpolation 
            for(int d = 0; d < target_dof; d++){
                interpolation_point_potentials(d, 0) = 
                                closest_surface_potential(d) + 
                                jump_correction_magnitude * jump_potential(d);
            }
            
            // using up target_dof*(L-1) doubles in buf, in addition to 
            // target_dof doubles already there from the on-surface corrected 
            // potential
            for(int j = 0; j < L-1; j++){
                for(int d = 0; d < target_dof; d++){
                    interpolation_point_potentials(d, 1 + j) = 
                        interp_point_potential_local(d, i*(L-1) + j);
                }
            }

            DblNumVec interpolation_weights(L);

            // compute the polynomial interpolation weights
            compute_interpolation_weights(
                    interpolation_nodes.size(),
                    &(interpolation_nodes[0]),
                    distance_to_closest_sample/h,
                    interpolation_weights.data());

            for(int ll = 0; ll < L; ll++){
                for(int d = 0; d < target_dof; d++){
                    final_potential_local(d, i) += 
                                interpolation_weights(ll) * 
                                interpolation_point_potentials(d, ll);

                }
            }
            clear(interpolation_point_potentials);
            clear(interpolation_weights);
        }
    }
  // MJM WARNING TODO this is untested; verify that this result
  // matches the pointer based loop commented out below
  DblNumMat intermediate_target_potential_local(target_dof, 
          num_local_intermed_targets, 
          false, 
          fvarr + num_local_near_targets*(L-1));

     // intermediate zone points, simply copy values
     // skipping points generated for interpolation 
     ///REFACTOR NOTE: this is where the assumption regarding intermediate
     // point ordering (i.e. interpolation points first) is used
    for(int i = 0; i < num_local_intermed_targets; i++){
        for(int d = 0; d < target_dof; d++){
            final_potential_local(d, num_local_near_targets + i) = 
                intermediate_target_potential_local(d,i);
        }
    }
/*
  int num_points_near = 0;
  int num_points_intermediate = 0;
  for(int ti=0; ti < num_local_near_targets; ti++) {
      Point3 tp(near_targets_ptr+ti*DIM);
      Point3 bp(closest_boundary_point_ptr+ti*DIM);
      Point3 dir(tp - bp);
      double distance_to_closest_sample = dir.length();
    if(distance_to_closest_sample<1e-8 ) { 
        // points on the surface DZ FIXME: replace 1e-8 with a parameter
		double* currv = rvarr + num_points_near*target_dof;
		double* curjv = jvarr + num_points_near*target_dof;
		double mag = ((int)(tiarr[ti])==1) ? 0.5 : -0.5;
		double* curvl = vlarr + ti*target_dof;

		// value = that produced by the on-the-surface evaluator +/- 1/2 of 
        // the jump (depending if it is inside or out)
		for(int d=0; d<target_dof; d++)
		  curvl[d] = currv[d] + mag * curjv[d];

		num_points_near ++;
	 } else if( distance_to_closest_sample<h ) { 
         // near zone points, interpolation 
         double buf[L*target_dof];
         double* currv = rvarr + num_points_near*target_dof;
         double* curjv = jvarr + num_points_near*target_dof;

         // as for on-the surface points, combine on-the-surface value with the jump
         // put it as a first point into the list of points for interpolation 
         double mag = ((int)(tiarr[ti])==1) ? 0.5 : -0.5;
         for(int d=0; d<target_dof; d++) 
             buf[d] = currv[d] + mag * curjv[d];
         
         // here we assume that in the list of the intermediate zone 
         // ("far" = Omega_1 or beyond) points, the points added for 
         // near zone interpolation
		 // come first, in groups of L-1
         
         // REFACTOR NOTE: this is where the assumption regarding intermediate 
         // interpolation point ordering is used.
		double* curfv = fvarr + num_points_intermediate*target_dof;

		// using up target_dof*(L-1) doubles in buf, in addition to target_dof doubles already there
		memcpy(buf+target_dof, curfv, target_dof*(L-1)*sizeof(double) );

        // MJM @BUG  only needs to be of size L, not L*target_dof
		double wgt[L*target_dof];

		// compute the polynomial interpolation weigts
        EvaluatorNear::compute_interpolation_weights(
                interpolation_nodes.size(),
                &(interpolation_nodes[0]),
                distance_to_closest_sample/h,
                wgt);

		double* curvl = vlarr + ti*target_dof;

		for(int l=0; l<L; l++)
		  for(int d=0; d<target_dof; d++)
			 curvl[d] += wgt[l] * buf[l*target_dof+d];

		num_points_near ++;
		num_points_intermediate += (L-1);
	 } else {
         // intermediate zone points, simply copy values
         // skip points generated for interpolation 
		double* curfv = fvarr + num_points_intermediate*target_dof;
		double* curvl = vlarr + ti*target_dof;
		for(int d=0; d<target_dof; d++)
		  curvl[d] = curfv[d];
		num_points_intermediate++;
	 }
  }
  */
  
  
  VecRestoreArray(jump_potential, &jvarr);
  VecRestoreArray(closest_surface_potential, &rvarr);
  VecRestoreArray(intermediate_potentials, &fvarr);
  VecRestoreArray(near_target_in_out, &tiarr);
  VecRestoreArray(near_targets, &near_targets_ptr);
  VecRestoreArray(closest_surface_targets, &closest_boundary_point_ptr);
  VecRestoreArray(val, &vlarr);
  cout << "interpolation near complete" << endl; 

}





END_EBI_NAMESPACE
