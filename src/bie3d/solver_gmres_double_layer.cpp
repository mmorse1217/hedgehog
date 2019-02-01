#include "solver_gmres_double_layer.hpp"
#include "evaluator_far.hpp"
#include "evaluator_near.hpp"
#include "common/mobo_blas.h"
#include "common/mobo_lapack.h" 
#include "common/vtk_writer.hpp" 
#include "bie3d/solver_utils.hpp"
#include "evaluator_qbkix.hpp"
#include "evaluator_qbkix_average.hpp"
#include "evaluator_near_interpolate.hpp"
#include "markgrid.hpp" 
#include "bdry3d/p4est_interface.hpp"
#include "bdry3d/p4est_refinement.hpp"
BEGIN_EBI_NAMESPACE

using std::cerr;

// ---------------------------------------------------------------------- 
SolverGMRESDoubleLayer::SolverGMRESDoubleLayer(const string& n, const string& p)
    : SolverInterface(n,p),
      _ksp(NULL),
      _mat(NULL),  
      _on_surface_evaluator(NULL),
      _pole_matrix(NULL),
      _constraint_matrix(NULL),
      _interpolation_directions(NULL),
      _target_in_out(NULL),
      precond_schur_complement_inverse(NULL),
      _evaluation_type(SINGULAR_EVAL), //legacy singular evaluation by default
      _compute_refined_surface(true),
      _nits(0)  {

    //int eqn_opt;

    //PetscBool flag;
    //PetscInt tmp;
    //PetscOptionsGetInt(NULL, "", "-kt", &tmp,  &flag);
    int eqn_opt = Options::get_int_from_petsc_opts("-kt");

    // Form enums to construct particular kernel
    // See common/kernel3d.hpp for details
    Kernel_type empty1; Kernel_variable empty2;
    Kernel3d::parse_kernel_enum(eqn_opt, _equation_type,  empty1, empty2 );
    Kernel3d temp_kernel(_equation_type + DOUBLE_LAYER + VAR_U, vector<double>());

    _is_incompressible = temp_kernel.is_incompressible();
    cout << "is incompressibile: " << _is_incompressible << endl;

}
SolverGMRESDoubleLayer::SolverGMRESDoubleLayer(PatchSurf* surface):
   SolverGMRESDoubleLayer("BIS3D_", "bis3d_"){
       _bdry = surface;
}

// ---------------------------------------------------------------------------- 
SolverGMRESDoubleLayer::~SolverGMRESDoubleLayer()
{
    if(_ksp != NULL) {
        KSPDestroy(&_ksp); 
        _ksp=NULL;
    }
    if(_mat != NULL) {
        MatDestroy(&_mat); 
        _mat=NULL;
    }
    if(_on_surface_evaluator != NULL) {
        delete _on_surface_evaluator; 
        _on_surface_evaluator=NULL;
    }
    if(_pole_matrix != NULL) {
        MatDestroy(&_pole_matrix); 
        _pole_matrix=NULL; 
    }
    if(_constraint_matrix != NULL) {
        MatDestroy(&_constraint_matrix);
        _constraint_matrix=NULL; 
    }
    if(precond_schur_complement_inverse != NULL) {
        MatDestroy(&precond_schur_complement_inverse);
        precond_schur_complement_inverse=NULL;
    }
    if(_reference_solution != NULL){
        delete _reference_solution;
    }
    if(_target_in_out != NULL){
        VecDestroy(&_target_in_out);
    }
    if(_interpolation_directions != NULL){
        VecDestroy(&_interpolation_directions);
    }
    /*if(_patch_samples!= NULL){
        delete _patch_samples;
    }
    if(_refined_patch_samples!= NULL){
        delete _refined_patch_samples;
    }*/
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::setFromOptions(){

    // make sure we have an initialized surface
  assert(_bdry);

  // make sure things aren't run with mpi
  int mpi_size;
  MPI_Comm_size(this->mpiComm(), &mpi_size);
  //assert(mpi_size == 1);

  // used to determine which patches are run on which process
  
  _patch_partition = vector<int>(_bdry->num_patches(), mpiRank());
  _dom = Options::get_int_from_petsc_opts("-dom"); 

    // specify which matvec to use inside GMRES
    //solver->set_evaluation_type(eval_type); 
    //solver->_compute_refined_surface = compute_refined_surface;

}
// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::setup()
{
  ebiFunctionBegin;
  //-----------------------
  // Validate that the input .wrl file is properly formatted with respect to the 
  // command-line arguments and the Petsc options file.
  // 
  // Convention: for bounded domains, the first surface in the .wrl file is the
  // bounding surface with orientation== 1, the remaining surfaces are the
  // internal bounding surfaces with orientation == -1. For unbounded surfaces, 
  // all surfaces have orientation == -1.
  
  vector<int>& boundary_orientation = this->bdry()->boundary_orientation();
  /*vector<int> boundary_orientation;
  boundary_orientation.reserve(_bdry->num_patches());
  for(auto p : _bdry){
      boundary_orientation.push_back(p->_orientation);
  }*/
  if(this->dom()==DOM_BND) { //verify dom and bdry match
	 ebiAssert(boundary_orientation[0]==1);
     
     for(int k=1; k<boundary_orientation.size(); k++)
         ebiAssert(boundary_orientation[k]==-1);
  } else {

	 for(int k=0; k<boundary_orientation.size(); k++)
         ebiAssert(boundary_orientation[k]==-1);
  }

  //---------------------------------------------------------------------------
  // 
  //---------------------------------------------------------------------------
  // Initialize the relavent kernels required for the solver and evaluators
  // All equations need single-layer, double-layer and the identity kernel.
  // The Stokes equation requires a Rotlet kernel and "Levi-Cevita kernel"
  // The identity kernel for n sources and m targets, fills out m x n
  // block matrix with identity matrices (assuming sdof = tdof, only used
  // as a part of constructing the matrix for the gmres solve, for which it
  // is always true
  // The L-C kernel fills an m x n block matrix with 3x3
  // cross-product  matrices corresponding to multiplication by trgPos - srcPos

  
  Kernel3d single_layer_kernel, double_layer_kernel, rotlet_kernel;
  single_layer_kernel = Kernel3d(_equation_type + SINGLE_LAYER + VAR_U, this->eqcoefs());
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + VAR_U, this->eqcoefs());
  
  _problem_kernel = single_layer_kernel;

  cout << "equation type: " << (int) _equation_type << endl;

  
  // Access kernel dependent degrees of freedom
  // Source degree of freedom
  //int sdof = single_layer_kernel.get_sdof();

  // Target degree of freedom
  //int tdof = single_layer_kernel.get_tdof();
  
  // Rotational degrees of freedom 
  //int rdof = single_layer_kernel.get_rdof();
  
  // Degrees of freedom contributed from multiply connected boundaries
  //int gdof = sdof; 

  // Total degrees of freedom of the current problem
  //int pdof = rdof + gdof;

  
  if (_problem_kernel.get_rdof() > 0) {
      // Rotlet kernel
      rotlet_kernel = Kernel3d(_equation_type + ROTLET + VAR_U, this->eqcoefs());
  }
  //---------------------------------------------------------------------------
  
  
  // Initialize two discretizations of the surface, one normal  and another 
  // refined, and populate them accordingly
  this->_patch_samples = new PatchSamples("", "");
  PatchSamples* patch_samples = this->_patch_samples;

  patch_samples->bdry()                = this->bdry();
  patch_samples->patch_partition()     = this->patch_partition();  

  patch_samples->setup();
  cerr<<"bis patch_samples setup"<<endl;
 PatchSamples* refined_patch_samples;
  if(Options::get_int_from_petsc_opts("-bdtype") == 2 && _compute_refined_surface){
      // TODO possibly remove this. could copy the patches from _bdry instead of running an
      // entirely new setup...?
      auto coarse_face_map = dynamic_cast<PatchSurfFaceMap*>(this->bdry());
      auto refined_face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");

      refined_face_map->_surface_type = coarse_face_map->_surface_type;
      refined_face_map->setup_from_existing_face_map(coarse_face_map);
      //refined_face_map->initialize_with_existing_p4est(coarse_face_map->_p4est);
      //refined_face_map->_quadrature_weights = coarse_face_map->_quadrature_weights;
      //refined_face_map->setFromOptions();
      //refined_face_map->setup();
      string upsampling_type = Options::get_string_from_petsc_opts("-upsampling_type");
      if(upsampling_type == "uniform"){
          int num_levels_upsampling = Options::get_int_from_petsc_opts("-uniform_upsampling_num_levels");
          refined_face_map->refine_uniform(num_levels_upsampling);
      } else if (upsampling_type == "adaptive"){

          refine_patches_for_fixed_qbkix_points(refined_face_map->_p4est, refined_face_map);
          refined_face_map->patches() = p4est_to_face_map_subpatches(refined_face_map->_p4est, refined_face_map);
      }
      /*vector<int> pids;
      for(int i =0; i < refined_face_map->num_patches(); i++)
          pids.push_back(i);
      write_face_map_patches_to_vtk(DblNumMat(0,0), pids, refined_face_map, 10, stats._file_prefix+"_upsampled_");*/
      stats.add_result("num upsampled_patches", refined_face_map->num_patches());
        /*
        // Refine until points are inside
        p4est_connectivity_t* connectivity = build_connectivity_from_face_map(refined_face_map);
        p4est_t* p4est = p4est_new(MPI_COMM_WORLD, connectivity, sizeof(RefinementData<FaceMapPatch>), NULL, NULL);


        refine_patches_for_fixed_qbkix_points(p4est, refined_face_map);
        //refine_patches_uniform(2, p4est, refined_face_map);
        
        vector<Patch*> subpatches = p4est_to_face_map_subpatches(p4est, refined_face_map);
        refined_face_map->patches() = subpatches;

        */
        _refined_surface = refined_face_map; 
        this->_refined_patch_samples = new PatchSamples("", "");
        refined_patch_samples = this->_refined_patch_samples;

        refined_patch_samples->bdry()            = _refined_surface;
        vector<int> refined_patch_partition(refined_face_map->patches().size(),this->mpiRank());
        refined_patch_samples->patch_partition() = refined_patch_partition;
        //refined_patch_samples->patch_partition() = this->patch_partition();

        refined_patch_samples->setup(); //refined=true for setup

  } else {

      // DZ the LL here is used in PatchSamples interpolate_data for density interpolation
      // the same LL is used for near evaluation interpolation -- these two
      // should be separate

      this->_refined_patch_samples = new PatchSamples("", "");
      refined_patch_samples = this->_refined_patch_samples;

      refined_patch_samples->bdry()            = this->bdry();
      refined_patch_samples->patch_partition() = this->patch_partition();

      refined_patch_samples->setup(true); //refined=true for setup
  }
  cerr<<"bis refined_patch_samples setup"<<endl;

  //---------------------------------------------------------------------------
  // Set up Petsc items required for GMRES. This includes a Matrix-vector
  // multiply operator (_mat) a KSP object to perform the GMRES (_ksp), and a
  // preconditioner wrapper (pc)
  Vec empty;
  VecCreateMPI(this->mpiComm(), 1, PETSC_DETERMINE, &empty);


  _reference_solution = new SolutionDensity(
          empty,                                        /* current solution vector*/
          patch_samples->local_num_sample_points(),     /* # of points on this processes */
          patch_samples->global_num_sample_points(),    /* # of points on all processes */
          this->bdry()->boundary_component_center().size(), /* # of poles/multiply connected comps of boundary */
          this->dom(),                                  /* 0 = bounded domain; 1 = unbounded domain */
          this->_problem_kernel                         /* kernel of the problem of interest */
      );
  // MEMORY LEAK
  //VecDestroy(&empty);


  // local_density_total_dofs = local total number of dof for sample points 
  // (sdof*num_local_sample_points)
  // local_pole_total_dofs = local total number of dof corresponding to 
  // connected components
  // (alpha_i, beta_i in the paper) 
  // global_density_total_dofs, global_pole_total_dofs -- corresponding global 
  // numbers, local_total_dofs = local_density_total_dofs+local_pole_total_dofs
  //
  int local_density_total_dofs, local_pole_total_dofs, local_total_dofs;
  int global_density_total_dofs, global_pole_total_dofs, global_total_dofs;

  _reference_solution->localSize(
          local_density_total_dofs,
          local_pole_total_dofs,
          local_total_dofs);

  _reference_solution->globalSize(
          global_density_total_dofs,
          global_pole_total_dofs,
          global_total_dofs);

  PetscBLASInt lablas = local_density_total_dofs;
  PetscBLASInt gbblas = global_pole_total_dofs; 

  MatCreateShell(this->mpiComm(),
          local_total_dofs, 
          local_total_dofs, 
          global_total_dofs, 
          global_total_dofs,
          (void*)this,
          &_mat);

  MatShellSetOperation(_mat, MATOP_MULT,
              (void(*)(void))SolverGMRESDoubleLayer::mmultWrapper); 

  KSPCreate(this->mpiComm(), &_ksp);
  KSPSetOperators(_ksp,  _mat,  _mat);

  char tmp[100];
  strcpy(tmp, this->prefix().c_str());

  KSPSetOptionsPrefix(_ksp, tmp);
  KSPSetFromOptions(_ksp);

  PC pc;
  KSPGetPC(_ksp, &pc);
  PCType pctype;

  PCGetType(pc, &pctype);
  ebiAssert( strcmp(pctype, "shell")==0 );
  void* ctx =  this;

  PCSetType(pc, PCSHELL);
  PCShellSetApply(pc, SolverGMRESDoubleLayer::pcmultWrapper);
  PCShellSetContext(pc, ctx);

  // list of poles in each domain (if multiply connected)
  // boundary_component_center is the list of interior points {z_m} at the
  // center of the volumes enclosed by the boundary componenets (see derivation
  // of equation (7) in paper (top of page 8))
  
  vector<Point3>& boundary_component_center =
         (this->bdry())->boundary_component_center();

  _poles.resize(DIM, boundary_component_center.size());

  for(int k=0; k<boundary_component_center.size(); k++)
	 for(int d=0; d<DIM; d++)
		_poles(d,k) = boundary_component_center[k](d);

  // We need to evaluate on surface at each step of the solver, so we need an
  // instance of EvaluatorOnSurface in order to perform the singular quadrature inside
  // GMRES
  if (this->_evaluation_type == SINGULAR_EVAL){
      
      assert(patch_samples != NULL);
      assert(patch_samples->sample_point_3d_position() != NULL);
      assert(patch_samples->sample_as_face_point() != NULL);

      EvaluatorOnSurface* evaluator = new EvaluatorOnSurface("", "");

      evaluator->knl() = Kernel3d(_equation_type + DOUBLE_LAYER + VAR_U, this->eqcoefs());
      evaluator->target_3d_position()  = patch_samples->sample_point_3d_position();
      evaluator->target_as_face_point() = patch_samples->sample_as_face_point();
      evaluator->set_surface_discretization(patch_samples);
      
      evaluator->setFromOptions();

      cout << "before solver on-surface evaluator setup" << endl;
      evaluator->setup();
      _on_surface_evaluator = (Evaluator*) evaluator;

  } else if (this->_evaluation_type == INTERIOR_EXTRAPOLATION){

      VecCreateMPI(this->mpiComm(),
              patch_samples->local_num_sample_points(),
              PETSC_DETERMINE,
              &_target_in_out);


      VecDuplicate(patch_samples->sample_point_normal(), &_interpolation_directions);
      VecCopy(patch_samples->sample_point_normal(), _interpolation_directions);
      double minus_one = -1.;
      VecScale(_interpolation_directions, minus_one);

      double one = 1.;
      //VecSet(_target_in_out, 3*one);
    VecSet(_target_in_out, 4);

      Vec on_surface_samples = patch_samples->sample_point_3d_position();
      Vec on_surface_face_points = patch_samples->sample_as_face_point();
      EvaluatorQBKIX* evaluator = 
      new EvaluatorQBKIX(
              Kernel3d(_equation_type + DOUBLE_LAYER + VAR_U, this->eqcoefs()),
              patch_samples,
              refined_patch_samples,
              on_surface_samples,           //targets
              on_surface_samples,           //closest samples == targets in solve
              on_surface_face_points,       //targets as face points
              _target_in_out,
              _interpolation_directions,
              patch_samples->sample_point_far_field(),
              patch_samples->sample_point_interpolant_spacing()
              );
      cout << "before qbkix interior only evaluator setup" << endl;
      evaluator->setup();
      _on_surface_evaluator = (Evaluator*) evaluator;

  } else if (this->_evaluation_type == EXTRAPOLATION_AVERAGE){

      VecCreateMPI(this->mpiComm(),
              patch_samples->local_num_sample_points(),
              PETSC_DETERMINE,
              &_target_in_out);


      VecDuplicate(patch_samples->sample_point_normal(), &_interpolation_directions);
      VecCopy(patch_samples->sample_point_normal(), _interpolation_directions);
      double minus_one = -1.;
      VecScale(_interpolation_directions, minus_one);

      double one = 1.;
      VecSet(_target_in_out, 3*one);

      Vec on_surface_samples = patch_samples->sample_point_3d_position();
      Vec on_surface_face_points = patch_samples->sample_as_face_point();
      EvaluatorQBKIX* evaluator = 
      new EvaluatorQBKIXAverage(
              Kernel3d(_equation_type + DOUBLE_LAYER + VAR_U, this->eqcoefs()),
              patch_samples,
              refined_patch_samples,
              on_surface_samples,           //targets
              on_surface_samples,           //closest samples == targets in solve
              on_surface_face_points,       //targets as face points
              _target_in_out,
              _interpolation_directions,
              _patch_samples->sample_point_far_field(),
              _patch_samples->sample_point_interpolant_spacing()
              );
      cout << "before qbkix averaged evaluator setup" << endl;
      evaluator->setup();
      _on_surface_evaluator = (Evaluator*) evaluator;

  } else if (this->_evaluation_type == INTERPOLATION_ACROSS_SURFACE){

      VecCreateMPI(this->mpiComm(),
              patch_samples->local_num_sample_points(),
              PETSC_DETERMINE,
              &_target_in_out);


      VecDuplicate(patch_samples->sample_point_normal(), &_interpolation_directions);
      VecCopy(patch_samples->sample_point_normal(), _interpolation_directions);
      double minus_one = -1.;
      VecScale(_interpolation_directions, minus_one);

      VecSet(_target_in_out, 3.);

      Vec on_surface_samples = patch_samples->sample_point_3d_position();
      Vec on_surface_face_points = patch_samples->sample_as_face_point();
      EvaluatorQBKIX* evaluator = 
      new EvaluatorNearInterpolate(
              Kernel3d(_equation_type + DOUBLE_LAYER + VAR_U, this->eqcoefs()),
              patch_samples,
              refined_patch_samples,
              on_surface_samples,           //targets
              on_surface_samples,           //closest samples == targets in solve
              on_surface_face_points,       //targets as face points
              _target_in_out,
              _interpolation_directions);

      cout << "before extrapolation evaluator setup" << endl;
      evaluator->setup();
      _on_surface_evaluator = (Evaluator*) evaluator;

  }
 
  //-----------------------
  // local_density_total_dofs = local density total number of dof, local_pole_total_dofs = additional dof from terms
  // corresponding to multiple boundary components (alphas and betas in 7-9
  // on p 253; local_total_dofs = local_density_total_dofs+local_pole_total_dofs
  // The complete system  matrix has the form: Ax, where 
  // x =  [ phi alpha_1 .. alpha_m beta_1 .. beta_m]^T
  // and
  // A = [(1/2)I + D (+ N for unbounded), S_1, ..., S_m, R_1, .., R_m]
  //     [ Tc                             0,   ...		          0  ]  
  //     [ Rc                             0,   ...                0  ]
  // where the upper-left block is evaluated using FMM,  S_1... R_m  block is stored _pole_matrix
  // and  [Tc Rc] block, corresponding to translation and rotation constraints in 7, 9 is stored in _constraint_matrix
  
  int local_num_sample_points = patch_samples->local_num_sample_points();

  if( this->dom()==DOM_UNBND || (this->dom()==DOM_BND && _poles.n()>1) ) {
      //use extra stuff 
     
     // if the domain is unbounded, then each interior point is inside the
     // volumes enclosed by the boundary components of the 
     // solution domain (Bottom of Page 7-top of page 8 of A High-order 3D 
     // Boundary solver... Ying, Biros & Zorin). If the domain is bounded, then 
     // by convention, the first point (0th point in _poles) is contained 
     // inside the solution domain, and the each of the remaining points are 
     // again contained in the volumes of the boundary components as before.
	 
     // Starting index of the interior points inside the boundary components

    int tdof = _problem_kernel.get_tdof();
    int pdof = _problem_kernel.get_pdof();
    int gdof = _problem_kernel.get_gdof();
    int rdof = _problem_kernel.get_rdof();

    // FIXME this is silly
    int bnd_domain_offset = (this->dom()==DOM_UNBND) ? 0 : 1;
    int num_components = _poles.n();  

     // number of boundary components, one interior point at the center of
     // each
     int num_poles = num_components-bnd_domain_offset;

	 double* parr; iC( VecGetArray(patch_samples->sample_point_3d_position(), &parr) );
	 double* narr; iC( VecGetArray(patch_samples->sample_point_normal(), &narr) );
	 double* warr; iC( VecGetArray(patch_samples->sample_point_combined_weight(), &warr) );
	 double* tarr; iC( VecGetArray(patch_samples->sample_as_face_point(), &tarr) );
	 double* iarr; iC( VecGetArray(patch_samples->sample_point_props(), &iarr) ); 

	 // _pole_matrix = [g r] =  [S_1.. S_m, R_1 .. R_m]
	 // pole_matrix 
	 MatCreateDense(this->mpiComm(),
             local_density_total_dofs,
             local_pole_total_dofs,
             global_density_total_dofs,
             global_pole_total_dofs,
             NULL,
             &_pole_matrix);
	 { 
        
		double* barr; 
        MatDenseGetArray(_pole_matrix, &barr);
		DblNumMat b(tdof*local_num_sample_points, pdof*num_poles, false, barr);
		DblNumMat current_pole_position(DIM,num_poles, false, _poles.clmdata(bnd_domain_offset));
        DblNumMat sample_point_3d_position(DIM,local_num_sample_points, false, parr);
        DblNumMat S_matrix, R_matrix;
        // S_matrix = [S_1, .. S_m]; R_matrix = [R_1, .. R_m]
	
        S_matrix = DblNumMat(tdof*local_num_sample_points, gdof*num_poles, false, b.data());
        iC( single_layer_kernel.kernel(current_pole_position, current_pole_position, 
                    sample_point_3d_position, S_matrix) );
        if (rdof > 0){
            // Indicates equation_type requires a Rotlet component
            R_matrix = DblNumMat(tdof*local_num_sample_points, rdof*num_poles, 
                    false, b.clmdata(num_poles*gdof) );

            iC( rotlet_kernel.kernel(current_pole_position, current_pole_position,
                        sample_point_3d_position, R_matrix) );

        } 

		iC( MatDenseRestoreArray(_pole_matrix, &barr) );
		iC( MatAssemblyBegin(_pole_matrix, MAT_FINAL_ASSEMBLY) );
		iC( MatAssemblyEnd(  _pole_matrix, MAT_FINAL_ASSEMBLY) );
	 }

     // encodes constraints in equations 7 and 9
     MatCreateDense(this->mpiComm(),
             local_density_total_dofs,
             local_pole_total_dofs,
             global_density_total_dofs,
             global_pole_total_dofs,
             NULL,
             &_constraint_matrix);
     {
         double* carr; 
         MatDenseGetArray(_constraint_matrix, &carr);
         DblNumMat c(tdof*local_num_sample_points, pdof*num_poles, false, carr);

#pragma omp parallel for
         for(int ui=0; ui<local_num_sample_points; ui++) {

             PatchSamples::Tag* tag = (PatchSamples::Tag*)(iarr+ui);
             // DZ RENAME component_id
             int cid = int(tag->_gid);
             //cout <<"cid: " << cid << endl;

             if(cid>=bnd_domain_offset) { //good
                 int p = cid - bnd_domain_offset;
                 // current_pole_position is the position of cid-th pole
                 DblNumVec current_pole_position(DIM, false, _poles.clmdata(cid));

                 // interface to a vector of sdof values corresponding to ui-th
                 // sample point in d->sample_point_3d_position			 
                 DblNumVec current_sample_point_3d_position(DIM, false, parr + ui*DIM);

                 ebiAssert(tdof == gdof);
                 DblNumMat id_matrix = Kernel3d::get_identity_matrix(tdof);
                 // note id_matrix is square by assertion

                 //\phi_k(c_k) = w_k*\varphi(g_k(c_k))*J_k
                 // warr is a pointer to wcb and 
                 // wcb[i] = jacobian * alpha *weight
                 // Note alpha = value of POU

                 for(int i=0; i<tdof; i++)
                     for(int j=0; j<gdof; j++){
                         c( ui*tdof+i, p*gdof+j ) = id_matrix(i,j) * warr[ui];
                     }

                 // this is 3x3 cross-product matrix only used for
                 // vector kernels which have the additional constraint
                 // in eqs 7 and 9
                 //ebiAssert(tdof==3 && rdof==3); // <--This might be redundant
                 DblNumMat cross_prod_matrix = 
                     Kernel3d::get_cross_product_matrix(current_pole_position, 
                             current_sample_point_3d_position);

                 // second block row of c is the row of cross-product
                 // matrices used to compute R \int (y-z_m) x phi(y)
                 // ds(y) in eqs 7 and 9
                 //
                 // Note that in the case of kernels with no rotlet part, rdof = 0 and
                 // so this loop never executes.
                 for(int i=0; i<tdof; i++)			 
                     for(int j=0; j<rdof; j++)				
                         c( ui*tdof+i, num_poles*gdof + p*rdof +j) =
                             cross_prod_matrix(i,j) * warr[ui];
             }
         }
         iC( MatDenseRestoreArray(_constraint_matrix, &carr) );
         iC( MatAssemblyBegin(_constraint_matrix, MAT_FINAL_ASSEMBLY) );
         iC( MatAssemblyEnd(  _constraint_matrix, MAT_FINAL_ASSEMBLY) );
     }
     // (constraint_matrix^T*pole_matrix)^{-1}

     // DZ RENAME ictb  precond_Schur_complement_inverse
     MatCreateDense(this->mpiComm(), 
             local_pole_total_dofs, 
             local_pole_total_dofs, 
             global_pole_total_dofs, 
             global_pole_total_dofs, 
             NULL, 
             &precond_schur_complement_inverse);
     {
         DblNumMat buf(global_pole_total_dofs, global_pole_total_dofs);
         double* pole_matrix_ptr; 
         MatDenseGetArray(_pole_matrix, &pole_matrix_ptr);
         double* constraint_matrix_ptr; 
         MatDenseGetArray(_constraint_matrix, &constraint_matrix_ptr);

		DblNumMat b(local_density_total_dofs, global_pole_total_dofs, false, pole_matrix_ptr);
		DblNumMat c(local_density_total_dofs, global_pole_total_dofs, false, constraint_matrix_ptr);

		{
		  char transa = 't';
		  char transb = 'n';
		  double alpha = 1.0;
		  double beta = 0.0;

          // Compute precond_schur_complement_inverse = _pole_matrix^T*_constraint_matrix
          // Note: lablas = local_density_total_dofs := number of locally held sample points
          //       gbblas = global_pole_total_dofs := number of global sample points (???)
		  DGEMM(&transa, &transb, &gbblas, &gbblas, &lablas, &alpha, b.data(), 
                  &lablas, c.data(), &lablas, &beta, buf.data(), &gbblas);
		}
		double* precond_ptr = NULL;
		// gather the matrix on a single processor:
		// it is a product of a wide by tall distributed matrices,
		// of size num_poles^2*sdof^2, i.e. small, so stored and
		// processed on the 0th processor  
		if(this->mpiRank()==0) {
		  MatDenseGetArray(precond_schur_complement_inverse, &precond_ptr);
		}

		MPI_Reduce(buf.data(),
                precond_ptr,
                global_pole_total_dofs*global_pole_total_dofs,
                MPI_DOUBLE,
                MPI_SUM,
			    0,
                this->mpiComm());

		// compute LU and invert  preconditioner Schur compl. matrix
		if(this->mpiRank()==0) {
		  // used by lapack for pivot elements
		  PetscBLASInt* ipiv = new PetscBLASInt[global_pole_total_dofs];
		  PetscBLASInt info;

		  //Compute LU Factorization of precond_schur_complement_inverse 
          //(L + U stored in precond_schur_complement_inverse on termination)
		  DGETRF(&gbblas, &gbblas, precond_ptr, &gbblas, ipiv, &info);
		  //ebiAssert(info==0); //MJM: FIXME
          
		  PetscBLASInt lwork = global_pole_total_dofs;
		  double* work = new double[lwork];

		  // Compute the inverse of precond_schur_complement_inverse 
          // (inverse restored in precond_schur_complement_inverse on exit)
		  DGETRI(&gbblas, precond_ptr, &gbblas, ipiv, work, &lwork, &info);

          cout << "info: " << info << endl;
	  // DZ check -- (_constraint_matrix^T*_pole_matrix)^{-1} or as below? compare to the
          // Greengard et al paper
          // In summary precond_schur_complement_inverse = 
          // (_pole_matrix^T*_constraint_matrix)^{-1} upon completion of set up
		  // ebiAssert(info==0);//MJM: FIXME
          
		  delete [] ipiv;
		  delete [] work;
		}
		if(this->mpiRank()==0) {
		  MatDenseRestoreArray(precond_schur_complement_inverse, &precond_ptr);

		}
		MatAssemblyBegin(precond_schur_complement_inverse, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(  precond_schur_complement_inverse, MAT_FINAL_ASSEMBLY);
		MatDenseRestoreArray(_pole_matrix, &pole_matrix_ptr);
		MatDenseRestoreArray(_constraint_matrix, &constraint_matrix_ptr);
	 }
	 VecRestoreArray(patch_samples->sample_point_3d_position(), &parr);
	 VecRestoreArray(patch_samples->sample_point_normal(),&narr); 
	 VecRestoreArray(patch_samples->sample_point_combined_weight(), &warr);
	 VecRestoreArray(patch_samples->sample_as_face_point(), &tarr);
	 VecRestoreArray(patch_samples->sample_point_props(), &iarr);
  }
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::solve(Vec b0, Vec x)
{
  // Copies around arrays and calls Petsc's KSPSolve function
  ebiFunctionBegin;
  
  int local_density_total_dofs, local_pole_total_dofs, local_total_dofs;
  
  _reference_solution->localSize(
          local_density_total_dofs,
          local_pole_total_dofs,
          local_total_dofs);

  int64_t tmp;
  iC( VecGetLocalSize(b0, &tmp) );
  cout << "rhs dofs: " << tmp << endl;
  cout << "required dofs: " << local_density_total_dofs << endl;
  ebiAssert(tmp==local_total_dofs);

  Vec b;
  iC( VecCreateMPI(this->mpiComm(), local_total_dofs, PETSC_DETERMINE, &b) );


  double zero = 0;
  iC( VecSet( b, zero) );

  double* barr;
  double* b0arr;
  
  iC( VecGetArray(b, &barr) );
  iC( VecGetArray(b0, &b0arr) );

  for(int i=0; i<local_density_total_dofs; i++)
      barr[i] = b0arr[i]; //copy b0 to b

  iC( VecRestoreArray(b, &barr) );
  iC( VecRestoreArray(b0, &b0arr) );

    double start = omp_get_wtime();
  iC( KSPSolve(_ksp, b, x) );
  stats.add_result("solve time", (omp_get_wtime() - start) );
  iC( VecDestroy(&b) );

  stats.add_result("number of GMRES iterations", this->nits());

  ebiFunctionReturn(0);
}


// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::mmult(SolutionDensity& solution_in, SolutionDensity& solution_out)
{
  // Copies around arrays and calls the built-in matrix-multiply routine.
  ebiFunctionBegin;
  
  double zero = 0; iC( VecSet(solution_out.solution(), zero) );

  int local_density_total_dofs, local_pole_total_dofs, local_total_dofs; 
  int global_density_total_dofs, global_pole_total_dofs, global_total_dofs; 

  solution_in.localSize( 
          local_density_total_dofs,
          local_pole_total_dofs,
          local_total_dofs);
  
  solution_in.globalSize(
          global_density_total_dofs,
          global_pole_total_dofs,
          global_total_dofs);

  Vec input_density, input_pole_coeffs, result_density, result_pole_coeffs;
  
  input_density     = solution_in.density();
  input_pole_coeffs = solution_in.pole_coeffs();
  result_density = solution_out.density();
  result_pole_coeffs = solution_out.pole_coeffs();

  //DL
  ebiAssert(_on_surface_evaluator!=NULL); // this assertion can probably be deleted
  cout << "before on-surface solver eval" << endl;
  iC( _on_surface_evaluator->eval(input_density, result_density) );
  cout << "after on-surface solver eval" << endl;

  if (this->_is_incompressible){
      PatchSamples* d = (this->patch_samples());

      if(this->dom()==DOM_BND) {
          int sdof = _problem_kernel.get_sdof();
          Vec density_intermed;
           
          VecCreateMPI(this->mpiComm(),
                   local_density_total_dofs,
                   global_density_total_dofs,
                   &density_intermed);
          
           denscale(sdof, d->sample_point_combined_weight(), input_density, density_intermed);

          double val;
          VecDot(d->sample_point_normal(), density_intermed, &val);
          VecAXPY( result_density, val,  d->sample_point_normal());
          VecDestroy(&density_intermed);
      }
  }
  //HALF
  double half = 0.5;
  iC( VecAXPY( result_density, half,  input_density) );

  //EXTRA
  if( this->dom()==DOM_UNBND || (this->dom()==DOM_BND && _poles.n()>1) ) {
     
      // Compute result_density := result_density + _b*input_pole_coeffs, where result_density is initialized to zeros
      // Recall that _pole_matrix is the interaction matrix of the poles (interior points of 
      // boundary components) and the sample points
	 iC( MatMultAdd(_pole_matrix, input_pole_coeffs, result_density, result_density) );

     // Compute result_pole_coeffs := _constraint_matrix^T*input_density
	 iC( MatMultTranspose(_constraint_matrix, input_density, result_pole_coeffs) );
  }
  
  this->nits()++;
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// Applying preconditioner from GREENGARD, KROPINSKI, AND MAYO
// "Integral Equation Methods for Stokes Flow and Isotropic Elasticity in the Plane"
//  

int SolverGMRESDoubleLayer::pcmult(SolutionDensity& solution_in, SolutionDensity& solution_out)
{
  ebiFunctionBegin;
  double zero = 0;
  iC( VecSet( solution_out.solution(), zero) );
  
  if( this->dom()==DOM_UNBND || (this->dom()==DOM_BND && _poles.n()>1) ) { 
	 // system matrix:
	 // x = in, y = out; matvec:  y = Ax
	 // y = [ya yb], x = [xa, xb];
	 // A = [A_aa, A_ab; A_ba Abb]
	 // A_aa = potential computation, implemented in FMM, A_ab is pole
	 // summation, A_ba = constraints, A_bb = 0

      int local_density_total_dofs, local_pole_total_dofs, local_total_dofs;
      int global_density_total_dofs, global_pole_total_dofs, global_total_dofs;
      
      solution_in.globalSize(
              global_density_total_dofs,
              global_pole_total_dofs,
              global_total_dofs);

      solution_in.localSize(
              local_density_total_dofs,
              local_pole_total_dofs,
              local_total_dofs);

	 
      Vec input_density, input_pole_coeffs, result_density, result_pole_coeffs;
      input_density = solution_in.density();
      input_pole_coeffs = solution_in.pole_coeffs();
      result_density = solution_out.density();
      result_pole_coeffs = solution_out.pole_coeffs();
     
      Vec density_intermed, pole_coeffs_intermed;
	 // vectors for densities and pole variables
	 VecCreateMPI(this->mpiComm(),
             local_density_total_dofs,
             global_density_total_dofs,
             &density_intermed); 

	 VecCreateMPI(this->mpiComm(), 
             local_pole_total_dofs, 
             global_pole_total_dofs, 
             &pole_coeffs_intermed);


	 // y_b = -1/2*Z*(x_b - 2 _constraint_matrix^T x_a)
	 // y_a = 2*x_a + _b*y_b
	 
	 //1...
	 // pole_coeffs_intermed = x_b - 2 _constraint_matrix^T x_a
	 // density_intermed = x_a;
	 
	 // pole_coeffs_intermed =  _constraint_matrix^T x_a	 
	 double half = 0.5;
	  VecCopy(input_density, density_intermed);
	 MatMultTranspose(_constraint_matrix, input_density, pole_coeffs_intermed);
	 
	 double movh = -1.0/half;

	 // pole_coeffs_intermed = -2*pole_coeffs_intermed
	 VecScale( pole_coeffs_intermed, movh);

	 // pole_coeffs_intermed = x_b + fbb
	 double one = 1.0;
	 VecAXPY( pole_coeffs_intermed, one,  input_pole_coeffs);

	 //2...
	 // y_b = -1/2*Z*pole_coeffs_intermed
	 // y_a = 2*density_intermed + _b*y_b
	 
	 MatMult(precond_schur_complement_inverse, pole_coeffs_intermed, result_pole_coeffs);
	 // y_a = density_intermed
	 VecCopy(density_intermed, result_density);

	 double ovh = 1.0/half;
	 // y_a = 2*y_a
	 VecScale( result_density, ovh);

	 // y_a = y_a + _b * y_b
	 MatMultAdd(_pole_matrix, result_pole_coeffs, result_density, result_density);

	 double mh = - half;
	 iC( VecScale( result_pole_coeffs, mh) );
	 // y_b = -1/2*y_b
	 
     iC( VecDestroy(&density_intermed) );
	 iC( VecDestroy(&pole_coeffs_intermed) );
  } else {
	 iC( VecCopy(solution_in.solution(), solution_out.solution()) );
  }
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::mmultWrapper(Mat M, Vec in, Vec ou)
{
  ebiFunctionBegin;
  SolverGMRESDoubleLayer* bis = NULL;

  iC( MatShellGetContext(M, (void**)&bis) );
  ebiAssert(bis!=NULL);

  // NOTE TO FUTURE SELF: Petsc Vec's in and ou are restored and sub-Vec's are 
  // appropriately freed in ~SolutionDensity() and is implicitly called when 
  // solution_in/solution_out go out of scope (i.e. at the end of this function body).
  // Be aware of this in the case of shuffling around object instantiations; 
  // these Vec's must be properly restored and freed at the end of the wrapper
  
  SolutionDensity solution_in(
          in,                                              /* current solution vector*/
          bis->patch_samples()->local_num_sample_points(),          /* # of points on this processes */
          bis->patch_samples()->global_num_sample_points(),         /* # of points on all processes */
          bis->bdry()->boundary_component_center().size(), /* # of poles/multiply connected comps of boundary */
          bis->dom(),                                      /* 0 = bounded domain; 1 = unbounded domain */
          bis->_problem_kernel                             /* kernel of the problem of interest */
      );

  SolutionDensity solution_out(
          ou,                                              /* current solution vector*/
          bis->patch_samples()->local_num_sample_points(),          /* # of points on this processes */
          bis->patch_samples()->global_num_sample_points(),         /* # of points on all processes */
          bis->bdry()->boundary_component_center().size(), /* # of poles/multiply connected comps of boundary */
          bis->dom(),                                      /* 0 = bounded domain; 1 = unbounded domain */
          bis->_problem_kernel                             /* kernel of the problem of interest */
      );

  iC( bis->mmult(solution_in, solution_out) );
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::pcmultWrapper(PC pc, Vec in, Vec ou) {
  ebiFunctionBegin;
  // NOTE TO FUTURE SELF: Petsc Vec's in and ou are restored and sub-Vec's are 
  // appropriately freed in ~SolutionDensity() and is implicitly called when 
  // solution_in/solution_out go out of scope (i.e. at the end of this function body).
  // Be aware of this in the case of shuffling around object instantiations; 
  // these Vec's must be properly restored and freed at the end of the wrapper
  
  SolverGMRESDoubleLayer* bis = NULL;
  PCShellGetContext(pc, (void**)&bis);
  ebiAssert(bis!=NULL);

  SolutionDensity solution_in(
          in,                                              /* current solution vector*/
          bis->patch_samples()->local_num_sample_points(),          /* # of points on this processes */
          bis->patch_samples()->global_num_sample_points(),         /* # of points on all processes */
          bis->bdry()->boundary_component_center().size(), /* # of poles/multiply connected comps of boundary */
          bis->dom(),                                      /* 0 = bounded domain; 1 = unbounded domain */
          bis->_problem_kernel                             /* kernel of the problem of interest */
      );

  SolutionDensity solution_out(
          ou,                                              /* current solution vector*/
          bis->patch_samples()->local_num_sample_points(),          /* # of points on this processes */
          bis->patch_samples()->global_num_sample_points(),         /* # of points on all processes */
          bis->bdry()->boundary_component_center().size(), /* # of poles/multiply connected comps of boundary */
          bis->dom(),                                      /* 0 = bounded domain; 1 = unbounded domain */
          bis->_problem_kernel                             /* kernel of the problem of interest */
      );


  iC( bis->pcmult(solution_in, solution_out) );
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 


void SolverGMRESDoubleLayer::evaluate(Vec target_points, 
        Vec density,
        Vec& potential,
        NumVec<OnSurfacePoint>& on_surface_points,
        Kernel_variable kernel_variable,
        bool precomputed_closest_points){
    int num_local_targets = num_local_points(target_points);
    int num_local_density_values = Petsc::get_vec_local_size(density)/_problem_kernel.srcDOF();
    int num_local_potential_values = Petsc::get_vec_local_size(potential)/_problem_kernel.trgDOF();
    cout << "num_local_targets: " << num_local_targets << endl; 
    cout << "num_local_density_values: " << num_local_density_values<< endl; 
    cout << "num_local_potential_values: " << num_local_potential_values<< endl; 
    assert(num_local_targets== on_surface_points.m());
    assert(num_local_density_values == _patch_samples->local_num_sample_points());
    assert(num_local_potential_values == num_local_targets);
    

    // mark points as near/far in/out
    DblNumMat target_points_local = get_local_vector(DIM, num_local_targets, target_points);
    cout << "marking" << endl;
    if(!precomputed_closest_points){
        on_surface_points = 
            Markgrid::mark_target_points(
                    target_points_local, 
                    dynamic_cast<PatchSurfFaceMap*>(_bdry));
    }
    cout << "marked" << endl;
    
    Vec target_in_out;
    Vec closest_on_surface_points_3d_position;
    Vec closest_on_surface_points_as_face_point;
    Vec target_far_field;
    Vec target_interpolant_spacing;
    Vec interpolation_directions;


    // make the necessary petsc Vec's to pass to the legacy evaluator interface
    Petsc::create_mpi_vec(this->mpiComm(), num_local_targets, target_in_out);
    Petsc::create_mpi_vec(this->mpiComm(), num_local_targets*DIM, 
            closest_on_surface_points_3d_position);
    Petsc::create_mpi_vec(this->mpiComm(), 
            num_local_targets*_patch_samples->bdry()->face_point_size_in_doubles(), 
              closest_on_surface_points_as_face_point);
    Petsc::create_mpi_vec(this->mpiComm(), num_local_targets*DIM, interpolation_directions);
    VecDuplicate( target_in_out, &target_far_field);
    VecDuplicate( target_in_out, &target_interpolant_spacing);
    cout << "created vecs" << endl;

    // set qbkix related junk TODO remove these
    VecSet(target_far_field, Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    VecSet(target_interpolant_spacing, Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    VecSet(interpolation_directions, 0.);
    


    DblNumMat target_in_out_local = get_local_vector(1, num_local_targets, target_in_out);
    DblNumMat target_far_field_local = get_local_vector(1, num_local_targets, target_far_field);
    DblNumMat target_interpolant_spacing_local = get_local_vector(1, num_local_targets, target_interpolant_spacing);
    DblNumMat interpolation_directions_local = get_local_vector(DIM, num_local_targets, interpolation_directions);

    DblNumMat closest_on_surface_points_3d_position_local = 
        get_local_vector(DIM, num_local_targets, closest_on_surface_points_3d_position);

    DblNumMat closest_on_surface_points_as_face_point_local = 
        get_local_vector(_patch_samples->bdry()->face_point_size_in_doubles(), 
                num_local_targets, closest_on_surface_points_as_face_point);

    cout << "get local vecs" << endl;
    // make index sets for near and far targets based on marking
    vector<int64_t> near_indices;
    vector<int64_t> far_indices;

    // compute the relevant data needed for each piece of evaluation
    string qbkix_convergence_type = Options::get_string_from_petsc_opts("-qbkix_convergence_type");
    bool qbkix_classical_conv = qbkix_convergence_type == "classic";
    bool qbkix_adaptive_conv = qbkix_convergence_type == "adaptive";
    
    cout << "setup up target data" << endl;
    for(int i =0; i < num_local_targets; i++){
        OnSurfacePoint p = on_surface_points(i);

        // this is free, I don't think we need it for far eval, but why not...
        if(p.inside_domain == OUTSIDE | p.inside_domain == OUTSIDE_SPATIAL_GRID){
            target_in_out_local(0,i) = 2; // dangerous... see common/NbrDefs.hpp to fix
        } else if (p.inside_domain == INSIDE){
            target_in_out_local(0,i) = 1;
        }

        //  keep track of indices associated with each near/far point. if the
        //  point is near, compute the needed data to pass to the legacy
        //  interface.
        if(p.region == NEAR){
            near_indices.push_back(i);

            FacePointFaceMap* face_point = 
                (FacePointFaceMap*) closest_on_surface_points_as_face_point_local.clmdata(i);
            *face_point = FacePointFaceMap::to_face_point(p);

            Point3 closest_on_surface_point(closest_on_surface_points_3d_position_local.clmdata(i));
            Point3 positions_and_derivs[3];
            FaceMapSubPatch* patch = 
                (FaceMapSubPatch*)_patch_samples->bdry()->patches()[p.parent_patch];

            patch->xy_to_patch_coords(p.parametric_coordinates.array(), 
                    PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, 
                    (double*)positions_and_derivs);
                    //closest_on_surface_point.array());
                    Point3 normal = -cross(positions_and_derivs[1], positions_and_derivs[2]).dir(); // interior normal!! note minus sign
            for(int d =0; d< DIM; d++){
                closest_on_surface_points_3d_position_local(d,i) = positions_and_derivs[0](d);
                interpolation_directions_local(d,i) = normal(d);
            }
            if(qbkix_classical_conv){
                target_far_field_local(0,i) *= sqrt(patch->characteristic_length());
                target_interpolant_spacing_local(0,i) *= sqrt(patch->characteristic_length());
            } else if(qbkix_adaptive_conv){
                target_far_field_local(0,i) *= patch->characteristic_length();
                target_interpolant_spacing_local(0,i) *= patch->characteristic_length();
            }

        } else if(p.region == FAR){
            far_indices.push_back(i);
        }

    }
    cout << "finished setup up target data" << endl;

    
    closest_on_surface_points_3d_position_local.restore_local_vector();
    closest_on_surface_points_as_face_point_local.restore_local_vector();
    target_in_out_local.restore_local_vector();
    target_points_local.restore_local_vector();
    target_far_field_local.restore_local_vector();
    target_interpolant_spacing_local.restore_local_vector();
    interpolation_directions_local.restore_local_vector();


    cout << "eval subvec" << endl;
    // initialize petsc index sets and get subvectors for points/potential for
    // near/far evaluation
    if( !near_indices.empty()){
        evaluate_subvector(
                // global arrays for the whole problem
                target_points, 
                density,
                potential,
                closest_on_surface_points_3d_position,
                closest_on_surface_points_as_face_point,
                interpolation_directions,
                target_far_field,
                target_interpolant_spacing,
                target_in_out,
                on_surface_points,
                // local indices for the current subvector.
                near_indices,
                kernel_variable,
                INTERIOR_EXTRAPOLATION);
    } 
    if(!far_indices.empty()){
        evaluate_subvector(
                // global arrays for the whole problem
                target_points, 
                density,
                potential,
                closest_on_surface_points_3d_position,
                closest_on_surface_points_as_face_point,
                interpolation_directions,
                target_far_field,
                target_interpolant_spacing,
                target_in_out,
                on_surface_points,
                // local indices for the current subvector.
                far_indices,
                kernel_variable,
                SMOOTH_QUADRATURE);
    }
    cout << "eval-ed subvec" << endl;
    // call near/far evaluation
    // restore vectors
    // finish
    Petsc::destroy_vec(target_in_out);
    Petsc::destroy_vec(closest_on_surface_points_3d_position);
    Petsc::destroy_vec(closest_on_surface_points_as_face_point);
    Petsc::destroy_vec(target_far_field);
    Petsc::destroy_vec(target_interpolant_spacing);
    Petsc::destroy_vec(interpolation_directions);


}


void SolverGMRESDoubleLayer::evaluate_subvector(
        // global arrays for the whole problem
        Vec target_points, 
        Vec density,
        Vec& potential,
        Vec closest_points_3d_position,
        Vec closest_points_as_face_point,
        Vec interpolation_directions,
        Vec target_far_field,
        Vec target_interpolant_spacing,
        Vec target_in_out,
        NumVec<OnSurfacePoint>& on_surface_points,
        // local indices for the current subvector.
        vector<int64_t>sub_vector_indices,
        Kernel_variable kernel_variable,
        EvaluationType evaluation_type
        ){

    
    Vec target_points_subvec;
    Vec potential_subvec;
    Vec closest_points_3d_position_subvec;
    Vec closest_points_as_face_point_subvec;
    Vec interpolation_directions_subvec;
    Vec target_far_field_subvec;
    Vec target_interpolant_spacing_subvec;
    Vec target_in_out_subvec;

    int64_t num_subvec_targets = sub_vector_indices.size();
    int64_t face_point_size =_patch_samples->bdry()->face_point_size_in_doubles();
    int64_t target_dof = _problem_kernel.trgDOF();
    
    // make new smaller Vec's
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets*DIM, target_points_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets*target_dof, potential_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets*DIM, closest_points_3d_position_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets*face_point_size, closest_points_as_face_point_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets*DIM, interpolation_directions_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets, target_far_field_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets, target_interpolant_spacing_subvec);
    Petsc::create_mpi_vec(this->mpiComm(), num_subvec_targets, target_in_out_subvec);

    copy_values_to_subvec(sub_vector_indices, DIM,              COPY_TO, target_points_subvec, target_points);
    copy_values_to_subvec(sub_vector_indices, target_dof,       COPY_TO, potential_subvec, potential);
    copy_values_to_subvec(sub_vector_indices, DIM,              COPY_TO, closest_points_3d_position_subvec, closest_points_3d_position);
    copy_values_to_subvec(sub_vector_indices, face_point_size,  COPY_TO, closest_points_as_face_point_subvec, closest_points_as_face_point);
    copy_values_to_subvec(sub_vector_indices, DIM,              COPY_TO, interpolation_directions_subvec, interpolation_directions);
    copy_values_to_subvec(sub_vector_indices, 1,                COPY_TO, target_far_field_subvec, target_far_field);
    copy_values_to_subvec(sub_vector_indices, 1,                COPY_TO, target_interpolant_spacing_subvec, target_interpolant_spacing);
    copy_values_to_subvec(sub_vector_indices, 1,                COPY_TO, target_in_out_subvec, target_in_out);


    Evaluator* evaluator;
    Kernel3d kernel(_equation_type + DOUBLE_LAYER + kernel_variable, this->eqcoefs());

    if(evaluation_type == INTERIOR_EXTRAPOLATION){
        evaluator = (Evaluator*)
            new EvaluatorQBKIX(
                    kernel,
                    _patch_samples,
                    _refined_patch_samples,
                    target_points_subvec,
                    closest_points_3d_position_subvec,
                    closest_points_as_face_point_subvec,
                    target_in_out_subvec,
                    interpolation_directions_subvec,
                    target_far_field_subvec,
                    target_interpolant_spacing_subvec);

    } else if (evaluation_type == SMOOTH_QUADRATURE){
        EvaluatorFar*  evaluator_far = new EvaluatorFar("", "");
        evaluator_far->knl() = kernel;
        evaluator_far->target_3d_position() = target_points_subvec;
        evaluator_far->set_surface_discretization(_patch_samples);

        evaluator = (Evaluator*) evaluator_far;
    }

    SolutionDensity solution(*_reference_solution);
    solution.set_petsc_vec(density);
    
    core_evaluation(solution, potential_subvec, kernel_variable, target_points_subvec, evaluator);

    delete evaluator;
    copy_values_to_subvec(sub_vector_indices, target_dof,       COPY_FROM, potential_subvec, potential);

    Petsc::destroy_vec(target_points_subvec);
    Petsc::destroy_vec(potential_subvec);
    Petsc::destroy_vec(closest_points_3d_position_subvec);
    Petsc::destroy_vec(closest_points_as_face_point_subvec);
    Petsc::destroy_vec(interpolation_directions_subvec);
    Petsc::destroy_vec(target_far_field_subvec);
    Petsc::destroy_vec(target_interpolant_spacing_subvec);
    Petsc::destroy_vec(target_in_out_subvec);


}

Vec SolverGMRESDoubleLayer::singular_correction(Vec sources, 
        Vec targets, Vec targets_as_face_points, Vec potential, Kernel3d kernel){
  
  Vec corrected_density;
  VecDuplicate(potential, &corrected_density);
  double zero = 0.;  
  VecSet(corrected_density, zero);

  EvaluatorOnSurface* on_surface_evaluator =
      new EvaluatorOnSurface(this->name()+"RE_", this->prefix()+"re_");
  
  on_surface_evaluator->knl()                   = kernel;
  on_surface_evaluator->target_3d_position()    = targets;
  on_surface_evaluator->target_as_face_point()  = targets_as_face_points;
  on_surface_evaluator->set_surface_discretization(this->patch_samples());

  on_surface_evaluator->setFromOptions();
  on_surface_evaluator->setup_no_fmm();
  

  
  on_surface_evaluator->singular_evaluation(potential, corrected_density);
  on_surface_evaluator->apply_singularity_cancellation(potential, corrected_density);
  
  //---------------------------------------------------------------------------
  delete on_surface_evaluator;
  return corrected_density;
}


int SolverGMRESDoubleLayer::core_evaluation(SolutionDensity& solution_in, Vec out,
        Kernel_variable qt, Vec tp, Evaluator* evaluator) {
  ebiFunctionBegin;

  Kernel3d single_layer_kernel, double_layer_kernel;
  single_layer_kernel = Kernel3d(_equation_type + SINGLE_LAYER + qt, this->eqcoefs());
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + qt, this->eqcoefs());
  
  int tdof = double_layer_kernel.trgDOF(); // MJM BUG for blended surface pressure eval, only trgDOF works instead of get_tdof(). Haven't figured out why.

  //1. separate
  double zero = 0.0;
  iC( VecSet( out, zero) );

  int local_density_total_dofs;
  int local_pole_total_dofs;
  int local_total_dofs;
  int global_density_total_dofs;
  int global_pole_total_dofs;
  int global_total_dofs;

  solution_in.localSize( 
          local_density_total_dofs, 
          local_pole_total_dofs, 
          local_total_dofs);
  
  solution_in.globalSize(
          global_density_total_dofs, 
          global_pole_total_dofs, 
          global_total_dofs);

  Vec input_density, input_pole_coeffs;
  input_density = solution_in.density();
  input_pole_coeffs = solution_in.pole_coeffs();

  iC( evaluator->setFromOptions() );
  iC( evaluator->setup() );
  iC( evaluator->eval(input_density, out) );
  
  //3. dense eval for pole
  if( this->dom()==DOM_UNBND || (this->dom()==DOM_BND && _poles.n()>1) ) {
	 int num_local_pts = num_local_points(tp);
	 int num_global_pts = num_global_points(tp);
    
     //int tdof = double_layer_kernel.get_tdof();
     int pdof = double_layer_kernel.get_pdof();
     int gdof = double_layer_kernel.get_gdof();
     int rdof = double_layer_kernel.get_rdof();

     // TODO make this less silly
	 int bnd_domain_offset = (this->dom()==DOM_UNBND) ? 0 : 1;
     int num_components = _poles.n();
     int num_poles = num_components-bnd_domain_offset;

	 double* parr;
     iC( VecGetArray(tp, &parr) );

     Mat _bn;
     MatCreateDense(this->mpiComm(), 
             tdof*num_local_pts, 
             local_pole_total_dofs, 
             tdof*num_global_pts, 
             global_pole_total_dofs, 
             NULL,
             &_bn);

     double* barr;
     iC( MatDenseGetArray(_bn, &barr) );

     DblNumMat b(tdof*num_local_pts, pdof*num_poles, false, barr);
     DblNumMat current_pole_position(DIM,num_poles,   false, _poles.clmdata(bnd_domain_offset));
     DblNumMat sample_point_3d_position(DIM,num_local_pts, false, parr);

     DblNumMat S_matrix, R_matrix;

     S_matrix = DblNumMat(tdof*num_local_pts, gdof*num_poles, false, b.data());
     iC( single_layer_kernel.kernel(current_pole_position, current_pole_position, 
                 sample_point_3d_position, S_matrix) );
     if (rdof > 0){
         // Equation has a kernel that requires a Rotlet part
         Kernel3d rotlet_kernel(_equation_type +  ROTLET + qt, this->eqcoefs());
         R_matrix = DblNumMat(tdof*num_local_pts, rdof*num_poles, false, b.clmdata(num_poles*gdof) );
         iC( rotlet_kernel.kernel(current_pole_position, current_pole_position, 
                     sample_point_3d_position, R_matrix) );

     } 
     iC( MatDenseRestoreArray(_bn, &barr) );
     iC( MatAssemblyBegin(_bn, MAT_FINAL_ASSEMBLY) );
	 iC( MatAssemblyEnd(  _bn, MAT_FINAL_ASSEMBLY) );
	 iC( MatMultAdd(_bn, input_pole_coeffs, out, out) );
	 iC( MatDestroy(&_bn) );
	 iC( VecRestoreArray(tp, &parr) );
  }
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::fareval(Vec tp, int qt, Vec in, Vec ou)
{
  ebiFunctionBegin;

  Kernel3d double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + qt, this->eqcoefs());

  EvaluatorFar* far_evaluator = new EvaluatorFar("", "");

  far_evaluator->knl()                  = double_layer_kernel;
  far_evaluator->target_3d_position()   = tp;
  far_evaluator->set_surface_discretization(this->patch_samples());
  
  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(in);

  core_evaluation(solution, ou, (Kernel_variable) qt, tp, far_evaluator);

  delete far_evaluator;
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::roneval(Vec tp, Vec tb, int qt, Vec in, Vec ou)
{
  ebiFunctionBegin;
  
  Kernel3d double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + qt, this->eqcoefs());

  EvaluatorOnSurface* on_surface_evaluator = new EvaluatorOnSurface("", "");
  
  on_surface_evaluator->knl()                   = double_layer_kernel;
  on_surface_evaluator->target_3d_position()    = tp;
  on_surface_evaluator->target_as_face_point()  = tb;
  on_surface_evaluator->set_surface_discretization(this->patch_samples());

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(in);

  core_evaluation(solution, ou, (Kernel_variable) qt, tp, on_surface_evaluator);
  delete on_surface_evaluator;
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::neaeval(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorNear* near_evaluator = new EvaluatorNear("", "");

  near_evaluator->knl()                            = double_layer_kernel;
  near_evaluator->target_3d_position()             = near_targets;
  near_evaluator->intermediate_target_3d_position()= intermediate_targets;
  near_evaluator->target_in_out()                  = target_in_out_near;
  near_evaluator->closest_sample_3d_position()     = closest_sample_3d_position;
  near_evaluator->closest_sample_as_face_point()   = closest_sample_as_face_point; 
  near_evaluator->set_surface_discretization(this->patch_samples());
  near_evaluator->set_refined_surface_discretization(this->refined_patch_samples());
  
  
  // MJM TODO figure what this does...
  int64_t near_interpolation_num_samples = (int64_t)(ceil((1.0/this->patch_samples()->spacing()))); 
  /*
  if (near_interpolation_num_samples > 8) {
      near_interpolation_num_samples = this->surface_interpolation_num_samples();
  }
  */

  //MJM TODO make sure this is ok to do
  //near_evaluator->interpolation_nodes() = get_interpolation_nodes();

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorQBKIX* near_evaluator = new EvaluatorQBKIX(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          NULL,
          _patch_samples->sample_point_far_field(),
          _patch_samples->sample_point_interpolant_spacing()
          );

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}

int SolverGMRESDoubleLayer::neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec dist_from_boundary,  //closest on surface sample for interpolation
        Vec interp_point_spacing,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorQBKIX* near_evaluator = new EvaluatorQBKIX(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          interpolation_directions,
            dist_from_boundary,
            interp_point_spacing
          );

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}



// ---------------------------------------------------------------------- 
int SolverGMRESDoubleLayer::neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorQBKIX* near_evaluator = new EvaluatorQBKIX(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          interpolation_directions,
          _patch_samples->sample_point_far_field(),
          _patch_samples->sample_point_interpolant_spacing()
          
          );

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}


int SolverGMRESDoubleLayer::neaeval_extrapolate_average(
        Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec dist_from_boundary,  //closest on surface sample for interpolation
        Vec interp_point_spacing,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorQBKIXAverage* near_evaluator = 
      new EvaluatorQBKIXAverage(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          interpolation_directions,
            dist_from_boundary,
            interp_point_spacing
          );

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}


int SolverGMRESDoubleLayer::neaeval_extrapolate_average(
        Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorQBKIXAverage* near_evaluator = 
      new EvaluatorQBKIXAverage(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          interpolation_directions,
          _patch_samples->sample_point_far_field(),
          _patch_samples->sample_point_interpolant_spacing()
          );

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}



int SolverGMRESDoubleLayer::neaeval_interpolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential) 
{
  ebiFunctionBegin;
  Kernel3d  double_layer_kernel;
  double_layer_kernel = Kernel3d(_equation_type + DOUBLE_LAYER + quantity, this->eqcoefs());
  
  EvaluatorNearInterpolate* near_evaluator = new EvaluatorNearInterpolate(
          double_layer_kernel,
          this->patch_samples(),
          this->refined_patch_samples(),
          near_targets,
          closest_sample_3d_position,
          closest_sample_as_face_point,
          target_in_out_near,
          interpolation_directions);

  SolutionDensity solution(*_reference_solution);
  solution.set_petsc_vec(density);

  core_evaluation(solution, potential, (Kernel_variable) quantity, near_targets, near_evaluator);

  delete near_evaluator;
  ebiFunctionReturn(0);
}


Vec greens_identity(
        MPI_Comm comm,
        Kernel3d problem_kernel, 
        Vec singularity_location,
        Vec singularity_strength,
        Vec target_points,
        SolverGMRESDoubleLayer* solver){
    // we need to evaluate 
    // u(x) = \int_\Gamma (du/dn)(y)*G(x,y) - u(y)*(dG/dn)(x,y) dy
    //
    // We use 
    //     u(y) = \int_\Gamma G(y,z)*f(z) dz
    // and 
    //     (du/dn)(y) = \int_\Gamma (dG/dn)(y,z)(y)*f(z) dz
    // We split this into four separate integrals: two to evaluate u(y) and
    // du/dn(y) for all sample points y \in \Gamma, and two more to evaluate 
    // \int_\Gamma (du/dn(y)*G(x,y) dy and \int_\Gamma u(y)*(dG/dn)(x,y) dy.
    // There is also another FMM call to evaluate the boundary data f at each
    // sample point.

    // single-layer singularities => only position and point strength

    // Create target potential vector
    Vec sample_point_3d_position = solver->patch_samples()->sample_point_3d_position();
    Vec sample_point_normal = solver->patch_samples()->sample_point_normal();
    
    // Kress says green's identity holds for interior faces normals
    //VecScale(sample_point_normal, -1);
    
    //Vec boundary_data;
    int num_local_samples = num_local_points(sample_point_3d_position);
    int target_dof = problem_kernel.get_tdof();
    int num_targets = Petsc::get_vec_local_size(target_points)/DIM;


    int num_singularities = Petsc::get_vec_local_size(singularity_location)/DIM;
    
    NumVec<OnSurfacePoint> sample_point_as_on_surface_point = 
            solver->patch_samples()->sample_point_as_on_surface_point();
    Vec potential = Test::compute_dirichlet_boundary_data(
            solver->mpiComm(),
            problem_kernel,
            singularity_location,
            //singularity_normals,
            singularity_strength,
            sample_point_3d_position,
            sample_point_normal);

    cout << "evaluated u(x)" << endl;
    VecView(potential, PETSC_VIEWER_STDOUT_SELF);
    Vec potential_dn = Test::compute_neumann_boundary_data(
            solver->mpiComm(),
            problem_kernel,
            singularity_location,
            //singularity_normals,
            singularity_strength,
            sample_point_3d_position,
            sample_point_normal);
    cout << "evaluated du/dn(x)" << endl;
    VecView(potential_dn, PETSC_VIEWER_STDOUT_SELF);
    
    Vec solution = greens_identity(comm, potential, potential_dn, target_points, solver);
    Petsc::destroy_vec(potential);
    Petsc::destroy_vec(potential_dn);

    return solution;
}

Vec greens_identity(
        MPI_Comm comm,
        Vec dirichlet_data,
        Vec neumann_data,
        Vec target_points,
        SolverGMRESDoubleLayer* solver){
    int num_local_targets = num_local_points(target_points);
    DblNumMat target_points_local = 
        get_local_vector(DIM, num_local_targets, target_points);

    NumVec<OnSurfacePoint> on_surface_points_targets = 
        Markgrid::mark_target_points(target_points_local, 
                dynamic_cast<PatchSurfFaceMap*>(solver->bdry()),false);
    target_points_local.restore_local_vector();
    cout << "marked target points"  << endl;
    /*for(int i = 0; i < on_surface_points_targets.m(); i++){
        cout << on_surface_points_targets(i).region << endl;
    }*/
    return greens_identity(comm, dirichlet_data, neumann_data, 
            on_surface_points_targets, target_points, solver);
}


Vec greens_identity(
        MPI_Comm comm,
        Vec dirichlet_data,
        Vec neumann_data,
        NumVec<OnSurfacePoint> on_surface_points_targets,
        Vec target_points,
        SolverGMRESDoubleLayer* solver){
    assert(solver->mpiComm() == MPI_COMM_WORLD);

    Vec first_boundary_integral;        // \int_\Gamma (du/dn_x)(y)*G(x,y)dy, x \in Gamma
    Vec second_boundary_integral;       //\int_\Gamma u(y)*(dG/dn_y)(x,y) dy, for x \in Gamma
    Vec final_potential;
    
    Kernel3d k = solver->problem_kernel();
    int target_dof = k.get_tdof();
    int num_targets = Petsc::get_vec_local_size(target_points)/DIM;
    Petsc::create_mpi_vec(comm, num_targets*target_dof, first_boundary_integral);
    Petsc::create_mpi_vec(comm, num_targets*target_dof, second_boundary_integral);
    Petsc::create_mpi_vec(comm, num_targets*target_dof, final_potential);
    // evaluate  \int_\Gamma G(x,y)*(du/dn)(y) dy
    cout << "evaluating sl" << endl;
    Options::set_value_petsc_opts("-dirichlet", "0");
    Options::set_value_petsc_opts("-neumann", "1");
    solver->evaluate(target_points,
            neumann_data, 
            first_boundary_integral,
            on_surface_points_targets,
            (Kernel_variable) -9, // single layer
            true); // don't re-mark target points
    cout << " evaluated  \\int_\\Gamma G(x,y)*(du/dn)(y) dy" <<  endl;

    // evaluate  \int_\Gamma (dG/dn)(x,y)(y)*u(y) dy
    cout << "evaluating dl" << endl;
    Options::set_value_petsc_opts("-dirichlet", "1");
    Options::set_value_petsc_opts("-neumann", "0");
    solver->evaluate(target_points,
            dirichlet_data, 
            second_boundary_integral,
            on_surface_points_targets,
            VAR_U, // double layer
            true); // don't re-mark target points
    cout << " evaluated  \\int_\\Gamma (dG/dn)(x,y)(y)*u(y) dy" <<  endl;


    // bug found. note that Kress's convention is the opposite sign of the
    // convention in this code (Kress says integral of constant density is -1,
    // we resolve +1), so we in fact need to negate his green's formula to
    // remain consistent
    VecCopy(second_boundary_integral, final_potential);

    VecAXPY(final_potential, -1., first_boundary_integral);
    cout << "evaluated  u(x) = \\int_\\Gamma (du/dn)(y)*G(x,y) - u(y)*(dG/dn)(x,y) dy" << endl;

    VecDestroy(&first_boundary_integral);
    VecDestroy(&second_boundary_integral);
    return final_potential;
}

void SolverGMRESDoubleLayer::populate_qbx_data(const NumVec<OnSurfacePoint> on_surface_points,
        Vec& target_in_out, 
        Vec& closest_on_surface_points_3d_position,
    Vec& closest_on_surface_points_as_face_point,
    Vec& interpolation_directions,
    Vec& target_far_field,
    Vec& target_interpolant_spacing
    ){

    int num_local_targets = on_surface_points.m();

    // make the necessary petsc Vec's to pass to the legacy evaluator interface
    Petsc::create_mpi_vec(this->mpiComm(), num_local_targets, target_in_out);
    Petsc::create_mpi_vec(this->mpiComm(), num_local_targets*DIM, closest_on_surface_points_3d_position);

    Petsc::create_mpi_vec(this->mpiComm(), 
            num_local_targets*_patch_samples->bdry()->face_point_size_in_doubles(), 
            closest_on_surface_points_as_face_point);

    // set qbkix related junk TODO remove these
    VecDuplicate( target_in_out, &target_far_field);
    VecDuplicate( target_in_out, &target_interpolant_spacing);
    VecSet(target_far_field, Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    VecSet(target_interpolant_spacing, Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    Petsc::create_mpi_vec(this->mpiComm(),num_local_targets*DIM,interpolation_directions);
    VecSet(interpolation_directions, 0.);


    DblNumMat target_in_out_local(1, target_in_out);
    DblNumMat target_far_field_local(1, target_far_field);
    DblNumMat target_interpolant_spacing_local(1, target_interpolant_spacing);
    DblNumMat interpolation_directions_local(DIM, interpolation_directions);

    DblNumMat closest_on_surface_points_3d_position_local(3,closest_on_surface_points_3d_position);
    DblNumMat closest_on_surface_points_as_face_point_local(
            _patch_samples->bdry()->face_point_size_in_doubles(),closest_on_surface_points_as_face_point);

    // compute the relevant data needed for each piece of evaluation
    string qbkix_convergence_type = Options::get_string_from_petsc_opts("-qbkix_convergence_type");
    bool qbkix_classical_conv = qbkix_convergence_type == "classic";
    bool qbkix_adaptive_conv = qbkix_convergence_type == "adaptive";
    bool is_face_map = Options::get_int_from_petsc_opts("-bdtype") == 2;
    for(int i =0; i < num_local_targets; i++){
        OnSurfacePoint p = on_surface_points(i);

        // this is free, I don't think we need it for far eval, but why not...
        if(p.inside_domain == OUTSIDE | p.inside_domain == OUTSIDE_SPATIAL_GRID){
            target_in_out_local(0,i) = 2; // dangerous... see common/NbrDefs.hpp to fix
        } else if (p.inside_domain == INSIDE){
            target_in_out_local(0,i) = 1;
        } else if (p.inside_domain == ON_SURFACE){
            target_in_out_local(0,i) = 3; // ????
        }


        //  keep track of indices associated with each near/far point. if the
        //  point is near, compute the needed data to pass to the legacy
        //  interface.
        //assert(p.region == NEAR);

        FacePointFaceMap* face_point = 
            (FacePointFaceMap*) closest_on_surface_points_as_face_point_local.clmdata(i);
        *face_point = FacePointFaceMap::to_face_point(p);

        Point3 closest_on_surface_point(closest_on_surface_points_3d_position_local.clmdata(i));
        Point3 positions_and_derivs[3];
        FaceMapSubPatch* patch = 
            (FaceMapSubPatch*)_patch_samples->bdry()->patches()[p.parent_patch];

        patch->xy_to_patch_coords(p.parametric_coordinates.array(), 
                PatchSamples::EVAL_VL|PatchSamples::EVAL_FD, 
                (double*)positions_and_derivs);
        //closest_on_surface_point.array());
        Point3 interior_normal = -cross(positions_and_derivs[1], positions_and_derivs[2]).dir(); // interior normal!! note minus sign
        for(int d =0; d< DIM; d++){
            closest_on_surface_points_3d_position_local(d,i) = positions_and_derivs[0](d);
            interpolation_directions_local(d,i) = interior_normal(d);
        }
        if(is_face_map){
            if(qbkix_classical_conv){
                target_far_field_local(0,i) *= sqrt(patch->characteristic_length());
                target_interpolant_spacing_local(0,i) *= sqrt(patch->characteristic_length());
            } else if(qbkix_adaptive_conv){
                target_far_field_local(0,i) *= patch->characteristic_length();
                target_interpolant_spacing_local(0,i) *= patch->characteristic_length();
            }
        }


    }
}

END_EBI_NAMESPACE

