#include "evaluator_on_surface.hpp"
#include "solver_utils.hpp"
#include "common/vecmatop.hpp"
#include <time.h>
#include "common/stats.hpp"
#include "common/vtk_writer.hpp"

BEGIN_EBI_NAMESPACE

using std::max;
using std::min;

EvaluatorOnSurface::EvaluatorOnSurface(Kernel3d kernel, Vec target_3d_position, 
        Vec target_as_face_point, PatchSamples* patch_samples): 
    EvaluatorOnSurface("",""){

    _target_3d_position = target_3d_position;
    _target_as_face_point = target_as_face_point;
    _knl = kernel;
    _patch_samples = patch_samples;

}


// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "EvaluatorOnSurface::setFromOptions"

int EvaluatorOnSurface::setFromOptions()
{
  _refinement_factor = Options::get_int_from_petsc_opts("-dnref");
  _surface_interpolation_num_samples = Options::get_int_from_petsc_opts("-LL");
  return 0;
}

// ----------------------------------------------------------------------
// should we use singularity cancellation (p. 259, sec 3.3.2) to evaluate the
// kernel (i.e. is it hypersingular). Only occurs in the case of Stokes
// pressure evaluation. 
//
#undef __FUNCT__
#define __FUNCT__ "EvaluatorOnSurface::use_singularity_cancellation"

bool EvaluatorOnSurface::use_singularity_cancellation(Kernel3d& knl){

    // Only use singularity cancellation for Stokes pressure.
    return knl.kernelType() == KNL_STK_D_P;
}

int EvaluatorOnSurface::setup_no_fmm(){

  ebiFunctionBegin;
  ebiAssert(_target_3d_position!=NULL && _target_as_face_point!=NULL);

  //do constant setup if necessary
  iC( _collocation_data->distribute_collocation_points(
              _target_3d_position, 
              _target_as_face_point,
              this->_patch_samples, 
              this->source_dof(),
              this->target_dof()) );

  // setup constant density solutions for hypersingular eval
  iC( setup_unit_vector_densities() );

  // copy distributed source positions, normals, and target positions to
  // non-distributed arrays to pass to FMM initialization;  it seems this only works
  // correctly for a single processor, as it goes over sposarr and other
  // arrays obtained from distributed vectors below from 0 to global vector DIMension
  
  ebiFunctionReturn(0);
}




// ---------------------------------------------------------------------- 
int EvaluatorOnSurface::setup()
{
  ebiFunctionBegin;
  ebiAssert(_target_3d_position!=NULL && _target_as_face_point!=NULL);
  PatchSamples* patch_samples = this->_patch_samples;

  //do constant setup if necessary
  iC( _collocation_data->distribute_collocation_points(
              _target_3d_position, 
              _target_as_face_point,
              patch_samples, 
              this->source_dof(),
              this->target_dof()) );

  // setup constant density solutions for hypersingular eval
  iC( setup_unit_vector_densities() );

  // copy distributed source positions, normals, and target positions to
  // non-distributed arrays to pass to FMM initialization;  it seems this only works
  // correctly for a single processor, as it goes over sposarr and other
  // arrays obtained from distributed vectors below from 0 to global vector DIMension
  
  fmm = new PvFMM(
          patch_samples->sample_point_3d_position(), /* source positions */
          patch_samples->sample_point_normal(),      /* source normals */
          _target_3d_position,                       /* target positions */
          this->knl()                                /* kernel to evaluate */
          );
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 

// precompute convolutions of the kernel with  constant densities
// corresponding to unit vectors (for Laplace - just 1, for Stokes [1,0,0],
// [0,1,0], [0,0,1]
#undef __FUNCT__
#define __FUNCT__ "EvaluatorOnSurface::setup_unit_vector_densities"

int EvaluatorOnSurface::setup_unit_vector_densities()
{
  ebiFunctionBegin;
  //1. untvecs
  
  // Check if evaluating Stokes pressure

  if(use_singularity_cancellation(this->knl())) {
	 PatchSamples* patch_samples = this->_patch_samples;
	 
     int source_dof = this->source_dof();  
     int target_dof = this->target_dof();

	 Vec source_3d_positions = this->_patch_samples->sample_point_3d_position();

     int num_local_targets = num_local_points(_target_3d_position);
	 int num_local_sources = num_local_points(source_3d_positions);

	 //space
	 // untvec = i-th Vec is the convolution of the kernel with i-th
	 // unit vector, i = 1.. source_dof
	 vector<Vec>& untvec = _singularity_cancellation_data.untvec();
	 untvec.resize(source_dof, NULL);
     
     // Initialize all unit vector densities to zero
     const double zero = 0.0;
	 for(int i=0; i<source_dof; i++) {
		
         VecCreateMPI(this->mpiComm(), 
                      num_local_targets*target_dof, 
                      PETSC_DETERMINE, 
                      &(untvec[i]));
        
         VecSet( untvec[i], zero);
	 }

	 //calculate untvec
	 Vec den;
     VecCreateMPI(this->mpiComm(), 
                 num_local_sources*source_dof, 
                 PETSC_DETERMINE, 
                 &den
                 );
     
     double* denarr;

	 for(int i=0; i<source_dof; i++) {
        iC( VecSet( den, zero) );
        iC( VecGetArray(den, &denarr) );
        
        // Set set up the unit vectors by setting the ith components to 1.0
		for(int k=0; k<num_local_sources; k++) {
            denarr[ k*source_dof+i ] = 1.0;
        }

		iC( VecRestoreArray(den, &denarr) );

		// apply kernel to const density with i-th component = 1,
		// rest = 0
		iC( singular_evaluation(den, untvec[i]) );
	 }
	 iC( VecDestroy(&den) );
  }
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "EvaluatorOnSurface::eval"

int EvaluatorOnSurface::eval(Vec density, Vec val)
{
  ebiFunctionBegin;  //ebiLogInfo( "eval.............");


  double matvec_start = omp_get_wtime(); 
  //1. fmm
  // scale densities by quadrature weights and jacobians
  Vec scaled_density;
  VecDuplicate(density, &scaled_density);
  denscale(this->source_dof(),                /* kernel source degrees of freedom */
           this->_patch_samples->sample_point_combined_weight(),  /* per point scaling factor */
           density,                           /* density to scale */
           scaled_density);                   /* result: scaled density */

  //  call FMM evaluation
  double start = omp_get_wtime(); 
  fmm->evaluate(scaled_density, val); 
    stats.result_plus_equals("total fmm time", omp_get_wtime() - start );
  VecDestroy(&scaled_density);

  //  usual singular kernel evaluation
  start = omp_get_wtime(); 
  singular_evaluation(density, val);
  stats.result_plus_equals("total singular eval time", omp_get_wtime() - start );

  // if hypersingular, do the singularity cancellation modification
  // this function does nothing unless use_singularity_cancellation
  // returns true
  apply_singularity_cancellation(density, val);
  stats.result_plus_equals("total matvec time", omp_get_wtime() - matvec_start );

  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// Implementation of Algorithm 1: Singular integral evaluation. The first sum is computed
// by FMM (First ~50 lines of singular_evaluation) to evaluate the kernel at 
// the quadrature points and the "prep" portion of code (to actually compute
// the integral. The second sum is computed by the "sub" and "add" pieces of
// code, which subtract the inaccurate porition compute by the FMM and add back
// the correctly computed singular integral component via polar integration,
// respectively.
int EvaluatorOnSurface::singular_evaluation(Vec den, Vec val)
{
  ebiFunctionBegin;
  
  // Overlapping patch representation
  PatchSamples* patch_samples = this->_patch_samples;
  PatchSurf* bd = patch_samples->bdry();

  // Regular sampling of overlapping patch representation
  
  //-------------------------
  
  //-------------------------
  //2. bruno-kunyansky: 
  int source_dof = this->source_dof();
  int target_dof = this->target_dof();
  double zero  = 0.0;

  //  allocate mod vecs
  Vec _colloc_point_3d_position = _collocation_data->colloc_point_3d_position();
  int num_local_colloc_points = num_local_points(_colloc_point_3d_position);
  int num_local_targets = num_local_points(_target_3d_position);

  // Vec storing modified densities of collocation points
  Vec density_at_collocs;
  VecCreateMPI(this->mpiComm(), 
              num_local_colloc_points*target_dof, 
              PETSC_DETERMINE, 
              &density_at_collocs
              );
  VecSet( density_at_collocs, zero);

  // Vec storing modified densities of target points
  Vec density_at_targets;
  VecCreateMPI(this->mpiComm(), 
               num_local_targets*target_dof, 
              PETSC_DETERMINE,
              &density_at_targets
              );  
  VecSet( density_at_targets, zero);

  // multiply density by blend function -- preparing for summation
  Vec dat;
  iC( VecDuplicate(den,&dat) );
   denscale(source_dof,                     /* kernel source degree of freedom */
              patch_samples->sample_point_blend_func_value(),  /* per point scaling factor */
              den,                          /* density to scale */
              dat); //dat = den * alf;      /* output: scaled density */

  // only needed for hypersingular?? not sure 

  double start = omp_get_wtime(); 
  iC( patch_samples->refine_data(source_dof, _refinement_factor, dat, _refined_datvec));
  stats.result_plus_equals("total fft upsample time", omp_get_wtime() - start );
  _dat = dat;
  cout << "prep " << endl;
  
  start = omp_get_wtime(); 
  //  sub + add

  // radius of the floating partition of unity function
  double floating_POU_radius_scale;
  radmult_spacing(patch_samples->spacing(), floating_POU_radius_scale);
  int face_point_sze_dbl = patch_samples->face_point_size_in_doubles();
  
  CollocationPatchSamples* colloc_data = _collocation_data;

  // summation over patches belonging to current partition
  // the loop is over all patches, but the for(int ci..) loop below
  // will not execute unless the point is owned by the current partition
  for(int pi=0; pi<bd->patches().size(); pi++) {
	 //-----------------
	 Patch* patch = bd->patches()[pi];

     // Collect 3d positions, face points, and to-be-computed density arrays for 
     // collocation points for the pi-th patch
	 DblNumMat colloc_point_as_face_point = 
             colloc_data->colloc_point_as_face_point(pi, face_point_sze_dbl);

	 DblNumMat colloc_point_3d_position = 
             colloc_data->colloc_point_3d_position(pi);

	 DblNumMat colloc_point_density = 
             colloc_data->coldat(pi, density_at_collocs, target_dof); 

     double* colloc_point_density_ptr = colloc_point_density.data();
	 //----------------
     
     // Number of collocation points on the pi-th patch
     int num_colloc_points = colloc_data->num_collocation_points(pi);
	 // compute the second integral in (19), p257
#pragma omp parallel for
	 for(int ci = 0; ci < num_colloc_points; ci++) {
		FacePointOverlapping* cur_face_point =
            (FacePointOverlapping*)(colloc_point_as_face_point.clmdata(ci));

		double xy[2];
        iC( patch->face_point_to_xy(cur_face_point,xy) );

		//current target
		DblNumMat cith_colloc_point(DIM, 1, false, colloc_point_3d_position.clmdata(ci));
		DblNumVec curval(target_dof, false, colloc_point_density_ptr+ci*target_dof);

		{
		  //source
		  DblNumMat srcpos, srcnor;
          DblNumVec source_quadrature_weight;
          DblNumMat srcdat;
		  subtract_inaccurate_part(
                      pi,                       /* patch number */
                      source_dof,               /* kernel source degree of freedom */
                      xy,                       /* xy position of 3d-point in param. space */
                      floating_POU_radius_scale, /* radius of floating POU for singular evaluation */
                      srcpos, 
                      srcnor, 
                      source_quadrature_weight, 
                      srcdat);

          //build the dense matrix evaluating the kernel at
          //the target point
		  for(int k=0; k<srcdat.n(); k++)
			 for(int d=0; d<srcdat.m(); d++)
				srcdat(d,k) *= source_quadrature_weight(k);

		  DblNumMat submat;
          submat.resize(target_dof, source_dof*srcpos.n());
          fmm->interaction_matrix(srcpos, srcnor, 
                  cith_colloc_point, submat);
		  //(this->knl()).kernel(srcpos,              /* source position */
          //                     srcnor,              /* source normal */
          //                     cith_colloc_point,   /* target position */
          //                     submat);             /* potential */

		  DblNumVec tmpdat(srcdat.m()*srcdat.n(), false, srcdat.data());

		  // compute the integral  
		  iC( dgemv(-1.0, submat, tmpdat, 1.0, curval) );
		}
	 }
  }
  cout << "sub " << endl;
  cout << "innaccurate local part" << endl;

stats.result_plus_equals("total inaccurate subtract time", omp_get_wtime() - start );
  start = omp_get_wtime();

  // compute discretization of integral (18) on p 256
  // (polar coordinates trapezoid rule)
  // the loops here are identical to "sub" part", the difference is in
  // add_singular_quadrature_part, and the integral is added rather than subtracted 
 Kernel3d true_kernel(this->knl().kernelType(),this->knl().coefs(), 0.);
  for(int pi=0; pi<bd->patches().size(); pi++) {
		//-----------------
		Patch* patch = bd->patches()[pi];

         // Collect 3d positions, face points, and to-be-computed density arrays for 
         // collocation points for the pi-th patch
        int num_colloc_points = colloc_data->num_collocation_points(pi);

		DblNumMat colloc_point_as_face_point(
                colloc_data->colloc_point_as_face_point(pi, face_point_sze_dbl));
		DblNumMat colloc_point_3d_position(
                colloc_data->colloc_point_3d_position(pi));
		DblNumMat colloc_point_density(
                colloc_data->coldat(pi, density_at_collocs, target_dof));
		//----------------

#pragma omp parallel for
		for(int ci=0; ci<num_colloc_points; ci++) {
		  FacePointOverlapping* cur_face_point =
              (FacePointOverlapping*)(colloc_point_as_face_point.clmdata(ci)); 

		  double xy[2];
          iC( patch->face_point_to_xy(cur_face_point,xy) );

		  //current target
		  DblNumMat cith_colloc_point(DIM,1,false,
                  colloc_point_3d_position.clmdata(ci));
		  DblNumVec curval(target_dof, false, colloc_point_density.data()+ci*target_dof);
		  {
			 //source
			 DblNumMat srcpos, srcnor;
             DblNumVec source_quadrature_weight;
             DblNumMat srcdat;
			 add_singular_quadrature_part(
                         pi,                            /* patch number */
                         source_dof,                    /* kernel source degree of freedom */
                         xy,                            /* xy position of 3d-point in param. space */
                         floating_POU_radius_scale,     /* radius of floating POU for singular evaluation */
                         srcpos, 
                         srcnor,
                         source_quadrature_weight,
                         srcdat);

			 for(int k=0; k<srcdat.n(); k++)
				for(int d=0; d<srcdat.m(); d++)
				  srcdat(d,k) *= source_quadrature_weight(k);
			 //matrix
			 DblNumMat addmat;
             addmat.resize(target_dof, source_dof*srcpos.n());

          //fmm->interaction_matrix(srcpos, srcnor, 
                  //cith_colloc_point, addmat);
			 //(this->knl()).kernel(srcpos,               /* source position */
             //                    srcnor,               /* source normal */
             //                    cith_colloc_point,    /* target position */
             //                     addmat);              /* potential */
			 true_kernel.kernel(srcpos,               /* source position */
                                 srcnor,               /* source normal */
                                 cith_colloc_point,    /* target position */
                                  addmat);              /* potential */

			 DblNumVec tmpdat(srcdat.m()*srcdat.n(), false, srcdat.data());
            
             // compute the integral
			 iC( dgemv( 1.0, addmat, tmpdat, 1.0, curval) );
		  }
		}
  }
  stats.result_plus_equals("total polar quad time", omp_get_wtime() - start );
 
  iC( VecDestroy(&dat) );
  _dat = NULL;

  iC( VecScatterBegin(colloc_data->colloc_to_target_value_scatter(),
                      density_at_collocs,
                      density_at_targets,
                      ADD_VALUES,
                      SCATTER_FORWARD) );

  iC( VecScatterEnd(colloc_data->colloc_to_target_value_scatter(),
                    density_at_collocs,
                    density_at_targets,
                    ADD_VALUES,
                    SCATTER_FORWARD) );

  cout << "add" << endl;
  
  //3. put modification back and deallocate mod vecs
  PetscScalar one = 1.0;

  // just testing
  int64_t tvmsze,valsze;
  iC( VecGetSize(density_at_targets,&tvmsze));
  iC( VecGetSize(val,&valsze));
  ebiAssert( tvmsze == valsze );

  // add discretization of (18) + (19)second part to the first part of (19) 
  iC( VecAXPY( val, one,  density_at_targets) );
  iC( VecDestroy(&density_at_collocs) );
  iC( VecDestroy(&density_at_targets) );
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// 
#undef __FUNCT__
#define __FUNCT__ "EvaluatorOnSurface::apply_singularity_cancellation"

int EvaluatorOnSurface::apply_singularity_cancellation(Vec den, Vec val)
{
  ebiFunctionBegin;
  
  // hypersingular only 
  if( use_singularity_cancellation(this->knl()) ) {
	 PatchSamples* patch_samples = this->_patch_samples;
	 PatchSurf* bd = patch_samples->bdry();
     CollocationPatchSamples* colloc_data = this->_collocation_data;
     int face_point_sze_dbl = patch_samples->face_point_size_in_doubles();

	 int source_dof = this->source_dof();
     int target_dof = this->target_dof();

	 Vec _colloc_point_3d_position = colloc_data->colloc_point_3d_position();	 
	 
     int num_local_colloc_points = 
         num_local_points(_colloc_point_3d_position);
	 int num_local_targets = num_local_points(_target_3d_position);
	 double zero = 0.0;
	 
     Vec density_at_collocs;
     iC( VecCreateMPI(this->mpiComm(), num_local_colloc_points*source_dof, 
                 PETSC_DETERMINE, &density_at_collocs) );
     iC( VecSet( density_at_collocs, zero) );

	 Vec density_at_targets;
     iC( VecCreateMPI(this->mpiComm(), num_local_targets*source_dof, 
                 PETSC_DETERMINE, &density_at_targets) );
     iC( VecSet( density_at_targets, zero) );

	 Vec dat;
     VecDuplicate(den,&dat); 
     denscale(source_dof, patch_samples->sample_point_blend_func_value(), den, dat);

	 // mystery refinement 
	 patch_samples->refine_data(source_dof, _refinement_factor, dat, _refined_datvec);
     _dat = dat;

	 for(int pi=0; pi<bd->patches().size(); pi++) {
		Patch* patch = bd->patches()[pi];
        int num_colloc_points = colloc_data->num_collocation_points(pi);

		DblNumMat colloc_point_as_face_point = 
                colloc_data->colloc_point_as_face_point(pi, face_point_sze_dbl);
		DblNumMat colden = colloc_data->coldat(pi, density_at_collocs, source_dof);

#pragma omp parallel for
		for(int ci = 0; ci < num_colloc_points; ci++) {

		  FacePointOverlapping* cur_face_point =
              (FacePointOverlapping*)(colloc_point_as_face_point.clmdata(ci));
		  
		  bool is_valid;
		  iC( patch->is_face_point_valid(cur_face_point, is_valid) ); 
		  ebiAssert(is_valid==true); // cur_face_point is inside the (pi)th patch
		  
		  double xy[2]; 
          iC( patch->face_point_to_xy(cur_face_point,xy) );

		  iC( patch_samples->interpolate_data(pi, xy, _refinement_factor, 
                      _surface_interpolation_num_samples,
                      PatchSamples::EVAL_VL, _refined_datvec, colden.clmdata(ci)) );
		}
	 }
	 VecScatterBegin(colloc_data->colloc_to_target_density_scatter(), 
                     density_at_collocs,  
                     density_at_targets,  
                     ADD_VALUES, 
                     SCATTER_FORWARD);
	VecScatterEnd(colloc_data->colloc_to_target_density_scatter(),   
                  density_at_collocs,
                  density_at_targets,  
                  ADD_VALUES,  
                  SCATTER_FORWARD);

	 iC( VecDestroy(&dat) );
	 _dat = NULL;
	 
	 //2. remove const part
	 // precomputed integrals for constant density
	 vector<Vec>& untvec = _singularity_cancellation_data.untvec();
	 double* denarr; 
	 double* valarr; 

     iC( VecGetArray(density_at_targets, &denarr) );
     iC( VecGetArray(val,    &valarr) );

	 vector<double*> untarrvec; 
     untarrvec.resize( source_dof, NULL );

	 for(int i=0; i<source_dof; i++) {
		iC( VecGetArray(untvec[i], &(untarrvec[i])) );
	 }

#pragma omp parallel for
	 // subtract the sum of precomputed integrals
	 for(int k=0; k<num_local_targets; k++) {
		for(int i=0; i<source_dof; i++) {
		  int soff = k*source_dof + i;
		  for(int j=0; j<target_dof; j++) {
			 int toff = k*target_dof + j;
			 valarr[toff] -= denarr[soff] * untarrvec[i][toff];
		  }
		}
	 }

	 for(int i=0; i<source_dof; i++) {
		iC( VecRestoreArray(untvec[i], &(untarrvec[i])) );
	 }

	 iC( VecRestoreArray(density_at_targets, &denarr) );
	 iC( VecRestoreArray(val,    &valarr) );
	 
	 iC( VecDestroy(&density_at_collocs) );
	 iC( VecDestroy(&density_at_targets) );
  }
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
// go over a regular grid of samples in patch pi inside the bounding box of
// the the support of the floating partition of unity function centered at
// xy
// retrieve face points with positions on this grid and corresponding 3d
// positions. normals jacobians and density return these in subpos etc
// TODO change the return type to void
int EvaluatorOnSurface::subtract_inaccurate_part(int pi, int source_dof, double* xy,
        int floating_POU_radius_scale,
        DblNumMat& subpos, DblNumMat& subnor, DblNumVec& subjaw, DblNumMat& subdat)
{
  ebiFunctionBegin;
  PatchSamples* patch_samples = this->_patch_samples;
  PatchSurf* bd = patch_samples->bdry();
  Patch* patch = bd->patches()[pi];  
  int num_samples = patch_samples->num_sample_points()[pi];
  double step = patch_samples->step_size()[pi];
  double init = -patch->bnd();

  // radius of the floating partition of unity 
  double RAD = step * floating_POU_radius_scale;  

  // vector of indices of sample points inside the square of halfwidth RAD
  // centered at xy 
  vector<int> ijvec;
  double xyll[2];
  xyll[0]=xy[0]-RAD;
  xyll[1]=xy[1]-RAD;

  int idll[2];
  idll[0]=(int)ceil( (xyll[0]-init)/step );
  idll[1]=(int)ceil( (xyll[1]-init)/step );

  idll[0] = max(idll[0], 0);
  idll[1] = max(idll[1], 0);

  double xyur[2];
  xyur[0]=xy[0]+RAD;
  xyur[1]=xy[1]+RAD;

  int idur[2];
  idur[0]=(int)ceil( (xyur[0]-init)/step );
  idur[1]=(int)ceil( (xyur[1]-init)/step );

  idur[0] = min(idur[0], num_samples);
  idur[1] = min(idur[1], num_samples);
  for(int j=idll[1]; j<idur[1]; j++) {
	 for(int i=idll[0]; i<idur[0]; i++) {
		int ij[2]; ij[0]=i; ij[1]=j;

		bool is_valid; 
		iC( patch_samples->is_sample_point_valid(pi, ij, is_valid) );

		if(is_valid==true) {	
		  ijvec.push_back(i);
          ijvec.push_back(j);
		}
	 }
  }
  //2. collect stuff
  int subcnt = ijvec.size() / 2;

  subpos.resize(DIM,subcnt);
  subnor.resize(DIM,subcnt);
  subjaw.resize(subcnt);
  subdat.resize(source_dof, subcnt);

  for(int k=0; k<subcnt; k++) {
      int ij[2];
      ij[0]=ijvec[2*k];
      ij[1]=ijvec[2*k+1];

      double polar_quad_point[2];
      polar_quad_point[0] = init + ij[0]*step;
      polar_quad_point[1] = init + ij[1]*step;

      //new xy
      double off[2];
      off[0] = polar_quad_point[0] - xy[0];
      off[1] = polar_quad_point[1] - xy[1];

	 // distance to the center of the sampled area for evaluting eta
	 double dist = sqrt(off[0]*off[0] + off[1]*off[1]);
	 //geo
	 iC( patch_samples->get_sample_on_patch(pi, ij, subpos.clmdata(k), subnor.clmdata(k), 
                 subjaw.data()+k) );
	 // trapezoidal quadrature weight
	 subjaw(k) *= step*step;

	 //dat
	 iC( patch_samples->get_sample_point(pi, ij, source_dof, _dat, subdat.clmdata(k)) );
	 double et = eta(dist/RAD);

	 // multiply density by eta 
	 for(int d=0; d<source_dof; d++) 
         subdat(d,k) *= et;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
// TODO change return type to void
int EvaluatorOnSurface::add_singular_quadrature_part(int pi, int source_dof, double* xy,
        int floating_POU_radius_scale, DblNumMat& polar_quad_point_positions, 
        DblNumMat& polar_quad_point_normals, DblNumVec& polar_quad_point_quad_weights, 
        DblNumMat& corrected_density) { 
    
    PatchSamples* patch_samples = this->_patch_samples;

    double step = patch_samples->step_size()[pi]; 
    double RAD = step * floating_POU_radius_scale;  
    int64_t dnref = Options::get_int_from_petsc_opts("-dnref");
    int64_t LL= Options::get_int_from_petsc_opts("-LL");

    //vector of x,y positions of uniform samples in radial coordinates
    vector<double> xyvec;

    // resolution in radial direcion ceil( RAD/step -0.5) = round(RAD/step)
    int num_radial_quad_points = (int)ceil((RAD-step/2.0)/step);
    double radial_quad_step = step;

    // twice as many points in the angular direction
    int num_angular_quad_points = 2*num_radial_quad_points;
    double angular_quad_step = PI/num_angular_quad_points;

    // compute patch domain xy positions of polar samples
    for(int ai=0; ai<num_angular_quad_points; ai++) {

        double ang = ai*angular_quad_step;
        for(int ri=-num_radial_quad_points; ri<num_radial_quad_points; ri++) {

            double rad = ri*radial_quad_step + radial_quad_step/2.0;

            // convert to polar coordinates
            double polar_quad_point[2];
            polar_quad_point[0] = rad*sin(ang) + xy[0];
            polar_quad_point[1] = rad*cos(ang) + xy[1];

            // Check that this quadrature point is in the (pi)th patch
            bool is_valid;
            patch_samples->is_interp_xy_valid(pi, polar_quad_point, is_valid);

            if(is_valid) {
                xyvec.push_back(polar_quad_point[0]);
                xyvec.push_back(polar_quad_point[1]);

            }
        }
    }
    //2. collect stuff
    int num_polar_quad_points = xyvec.size() / 2;

    polar_quad_point_positions.resize(DIM,num_polar_quad_points);
    polar_quad_point_normals.resize(DIM,num_polar_quad_points);
    polar_quad_point_quad_weights.resize(num_polar_quad_points);

    corrected_density.resize(source_dof,num_polar_quad_points);
    for(int k=0; k<num_polar_quad_points; k++) {
        double polar_quad_point[2];
        polar_quad_point[0] = xyvec[2*k];
        polar_quad_point[1] = xyvec[2*k+1];

        double off[2];
        off[0] = polar_quad_point[0] - xy[0];
        off[1] = polar_quad_point[1] - xy[1];

        double dist = sqrt(off[0]*off[0] + off[1]*off[1]);
        //geo
        Point3 point_position_and_derivs[3];
        // compute surface position and derivatives (tangents)
        patch_samples->interpolated_position_and_derivatives(
                pi,
                polar_quad_point, 
                PatchSamples::EVAL_VL | PatchSamples::EVAL_FD,
                (double*)point_position_and_derivs) ;
        
        Point3 cp = point_position_and_derivs[0]; // positions
        Point3 cn = cross(point_position_and_derivs[1],
                point_position_and_derivs[2]); // n = tangent1 x tangent2 	 
        double jac = cn.l2();  // Jacobian 
        cn = cn/jac; // unit normal

        for(int d=0; d<DIM; d++)
            polar_quad_point_positions(d,k) = cp(d);
        for(int d=0; d<DIM; d++)
            polar_quad_point_normals(d,k) = cn(d);

        polar_quad_point_quad_weights(k) = jac * radial_quad_step*angular_quad_step*dist;

        // interpolate data to position polar_quad_point	 
        patch_samples->interpolate_data(
                pi,
                polar_quad_point, 
                _refinement_factor,
                _surface_interpolation_num_samples,
                PatchSamples::EVAL_VL, 
                _refined_datvec, 
                corrected_density.clmdata(k));

        // floating POU function evaluation
        double et = eta(dist/RAD);
        for(int d=0; d<source_dof; d++) 
            corrected_density(d,k) *= et;
    }
    return 0;
}

// ---------------------------------------------------------------------- 
double EvaluatorOnSurface::eta(double t)
{  
    double res[3];

    double eps = 1e-7;
    double a;  
    pou1d(PatchSamples::EVAL_VL, t, eps, 1.0-eps, res, 
            this->_patch_samples->bdry()->pouctrl(),
            this->_patch_samples->bdry()->pdeg());
   a = res[0];
    /*double a;  
    pou1d(PatchSamples::EVAL_VL, t, 0., 1.0, res, 
            this->_patch_samples->bdry()->pouctrl(),
            this->_patch_samples->bdry()->pdeg());
   a = res[0];
   */
    //a = exp(1.-pow(2*t,8))/exp(1.);
    //a = exp(1.-pow(1.5*t,4))/exp(1.);
    //a = exp(-20.*pow(t,8));
    
    // good convergence rate for rad=2/sqrt(h)
    //a = exp(-40.*pow(t,4));
    //a = exp(-40.*pow(t,5));
    /*if(t <= 1e-7) {
        a = 1.0;
    } else if (t >= 1-1e-7){
        a = 0.0;
    } else {
    a = exp(-25.*pow(t,6));
    //a = exp(-36.*pow(t,10));
    }*/
    //a = exp(-25.*pow(t,6));
    /*double eps = 1e-5
        double LB = eps;
        double UB = 1.-eps;
    if(t <= LB) {
        a = 1.0;
    } else if (t >=UB){
        a = 0.0;
    } else {
        t = (t-LB)/(UB-LB);
        a = exp(1-1.0/(1-exp(1-1.0/t)));
        //a =a/(UB-LB);
    }*/

    return a;
}

END_EBI_NAMESPACE
