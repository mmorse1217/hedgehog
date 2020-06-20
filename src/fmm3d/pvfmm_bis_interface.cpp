#include "pvfmm_bis_interface.hpp" //interface in fmm3d
#include <memory>
#include <vec3t.hpp>
#include "common/utils.hpp"


// from pvfmm_interface.h
#include <pvfmm.hpp>
#include <cstddef>  
#include <vector>
#include "common/kernel3d.hpp"
#include <cassert>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <time.h>
#include "bie3d/solver_utils.hpp"
typedef double real_t;
#include "pvfmm/pvfmm_stks_kernel.hpp"
#include "pvfmm/pvfmm_laplace_kernel.hpp"
#include "pvfmm/pvfmm_mod_helmholtz_kernel.hpp"
#include "pvfmm/pvfmm_navier_kernel.hpp"
#define COUT(str) (std::cout<<str<<std::endl)
#define CERR(str,action) (std::cerr<<"[ERROR]["<< __FUNCTION__ <<"] "<<str<<std::endl, action)
#define ASSERT(expr, msg) ( (expr) ? assert(true) : CERR(msg,abort()))
#define COUTDEBUG(str) (std::cout<<"[DEBUG] "<<str<<std::endl)


double POISS;
double SHEAR;

using std::vector;

typedef std::vector<real_t> vec;
typedef pvfmm::FMM_Node<pvfmm::MPI_Node<real_t> > Node_t;

class Ebi::PvFMM::PVFMMImpl{
  public:
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  PVFMMImpl() : tree(NULL), mat(NULL) {

    max_depth = Options::get_int_from_petsc_opts("-bis3d_maxlevel");
    mult_order = Options::get_int_from_petsc_opts("-bis3d_np");
    max_pts = Options::get_int_from_petsc_opts("-bis3d_ptsmax");
    periodic = false;
    dense_eval = Options::get_int_from_petsc_opts("-direct_eval");
    
    // TODO factor this out when user specified communicators are added
    comm = MPI_COMM_WORLD;
  };

  ~PVFMMImpl(){
    COUTDEBUG("destructing context object");
    if (this->tree != NULL){
        delete this->tree;
    }
    if (this->mat != NULL){
        delete this->mat;
    }
  
    COUTDEBUG("deleted all pointers");
  };

  int mult_order;
  int max_pts;
  int max_depth;
  int sl;
  int dl;
  int periodic;
  int source_dof;
  int target_dof;
  int dense_eval;
  MPI_Comm comm;

  pvfmm::BoundaryType boundary;
  const pvfmm::Kernel<real_t>* ker;

  Tree_t* tree;
  Mat_t* mat;
};

const pvfmm::Kernel<double>* get_pvfmm_kernel(int _kernel_type,
                                              Ebi::KernelOptions _kernel_opts) {
  const pvfmm::Kernel<double>* kernel;
  switch (_kernel_type) {
  case Ebi::KNL_STK_S_U:
    kernel = &pvfmm::StokesKernel<double>::velocity();
    break;
  case Ebi::KNL_STK_S_P:
    kernel = &pvfmm::StokesKernel<double>::pressure();
    break;
  case Ebi::KNL_STK_D_U:
    kernel = &ker_stokes_dl;
    break;
  case Ebi::KNL_STK_D_P:
    // kernel = &pvfmm::StokesKernel<double>::pressure();
    kernel = &ker_stokes_pressure_dl;
    break;
  case Ebi::KNL_LAP_D_U:
    // kernel = &pvfmm::LaplaceKernel<double>::potential();
    kernel = &ker_laplace_dl;
    break;
  case Ebi::KNL_LAP_S_U:
    // kernel = &pvfmm::LaplaceKernel<double>::potential();
    kernel = &ker_laplace_sl;
    break;
  case Ebi::KNL_MODHEL_D_U:
    assert(_kernel_opts.initialized);

    helmholtz_frequency(_kernel_opts.helmholtz_frequency);
    kernel = &ker_mod_helmholtz_dl;
    const_cast<std::string &>(kernel->ker_name) = helmholtz_name(false);
    break;
  case Ebi::KNL_MODHEL_S_U:
    assert(_kernel_opts.initialized);

    helmholtz_frequency(_kernel_opts.helmholtz_frequency);
    kernel = &ker_mod_helmholtz_sl;
    const_cast<std::string &>(kernel->ker_name) = helmholtz_name(true);

    break;
  case Ebi::KNL_NAV_D_U:
    assert(_kernel_opts.initialized);

    navier_mu(_kernel_opts.navier_mu);
    navier_nu(_kernel_opts.navier_nu);
    kernel = &ker_navier_dl;
    POISS = navier_nu();
    SHEAR = navier_mu();
    // kernel = &ElasticKernel<double>::Disp();
    const_cast<std::string &>(kernel->ker_name) = navier_name(false);
    // kernel->ker_poten = navier_sl;
    break;
  case Ebi::KNL_NAV_S_U:
    assert(_kernel_opts.initialized);
    navier_mu(_kernel_opts.navier_mu);
    navier_nu(_kernel_opts.navier_nu);
    POISS = navier_nu();
    SHEAR = navier_mu();
    // kernel = &ElasticKernel<double>::Disp();
    kernel = &ker_navier_dl;
    const_cast<std::string &>(kernel->ker_name) = navier_name(false);
    break;
  default:
    ASSERT(false, "KernelNotImplementedError");
  }
  return kernel;
}

void Ebi::PvFMM::_unscale_potential(int ntrg, int target_dof, real_t *pot) {
  int omp_p = omp_get_max_threads();

  Kernel3d temp(_kernel_type, vector<double>());
  double density_unscale_factor =  temp.density_unscaling_value(_scale_factor);

#pragma omp parallel for
  for (int i = 0; i < omp_p; i++) {
    size_t a = (i * ntrg * target_dof) / omp_p;
    size_t b = ((i + 1) * ntrg * target_dof) / omp_p;
    for (size_t iT = a; iT < b; iT++) {
      pot[iT] = density_unscale_factor*pot[iT];
      if (std::isnan(pot[iT])) {
        CERR("fmm potential is NaN", abort());
      }
    }
  }
}

void Ebi::PvFMM::_copy_potential(int ntrg, int target_dof, std::vector<real_t> potv,
                    real_t *pot, bool rescale) {
  int omp_p = omp_get_max_threads();

  Kernel3d temp(_kernel_type, vector<double>());
  double density_unscale_factor = rescale ? temp.density_unscaling_value(_scale_factor) : 1.;

#pragma omp parallel for
  for (int i = 0; i < omp_p; i++) {
    size_t a = (i * ntrg * target_dof) / omp_p;
    size_t b = ((i + 1) * ntrg * target_dof) / omp_p;
    for (size_t iT = a; iT < b; iT++) {
      pot[iT] = density_unscale_factor*potv[iT];
      if (std::isnan(pot[iT])) {
        CERR("fmm potential is NaN", abort());
      }
    }
  }
}




/*****************************************************************************/
// end pvfmm_interface.cc
/*****************************************************************************/
BEGIN_EBI_NAMESPACE

using namespace std;
using namespace Ebi;


PvFMM::PvFMM() {
}
int PvFMM::setup(){

}
void PvFMM::_setup_context(){
    ebiAssert(_kernel_type != 0);
    
    ebiAssert(_equation_type == LAPLACE || _equation_type == STOKES 
            || _equation_type == MOD_HELMHOLTZ || _equation_type == NAVIER);
    _rebuild_tree = true; 


    int periodic = false; 
  pvfmm::BoundaryType bndry(periodic ? pvfmm::Periodic : pvfmm::FreeSpace);

  // Create new context.
  unique_ptr<PVFMMImpl> impl(new PVFMMImpl);

  impl->comm = MPI_COMM_WORLD;
  impl->source_dof = _source_dof;
  impl->target_dof = _target_dof;
  impl->sl = _kernel_layer == SINGLE_LAYER;
  impl->dl = _kernel_layer == DOUBLE_LAYER;

  impl->boundary = bndry;
  impl->ker = get_pvfmm_kernel(_kernel_type, _kernel_opts);

    COUTDEBUG("Initializing a pvfmm context object"
        <<", kernel=" << _kernel_type 
	    <<", sl="<<impl->sl
	    <<", dl="<<impl->dl
	    <<", sdof="<<impl->source_dof
	    <<", tdof="<<impl->target_dof
	    <<", mult_order="<<impl->mult_order
	    <<", max_pts="<<impl->max_pts
	    <<", max_depth="<<impl->max_depth
	    <<", dense_eval="<<impl->dense_eval
	    <<", periodic="<<impl->periodic
	    );


  //pvfmm::mem::MemoryManager mem_mgr(10000000);
  //impl->mat        = new PVFMMImpl::Mat_t(&mem_mgr);
  if (!impl->dense_eval) {
    impl->mat = new PVFMMImpl::Mat_t();

    impl->mat->Initialize(impl->mult_order, impl->comm, impl->ker);
  }
  COUTDEBUG("Finished");

  _impl = std::move(impl);
}
void PvFMM::_store_problem_data(Vec source_positions, Vec source_normals,
                               Vec target_positions, Kernel3d kernel) {

  DblNumMat srcPos = get_local_vector(DIM, num_local_points(source_positions), source_positions);
  DblNumMat srcNor = get_local_vector(DIM, num_local_points(source_normals), source_normals);
  DblNumMat trgPos = get_local_vector(DIM, num_local_points(target_positions), target_positions);

  DblNumMat *srcPos_ptr = &srcPos;
  DblNumMat *srcNor_ptr = &srcNor;
  DblNumMat *trgPos_ptr = &trgPos;

  set_kernel(kernel);
  set_src_positions(srcPos_ptr);
  set_src_normals(srcNor_ptr);
  set_trg_positions(trgPos_ptr);

  srcPos.restore_local_vector();
  srcNor.restore_local_vector();
  trgPos.restore_local_vector();
}

PvFMM::PvFMM(Vec source_positions, Vec source_normals,
        Vec target_positions, Kernel3d kernel){

    
    _store_problem_data(source_positions, source_normals, target_positions, kernel);
  
    _setup_context();
  
    // Scale sources and targets to [0,1] x [0,1] x [0,1]
    //Note: resulting scaled positions are restored in _src_positions and
    //_trg_positions in order to provide access to points in wrapper lib
    _scale_factor = scale_sources_and_targets(_src_positions, _trg_positions);
    
    _kernel_opts.scale_factor = _scale_factor;
    
}


PvFMM::~PvFMM(){

}

double PvFMM::scale_sources_and_targets(vec sources, vec targets){
    double bbox_half_width = 1.0;
    //double eps = 0.05;
    double eps = 0.2;
    vec srcs(sources);
    vec trgs(targets);
    int m = DIM;
    double x;
    bool need_to_scale = false;
    
    bbox_half_width = 0.0;
    cout << "begin scaling" << endl; 
    // TODO add openmp
    //Compute the half-width of the bounding box of size 
    // [-bbox_half_width, bbox_half_width]^3
    // that encloses all source and targets
    for(int i = 0; i < m; i++){
        for(int j = 0; j <_num_sources; j++){
            x = srcs[m*j + i];
            bbox_half_width = fabs(x) > bbox_half_width ? fabs(x) : bbox_half_width;
        }
        // Same procedure for target points
        for(int j = 0; j <_num_targets; j++){
            x = trgs[m*j + i];
            bbox_half_width = fabs(x) > bbox_half_width ? fabs(x) : bbox_half_width;
        }
   }
   double buffer;
   MPI_Allreduce(&bbox_half_width, // send
                 &buffer,          // recieve
                 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   int r;
   MPI_Comm_rank(PETSC_COMM_WORLD, &r);
   cout << "bbox_half_width rank " << r << ", " << buffer << endl;
   // increase scaling by a small amount to ensure points are strictly
   // contained in [0,1)^3 instead of [0,1]^3
   bbox_half_width = buffer + eps;
   // Now rescale all sources and targets to by  twice this value (half-width
   // *2 = width) and translate by (min_x, min_y, min_z) to bring points
   // within in the box [-.5,.5]^3, and translate by (.5, .5, .5) to move unit
   // box to [0,1]^3
   for (int i = 0; i < m; i++) {
     for (int j = 0; j < _num_sources; j++) {
       x = srcs[m * j + i];
       srcs[m * j + i] = x / (2 * bbox_half_width) + .5;
     }
     // Same procedure for target points
     for (int j = 0; j < _num_targets; j++) {
       x = trgs[m * j + i];
       trgs[m * j + i] = x / (2 * bbox_half_width) + .5;
     }
   }
    bbox_half_width *= 2;

    _src_positions = srcs;
    _trg_positions = trgs;
    return bbox_half_width;
}


// Call PvFMM's setFromOptions() functions
int PvFMM::setFromOptions(){
    return 0; //change to ebiFunctionReturn(0);
}


void PvFMM::evaluate(const Vec& den, Vec& val){
//
  ebiFunctionBegin;
  double zero = 0.0;
  iC( VecSet( val, zero) );

  // density and value -- put DblNumVec interface on arrays obtained from Vec
  // @BUG also looks like it cannot work unless it is a single processor
  int64_t dl;
  VecGetLocalSize(den, &dl);

  double* den_arr;
  iC( VecGetArray(den, &den_arr) );
  double* dbuf=den_arr;
  VecRestoreArray(den, &den_arr);
  DblNumVec srcDen(dl, false, dbuf);

  int64_t vl;
  VecGetLocalSize(val, &vl);

  double* valarr;
  iC( VecGetArray(val, &valarr) );
  double* vbuf=valarr;
  VecRestoreArray(val, &valarr);
  DblNumVec trgVal(vl, false, vbuf);

  evaluate(srcDen, trgVal);


  iC( VecGetArray(val, &valarr) ); 
  for (int i = 0; i < vl; i++){                                                    
      // copy into Petsc vector and unscale
      valarr[i] = trgVal(i);
  }                                                                                

  VecRestoreArray(val, &valarr);
}

void PvFMM::interaction_matrix(DblNumMat source_positions,
                              DblNumMat source_normals,
                              DblNumMat target_positions,
                              DblNumMat &interaction_matrix) {
  int num_sources = source_positions.n();
  int num_targets = target_positions.n();
  // context object
  COUTDEBUG("pvfmm_interaction_matrix with nsrc=" << num_sources
                                                  << ", ntrg=" << num_targets);

  COUTDEBUG("Using context with"
            << " mult_order=" << _impl->mult_order
            << ", max_pts=" << _impl->max_pts
            << ", max_depth=" << _impl->max_depth << ", sl=" << _impl->sl
            << ", dl=" << _impl->dl << ", periodic=" << _impl->periodic);
  assert(num_sources == source_normals.n());
  assert(DIM == source_normals.m());
  assert(DIM == source_positions.m());
  assert(DIM == target_positions.m());
  assert(interaction_matrix.m() == num_targets * _target_dof);
  assert(interaction_matrix.n() == num_sources * _source_dof);
  vector<double> normals((DIM + _source_dof) * num_sources, 0.);

  for (int i = 0; i < num_sources; i++) {
    for (int j = 0; j < DIM; j++) {
      normals[j + i * (DIM + _source_dof)] = source_normals(j, i);
    }
  }

  // yanked from Kernel<T>::BuildMatrix in kernel.txx
  //(T* r_src, int src_cnt, T* r_trg, int trg_cnt, T* k_out) const
  if (_impl->sl && !_impl->dl) {

    for (int i = 0; i < num_sources; i++) {
      for (int j = 0; j < _source_dof; j++) {
        std::vector<real_t> v_src(_source_dof, 0.);
        v_src[j] = 1.0;
        _impl->ker->ker_poten(
                source_positions.data() + DIM * i, 1, 
                &v_src[0], 1, 
                target_positions.data(), num_targets,
                interaction_matrix.data() + (i * _source_dof + j) * num_targets * _target_dof,
                NULL);
      }
    }
  } else if (_impl->dl && !_impl->sl) {

    for (int i = 0; i < num_sources; i++) {
      for (int j = 0; j < _source_dof; j++) {
        std::vector<real_t> v_src(_source_dof + DIM, 0.);
        for (int d = 0; d < DIM; d++) {
          v_src[d] = normals[(DIM + _source_dof) * i + d];
        }
        v_src[DIM + j] = 1.0;
        _impl->ker->dbl_layer_poten(
            source_positions.data() + DIM * i, 1, 
            &v_src[0], 1,
            target_positions.data(), num_targets,
            interaction_matrix.data() + (i * _source_dof + j) * num_targets * _target_dof,
            NULL);
      }
    }
  }

}
void PvFMM::evaluate_direct(const Vec& srcDen, Vec& trgVal){
    DblNumVec density(srcDen);
    DblNumVec potential(trgVal);
    evaluate_direct(density, potential);
}

void PvFMM::_interleave_density(DblNumVec srcDen){
    for (int i = 0; i < _num_sources; i++) {
      for (int j = 0; j < _source_dof; j++) {
        _src_normals[DIM + j + i * (DIM + _source_dof)] = srcDen(i * _source_dof + j);
      }
    }
}
void PvFMM::evaluate_direct(const DblNumVec &srcDen, DblNumVec &trgVal) {
  vec source_density;
  if (_kernel_layer == DOUBLE_LAYER) {
    // Copy source density into source_normal array
    
    _interleave_density(srcDen);
    source_density = _src_normals;

  } else if (_kernel_layer == SINGLE_LAYER) {
    source_density =
        vec(srcDen.data(), srcDen.data() + _source_dof * _num_sources);

  } else {
    ebiAssert(false);
  }
    vec target_potential(_target_dof*_num_targets, 0.);
    const size_t num_sources = _src_positions.size()/DIM;
    const size_t ntrg = _trg_positions.size()/DIM;
  ASSERT(_impl, "was called with null context. "
                "Generate a context by calling make_pvfmm_context.");

  // context object
  COUTDEBUG("pvfmm_direct_summation with num_source_points="
            << num_sources << ", ntrg=" << ntrg);

  COUTDEBUG("Using context with"
            << " mult_order=" << _impl->mult_order
            << ", max_pts=" << _impl->max_pts
            << ", max_depth=" << _impl->max_depth << ", sl=" << _impl->sl
            << ", dl=" << _impl->dl << ", periodic=" << _impl->periodic);

  int source_dof = _impl->source_dof;
  int target_dof = _impl->target_dof;
  // copy data
  
  assert(_impl->dense_eval);

  pvfmm::Kernel<real_t>::Ker_t kernel_function;

  if (_impl->sl && !_impl->dl) {
    kernel_function = _impl->ker->ker_poten;

  } else if (!_impl->sl && _impl->dl) {
    kernel_function = _impl->ker->dbl_layer_poten;
  } else {
    abort();
  }
  int omp_p = omp_get_max_threads();
#pragma omp parallel for
  for (int i = 0; i < omp_p; i++) {
    size_t a = (i * ntrg) / omp_p;
    size_t b = ((i + 1) * ntrg) / omp_p;

    kernel_function(_src_positions.data(), num_sources, source_density.data(), 1,
                    _trg_positions.data() + a * 3, b - a,
                    trgVal.data() + a * target_dof, NULL);
  }
  // remove
  //_copy_potential(_num_targets, _target_dof, target_potential, trgVal.data(), true);
  _unscale_potential(_num_targets, _target_dof, trgVal.data());
  
  _rebuild_tree = false;
 

  cout << "fmm return" << endl;
}

void PvFMM::_setup_tree(vec source_density) {
  COUTDEBUG("(Re)building sldl pvfmm tree");
  ASSERT(_impl->mat != NULL, "badly initialized context");

  // FMM Setup
  if (_rebuild_tree && _impl->tree != NULL) {
    delete _impl->tree;
  }
  vec empty;
      // different tree setup arguments passed for single- vs. double-layer.
  if (_impl->dl) {

    _impl->tree = pvfmm::PtFMM_CreateTree<double>(
        empty, empty, _src_positions, source_density, _trg_positions,
        _impl->comm, _impl->max_pts, _impl->boundary);
  } else {

    _impl->tree = pvfmm::PtFMM_CreateTree<double>(
        _src_positions, source_density, empty, empty, _trg_positions,
        _impl->comm, _impl->max_pts, _impl->boundary);
  }

  _impl->tree->SetupFMM(_impl->mat);
}

void PvFMM::_interface_evaluate(vec source_density, vec &potv) {
  int ntrg = _trg_positions.size() / DIM;
  ASSERT(_impl, "was called with null context. "
                "Generate a context by calling make_pvfmm_context.");

  // context object
  int nsrc = _impl->dl?source_density.size()/(_source_dof+DIM): source_density.size()/(_source_dof);
  COUTDEBUG("pvfmm evaluation with nsrc="
            <<  nsrc << ", ntrg=" << ntrg
            << ", rebuild_tree=" << _rebuild_tree);

  COUTDEBUG("Using context with"
            << " mult_order=" << _impl->mult_order
            << ", max_pts=" << _impl->max_pts
            << ", max_depth=" << _impl->max_depth << ", sl=" << _impl->sl
            << ", dl=" << _impl->dl << ", periodic=" << _impl->periodic);
  if (_rebuild_tree) {
    // Run FMM
    pvfmm::PtFMM_Evaluate<double>(_impl->tree, potv, ntrg);
  } else {
      // different evaluation arguments passed for single- vs. double-layer.
    if (_impl->dl) {
      _impl->tree->ClearFMMData();
      pvfmm::PtFMM_Evaluate<double>(_impl->tree, potv, ntrg, NULL,
                                    &source_density);
    } else {
      pvfmm::PtFMM_Evaluate<double>(_impl->tree, potv, ntrg, &source_density,
                                    NULL);
    }
  }
}

void PvFMM::evaluate(const DblNumVec &srcDen, DblNumVec &trgVal) {
  int ntrg = _trg_positions.size() / DIM;
  vec potv(_target_dof * ntrg);

  vec source_density;
  if (_kernel_layer == DOUBLE_LAYER) {
    // Copy source density into source_normal array

    _interleave_density(srcDen);
    source_density = _src_normals;

  } else if (_kernel_layer == SINGLE_LAYER) {
    source_density =
        vec(srcDen.data(), srcDen.data() + _source_dof * _num_sources);
  } else {
    abort();
  }
  _setup_tree(source_density);
  _interface_evaluate(source_density, potv);

  _copy_potential(ntrg, _target_dof, potv, trgVal.data(), true);

  _rebuild_tree = false;

  cout << "fmm return" << endl;
}

void PvFMM::set_src_positions(DblNumMat*& src_positions) { 
    _num_sources = src_positions->n();
    vec temp(src_positions->data(),
                src_positions->data() + 
                src_positions->m()*src_positions->n());
    _src_positions = temp;
}

 
void PvFMM::set_src_normals(DblNumMat*& src_normals) {
    ebiAssert(src_normals->n() == _num_sources);
    vec source_normals_padded(_num_sources*(DIM+_source_dof), 0.);

    for(int j = 0; j < _num_sources; j++){
        for(int i =0; i < DIM; i++){
            source_normals_padded[j*(DIM+_source_dof)+i] = (*src_normals)(i,j);
        }
    }
    _src_normals = source_normals_padded;
}

 
void PvFMM::set_trg_positions(DblNumMat*& trg_positions) {
    _num_targets = trg_positions->n();
    _markgrid_targets = trg_positions;
    vec temp(trg_positions->data(),
            trg_positions->data() + trg_positions->m()*trg_positions->n());
    _trg_positions = temp;
}

 
void PvFMM::set_kernel(Kernel3d& kernel){
    _kernel_type = kernel.kernelType(); //make this safer
    
    Kernel3d::parse_kernel_enum((int)_kernel_type, _equation_type, 
            _kernel_layer, _kernel_variable);

    if (_equation_type == MOD_HELMHOLTZ){
        _kernel_opts.helmholtz_frequency = kernel.coefs(1);
        _kernel_opts.initialized = true;
    } else if (_equation_type == NAVIER){
        _kernel_opts.navier_mu = kernel.coefs(0);
        _kernel_opts.navier_nu = kernel.coefs(1);
        _kernel_opts.initialized = true;
    } else {
        _kernel_opts.initialized = false;
    }

    _source_dof = kernel.srcDOF();
    _target_dof = kernel.trgDOF();

}

END_EBI_NAMESPACE
