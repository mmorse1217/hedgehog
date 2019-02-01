#include "pvfmm_bis_interface.hpp" //interface in fmm3d
#include <pvfmm_interface.h> //interface in pvfmm/pvfmm_utils/
#include <vec3t.hpp>
#include "legacy/syms.hpp"
#include "common/utils.hpp"
BEGIN_EBI_NAMESPACE

using namespace std;
using namespace Ebi;

 
PvFMM::PvFMM(){
        _num_points = 0;
}
/*
PvFMM::PvFMM(Vec source_positions, Vec source_normals,
                Vec target_positions, Kernel3d kernel){
    FMM::FMM(source_positions, source_normals, target_positions, kernel);

}
*/

// Delete internal FMM object
PvFMM::~PvFMM(){
    clear_pvfmm_context(&_context);
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

    MPI_Allreduce(&bbox_half_width, //send
            &bbox_half_width, //recieve
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD);
    // increase scaling by a small amount to ensure points are strictly
    // contained in [0,1)^3 instead of [0,1]^3
    bbox_half_width = bbox_half_width +  eps;
    // Now rescale all sources and targets to by  twice this value (half-width
    // *2 = width) and translate by (min_x, min_y, min_z) to bring points 
    // within in the box [-.5,.5]^3, and translate by (.5, .5, .5) to move unit
    // box to [0,1]^3
    for(int i = 0; i < m; i++){
        for(int j = 0; j <_num_sources; j++){
            x = srcs[m*j + i];
            srcs[m*j + i]   = x/(2*bbox_half_width) + .5;
        }
        // Same procedure for target points
        for(int j = 0; j <_num_targets; j++){
            x = trgs[m*j + i];
            trgs[m*j + i]   = x/(2*bbox_half_width) + .5;
        }
   }
    bbox_half_width *= 2;

    _src_positions = srcs;
    _trg_positions = trgs;
    return bbox_half_width;
}

int PvFMM::setup(){
    
    int sl, dl, dim, periodic; 

    
    periodic = false; 

    ebiAssert(_kernel_type != 0);
    
    ebiAssert(_equation_type == LAPLACE || _equation_type == STOKES 
            || _equation_type == MOD_HELMHOLTZ || _equation_type == NAVIER);
    if (_kernel_layer == SINGLE_LAYER){
        sl = 1; 
        dl = 0;

    } else if (_kernel_layer == DOUBLE_LAYER){
        sl = 0;
        dl = 1;
    }
    _rebuild_tree = true; 

    // _multipole_order: square root of the number of points on one face of check
    // surface (TODO double check)
    cout << "kernel: " << _kernel_type << endl;
    cout << "m:" << _multipole_order<< endl;
    make_pvfmm_context(_multipole_order, _max_points_per_box, _max_level, sl, dl,
                        periodic, _kernel_type, _source_dof, _target_dof, _direct_eval,
                        &_context, &_kernel_opts);
    return 0; 
}

// Call PvFMM's setFromOptions() functions
int PvFMM::setFromOptions(){
    // Scale sources and targets to [0,1] x [0,1] x [0,1]
    //Note: resulting scaled positions are restored in _src_positions and
    //_trg_positions in order to provide access to points in wrapper lib
    _scale_factor = scale_sources_and_targets(_src_positions, _trg_positions);
    
    _direct_eval = Options::get_int_from_petsc_opts("-direct_eval");

    _kernel_opts.scale_factor = _scale_factor;
    return 0; //change to ebiFunctionReturn(0);
}


int PvFMM::evaluate(const Vec& den, Vec& val){
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
    ebiFunctionReturn(0);
}

int PvFMM::interaction_matrix(DblNumMat source_positions, DblNumMat source_normals, 
        DblNumMat target_positions, DblNumMat& interaction_matrix){
    int num_sources = source_positions.n();
    int num_targets = target_positions.n();
    assert(num_sources == source_normals.n());
    assert(DIM == source_normals.m());
    assert(DIM == source_positions.m());
    assert(DIM == target_positions.m());
    assert(interaction_matrix.m() == num_targets*_target_dof);
    assert(interaction_matrix.n() == num_sources*_source_dof);
    vector<double> normals((DIM+_source_dof)*num_sources,0.);
        for(int i = 0; i < num_sources; i++){
            for (int j = 0; j < DIM; j++) {
                normals[j + i*(DIM+_source_dof)] = source_normals(j,i);
            }
        }
        pvfmm_interaction_matrix( //Single layer source data
                            num_sources, source_positions.data(),
                           NULL, 
                           // Double layer source data
                           num_sources, source_positions.data(),
                           normals.data(),
                           // Single/Double layer target data
                            num_targets, target_positions.data(),
                            interaction_matrix.data(),
                            //PvFMM tree data
                            _rebuild_tree, &_context);
    return 0;
}
int PvFMM::evaluate_direct(const Vec& srcDen, Vec& trgVal){
    DblNumVec density(srcDen);
    DblNumVec potential(trgVal);
    evaluate_direct(density, potential);
    return 0;
}
int PvFMM::evaluate_direct(const DblNumVec& srcDen, DblNumVec& trgVal){
    if (_kernel_layer == DOUBLE_LAYER){
        //Copy source density into source_normal array
        for(int i = 0; i < _num_sources; i++){
          for(int j =0; j < _source_dof; j++){
            _src_normals[DIM +  j + i*(DIM+_source_dof)] = srcDen(i*_source_dof + j);
          }
        }
        pvfmm_direct_summation( //Single layer source data
                            _num_sources, _src_positions.data(),
                           srcDen.data(), 
                           // Double layer source data
                           _num_sources, _src_positions.data(),
                           _src_normals.data(),
                           // Single/Double layer target data
                            _num_targets, _trg_positions.data(),
                            trgVal.data(),
                            //PvFMM tree data
                            _rebuild_tree, &_context);
    } else if ( _kernel_layer == SINGLE_LAYER){
        pvfmm_direct_summation( //Single layer source data
                            _num_sources, _src_positions.data(),
                           srcDen.data(), 
                           // Double layer source data
                           NULL, NULL,
                           NULL,
                           // Single/Double layer target data
                            _num_targets, _trg_positions.data(),
                            trgVal.data(),
                            //PvFMM tree data
                            _rebuild_tree, &_context);
    } else {
        ebiAssert(false);
    }
    _rebuild_tree = false;
  Kernel3d temp(_kernel_type, vector<double>());
  double density_unscale_factor = temp.density_unscaling_value(_scale_factor);

  cout << "density unscale: " << density_unscale_factor << endl;
    // Scale density from PvFMM to old KIFMM assumption
    for (int i= 0; i < _num_targets; i++){
        for (int j= 0; j < _target_dof; j++){
            trgVal(j + i*_target_dof) = trgVal(j + i*_target_dof)*density_unscale_factor;
        }
    }


  cout << "fmm return"  << endl;
    return 0;
}


int PvFMM::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal){
    if (_kernel_layer == DOUBLE_LAYER){
        //Copy source density into source_normal array
        for(int i = 0; i < _num_sources; i++){
          for(int j =0; j < _source_dof; j++){
            _src_normals[DIM +  j + i*(DIM+_source_dof)] = srcDen(i*_source_dof + j);
          }
        }
        laplace_sldl_pvfmm( //Single layer source data
                            _num_sources, _src_positions.data(),
                           srcDen.data(), 
                           // Double layer source data
                           _num_sources, _src_positions.data(),
                           _src_normals.data(),
                           // Single/Double layer target data
                            _num_targets, _trg_positions.data(),
                            trgVal.data(),
                            //PvFMM tree data
                            _rebuild_tree, &_context);
    } else if ( _kernel_layer == SINGLE_LAYER){
        laplace_sldl_pvfmm( //Single layer source data
                            _num_sources, _src_positions.data(),
                           srcDen.data(), 
                           // Double layer source data
                           NULL, NULL,
                           NULL,
                           // Single/Double layer target data
                            _num_targets, _trg_positions.data(),
                            trgVal.data(),
                            //PvFMM tree data
                            _rebuild_tree, &_context);
    } else {
        ebiAssert(false);
    }
    _rebuild_tree = false;
  Kernel3d temp(_kernel_type, vector<double>());
  double density_unscale_factor = temp.density_unscaling_value(_scale_factor);

  cout << "density unscale: " << density_unscale_factor << endl;
    // Scale density from PvFMM to old KIFMM scaling 
    for (int i= 0; i < _num_targets; i++){
        for (int j= 0; j < _target_dof; j++){
            trgVal(j + i*_target_dof) = trgVal(j + i*_target_dof)*density_unscale_factor;
        }
    }


  cout << "fmm return"  << endl;
    return 0;
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
    DblNumMat src_normals_transpose(DIM + _source_dof, _num_sources );

    // TOOD remove
    for(int i =0; i < DIM; i++){
        for(int j = 0; j < _num_sources; j++){
            src_normals_transpose(i,j) = (*src_normals)(i,j);
        }
    }
    vec temp(src_normals_transpose.data(),
              src_normals_transpose.data() + 
              src_normals_transpose.m()*src_normals_transpose.n());
    _src_normals = temp;
}

 
void PvFMM::set_trg_positions(DblNumMat*& trg_positions) {
    _num_targets = trg_positions->n();
    _markgrid_targets = trg_positions;
    vec temp(trg_positions->data(),
            trg_positions->data() + trg_positions->m()*trg_positions->n());
    _trg_positions = temp;
}

 
int& PvFMM::set_num_points(int64_t& num_points) {
    _num_points = num_points;
}
 
void PvFMM::set_kernel(Kernel3d& kernel){
    _kernel_type = (Kernel) kernel.kernelType(); //make this safer
    
    Kernel3d::parse_kernel_enum((int)_kernel_type, _equation_type, 
            _kernel_layer, _kernel_variable);

    if (_equation_type == MOD_HELMHOLTZ){
        _kernel_opts.helmholtz_frequency = kernel.coefs(1);
    } else if (_equation_type == NAVIER){
        _kernel_opts.navier_mu = kernel.coefs(0);
        _kernel_opts.navier_nu = kernel.coefs(1);
    }

    _source_dof = kernel.srcDOF();
    _target_dof = kernel.trgDOF();

}
END_EBI_NAMESPACE
