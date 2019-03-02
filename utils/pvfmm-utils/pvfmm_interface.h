#ifndef _PVFMMINTERFACE_H_
#define _PVFMMINTERFACE_H_

#include <cstddef>  //for size_t
#include <vector>
#include <map>


typedef double real_t;


// MJM 2/18: TODO clean up this nightmare

enum Kernel{ // from src/ebi/common/kernel3d.hpp
  /*! Laplace kernel - Single Layer */
  KNL_LAP_S_U = 111,
  KNL_LAP_E_U = 112,
  /*! Laplace kernel - Double Layer */
  KNL_LAP_D_U = 121,
  /*! Laplace kernel - Identity Tensor */
  KNL_LAP_I   = 191,
  /*! Modified La Single Layer */
  KNL_MODHEL_S_U = 211,
  /*! Modified Lap Double Layer */
  KNL_MODHEL_D_U = 221,
  /*! Stokes kernel - F Velocity */
  ///kernel used by FMM3d algorithm for stokes equation
  KNL_STK_F_U = 301,
  /*! Stokes kernel - Single Layer Velocity */
  KNL_STK_S_U = 311,
  /*! Stokes kernel - Single Layer Pressure */
  KNL_STK_S_P = 312,
  /*! Stokes kernel - Double Layer Velocity */
  KNL_STK_D_U = 321,
  /*! Stokes kernel - Double Layer Pressure */
  KNL_STK_D_P = 322,
  /*! Stokes kernel - R Velocity */
  KNL_STK_R_U = 331,
  /*! Stokes kernel - R Pressure */
  KNL_STK_R_P = 332,
  /*! Stokes kernel - Identity Tensor */
  KNL_STK_I   = 391,
  /*! Stokes kernel - Levi-Civita Tensor */
  KNL_STK_E   = 392,

  /*! Unsteady Stokes */
  KNL_UNSTK_F_U = 401,
  KNL_UNSTK_S_U = 411,
  /*! Unsteady Stokes kernel - Single Layer Pressure */
  KNL_UNSTK_S_P = 412,
  /*! Unsteady Stokes kernel - Double Layer Velocity */
  KNL_UNSTK_D_U = 421,
  /*! Unsteady Stokes kernel - Double Layer Pressure */
  KNL_UNSTK_D_P = 422,

  //navier kernels  //KNL_NAV_F_U = 501, //used for fmm
  /*! Stokes kernel - Levi-Civita Tensor */
  KNL_NAV_S_U = 511, //single displacement
  KNL_NAV_D_U = 521, //double displacement
  KNL_NAV_R_U = 531,
  KNL_NAV_I   = 591, //identity tensor
  KNL_NAV_E   = 592, //levi-civita tensor
  //other kernels
  KNL_SQRTLAP = 901,
  KNL_EXP     = 902,
  //error
  KNL_ERR = -1
};
extern "C" {

//-----------------------------------------------------------------------------
// Pvfmm Tree functions 
//-----------------------------------------------------------------------------

bool is_leaf(int node_id, void **context);
int depth(int node_id, void **context);
std::vector<int> child_indices(int node_id, void **context);
std::vector<double> box_center(int node_id, void **context);

std::vector<int> top_down_level_order(void **context);
std::vector<int> bottom_up_level_order(void **context);
std::map<int, std::vector<int> > build_level_to_box_ids_map(void **context);
std::vector<int> target_indices_in_box(int node_id, void **context);
std::vector<int> source_indices_in_box(int node_id, void **context);
void initialize_pvfmm_tree(std::vector<double> sl_srcv, 
                           std::vector<double> sl_denv,
                           std::vector<double> dl_srcv,
                           std::vector<double> dl_denv,
                           std::vector<double> trgv,
                           void **context);
void clear_pvfmm_tree(void **context);

struct KernelOptions{
    double scale_factor;
    double navier_mu;
    double navier_nu;
    double helmholtz_frequency;
};
// MPI interface
void mpi_init_pvfmm(int &rank,int &size);
void mpi_finalize_pvfmm();

//-----------------------------------------------------------------------------
// Pvfmm Wrapper functions 
//-----------------------------------------------------------------------------
void make_pvfmm_context(const int &mult_order, const int &max_pts,
          const int &max_depth, const int  &sl,  const int &dl,
          const int &periodic, const Kernel &kernel, const int &source_dof,
          const int &target_dof, const int &dense_eval,
          void **context, const KernelOptions* kernel_opts=NULL);

void clear_pvfmm_context(void **context);

void stokes_sl_pvfmm(const size_t &nsrc, const real_t *src, const real_t *den,
                       const size_t &ntrg, const real_t *trg, real_t *pot,
                       const int &rebuild_tree, void **context);

void laplace_sl_pvfmm(const size_t &nsrc, const real_t *src, 
                    const real_t *den, const size_t &ntrg, const real_t *trg,
                    real_t *pot, const int &rebuild_tree, void **context);


void laplace_sldl_pvfmm(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);
void pvfmm_direct_summation(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);
void pvfmm_interaction_matrix(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);
void stokes_sldl_pvfmm(const size_t &sl_nsrc, const real_t *sl_src, 
                    const real_t *sl_den, const size_t &dl_nsrc, 
                    const real_t *dl_src, const real_t *dl_den_nor,
                    const size_t &ntrg, const real_t *trg, real_t *pot,
                    const int &rebuild_tree, void **context);

void pvfmm_mpi_gather_(void* send_data, int &send_count,
                        void* recv_data, int &recv_count,
                        int &root);

void pvfmm_mpi_scatter_(void* send_data, int &send_count,
                          void* recv_data, int &recv_count,
                          int &root);


void pvfmm_mpi_reduce_(void* send_data, void* recv_data, 
                    int &send_count, int &root);


void pvfmm_mpi_barrier_();

void pvfmm_mpi_bcast_(void* send_data, int &send_count, int &root);
}

#endif //_PVFMMINTERFACE_H_
