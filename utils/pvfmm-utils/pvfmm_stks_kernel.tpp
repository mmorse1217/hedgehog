//#include <mpi.h>
//#include <omp.h>
//#include <iostream>

#include <pvfmm.hpp>
#include <math_utils.hpp>
#include <kernel.hpp>

template <class T>
void stokes_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T invR3=invR2*invR;
        T* f=&v_src[s*4];

        T inner_prod=(f[0]*dR[0] +
                      f[1]*dR[1] +
                      f[2]*dR[2])* invR3;
        
        /*
         * The stress tensor for incompressible Stokes flow is given by
            T_{ijk} = - P_{j}* g_j + 
                        \mu ( \frac{\partial G_{ij} g_j}}{\partial x_k} + 
                              \frac{\partial G_{ji} g_j}}{\partial x_k})
            where P_j = r_j/|r|^3.

            Since the double layer kernel is D_{ij} = T_{ijk} n_k for normal
            vector n_k, and we're using single layer kernel G_{ij}  for
            translation, we need to make sure that the double layer solution 
            can be written as a sum of derivatives of single layer kernels.
            Since it can't, we need we need to add the 
            a P_j*g_j term to the solution below: to add a dimensional dependence
            on pressure to the linear combination of single layer solutions.
            Otherwise, we can't resolve double layer solutions from just single layer
            evaluations.
         */

        T inner_prod_plus_f3_invR3=inner_prod+f[3]*invR3;

        p[0] += f[0]*invR + dR[0]*inner_prod_plus_f3_invR3;
        p[1] += f[1]*invR + dR[1]*inner_prod_plus_f3_invR3;
        p[2] += f[2]*invR + dR[2]*inner_prod_plus_f3_invR3;
      }
    }
    k_out[t*3+0] += p[0];
    k_out[t*3+1] += p[1];
    k_out[t*3+2] += p[2];
  }
}

const pvfmm::Kernel<real_t> ker_stokes_sl_m2l=pvfmm::BuildKernel<real_t, stokes_sl_m2l>("stokes_sl_m2l", 3, std::pair<int,int>(4,3));

template <class T>
void stokes_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T mu=1.0;
  const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*mu);
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

      if (R!=0){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T invR3=invR2*invR;
        T* f=&v_src[s*3];

        T inner_prod=(f[0]*dR[0] +
                      f[1]*dR[1] +
                      f[2]*dR[2])* invR3;

        p[0] += f[0]*invR + dR[0]*inner_prod;
        p[1] += f[1]*invR + dR[1]*inner_prod;
        p[2] += f[2]*invR + dR[2]*inner_prod;
      }
    }
    k_out[t*3+0] += p[0]*SCAL_CONST;
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
  }
}

template <class T>
void stokes_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
  const T mu=1.0;
  const T TOEPMU = -6.0/(8.0*pvfmm::const_pi<T>()*mu);
  for(int t=0;t<trg_cnt;t++){
    for(int i=0;i<dof;i++){
      T p[3]={0,0,0};
      for(int s=0;s<src_cnt;s++){
        T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                 r_trg[3*t+1]-r_src[3*s+1],
                 r_trg[3*t+2]-r_src[3*s+2]};
        T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

        if (R!=0){
          T invR2=1.0/R;
          T invR=sqrt(invR2);
          T invR3=invR2*invR;
          T invR5=invR2*invR3;

          //T* f=&v_src[(s*dof+i)*6+0];
          //T* n=&v_src[(s*dof+i)*6+3];
          T* n=&v_src[(s*dof+i)*6+0];
          T* f=&v_src[(s*dof+i)*6+3];

          T r_dot_n=(n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2]);
          T r_dot_f=(f[0]*dR[0]+f[1]*dR[1]+f[2]*dR[2]);
          T p_=r_dot_n*r_dot_f*invR5;

          p[0] += dR[0]*p_;
          p[1] += dR[1]*p_;
          p[2] += dR[2]*p_;
        }
      }
      k_out[(t*dof+i)*3+0] += p[0]*TOEPMU;
      k_out[(t*dof+i)*3+1] += p[1]*TOEPMU;
      k_out[(t*dof+i)*3+2] += p[2]*TOEPMU;
    }
  }
}

//extern const real_t default_eps2(1e-6);

template <class T>//, const T* regularization_eps=1e-6>
void stokes_sl_s2t(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T mu=1.0;
  const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*mu);
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

      if (R!=0){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T invR3=invR2*invR;
        T* f=&v_src[s*3];

        T inner_prod=(f[0]*dR[0] +
                      f[1]*dR[1] +
                      f[2]*dR[2])* invR3;

        p[0] += f[0]*invR + dR[0]*inner_prod;
        p[1] += f[1]*invR + dR[1]*inner_prod;
        p[2] += f[2]*invR + dR[2]*inner_prod;
      }
    }
    k_out[t*3+0] += p[0]*SCAL_CONST;
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
  }
}


template <class T>
void stokes_press_dl(T* r_src, int src_cnt, T* v_src_, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(17*dof));
#endif

  const T mu=1.0;
  const T OOFP = mu/(2.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
      T p = 0;
      for(int s=0;s<src_cnt;s++){
        T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                 r_trg[3*t+1]-r_src[3*s+1],
                 r_trg[3*t+2]-r_src[3*s+2]};
        T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
        if (R!=0){
          T invR2=1.0/R;
          T invR=pvfmm::sqrt<T>(invR2);
          T invR3=invR2*invR;
          T invR5=invR2*invR3;

          //T* f = &v_src_[6*s + 0];
          //T* n = &v_src_[6*s + 3];
          T* f = &v_src_[6*s + 3];
          T* n = &v_src_[6*s + 0];

          T f_dot_r =  (f[0]*dR[0] +
                        f[1]*dR[1] +
                        f[2]*dR[2]);
          T n_dot_r =  (n[0]*dR[0] +
                        n[1]*dR[1] +
                        n[2]*dR[2]);
          T n_dot_f =  (n[0]*f[0] +
                        n[1]*f[1] +
                        n[2]*f[2]);

          p +=  OOFP*(n_dot_f*invR3 - 3*(n_dot_r*f_dot_r)*invR5);
        }
      }
      k_out[t] += p*OOFP;
    }
}




// This is a hack , don't use it
const pvfmm::Kernel<real_t> ker_stokes_sl_s2t=pvfmm::BuildKernel<real_t, stokes_sl_s2t>(
        "stokes_sl_s2t",
        3,
        std::pair<int,int>(3,3));

const pvfmm::Kernel<real_t> ker_stokes_sl=pvfmm::BuildKernel<real_t, stokes_sl>(
        "stokes_sl"		/* name  */,
	3			/* dim   */,
	std::pair<int,int>(3,3)	/* k_dim */,
	NULL			/* s2m   */,
	NULL			/* s2l   */,
	NULL        	/* s2t   */,
	&ker_stokes_sl_m2l	/* m2m   */,
	&ker_stokes_sl_m2l	/* m2l   */,
	&ker_stokes_sl_m2l	/* m2t   */,
	NULL			/* l2l   */,
	NULL			/* l2t   */);


// This one works
const pvfmm::Kernel<real_t> ker_stokes_dl=pvfmm::BuildKernel<real_t, stokes_sl, stokes_dl>(
        "stokes_sldl"		/* name  */,
	3			/* dim   */,
	std::pair<int,int>(3,3)	/* k_dim */,
	NULL			/* s2m   */,
	NULL			/* s2l   */,
	NULL	                /* s2t   */,
	&ker_stokes_sl_m2l	/* m2m   */,
	&ker_stokes_sl_m2l	/* m2l   */,
	&ker_stokes_sl_m2l	/* m2t   */,
	NULL			/* l2l   */,
	NULL			/* l2t   */);


const pvfmm::Kernel<real_t> ker_stokes_pressure_dl=pvfmm::BuildKernel<real_t, pvfmm::stokes_press, stokes_press_dl>(
        "stokes_press_sldl"		/* name  */,
	3			/* dim   */,
	std::pair<int,int>(3,1));










