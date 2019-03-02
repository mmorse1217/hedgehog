//#include <mpi.h>
//#include <omp.h>
//#include <iostream>

#include <pvfmm.hpp>
#include <math_utils.hpp>

template <class T>
void laplace_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif
  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
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
        T* f=&v_src[s];


        p[0] += f[0]*invR * SCAL_CONST;
      }
    }
    k_out[t] += p[0];
  }
}

const pvfmm::Kernel<real_t> ker_laplace_sl_m2l=pvfmm::BuildKernel<real_t, laplace_sl_m2l>("laplace_sl_m2l", 3, std::pair<int,int>(4,3));

template <class T>
void laplace_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
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
        T* f=&v_src[s];


        p[0] += f[0]*invR;
      }
    }
    k_out[t] += p[0]*SCAL_CONST;
  }
}
template <class T>
void laplace_sl_eps(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

      //if (R!=0){
      //if (R>=1e-6){
        if (R>=1e-8){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T* f=&v_src[s];


        p[0] += f[0]*invR;
      }
      /*if(s == 68245 || t == 68245){
          std::cout << "source: " << r_src[3*s] << ", " << r_src[3*s+1] << ", "<< r_src[3*s+2] << std::endl;
          std::cout << "target: " << r_trg[3*t] << ", " << r_trg[3*t+1] << ", "<< r_trg[3*t+2] << std::endl;
          std::cout << dR[0] << ", " << dR[1] << ", "<< dR[2] << std::endl;
          std::cout << R << ", "  << p[0] << std::endl;
      }*/
    }
    k_out[t] += p[0]*SCAL_CONST;
  }
}

template <class T>
void laplace_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
  const T SCAL_CONST = -1.0/(4.0*pvfmm::const_pi<T>());
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

          T* n=&v_src[(s*dof+i)*4+0];
          T* f=&v_src[(s*dof+i)*4+3];

          T r_dot_n=(n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2]);

          p[0] += invR3*f[0]*r_dot_n;
        }
      
      }
      k_out[(t*dof+i)] += p[0]*SCAL_CONST;
    }
  }
}


template <class T>
void laplace_dl_eps(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
  const T SCAL_CONST = -1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    for(int i=0;i<dof;i++){
      T p[3]={0,0,0};
      for(int s=0;s<src_cnt;s++){
        T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                 r_trg[3*t+1]-r_src[3*s+1],
                 r_trg[3*t+2]-r_src[3*s+2]};
        T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

        /*if (sqrt(R)<= 1e-3){
          std::cout.precision(16);
          std::cout << "source: " << r_src[3*s] << ", " << r_src[3*s+1] << ", "<< r_src[3*s+2] << std::endl;
          std::cout << "target: " << r_trg[3*t] << ", " << r_trg[3*t+1] << ", "<< r_trg[3*t+2] << std::endl;
          std::cout << sqrt(R) << std::endl;

        }*/
        //if (R!=0){
        //if (R>=1e-8){
        if (R>=1e-12){
          T invR2=1.0/R;
          T invR=sqrt(invR2);
          T invR3=invR2*invR;

          T* n=&v_src[(s*dof+i)*4+0];
          T* f=&v_src[(s*dof+i)*4+3];

          T r_dot_n=(n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2]);

          p[0] += invR3*f[0]*r_dot_n;
        }
      /*if(s == 68245 || t == 68245){
          std::cout.precision(16);
          std::cout << "source: " << r_src[3*s] << ", " << r_src[3*s+1] << ", "<< r_src[3*s+2] << std::endl;
          std::cout << "target: " << r_trg[3*t] << ", " << r_trg[3*t+1] << ", "<< r_trg[3*t+2] << std::endl;
          std::cout << dR[0] << ", " << dR[1] << ", "<< dR[2] << std::endl;
          std::cout << R << ", "  << p[0] << std::endl;
      }*/
      }
      k_out[(t*dof+i)] += p[0]*SCAL_CONST;
    }
          //std::cout << "potential: " << k_out[t] << std::endl;
  }
}

const real_t default_eps2(1e-3);

template <class T, const T* regularization_eps=&default_eps2>
void laplace_sl_s2t(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      //R  += *regularization_eps;

      if (R>=*regularization_eps){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T invR3=invR2*invR;
        T* f=&v_src[s];

        /*T inner_prod=(f[0]*dR[0] +
                      f[1]*dR[1] +
                      f[2]*dR[2])* invR3;*/

        p[0] += f[0]*invR * SCAL_CONST;
      }
    }
    k_out[t] += p[0];
    /*
        p[0] += f[0]*invR + dR[0]*inner_prod;
        p[1] += f[1]*invR + dR[1]*inner_prod;
        p[2] += f[2]*invR + dR[2]*inner_prod;
      }
    }
    k_out[t*3+0] += p[0]*SCAL_CONST;
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
    */
  }
}

const pvfmm::Kernel<real_t> ker_laplace_sl_eps=pvfmm::BuildKernel<real_t, laplace_sl_eps>(
        "laplace_sl_eps",
        3,
        std::pair<int,int>(1,1));
const pvfmm::Kernel<real_t> ker_laplace_dl_eps=pvfmm::BuildKernel<real_t, laplace_dl_eps>(
        "laplace_dl_eps",
        3,
        std::pair<int,int>(1,1));
const pvfmm::Kernel<real_t> ker_laplace_dl_temp=pvfmm::BuildKernel<real_t, laplace_dl>(
        "laplace_dl_temp",
        3,
        std::pair<int,int>(1,1));
// This is a hack for regularization of the kernel
const pvfmm::Kernel<real_t> ker_laplace_sl_s2t=pvfmm::BuildKernel<real_t, laplace_sl_s2t>(
        "laplace_sl_s2t",
        3,
        std::pair<int,int>(1,1));

const pvfmm::Kernel<real_t> ker_laplace_sl=pvfmm::BuildKernel<real_t, laplace_sl>(
        "laplace_sl"		/* name  */,
	3			/* dim   */,
	std::pair<int,int>(1,1)	/* k_dim */,
	NULL			/* s2m   */,
	NULL			/* s2l   */,
	NULL        	/* s2t   */,
	NULL	/* m2m   */,
	NULL	/* m2l   */,
	NULL	/* m2t   */,
	NULL			/* l2l   */,
	NULL			/* l2t   */);


// This one works
const pvfmm::Kernel<real_t> ker_laplace_dl=pvfmm::BuildKernel<real_t, laplace_sl, laplace_dl_eps>(
        "laplace_sldl"		/* name  */,
	3			/* dim   */,
	std::pair<int,int>(1,1)	/* k_dim */,
	NULL 			/* s2m   */,
	NULL 		/* s2l   */,
	NULL /* s2t   */,
	NULL	           /* m2m   */,
	NULL               /* m2l   */,
	NULL   /* m2t   */,
	NULL   /* l2l   */,
	NULL   /* l2t   */);













