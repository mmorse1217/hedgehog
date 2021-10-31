#include <pvfmm.hpp>
#include <math_utils.hpp>


inline double helmholtz_frequency(double in_frequency=0){
    static double frequency;
    if (in_frequency != 0){
        frequency = in_frequency;
    }
    assert(frequency != 0);
    return frequency;
}

inline const char* helmholtz_name(bool single_layer, bool m2l=false){
    static std::string name;
    std::ostringstream os;
    os << "mod_helmholtz_";
    if(single_layer){
        os << "sl_";
    } else {
        os << "sldl_";
    }
    if (m2l){
        os << "m2l_";
    }
    os << "freq_" << helmholtz_frequency();
    name = os.str();
    assert(!name.empty());
    return name.c_str();
}


template <class T>
void mod_helmholtz_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif
  T mu = helmholtz_frequency();
  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[1]={0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR = 1.0/R;
        T e_mu_R = exp(-mu * R);

        T f = v_src[s];

        p[0] += SCAL_CONST * invR * e_mu_R * f;
      }
    }
    k_out[t] += p[0];
  }
}


template <class T>
void mod_helmholtz_dl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif

  T mu = helmholtz_frequency();
  const T SCAL_CONST = -1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[1]={0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR = 1.0/R;
        T invR2 = invR*invR;
        T invR3 = invR*invR2;
        T e_mu_R = exp(-mu * R);
        
        T* n= &v_src[4*s]; //normal vector \in \mathbb{R}^3
        T f = v_src[4*s+3]; // source density is a scalar
        T r_dot_n = n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2];
        p[0] += SCAL_CONST * (invR3 + mu*invR2) * e_mu_R * r_dot_n * f;
      }
    }
    k_out[t] += p[0];
  }
}


const pvfmm::Kernel<real_t> ker_mod_helmholtz_sl_m2l=pvfmm::BuildKernel<real_t, mod_helmholtz_sl_m2l>("mod_helmholtz_sl_m2l", 3, std::pair<int,int>(1,1));

const pvfmm::Kernel<real_t> ker_mod_helmholtz_dl_m2l=pvfmm::BuildKernel<real_t, mod_helmholtz_dl_m2l>("mod_helmholtz_dl_m2l", 3, std::pair<int,int>(1,1));

template <class T>
void mod_helmholtz_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  T mu = helmholtz_frequency();

  const T SCAL_CONST = 1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[1]={0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR = 1.0/R;
        T e_mu_R = exp(-mu * R);

        T f = v_src[s];

        p[0] += SCAL_CONST * invR * e_mu_R * f;
      }
    }
    k_out[t] += p[0];
  }
}

template <class T>
void mod_helmholtz_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif

  T mu = helmholtz_frequency();
  const T SCAL_CONST = -1.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    T p[1]={0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR = 1.0/R;
        T invR2 = invR*invR;
        T invR3 = invR*invR2;
        T e_mu_R = exp(-mu * R);
        
        T* n= &v_src[4*s]; //normal vector \in \mathbb{R}^3
        T f = v_src[4*s+3]; // source density is a scalar
        T r_dot_n = n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2];
        p[0] += SCAL_CONST * (invR3 + mu*invR2) * e_mu_R * r_dot_n * f;
      }
    }
    k_out[t] += p[0];
  }
}



const pvfmm::Kernel<real_t> ker_mod_helmholtz_sl=pvfmm::BuildKernel<real_t, mod_helmholtz_sl>(
        "mod_helmholtz_sl"		/* name  */,
	3			                /* dim   */,
	std::pair<int,int>(1,1)	    /* k_dim */,
	NULL			            /* s2m   */,
	NULL			            /* s2l   */,
	NULL        	            /* s2t   */,
	&ker_mod_helmholtz_sl_m2l	/* m2m   */,
	&ker_mod_helmholtz_sl_m2l	/* m2l   */,
	&ker_mod_helmholtz_sl_m2l	/* m2t   */,
	NULL			            /* l2l   */,
	NULL			            /* l2t   */);


const pvfmm::Kernel<real_t> ker_mod_helmholtz_dl=pvfmm::BuildKernel<real_t, mod_helmholtz_sl, mod_helmholtz_dl>(
    "mod_helmholtz_sldl"		/* name  */,
	3               			/* dim   */,
	std::pair<int,int>(1,1)	    /* k_dim */,
	NULL		            	/* s2m   */,
	NULL			            /* s2l   */,
	NULL	                    /* s2t   */,
	&ker_mod_helmholtz_sl_m2l	/* m2m   */,
	&ker_mod_helmholtz_sl_m2l	/* m2l   */,
	&ker_mod_helmholtz_sl_m2l	/* m2t   */,
	NULL		            	/* l2l   */,
	NULL			            /* l2t   */);

