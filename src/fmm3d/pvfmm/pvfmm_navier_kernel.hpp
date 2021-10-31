#include <pvfmm.hpp>
#include <math_utils.hpp>
#include <string>

inline double navier_nu(double in_nu=0){
    static double nu;
    if (in_nu != 0){
        nu = in_nu;
    }
    assert(nu != 0);
    return nu;
}

inline double navier_mu(double in_mu=0){
    static double mu;
    if (in_mu != 0){
        mu = in_mu;
    }
    assert(mu != 0);
    return mu;
}

inline const char* navier_name(bool single_layer, bool m2l=false){
    static std::string name;
        std::ostringstream os;
        os << "navier_";
        if(single_layer){
            os << "sl_";
        } else {
            os << "sldl_";
        }
        if (m2l){
            os << "m2l_";
        }
        os << "nu_" << navier_nu() << "_mu_" << navier_mu();
        name = os.str();
    assert(!name.empty());
    return name.c_str();
}

template <class T>
void navier_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif
  T mu = navier_mu();
  T nu = navier_nu();
  const T SCAL_CONST = 1.0/(16.0*pvfmm::const_pi<T>()*mu*(1.-nu));
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0, 0, 0};
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
                                                                                   
        T f_dot_r=(f[0]*dR[0] +                                                 
                   f[1]*dR[1] +                                                 
                   f[2]*dR[2])* invR3;                                          
        
        /*
         * The stress tensor for linear elasticity is given by
           \sigma_{ik} = - P_{j}*g_j + 
                        \mu ( \frac{\partial G_{ij} g_j}}{\partial x_k} + 
                              \frac{\partial G_{ji} g_j}}{\partial x_k})
                          + T_{ijk}'
                      = T_ijk g_j
            where P_j = r_j/|r|^3, T_{ijk} is a mystery tensor that gives 
            stress when contracted with density and DL operator when 
            contracted with normal. T_{ijk}' is some left over terms made up
            of derivatives of single layer kernel (?).
            
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
        f_dot_r += f[3]*invR3;
                                                                                   
        p[0] += (3. - 4.*nu)*f[0]*invR + dR[0]*f_dot_r;
        p[1] += (3. - 4.*nu)*f[1]*invR + dR[1]*f_dot_r;
        p[2] += (3. - 4.*nu)*f[2]*invR + dR[2]*f_dot_r;
      }
    } 
    k_out[t*3+0] += p[0]*SCAL_CONST;                                               
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
        
  }
}


template <class T>
void navier_dl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
  T mu = navier_mu();
  T nu = navier_nu();
  const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*(1.-nu));

  for(int t=0;t<trg_cnt;t++){
    T p[3]={0, 0, 0};
    for(int s=0;s<src_cnt;s++){
      
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      
      if (R!=0){                                                                   
        T invR2 = 1.0/R;
        T invR  = sqrt(invR2);
        T invR3 = invR2*invR;
        T invR5 = invR2*invR3;

        T* n = &v_src[6*s +3   ]; 
        T* f = &v_src[6*s ];

        T f_dot_r = (f[0]*dR[0] + f[1]*dR[1] + f[2]*dR[2]);
        T n_dot_f =  (f[0]*n[0] +  f[1]*n[1]  + f[2]*n[2]);
        T r_dot_n = (dR[0]*n[0] + dR[1]*n[1] + dR[2]*n[2]);


        p[0] +=  -(r_dot_n*f[0] + f_dot_r*n[0])*invR3*(1.0 - 2*nu)
                 + n_dot_f*dR[0]*invR3*(1.0 - 2*nu)
                 -3.*(r_dot_n*f_dot_r*dR[0])*invR5;

        p[1] +=  -(r_dot_n*f[1] + f_dot_r*n[1])*invR3*(1.0 - 2*nu)
                 + n_dot_f*dR[1]*invR3*(1.0 - 2*nu)
                 -3.*(r_dot_n*f_dot_r*dR[1])*invR5;

        p[2] +=  -(r_dot_n*f[2] + f_dot_r*n[2])*invR3*(1.0 - 2*nu)
                 + n_dot_f*dR[2]*invR3*(1.0 - 2*nu)
                 -3.*(r_dot_n*f_dot_r*dR[2])*invR5;
      }
    } 
    k_out[3*t+0] += p[0]*SCAL_CONST;
    k_out[3*t+1] += p[1]*SCAL_CONST;
    k_out[3*t+2] += p[2]*SCAL_CONST;

        
  }
}

const pvfmm::Kernel<real_t> ker_navier_sl_m2l=pvfmm::BuildKernel<real_t, navier_sl_m2l>("navier_sl_m2l", 3, std::pair<int,int>(4,3));
const pvfmm::Kernel<real_t> ker_navier_dl_m2l=pvfmm::BuildKernel<real_t, navier_dl_m2l>("navier_dl_m2l", 3, std::pair<int,int>(4,3));

template <class T>
void navier_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
    pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
    T mu = navier_mu();
    T nu = navier_nu();
    const T SCAL_CONST = 1.0/(16.0*pvfmm::const_pi<T>()*mu*(1.-nu));
    for(int t=0;t<trg_cnt;t++){
        T p[3]={0, 0, 0};
        for(int s=0;s<src_cnt;s++){
            T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                r_trg[3*t+1]-r_src[3*s+1],
                r_trg[3*t+2]-r_src[3*s+2]};
            T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

            if (R!=0){                                                                   
                T invR2=1.0/R;                                                             
                T invR=sqrt(invR2);                                                        
                T invR3=invR2*invR;                                                        
                T* f=&v_src[3*s];                                                          

                T f_dot_r=(f[0]*dR[0] +
                        f[1]*dR[1] +
                        f[2]*dR[2])* invR3;

                p[0] += (3. - 4.*nu)*f[0]*invR + dR[0]*f_dot_r;
                p[1] += (3. - 4.*nu)*f[1]*invR + dR[1]*f_dot_r;
                p[2] += (3. - 4.*nu)*f[2]*invR + dR[2]*f_dot_r;
            }
        } 
        k_out[t*3+0] += p[0]*SCAL_CONST;                                               
        k_out[t*3+1] += p[1]*SCAL_CONST;
        k_out[t*3+2] += p[2]*SCAL_CONST;

    }
}

template <class T>
void navier_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
    pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
    T mu = navier_mu();
    T nu = navier_nu();
    const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*(1.-nu));

    for(int t=0;t<trg_cnt;t++){
        T p[3]={0, 0, 0};
        for(int s=0;s<src_cnt;s++){
            T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                r_trg[3*t+1]-r_src[3*s+1],
                r_trg[3*t+2]-r_src[3*s+2]};

            T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

            //if (R!=0){ 
            if (R>=1e-12){ 
                T invR2 = 1.0/R;
                T invR  = sqrt(invR2);
                T invR3 = invR2*invR;
                T invR5 = invR2*invR3;

                T* n = &v_src[6*s];
                T* f = &v_src[6*s  +3]; 

                T f_dot_r = (f[0]*dR[0] + f[1]*dR[1] + f[2]*dR[2]);
                T n_dot_f =  (f[0]*n[0] +  f[1]*n[1]  + f[2]*n[2]);
                T r_dot_n = (dR[0]*n[0] + dR[1]*n[1] + dR[2]*n[2]);


                p[0] +=  - (1.0 - 2*nu) * (r_dot_n * f[0] + f_dot_r * n[0]) * invR3
                    + (1.0 - 2*nu) * n_dot_f * dR[0] * invR3
                    - 3. * (r_dot_n * f_dot_r * dR[0]) * invR5;

                p[1] +=  - (1.0 - 2*nu) * (r_dot_n * f[1] + f_dot_r * n[1]) * invR3
                    + (1.0 - 2*nu) * n_dot_f * dR[1] * invR3
                    - 3. * (r_dot_n * f_dot_r * dR[1]) * invR5;

                p[2] +=  - (1.0 - 2*nu) * (r_dot_n * f[2] + f_dot_r * n[2]) * invR3
                    + (1.0 - 2*nu) * n_dot_f * dR[2] * invR3
                    - 3. * (r_dot_n * f_dot_r * dR[2]) * invR5;
            }
        }  
        k_out[3*t+0] += p[0]*SCAL_CONST;
        k_out[3*t+1] += p[1]*SCAL_CONST;
        k_out[3*t+2] += p[2]*SCAL_CONST;
    }    
}



const pvfmm::Kernel<real_t> ker_navier_sl=pvfmm::BuildKernel<real_t, navier_sl>(
    "navier_sl"		/* name  */,
	3			                /* dim   */,
	std::pair<int,int>(3,3)	    /* k_dim */,
	NULL			            /* s2m   */,
	NULL			            /* s2l   */,
	NULL        	            /* s2t   */,
    NULL			            /* m2m   */,
    NULL			            /* m2l   */,
    NULL        	            /* m2t   */,
	NULL			            /* l2l   */,
	NULL			            /* l2t   */);


const pvfmm::Kernel<real_t> ker_navier_dl=pvfmm::BuildKernel<real_t, navier_sl, navier_dl>(
    "navier_sldl"		/* name  */,
	3               			/* dim   */,
	std::pair<int,int>(3,3)	    /* k_dim */,
	NULL		            	/* s2m   */,
	NULL			            /* s2l   */,
	NULL	                    /* s2t   */,
	&ker_navier_sl_m2l	/* m2m   */,
	&ker_navier_sl_m2l	/* m2l   */,
	&ker_navier_sl_m2l	/* m2t   */,
	NULL		            	/* l2l   */,
	NULL			            /* l2t   */);

