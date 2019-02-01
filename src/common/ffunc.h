#ifndef __FUNC_H__

extern double glb_xtarg, glb_ytarg, glb_ztarg;
extern int glb_kt;
extern double glb_lambda;
extern int glb_nk;
extern int glb_k;
extern int glb_sdof;
extern int glb_tdof;


void func(long *ndim, double *pnt, long *Xnumfun, double *funvls);


#endif
