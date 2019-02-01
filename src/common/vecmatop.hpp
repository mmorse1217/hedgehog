/*! \file */
#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

//#include "ebi_namespace.hpp"
#include "nummat.hpp"
#include "vec3t.hpp"

BEGIN_EBI_NAMESPACE

//--------------------------------------------------
// matrix operation, call lapack

//x = a x
int dscal(double alpha, DblNumVec& X);
int dscal(int n, double alpha, double* X);
//y = a x + y
int daxpy(double a, const DblNumVec& X, DblNumVec& Y);
int daxpy(double a, const DblNumMat& X, DblNumMat& Y);
int daxpy(int n, double a, double* X, double* Y);
//y = x
int dcopy(int n, double *X, double *Y);
int dcopy(const DblNumVec& X, DblNumVec& Y);
int dcopy(const DblNumMat& X, DblNumMat& Y);
// c = alpha*a*b + beta*c
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C);
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C);
// a = alpha* x*y' + a
int dger(double alpha, const DblNumVec& X, const DblNumVec& Y, DblNumMat& A);
int dger(int m, int n, double alpha, double* X, double* Y, double* A);
// y <= alpha A x + beta y
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y);
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y);
int dgemv(const DblNumMat& A, const DblNumVec& X, DblNumVec& Y);

// R <= tran(M)
int tran(const DblNumMat& M, DblNumMat& R);
// R <= pinv(M, epsilon)
int pinv(const DblNumMat& M, double epsilon, DblNumMat& R, const int type=0);

int pinv_mult(const DblNumMat& M, const DblNumVec& f, DblNumVec& c);

// R <= inv(M);
int inv(const DblNumMat& M, DblNumMat& R);
//--------------------------------------------------
// interpolation, etc.
//evaluation flags
enum {  EVFLAG_VL = 1,  EVFLAG_FD = 2,  EVFLAG_SD = 4 };
//domain flag
enum {  DMFLAG_PERIOD = 0,  DMFLAG_CLOSED = 1 };

// cubic spline interpolation
//int spev1d( int evflag, int dmflag, const DblNumMat& M, int n,   int i,   double u,   DblNumMat& res);
//int spev2d( int evflag, int dmflag, const DblNumMat& M, int* mn, int* ij, double* uv, DblNumMat& res);
int spev1d( int evflag, int dmflag, int dof, double* M, int n,   double e,   int i,   double u,   double* res);
int spev2d( int evflag, int dmflag, int dof, double* M, int* mn, double* ef, int* ij, double* uv, double* res);
int spcoef( int evflag, double u, double* us);
// cubic lagrangian interpolation
//int lagev1d(int evflag, int pdflag, const DblNumMat& M, int n,   int i,   double u,   DblNumMat& res);
//int lagev2d(int evflag, int pdflag, const DblNumMat& M, int* mn, int* ij, double* uv, DblNumMat& res);
int lagev1d(int evflag, int pdflag, int dof, double* M, int n,   double e,   int i,   double u,   double* res);
int lagev2d(int evflag, int pdflag, int dof, double* M, int* mn, double* ef, int* ij, double* uv, double* res, int LL);
int lagev3d(int evflag, int pdflag, int dof, double* M, int* mno, double* efg, int* ijk, double* uvw, double* res);

int lagcoef_old(int evflag, double u, double* us);
int lagcoef(int num, int evflag, double u, double* us);
// fft-based periodic refinement
//int fftrf1d(const DblNumMat& M, int  m,  int  ref, DblNumMat& R);
//int fftrf2d(const DblNumMat& M, int* mn, int* ref, DblNumMat& R);
int fftrf1d(int dof, double* M, int  m,  int  ref, double* res);
int fftrf2d(int dof, double* M, int* mn, int* ref, double* res);

// 3d shephard's interpolation from scattered points SPOS with values SDEN to target points TPOS with nbew values TDEN
int shep3d(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN);


//int dsgridinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN, const int KVAL, const double R, const Point3 ctr);
int dsgridinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN, const int KVAL, const double R, const Point3 ctr);
int dspntinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec &SDEN, const DblNumMat& TPOS, DblNumVec& TDEN);
int csapntinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec &SDEN, const DblNumMat& TPOS, const Index3 knots, const double smth, DblNumVec& TDEN);

int csgridtest();

int pou1d(int flags, double a, double b, double c, double* d, int pouctrl, int pdge);


/**
 * Creates nth degree spline basis functions in one dimension.
 * @param double    u       point where the function is evaluated
 * @param int       k       grid point that marks the end of the segment pt falls in
 * @param int       n       degree of spline 
 * @param int       m       number of gridpts 
 * @param double*   N       array to store the values at each grid pt
 */
inline int SplineBasis(double u, int k, int n, int m,double* N);
//int createN(double , int , int , int , double , double, double, double*);

/**
 * Partition of unity calculation using spline basis functions
 * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
 * @param int       deg     degree of spline basis functions
 * @param double    t       position in cd coordinates normalized by actual 
 *                          length of basis function (ub - lb)
 * @param double    scale
 * @param double*   res     Coordinates after POU
 * 
 * The expected dimension of res depends on flags (1 | 2 | 3)
 * (val, 1st deriv, 2nd deriv)
 */
int splineBlendFunc1D(int flags, int deg, double t, double scale, double* res);


/**
 * Partition of unity calculation using the optimal blending function constructed using the lower bound proof.
 * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
 * @param double    t       position in cd coordinates normalized by actual 
 *                          length of basis function (ub - lb)
 * @param double    deg     deg of blending function
 * @param double*   res     position|derivatives after POU
 * 
 * The expected dimension of res depends on flags (1 | 2 | 3)
 * (val, 1st deriv, 2nd deriv)
 */
int  optimBlendFunc1D(int flags, int deg, double t, double scale, double* res);

END_EBI_NAMESPACE

#endif
