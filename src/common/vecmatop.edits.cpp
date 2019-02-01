/*! \file */
#include "mobo_blas.h"
#include "mobo_lapack.h"
#include "svdrep.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

#include "vecmatop.hpp"
#include "ebi_uti.hpp"

extern "C" { 	 	
#include <stddef.h>
//#include "ngmath/include/ngmath.h"

}

BEGIN_EBI_NAMESPACE

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dscal"
ebiEC dscal(double alpha, DblNumVec& X)
{
  ebiFunctionBegin;
  dscal( X.m(), alpha, X.data() );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dscal"
ebiEC dscal(int n, double alpha, double* X)
{
  ebiFunctionBegin;
  PetscBLASInt incx = 1;
  PetscBLASInt nblas = n; 
  DSCAL(&nblas, &alpha, X, &incx);
  iC( PetscLogFlops( n ) );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "daxpy"
ebiEC daxpy(double a, const DblNumVec& X, DblNumVec& Y)
{
  ebiFunctionBegin;
  ebiAssert( X.m() == Y.m() );
  daxpy(X.m(), a, X.data(), Y.data());
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "daxpy"
ebiEC daxpy(double a, const DblNumMat& X, DblNumMat& Y)
{
  ebiFunctionBegin;
  if (!(X.m() == Y.m()) || !(X.n() == Y.n())) { cerr << X.m() << " " << Y.m() << " " << X.n() << " "<< Y.n() << endl; }
  ebiAssert( X.m() == Y.m() );  ebiAssert( X.n() == Y.n() );
  iC( daxpy(X.m()*X.n(), a, X.data(), Y.data()) );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "daxpy"
ebiEC daxpy(int n, double a, double* X, double* Y)
{
  ebiFunctionBegin;
  PetscBLASInt incx = 1; PetscBLASInt incy = 1;
  PetscBLASInt nblas = n; 
  DAXPY(&nblas, &a, X, &incx, Y, &incy);
  iC( PetscLogFlops( 2*n ) );
  ebiFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dcopy"
ebiEC dcopy(const DblNumVec& X, DblNumVec& Y)
{
  ebiFunctionBegin
  assert( X.m() == Y.m() );
  dcopy(X.m(), X.data(), Y.data());
  ebiFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "dcopy"
// ---------------------------------------------------------------------- 
ebiEC dcopy(const DblNumMat& X, DblNumMat& Y)
{
  ebiFunctionBegin;
  assert( X.m() == Y.m() );  assert( X.n() == Y.n() );
  iC( dcopy(X.m()*X.n(), X.data(), Y.data()) );
  ebiFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "dcopy"
// ---------------------------------------------------------------------- 
ebiEC dcopy(int n, double* X, double* Y)
{
  ebiFunctionBegin;
  PetscBLASInt incx = 1; PetscBLASInt incy = 1;
  PetscBLASInt nblas = n; 
  DCOPY(&nblas, X, &incx, Y, &incy);
  ebiFunctionReturn(0);
}

//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dgemm"
ebiEC dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  ebiFunctionBegin;
  ebiAssert( A.m() == C.m() );  ebiAssert( A.n() == B.m() );  ebiAssert( B.n() == C.n() );
  iC( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dgemm"
ebiEC dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  ebiFunctionBegin;
  if (m != 0 && n!= 0 && k!= 0){
	 char transa = 'N';
	 char transb = 'N';
	 ebiAssert(m!=0 && n!=0 && k!=0);
	 PetscBLASInt mblas = m; PetscBLASInt nblas = n; PetscBLASInt kblas = k;   
	 DGEMM(&transa, &transb, &mblas, &nblas, &kblas,
			 &alpha, A, &mblas, B, &kblas, &beta, C, &mblas);
	 iC( PetscLogFlops( 2*m*n*k ) );
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dger"
ebiEC dger(double alpha, const DblNumVec& X, const DblNumVec& Y, DblNumMat& A)
{
  ebiFunctionBegin;	
  ebiAssert(X.m() == A.m());
  ebiAssert(Y.m() == A.n());
  iC( dger(A.m(), A.n(), alpha, X.data(), Y.data(), A.data()) );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dger"
ebiEC dger(int m, int n, double alpha, double* X, double* Y, double* A)
{
  ebiFunctionBegin;
  if (m != 0 && n!= 0){
	 ebiAssert(m!=0 && n!=0);
	 PetscBLASInt incx = 1; PetscBLASInt incy = 1;
	 PetscBLASInt mblas = m; PetscBLASInt nblas = n;   
	 DGER(&mblas, &nblas, &alpha, X, &incx, Y, &incy, A, &mblas);
	 iC( PetscLogFlops( 2*m*n ) );
  }
  ebiFunctionReturn(0);
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dgemv"
ebiEC dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  ebiFunctionBegin;
  ebiAssert(Y.m() == A.m());
  ebiAssert(A.n() == X.m());
  iC( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dgemv"
ebiEC dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y)
{
  ebiFunctionBegin;
  if (m != 0 && n!= 0){
	 char trans = 'N';
	 ebiAssert(m!=0 && n!=0);
	 PetscBLASInt incx = 1; PetscBLASInt incy = 1;
	 PetscBLASInt mblas = m; PetscBLASInt nblas = n;  
	 DGEMV(&trans, &mblas, &nblas, &alpha, A, &mblas, X, &incx, &beta, Y, &incy);
	 iC( PetscLogFlops( 2*m*n ) );
  }
  ebiFunctionReturn(0);
}
//Assume alhpha == 1, beta == 1
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "dgemv"
ebiEC dgemv(const DblNumMat& A, const DblNumVec& X, DblNumVec& Y)
{
  ebiFunctionBegin;
  double one = 1.0;
  iC( dgemv(A.m(), A.n(), one, A.data(), X.data(), one, Y.data()));
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "tran"
ebiEC tran(const DblNumMat& M, DblNumMat& R)
{
  ebiFunctionBegin;
  ebiAssert(R.m()==M.n() && R.n()==M.m());  //R.resize(M.n(), M.m());
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		R(j,i) = M(i,j);
  iC( PetscLogFlops( M.n()*M.m() ) );
  ebiFunctionReturn(0);
}
// ----------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "pinv"
ebiEC pinv(const DblNumMat& M, double epsilon, DblNumMat& R, const int type)
{
  ebiFunctionBegin;
  ebiAssert(M.m() == R.n());  ebiAssert(M.n() == R.m());
  SVDRep svd;
  iC( svd.construct(epsilon, M, type) );
  //invert Svd
  double cutoff = 0.0;
  if (type == 0){
	 cutoff = svd.S()(0)*epsilon;
  }
  else {
	 ebiAssert(type == 1);
	 pinv(M,epsilon,R,0);
	 double mx = max((double)(M.m()),(double)(M.n()));
	 cutoff = (svd.S())(0)*mx * 2e-15;
	 cutoff = mx * 2e-15;
	 //std::cout << cutoff << endl;
	 //if (cutoff > 1e-15) { cutoff = 1e-15; }
	 //cutoff = 1e-10;
	 //std::cout << cutoff << endl;
  }
  for(int i=0; i<svd.S().m(); i++) {
    if( svd.S()(i) >= cutoff) {
      svd.S()(i) = 1.0/(svd.S()(i));
	 } else {
		//ebiAssert(0);
		svd.S()(i) = 0.0;
	 }
  }
  DblNumMat UT(svd.U().n(),  svd.U().m());
  DblNumMat V( svd.VT().n(), svd.VT().m());
  iC( tran(svd.U(), UT) );
  iC( tran(svd.VT(), V) );
  for(int i=0; i<V.m(); i++)
    for(int j=0; j<V.n(); j++) {
      V(i,j) = V(i,j) * svd.S()(j);
	 }
  iC( PetscLogFlops(V.m()*V.n()) );
  char transa = 'N';
  char transb = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  PetscBLASInt m = V.m();
  PetscBLASInt n = UT.n();
  PetscBLASInt k = V.n();
  DGEMM(&transa, &transb, &m, &n, &k, &alpha,
		  V.data(), &m, UT.data(), &k, 
		  &beta, R.data(), &m);  
  iC( PetscLogFlops( 2*m*n*k ) );
  ebiFunctionReturn(0);
}

// ----------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "pinv_mult"
ebiEC pinv_mult(const DblNumMat& M, const DblNumVec& f, DblNumVec& c)
{
  ebiFunctionBegin;
  //ebiAssert(M.m() == R.n());  ebiAssert(M.n() == R.m());
  SVDRep svd;
  iC( svd.construct(0.0, M, 1) );
  //invert Svd
  double cutoff = 0.0;
  double mx = max((double)(M.m()),(double)(M.n()));
  cutoff = (svd.S())(0)*mx * 2e-15;
  cutoff = mx * 2e-15;
  //std::cout << cutoff << endl;
  if (cutoff > 1e-15) { cutoff = 1e-15; }
  //cutoff = 1e-12;
  //std::cout << cutoff << endl;
  
  
  DblNumMat UT(svd.U().n(),  svd.U().m());
  DblNumMat V( svd.VT().n(), svd.VT().m());
  iC( tran(svd.U(), UT) );
  iC( tran(svd.VT(), V) );

  //cout << UT << " " << f << endl;
  //cout << c << " " << c.m() << endl;
  DblNumVec resid(c.m());
  iC( dgemv(1.0, UT, f, 0.0, resid));
  cout << resid << endl;
  for(int i=0; i<svd.S().m(); i++) {
    if( svd.S()(i) >= cutoff) {
      resid(i) *= 1.0/(svd.S()(i));
	 }
	 else {
		resid(i) = 0.0;
	 }
  }
  //cout << svd.S() <<  " " << resid << endl;
  iC( dgemv(1.0, V, resid, 0.0, c));
  //std::cout << c << std::endl; exit(0);
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "inv"
ebiEC inv(const DblNumMat& M, DblNumMat& R) //Gaussian Elimination
{
  ebiFunctionBegin;  //OR pinv(M, 0.0, R);
  ebiAssert(M.m()==M.n() && R.m()==R.n() && M.m()==R.m());
  memcpy(R.data(), M.data(), M.m()*M.n()*sizeof(double));
  PetscBLASInt info;
  PetscBLASInt m = M.m();
  PetscBLASInt* ipiv = new PetscBLASInt[m];
  DGETRF(&m, &m, R.data(), &m, ipiv, &info); ebiAssert(info==0);
  PetscBLASInt lwork = m;
  double* work = new double[lwork];
  DGETRI(&m, R.data(), &m, ipiv, work, &lwork, &info);  ebiAssert(info==0);
  delete [] ipiv;
  delete [] work;
  ebiFunctionReturn(0);
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "spev1d"
int spev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  ebiFunctionBegin;
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 ebiAssert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //ebiAssert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( spcoef(EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( spcoef(EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( spcoef(EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "spev2d"
int spev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  ebiFunctionBegin;
  
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  
  int is[4]; int js[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 ebiAssert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++)	is[k]=(i+k-1);
	 ebiAssert(j>=1 && j<=n-3);
	 for(int k=0; k<4; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  double scl;
  double us[4], vs[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 scl = double(m)/e;
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 scl = double(m*m)/(e*e);
	 iC( spcoef(EVFLAG_SD, u, us) );
	 iC( spcoef(EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( spcoef(EVFLAG_FD, u, us) );
	 iC( spcoef(EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( spcoef(EVFLAG_VL, u, us) );
	 iC( spcoef(EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "spcoef"
int spcoef(int evflag, double u, double* us)
{
  ebiFunctionBegin;
  double u1 = u;
  double u2 = u*u;
  double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {
	 us[0] = (  1 - 3*u1 + 3*u2 -   u3)/6.0;
	 us[1] = (  4        - 6*u2 + 3*u3)/6.0;
	 us[2] = (  1 + 3*u1 + 3*u2 - 3*u3)/6.0;
	 us[3] = (                  +   u3)/6.0;
  } else if(evflag==EVFLAG_FD) {
	 us[0] = (- 3 + 6*u1 - 3*u2)/6.0;
	 us[1] = (    -12*u1 + 9*u2)/6.0;
	 us[2] = (  3 + 6*u1 - 9*u2)/6.0;
	 us[3] = (             3*u2)/6.0;
  } else if(evflag==EVFLAG_SD) {
	 us[0] = (  6 - 6*u1 ) / 6.0;
	 us[1] = (-12 +18*u1 ) / 6.0;
	 us[2] = (  6 -18*u1 ) / 6.0;
	 us[3] = (      6*u1 ) / 6.0;	 //ebiAssert(0); //TODO;
  }  
  ebiFunctionReturn(0);
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "lagev1d"
int lagev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  ebiFunctionBegin;
  //exit(-1);
  int is[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
  } else {
	 ebiAssert(i>=1 && i<=m-3);
	 for(int k=0; k<4; k++) is[k]=(i+k-1);
  }
  DblNumMat M(dof,m,false,data); //ebiAssert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 double scl = 1.0;
	 double us[4]; iC( lagcoef(4, EVFLAG_VL, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++) 
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 double scl = double(m) / e;
	 double us[4]; iC( lagcoef(4, EVFLAG_FD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 double scl = double(m*m)/(e*e);
	 double us[4]; iC( lagcoef(4, EVFLAG_SD, u, us) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int d=0; d<dof; d++)
		for(int a=0; a<4; a++) {
		  res[d] += us[a] * M(d,is[a]);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
	 res+=dof; //cnt ++;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "lagev2d"
int lagev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res, int LL)
{
  ebiFunctionBegin;
  ebiAssert(LL >= 3);
  int m = mn[0];  int n = mn[1];
  double e = ef[0];  double f = ef[1];
  int i = ij[0];  int j = ij[1];
  double u = uv[0];  double v = uv[1];
  //int LL = 6;
  int is[LL]; int js[LL];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int k=0; k<LL; k++) is[k]=(i+k-1 + m) % m;
	 for(int k=0; k<LL; k++) js[k]=(j+k-1 + n) % n;
  } else {
	 ebiAssert(i>=1 && i<=m-3);
	 for(int k=0; k<LL; k++)	is[k]=(i+k-1);
	 ebiAssert(j>=1 && j<=n-3);
	 for(int k=0; k<LL; k++) js[k]=(j+k-1);
  }
  DblNumMat M(dof,m*n,false,data);
  
  double scl; 
  double us[LL], vs[LL];

  //  Replace all 4's with LL?
  
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(LL, EVFLAG_VL, u, us) );
	 iC( lagcoef(LL, EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<LL; a++)
		for(int b=0; b<LL; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_FD) {
	 //exit(-1);
	 scl = double(m)/e;
	 iC( lagcoef(4, EVFLAG_FD, u, us) );
	 iC( lagcoef(4, EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b];
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n)/f;
	 iC( lagcoef(4, EVFLAG_VL, u, us) );
	 iC( lagcoef(4, EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  //---------------------------
  if(evflag & EVFLAG_SD) {
	 //exit(-1);
	 scl = double(m*m)/(e*e);
	 iC( lagcoef(4, EVFLAG_SD, u, us) );
	 iC( lagcoef(4, EVFLAG_VL, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(m*n)/(e*f);
	 iC( lagcoef(4, EVFLAG_FD, u, us) );
	 iC( lagcoef(4, EVFLAG_FD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
	 //...
	 scl = double(n*n)/(f*f);
	 iC( lagcoef(4, EVFLAG_VL, u, us) );
	 iC( lagcoef(4, EVFLAG_SD, v, vs) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++) {
		  double coef = us[a]*vs[b]; 
		  for(int d=0; d<dof; d++)
			 res[d] += coef * M(d, is[a]+js[b]*m);
		}
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "lagev3d"
int lagev3d(int evflag, int dmflag, int dof, double* data, int* mno, double* efg, int* ijk, double* uvw, double* res)
{
  ebiFunctionBegin;
  ebiAssert(evflag==EVFLAG_VL);
  //exit(-1);
  int m = mno[0];  int n = mno[1];  int o = mno[2];
  //double e = efg[0];  double f = efg[1];  double g = efg[2];
  int i = ijk[0];  int j = ijk[1];  int k = ijk[2];
  double u = uvw[0];  double v = uvw[1];  double w = uvw[2];
  
  int is[4]; int js[4];  int ks[4];
  if(dmflag==DMFLAG_PERIOD) {
	 for(int h=0; h<4; h++) is[h]=(i+h-1 + m) % m;
	 for(int h=0; h<4; h++) js[h]=(j+h-1 + n) % n;
	 for(int h=0; h<4; h++) ks[h]=(k+h-1 + o) % o;
  } else {
	 ebiAssert(i>=1 && i<=m-3);
	 for(int h=0; h<4; h++)	is[h]=(i+h-1);
	 ebiAssert(j>=1 && j<=n-3);
	 for(int h=0; h<4; h++) js[h]=(j+h-1);
	 ebiAssert(k>=1 && k<=o-3);
	 for(int h=0; h<4; h++) ks[h]=(k+h-1);
  }
  DblNumMat M(dof,m*n*o,false,data);
  double scl; 
  double us[4], vs[4], ws[4];
  //---------------------------
  if(evflag & EVFLAG_VL) {
	 scl = 1.0;
	 iC( lagcoef(4, EVFLAG_VL, u, us) );
	 iC( lagcoef(4, EVFLAG_VL, v, vs) );
	 iC( lagcoef(4, EVFLAG_VL, w, ws) );
	 for(int d=0; d<dof; d++)		res[d] = 0;
	 for(int a=0; a<4; a++)
		for(int b=0; b<4; b++)
		  for(int c=0; c<4; c++) {
			 double coef = us[a]*vs[b]*ws[c]; 
			 for(int d=0; d<dof; d++)
				res[d] += coef * M(d, is[a]+js[b]*m+ks[c]*m*n);
		  }
	 for(int d=0; d<dof; d++)		res[d] *= scl;
	 res+=dof;
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "lagcoef"
int lagcoef(int num, int evflag, double u, double* us)
{
  ebiFunctionBegin;
  //double u1 = u;
  double u2 = u*u;
  //double u3 = u*u*u;
  if(       evflag==EVFLAG_VL) {

	 int n = num;
	 double p[8];// = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
	 p[0] = -1.0; p[1] = 0.0; p[2] = 1.0;
	 p[3] = 2.0; p[4] = 3.0; p[5] = 4.0;
	 p[6] = 5.0; p[7] = 6.0;

	 double p01 = -1.0; double p02 = -0.5; double p03 = -0.333333333; double p04 = -0.25; double p05 = -5.0; double p06 = -6.0;  double p07 = -7.0; 
	 
	 //for (int cc = 0; cc < n; cc++){
	 //p[cc] = -1.0 + cc;
	 //}
	 
	 for(int j=0; j<n; j++) {
		us[j] = 1.0;
		for(int i=0; i<n; i++) {
		  if(i!=j) us[j] *= (u-p[i])/(p[j]-p[i]);
		}
	 }
	 

	 /*
	 us[0] = (u3-3*u2+2*u)/(-6.0);
	 us[1] = (u3-2*u2-u+2)/(2.0);
	 us[2] = (u3-u2-2*u)  /(-2.0);
	 us[3] = (u3-u)       /(6.0);
	 */

	 /*
	 us[0] = (u-3)*(u3-3*u2+2*u)/(-6.0*(-1-3));
	 us[1] = (u-3)*(u3-2*u2-u+2)/(2.0*(0-3));
	 us[2] = (u-3)*(u3-u2-2*u)  /(-2.0*(1-3));
	 us[3] = (u-3)*(u3-u)       /(6.0*(2-3));
	 us[4] = (u-2)*(u-1)*(u-0)*(u+1)/((3-2)*(3-1)*(3-0)*(3+1));
	 */
	 
	 /*
	 us[0] = (u-3)*(u-2)*(u-1)*(u-0)*(u+1)*(u+2)/((-3-3)*(-3-2)*(-3-1)*(-3+0)*(-3+1)*(-3+2));
	 us[1] = (u-3)*(u-2)*(u-1)*(u-0)*(u+1)*(u+3)/((-2-3)*(-2-2)*(-2-1)*(-2+0)*(-2+1)*(-2+3));
	 us[2] = (u-3)*(u-2)*(u-1)*(u-0)*(u+2)*(u+3)/((-1-3)*(-1-2)*(-1-1)*(-1+0)*(-1+2)*(-1+3));
	 us[3] = (u-3)*(u-2)*(u-1)*(u+1)*(u+2)*(u+3)/((0-3)*(0-2)*(0-1)*(0+1)*(0+2)*(0+3));
	 us[4] = (u-3)*(u-2)*(u-0)*(u+1)*(u+2)*(u+3)/((1-3)*(1-2)*(1-0)*(1+1)*(1+2)*(1+3));
	 us[5] = (u-3)*(u-1)*(u-0)*(u+1)*(u+2)*(u+3)/((2-3)*(2-1)*(2-0)*(2+1)*(2+2)*(2+3));
	 us[6] = (u-2)*(u-1)*(u-0)*(u+1)*(u+2)*(u+3)/((3-2)*(3-1)*(3-0)*(3+1)*(3+2)*(3+3));
	 */

	 
  } else if(evflag==EVFLAG_FD) {
	 //cerr << "NOT CALLED!" << endl; exit(0);
	 us[0] = (3*u2-6*u+2)/(-6.0);
	 us[1] = (3*u2-4*u-1)/(2.0);
	 us[2] = (3*u2-2*u-2)/(-2.0);
	 us[3] = (3*u2-1)    /(6.0);
  } else if(evflag==EVFLAG_SD) {
	 ebiAssert(0); //TODO
	 exit(-1);
  }
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "fftrf1d"
int fftrf1d(int, double*, int, int, double*)//int fftrf1d(int dof, double* data, int n, int ref, double* res)
{
  ebiFunctionBegin;
  //TODO
  ebiFunctionReturn(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "fftrf2d"
int fftrf2d(int dof, double* data, int* mm, int* rr, double* res)
{
  ebiFunctionBegin;
  //DblNumMat M(dof,m   * n,  false,data);
  //DblNumMat R(dof,m*r * n*s,false,data);  //ebiAssert(M.m()==R.m() && M.n()==m*n && R.n()==m*r*n*s);  //int dof = M.m();
  //int m = mn[0];  int n = mn[1];  //int r = rs[0];  int s = rs[1];  //int m1 = m;  int m2 = n;  //int n1 = m*r;  int n2 = n*s;
  int m1 = mm[0];        int m2 = mm[1]; //old matrix size
  int n1 = mm[0]*rr[0];  int n2 = mm[1]*rr[1]; //new matrix size
  int nn[2]; nn[0] = n1; nn[1] = n2;
  int mq1 = m1/2+1;  int mq2 = m2/2+1;
  int nq1 = n1/2+1;  //int nq2 = n2/2+1;
  //allocate space
  double* _sorg = data;//M.data();
  double* _sref = res; //R.data();
  double* _forg = new double[dof * 2*mq1 * m2];  ebiAssert(_forg!=NULL);
  double* _fref = new double[dof * 2*nq1 * n2];  ebiAssert(_fref!=NULL);
  memset( _forg, 0, dof*2*mq1*m2*sizeof(double) );
  memset( _fref, 0, dof*2*nq1*n2*sizeof(double) );
  //scale  DblNumVec sorgvec(dof*m1*m2, false, _sorg);  iC( dscal(1.0/double(m1*m2), sorgvec) );
#ifdef FFTW3
  fftw_plan _forplan = fftw_plan_many_dft_r2c(2,mm,dof,_sorg,                NULL, dof, 1, (fftw_complex*)_forg, NULL, dof, 1, FFTW_ESTIMATE);
#else
  rfftwnd_plan regfftplan = rfftw2d_create_plan(m1,m2,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
  rfftwnd_plan refined_fftplan = rfftw2d_create_plan(n1,n2,FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE);
#endif
  //
  //do fft
#ifdef FFTW3
  fftw_execute(_forplan);
#else
  rfftwnd_real_to_complex(regfftplan, dof, _sorg, dof, 1, (fftw_complex*)_forg, dof, 1);
#endif

  //rearrange
  NumMat<fftw_complex> forgmat( dof*mq1, m2, false, (fftw_complex*)_forg );
  NumMat<fftw_complex> frefmat( dof*nq1, n2, false, (fftw_complex*)_fref );
  //    j direction
  for(int i=0; i<dof*mq1; i++) {
	 for(int j=0; j<mq2; j++) {
#ifdef FFTW3
		frefmat(i,j)[0] = forgmat(i,j)[0];
		frefmat(i,j)[1] = forgmat(i,j)[1];
#else
		frefmat(i,j) = forgmat(i,j);
#endif
	 }
	 for(int j=0; j<m2-mq2; j++) {
#ifdef FFTW3
		frefmat(i,j+mq2+n2-m2)[0] = forgmat(i,j+mq2)[0];
		frefmat(i,j+mq2+n2-m2)[1] = forgmat(i,j+mq2)[1];
#else
		frefmat(i,j+mq2+n2-m2) = forgmat(i,j+mq2);
#endif
	 }
	 if(m2 % 2 == 0) {
#ifdef FFTW3
		frefmat(i,mq2-1)[0] /= 2.0;
		frefmat(i,mq2-1)[1] /= 2.0;
		frefmat(i,mq2+n2-m2-1)[0] = frefmat(i,mq2-1)[0];
		frefmat(i,mq2+n2-m2-1)[1] = frefmat(i,mq2-1)[1];
#else
		frefmat(i,mq2-1).re /= 2.0;
		frefmat(i,mq2-1).im /= 2.0;
		frefmat(i,mq2+n2-m2-1).re = frefmat(i,mq2-1).re;
		frefmat(i,mq2+n2-m2-1).im = frefmat(i,mq2-1).im;
#endif
	 }
  }
  //    i direction
  if(m1 % 2 == 0) {
	 for(int i=(mq1-1)*dof; i<mq1*dof; i++) {
		for(int j=0; j<n2; j++) {
#ifdef FFTW3
		  frefmat(i,j)[0] /= 2.0;
		  frefmat(i,j)[1] /= 2.0;
#else
		  frefmat(i,j).re /= 2.0;
		  frefmat(i,j).im /= 2.0;
#endif
		}
	 }
  }
  //do ifft
#ifdef FFTW3
  fftw_plan _invplan = fftw_plan_many_dft_c2r(2,nn,dof,(fftw_complex*)_fref, NULL, dof, 1, _sref,                NULL, dof, 1, FFTW_ESTIMATE);
  fftw_execute(_invplan);
#else
  rfftwnd_complex_to_real(refined_fftplan, dof, (fftw_complex*)_fref, dof, 1, _sref, dof, 1);
#endif	
  //scale
  DblNumVec srefvec(dof*n1*n2, false, _sref);  iC( dscal(1.0/double(m1*m2), srefvec) );
  //deallocate space
  delete [] _forg;
  delete [] _fref;
#ifdef FFTW3
  fftw_destroy_plan(_forplan); _forplan = NULL; 
  fftw_destroy_plan(_invplan); _invplan=NULL;
#else  
  rfftwnd_destroy_plan(regfftplan);
  rfftwnd_destroy_plan(refined_fftplan);
#endif
  
  ebiFunctionReturn(0);
}


// ----------------------------------------------------------------------
/*
int shep3d(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN){
  ebiFunctionBegin;
  long int n = SPOS.n();
  ebiAssert( srcdof*n == SDEN.m() );
  long int t = TPOS.n();
  ebiAssert( srcdof*t == TDEN.m() );

  //ebiAssert( n >= 10 );
  //std::cout << SPOS << endl;
  // Set these parameters to standard values - see qshep3d.f 
  long int nq = 17; long int nw = 32; long int nr = (long int)(pow(((double)n)/3.0,(1.0/3.0))); long int ier;
  if (nr < 3) nr = 3;
  if (n < 33) { nw = n - 1; }
  if (n < 18) { nq = n - 1; }
  long int lout = 6; long int lcell[nr*nr*nr]; long int lnext[n];
  float xyzmin[3]; float xyzdel[3];
  float rsq[n]; float rmax;
  
  float x[n]; float y[n]; float z[n]; float f[n*srcdof];
  float a[9*n];

  //cout <<
  for (int k = 0; k < n; k++) {
	 x[k] = (float)SPOS(0,k); y[k] = (float)SPOS(1,k); z[k] = (float)SPOS(2,k);
  }

  for (int d = 0; d < srcdof; d++){
	 for (int k = 0; k < n; k++) {
		f[k] = (float)SDEN(k*srcdof + d);
	 }
	 
	 QSHEP3(&n, x, y, z, f, &nq, &nw, &nr, lcell, lnext, xyzmin, xyzdel, &rmax, rsq, a, &ier);
	 //cout << xyzmin[0] << " " << xyzmin[1] << " " << xyzmin[2] << " " << xyzdel[0] << " " << xyzdel[1] << " " << xyzdel[2] << endl; 
	 if ((int)ier != 0){
		//cout << SPOS << endl;
		//
		//cout << "ier = " << (int)ier << " " << ier << " " << (short)ier << endl;
		//cout << "nw = " << nw << " nr = " << nr << " nq =  " << nq << " n = " << n << endl;
		//for (int k = 0; k < n; k++) {
		//cout << x[k] << " " << y[k] << " " << z[k] << " " << f[k] << endl;
		//}
		//cout << "SPOS = " << SPOS << endl;
		//cout << "TPOS = " << TPOS << endl;
		//cout << "SDEN = " << SDEN << endl;
		//ebiAssert( (int)ier == 0);
	 }

	 if ((int)ier == 0){
		for (int i = 0; i < t; i++){
		  float px = TPOS(0,i); float py = TPOS(1,i); float pz = TPOS(2,i);
		  //cout << px << " " << py << " " << pz << endl;
		  TDEN(i*srcdof + d) =  (double)QS3VAL(&px, &py, &pz, &n, x, y, z, f, &nr, lcell, lnext, xyzmin, xyzdel, &rmax, rsq, a);
		}
		for (int i = 0; i < nr*nr*nr; i++) lcell[i] = 0;
		for (int i = 0; i < n; i++) f[i] = 0;
	 }
	 
  }
  
  ebiFunctionReturn(ier);
}
*/

#undef __FUNCT__
#define __FUNCT__ "csapntinterp"
int csapntinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec &SDEN, const DblNumMat& TPOS, const Index3 knots, const double ssmth, DblNumVec& TDEN){
  ebiFunctionBegin;
    assert(0); // WARNING commented out key function call below
  long int SRCNUM = SPOS.n();  ebiAssert( srcdof*SRCNUM == SDEN.m() );
  long int TRGNUM = TPOS.n();  ebiAssert( srcdof*TRGNUM == TDEN.m() );

  double wts[SRCNUM], xi[SRCNUM], yi[SRCNUM], zi[SRCNUM], ui[SRCNUM], xo[TRGNUM], yo[TRGNUM], zo[TRGNUM];
  int i,j,dknots[3],nderiv[3]={0,0,0},ier;
  dknots[0] = knots(0); dknots[1] = knots(1); dknots[2] = knots(2);

  for (i=0; i < SRCNUM; i++){
	 xi[i] = SPOS(0,i); yi[i] = SPOS(1,i); zi[i] = SPOS(2,i);
  }
  for (i=0; i < TRGNUM; i++){
	 xo[i] = TPOS(0,i); yo[i] = TPOS(1,i); zo[i] = TPOS(2,i);
  }

  //Not messing with weights right now
  wts[0] = -1.0;

  //Do each dof separately - target and source dof must be the same obviously
  for (j = 0; j < srcdof; j++){
	 //Copy in source values
	 for (i=0; i < SRCNUM; i++){ ui[i] = SDEN(i*srcdof + j); }
	 double *output;
	 //output = c_csa3lxd(SRCNUM,xi,yi,zi,ui,wts,dknots,ssmth,nderiv,TRGNUM,xo,yo,zo,&ier);
	 if (ier != 0) { cerr << "IER = " << ier << endl; for(i=0;i<TRGNUM;i++) output[i] = 0.0; }

	 int errval = 0;
	 //Copy to target values
	 for (i=0; i < TRGNUM; i++){
		if (std::isnan(output[i]) || ((output[i] < output[i]) && (output[i] > output[i]))) {
		  TDEN(i*srcdof + j) = 0.0;
		  errval = 100;
		}
		else {
		  TDEN(i*srcdof + j) = output[i];
		  errval = 0;
		}
	 }
	 if (errval != 0){ ier = errval; }
  }

  ebiFunctionReturn(ier);
}

#undef __FUNCT__
#define __FUNCT__ "dsgridinterp"
int dspntinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec &SDEN, const DblNumMat& TPOS, DblNumVec& TDEN){
  ebiFunctionBegin;
    assert(0);// WARNING commented out key function call below
  long int SRCNUM = SPOS.n();  ebiAssert( srcdof*SRCNUM == SDEN.m() );
  long int TRGNUM = TPOS.n();  ebiAssert( srcdof*TRGNUM == TDEN.m() );

  double xi[SRCNUM], yi[SRCNUM], zi[SRCNUM], ui[SRCNUM], xo[TRGNUM], yo[TRGNUM], zo[TRGNUM];
  int i,j,k,ier;

  for (i=0; i < SRCNUM; i++){
	 xi[i] = SPOS(0,i); yi[i] = SPOS(1,i); zi[i] = SPOS(2,i);
  }
  for (i=0; i < TRGNUM; i++){
	 xo[i] = TPOS(0,i); yo[i] = TPOS(1,i); zo[i] = TPOS(2,i);
  }

  //Do each dof separately - target and source dof must be the same obviously
  for (j = 0; j < srcdof; j++){
	 //Copy in source values
	 for (i=0; i < SRCNUM; i++){ ui[i] = SDEN(i*srcdof + j); }
	 double output[TRGNUM];
	 for (k = 0; k < TRGNUM; k++){
		//c_dspnt3d(SRCNUM,xi,yi,zi,ui,1,&xo[k],&yo[k],&zo[k],&output[k],&ier);
		//cerr << k << " " << output[k] << endl;
		if (ier != 0) { cerr << "IER = " << ier << endl; for(i=0;i<TRGNUM;i++) output[i] = 0.0; exit(-1); }

		int errval = 0;
		//Copy to target values
		if (std::isnan(output[i]) || ((output[i] < output[i]) && (output[i] > output[i]))) {
		  TDEN(k*srcdof + j) = 0.0;
		  errval = 100;
		}
		else {
		  TDEN(k*srcdof + j) = output[i];
		  errval = 0;
		}
		if (errval != 0){ ier = errval; }
	 }
  }
  ebiFunctionReturn(ier);
}

#undef __FUNCT__
#define __FUNCT__ "csgridtest"
int csgridtest(){
  ebiFunctionBegin;
    assert(0);// WARNING commented out key function call below
  int XMIN = -1.0, YMIN = -1.0, ZMIN = -1.0, XMAX = 1.0, YMAX = 1.0, ZMAX = 1.0;
  int NDATA = 1000;
  int NX = 3, NY = 3, NZ = 3;
  int NO = NX*NY*NZ;
  int N1 = 4, N2 = 4, N3 = 4;
  double ssmth = 0.5;

  DblNumMat SPOS(3,NDATA);
  int sdof = 1;
  DblNumVec SDEN(sdof*NDATA);
  DblNumMat TPOS(3,NO);
  DblNumVec TDEN(sdof*NO);

  int i,j,k;
  for (i = 0; i < NDATA; i++) {
    SPOS(0,i) = XMIN+(XMAX-(XMIN))*(((double) rand()/ (double) RAND_MAX));
    SPOS(1,i) = YMIN+(YMAX-(YMIN))*(((double) rand()/ (double) RAND_MAX));
	 SPOS(2,i) = ZMIN+(ZMAX-(ZMIN))*(((double) rand()/ (double) RAND_MAX));
	 double x = SPOS(0,i); double y = SPOS(1,i); double z = SPOS(2,i);
    SDEN(i) = x + y +z;
    double t1 = 1.0/(pow(fabs(x-0.1),2.75) + pow(fabs(y),2.75) + pow(fabs(z),2.75) + 0.09);
    double t2 = 1.0/(pow(fabs(x+0.1),2.75) + pow(fabs(y),2.75) + pow(fabs(z),2.75) + 0.09);
	 double t3 = 1.0/(pow(fabs(x+0.0),2.75) + pow(fabs(y),2.75) + pow(fabs(z),2.75) + 0.09);
	 SDEN(i) = exp(0.3*(SDEN(i)+t3-(t1-t2)));
  }

  int indx = 0;
  for (k = 0; k < NZ; k++) {
	 for (j = 0; j < NY; j++) {
		for (i = 0; i < NX; i++) {
		  TPOS(0,indx) = -2.0+((double)i/(double)(NX-1))*(XMAX-(XMIN));
		  TPOS(1,indx) = -2.0+((double)j/(double)(NY-1))*(YMAX-(YMIN));
		  TPOS(2,indx) = -2.0+((double)k/(double)(NZ-1))*(ZMAX-(ZMIN));
		  indx++;
		}
    }
  }

  Index3 knots(N1,N2,N3);

  //csapntinterp(sdof, SPOS, SDEN, TPOS, knots, ssmth, TDEN);

  indx = 0;
  for (k = 0; k < NZ; k++){
	 for (j = 0; j < NY; j++){
		for (i = 0; i < NX; i++){
		  double val = TPOS(0,indx) + TPOS(1,indx) + TPOS(2,indx);
		  double t1 = 1.0/(pow(fabs(TPOS(0,indx)-0.1),2.75) + pow(fabs(TPOS(1,indx)),2.75) + pow(fabs(TPOS(2,indx)),2.75) + 0.09);
		  double t2 = 1.0/(pow(fabs(TPOS(0,indx)+0.1),2.75) + pow(fabs(TPOS(1,indx)),2.75) + pow(fabs(TPOS(2,indx)),2.75) + 0.09);
		  double t3 = 1.0/(pow(fabs(TPOS(0,indx)+0.0),2.75) + pow(fabs(TPOS(1,indx)),2.75) + pow(fabs(TPOS(2,indx)),2.75) + 0.09);
		  val = exp(0.3*(val+t3-(t1-t2)));
		  printf("%d (%f,%f,%f) [%f, %f]\n", indx, TPOS(0,indx),TPOS(1,indx),TPOS(2,indx),val,TDEN(indx));
		  indx++;
		}
	 }
  }
  
  ebiFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dspntinterp"
int dsgridinterp(const int srcdof, const DblNumMat& SPOS, const DblNumVec& SDEN, const DblNumMat& TPOS, DblNumVec& TDEN, const int KVAL, const double R, const Point3 ctr){
  ebiFunctionBegin;
  assert(0);// WARNING commented out key function call below
  //ebiAssert(0);
  long int n = SPOS.n();
  ebiAssert( srcdof*n == SDEN.m() );
  long int t = TPOS.n();
  ebiAssert( srcdof*t == TDEN.m() );
  int NUM = n;
  
  int ier;
  int grdnum = KVAL*KVAL*KVAL;
  ebiAssert(TPOS.n() == grdnum);
  double h = 2.0/((double)(KVAL));
  double hh = h/2.0;
  double init = -1.0+hh;

  double xi[NUM], yi[NUM], zi[NUM], fi[NUM];
  double xo[KVAL], yo[KVAL], zo[KVAL];

  //cerr << R << " " << ctr << endl;
  for (int i = 1; i <= KVAL; i++) {
	 double tmp = R*(init + (i-1)*h);
    xo[i] = tmp + ctr(0);
	 yo[i] = tmp + ctr(1);
	 zo[i] = tmp + ctr(2);	 
	 //printf("%d %f %f %f\n", i, xo[i], yo[i], zo[i]);
  }

  for (int i = 0; i < NUM; i++) {
	 xi[i] = SPOS(0,i);
	 yi[i] = SPOS(1,i);
	 zi[i] = SPOS(2,i);
  }

  /* Interpolate for each degree of freedom separately */
  for (int d = 0; d < srcdof; d++){
	 for (int i = 0; i < NUM; i++) {
		fi[i*srcdof] = SDEN(i*srcdof + d);
	 }

	 //double *output;
	 //output = c_dsgrid3d(NUM, xi, yi, zi, fi, KVAL, KVAL, KVAL, xo, yo, zo, &ier);

	 double output[KVAL*KVAL*KVAL];
	 //c_dsgrid3d_mod(NUM, xi, yi, zi, fi, KVAL, KVAL, KVAL, xo, yo, zo, output, &ier);

	 if (ier != 0) {
		printf(" Error %d returned from c_dsgrid3s\n",ier);
		exit(1);
	 }
	 for (int i = 0; i < KVAL; i++){
		for (int j = 0; j < KVAL; j++){
		  for (int k = 0; k < KVAL; k++){
			 TDEN((KVAL*KVAL*i + KVAL*j + k)*srcdof + d) = output[KVAL*KVAL*i + KVAL*j + k];
		  }
		}
	 }
  }

  ebiFunctionReturn(ier);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "pou1d"
int pou1d(int flags, double u, double LB, double UB, double* res, int pouctrl, int pdeg)
{
  ebiFunctionBegin;
  enum {	 EVAL_VL=1,	 EVAL_FD=2,	 EVAL_SD=4  };
  ebiAssert( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  double* val = res;
  double* ext = res+1;
  double* fnl = res+2;
  double t = (u-LB)/(UB-LB);
  double ERR = 1e-7;
  
  if(       t<=0+ERR) {
	 //--------------------------------------------------
	 if(flags & EVAL_VL) val[0]=1;
	 if(flags & EVAL_FD) ext[0]=0;
	 if(flags & EVAL_SD) fnl[0]=0;
  } else if(t>=1-ERR) {
	 //--------------------------------------------------
	 if(flags & EVAL_VL) val[0]=0;
	 if(flags & EVAL_FD) ext[0]=0;
	 if(flags & EVAL_SD) fnl[0]=0;
  } else {
	 //--------------------------------------------------
	 if(       pouctrl==0 ) {
		//--------------------------------
		double s = 1-t;
		double t2 = t*t;		double t3 = t2*t;		double t4 = t2*t2;
		double s2 = s*s;		double s3 = s2*s;		double s4 = s2*s2;
		double a = 2*exp(-1/t)/(t-1);
		double b = 2*exp(-1/s)/(s-1);
		double da =  a*(1/(t*t) - 1/(t-1));
		double db = -b*(1/(s*s) - 1/(s-1));
		double dda = a*(-4*t3+7*t2-4*t+1+2*t4)/((t-1)*(t-1))/t4;
		double ddb = b*(-4*s3+7*s2-4*s+1+2*s4)/((s-1)*(s-1))/s4;
		double ea = exp(a);
		double eb = exp(b);
		double f = ea;
		double g = ea+eb;
		double df = ea*da;
		double dg = ea*da + eb*db;
		double ddf = ea*da*da + ea*dda;
		double ddg = ea*da*da + ea*dda + eb*db*db + eb*ddb;
		if(flags & EVAL_VL) {
		  val[0] = f/g;		  
		}
		if(flags & EVAL_FD) {
		  ext[0] = (df/g - f/(g*g)*dg) / (UB-LB); //ext[0] = (dea*eab-ea*deab) / (eab*eab) * 1.0/(UB-LB); //scaling
		}
		if(flags & EVAL_SD) {
		  fnl[0] = (ddf/g - 2*df/(g*g)*dg + 2*f/(g*g*g)*dg*dg - f/(g*g)*ddg) / ((UB-LB)*(UB-LB));
		}
	 } else if(pouctrl==1) {
		//--------------------------------
		double poucp[6] = {1,1,1,0,0,0}; //ctrl points for pou
		double x = t*3 + 1.0;
		int i = min(max(int(floor(x)),1),3);
		double f = x - i;		//double tmp[2]; //iC( spev1d(flags, _poucp, i, f, tmp) );//res[0] = tmp[0];res[1] = tmp[1] * 3 / (UB-LB);
		iC( spev1d(flags, DMFLAG_CLOSED, 1, poucp, 6, 2*(UB-LB), i, f, res) ); //cerr<<res[0]<<endl;
	 }
	 else if(pouctrl==2) {
		splineBlendFunc1D(flags, pdeg, t, 2*(UB-LB), res);
		res[1] = res[1]/(UB-LB);
    }
	 else if(pouctrl==3) {
		optimBlendFunc1D(flags, pdeg, 1-t, (UB-LB),res);
		res[1] = res[1]/(UB-LB);
    }
  }
  ebiFunctionReturn(0);
}


// ----------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "optimBlendFunc1D"
int optimBlendFunc1D(int flags, int deg, double z, double scale,double* res){
  ebiFunctionBegin;
  enum {	 EVAL_VL=1,	 EVAL_FD=2,	 EVAL_SD=4  };
  assert(deg<=10);
  //assert(z>=0 && z<=1); 
  double z2,z3,z4,z5,z6,z7,z8,z9,z10;
  z2 = z*z;  z3 = z2*z;  z4 = z3*z;  z5 = z4*z;  z6=z5*z;   z7= z6*z;    z8 = z7*z;  z9 = z8*z;  z10 = z9*z;
  
  
  switch(deg){
  case 2:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= .50000000000000000000) res[0] =2.*z2;
	 else if(z <= 1.) res[0] = -1.+(4.-2.*z)*z;
	 else res[0] = 1;
	 break;
  case 3:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= .2500000000) res[0] = 5.333333333*z3;
	 else if (z <= .7500000000) res[0] = .1666666667-2.*z+8.*z2-5.333333333*z3;
	 else if (z <= 1.) res[0] = -4.333333333+16.*z-16.*z2+5.333333333*z3;
	 else res[0] = 1;
	 break;
  case 4:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= .14644660940672623780) res[0] = 16.*z4;
	 else if (z <= .50000000000000000000) res[0] = -0.14718625761429707190e-1+.40202025355333863355*z-4.1177490060914376575*z2+18.745166004060958438*z3-16.*z4;
	 else if (z <= .85355339059327376220) res[0] = 1.9852813742385702924-15.597979746446661366*z+43.882250993908562342*z2-45.254833995939041560*z3+16.*z4;
	 else if ( z <= 1.) res[0] = -14.999999999999999999+64.*z-96.*z2+64.*z3-16.*z4;
	 else res[0] = 1.0000000000000000010+0.53573911150189318880e-19*z-0.26786955575094659440e-19*z2;
	 break;
  case 5:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= 0.95491502812526287949e-1) res[0] = 51.200000000000000000*z5;
	 else if (z <= .34549150281252628795) res[0] = 0.81306187557833487499e-3-0.42572472504416375409e-1*z+.89164944001345942980*z2-9.3374741600201891447*z3+48.891649440013459430*z4-51.200000000000000000*z5;
	 else if (z <= .65450849718747371205) res[0] = -.50325224750231333941+7.2523292150113563936*z-41.337474160020189144*z2+112.89164944001345943*z3-128.*z4+51.200000000000000000*z5;
	 else if (z <= .90450849718747371205) res[0] = 11.795934690622108325-86.705098312484227236*z+245.77087639996635144*z2-325.77087639996635143*z3+207.10835055998654057*z4-51.200000000000000000*z5;
	 else if ( z <= 1.) res[0] = -50.199999999999999992+255.99999999999999999*z-512.*z2+512.*z3-256.*z4+51.200000000000000000*z5;
	 else res[0] = 1.0000000000000000005+0.50335829998180960801e-17*z-0.90945424859373685600e-17*z2+0.20480000000000000000e-17*z3-0.51200000000000000000e-18*z4;
	 break;
  case 6:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= 0.66987298107780676618e-1) res[0] = 170.66666666666666667*z6;
	 else if ( z <= .25000000000000000000) res[0] = -0.30841356309254049329e-4+0.27624362092913055274e-2*z-.10309552285743124926*z2+2.0520412231296636895*z3-22.974966311837064286*z4+137.18998652473482572*z5-170.66666666666666667*z6;
	 else if (z <= .50000000000000000000) res[0] = 0.83302491977024079238e-1-1.9972375637907086944*z+19.896904477142568750*z2-104.61462544353700297*z3+297.02503368816293572*z4-374.81001347526517428*z5+170.66666666666666667*z6;
	 else if ( z <= .75000000000000000000) res[0] = -5.2500308413563092643+62.002762436209291334*z-300.10309552285743126*z2+748.71870788979633037*z3-982.97496631183706432*z4+649.18998652473482570*z5-170.66666666666666667*z6;
	 else if (z <= .93301270189221932338) res[0] = 55.499969158643690793-423.99723756379070876*z+1319.8969044771425687*z2-2131.2812921102036695*z3+1897.0250336881629356*z4-886.81001347526517424*z5+170.66666666666666667*z6;
	 else if (z <= 1.) res[0] = -169.66666666666666681+1023.9999999999999999*z-2559.9999999999999998*z2+3413.3333333333333333*z3-2560.*z4+1024.*z5-170.66666666666666667*z6;
	 else res[0]= 1.0000000000000001158-0.16720452122920761250e-15*z+0.18244964361477047929e-15*z2-0.17930021879695548175e-15*z3+0.62756570302591872875e-16*z4-0.40960000000000000000e-17*z5;
	 break;
  case 7:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= 0.49515566048790436882e-1) res[0] = 585.14285714285714286*z7;
	 else if (z <= .18825509907063323474) res[0] =0.85405165909716607756e-6-0.12073701445297728894e-3*z+0.73150944695255066160e-2*z2-.24622204871620700105*z3+4.9726190845438532088*z4-60.255222525103182580*z5+405.63151707169125893*z6-585.14285714285714286*z7;
	 else if (z <= .38873953302184279786) res[0] = -0.98057613154314852871e-2+.36452440624340911689*z-5.8036054870927415800*z2+51.199224495620706333*z3-268.30256335837421735*z4+810.71778312098134436*z5-1136.5542545149362001*z6+585.14285714285714286*z7; 
	 else if ( z <= .61126046697815720214) res[0] = 1.5602058488570656905-27.906542304612471398*z+212.37127236306429860*z2-884.19521195470236570*z3+2137.9215455689892500*z4-2903.1686182275957000*z5+2048.*z6-585.14285714285714286*z7;
	 else if ( z <= .81174490092936676526) res[0] = -35.754160043984012320+399.40812200071269746*z-1884.8427013352945168*z2+4834.0817119732132800*z3-7216.9725345224895032*z4+6279.3922560313641436*z5-2959.4457454850638000*z6+585.14285714285714286*z7;
	 else if ( z <= .95048443395120956312) res[0] = 236.03297034893468923-1944.3206905514169708*z+6776.9751057200789195*z2-12950.277629527747442*z3+14691.830737465603176*z4-9914.4661200949556290*z5+3690.3684829283087410*z6-585.14285714285714286*z7;
	 else if(z <= 1.) res[0] = -584.14285714285714216+4095.9999999999999974*z-12287.999999999999996*z2+20479.999999999999998*z3-20480.*z4+12288.*z5-4096.*z6+585.14285714285714286*z7;
	 else res[0]= 1.0000000000000000484-0.36175677281249441305e-15*z+0.63749634654633829355e-15*z2-0.26517790529631145793e-16*z3-0.42372151766370038882e-15*z4+0.36623100231065481024e-16*z5+0.16384000000000000000e-16*z6;
	 break;
  case 8:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= 0.38060233744356621936e-1) res[0] = 2048.*z8;
	 else if (z <= .14644660940672623780) res[0] = -0.18035639965439994271e-7+0.37909677773566962999e-5*z-0.34861549484613467622e-3*z2+0.18319146287314938525e-1*z3-.60164982204133232780*z4+12.646266464520426781*z5-166.13490276311756640*z6+1247.1577393350777876*z7-2048.*z8;
	 else if (z <= .30865828381745511414) res[0] = 0.86653374158012347544e-3-0.47333694316756025668e-1*z+1.1309933743747886546*z2-15.432253515692972281*z3+131.27723700036884408*z4-707.77402790306240450*z5+2293.5338368755011943*z6-3551.6047577045275726*z7+2048.*z8;
	 else if ( z <= .50000000000000000000) res[0] = -.33656327434686270681+8.6983848113431433435*z-98.040214315981076030*z2+627.16321785035092847*z3-2471.0969137047681297*z4+6037.2233679723528780*z5-8632.7846024904456260*z6+6562.5098864258416076*z7-2048.*z8;
	 else if (z <= .69134171618254488586) res[0] = 15.663436725653137396-247.30161518865685704*z+1693.9597856840189244*z2-6540.8367821496490713*z3+15448.903086295231870*z4-22634.776632027647122*z5+20039.215397509554373*z6-9821.4901135741583924*z7+2048.*z8;
	 else if(z <= .85355339059327376220) res[0] = -198.08456096638670317+2226.1264183148314156*z-10828.065017136323776*z2+29684.446633306717202*z3-50049.248130959109695*z4+53157.729081554865736*z5-34776.300532943808187*z6+12832.395242295472428*z7-2048.*z8;
	 else if(z <= .96193976625564337806) res[0]= 955.91457248183607910-8589.8262441998840622*z+33522.803640873806616*z2-74236.102794031302533*z3+102138.87298221848013*z4-89481.850624077551436*z5+48780.030727417573053*z6-15136.842260664922211*z7+2048.*z8;
	 else if(z <= 1.) res[0] = -2047.0000000000000043+16384.000000000000064*z-57344.000000000000130*z2+114688.00000000000009*z3-143360.00000000000002*z4+114688.*z5-57344.*z6+16384.*z7-2048.*z8;
	 else res[0] = 1.0000000000000049597-0.54820642047180343368e-14*z+0.39224496026381535292e-14*z2-0.11250317116645488317e-13*z3+0.10613072687253747715e-13*z4-0.51840487459490730284e-14*z5+0.13760453616361108574e-14*z6-0.13107200000000000000e-15*z7;
	 break;
  case 9:
	 if      (z<= 0) res[0] = 0.0;
	 else if (z <= 0.30153689607045807973e-1) res[0] = 7281.7777777777777778*z9;
	 else if (z <= .11697777844051098240) res[0] = 0.30014530688699986480e-9-0.89584651072079832882e-7*z+0.11883739899099736676e-4*z2-0.91957988533361885860e-3*z3+0.45744645049278490318e-1*z4-1.5170496760233828764*z5+33.540388054068717260*z6-476.70624094417708557*z7+3952.3044041747081425*z8-7281.7777777777777778*z9;
	 else if (z <= .25000000000000000000) res[0] = -0.59730418570177560056e-4+0.45954539075066661977e-2*z-.15713056001046843486*z2+3.1335706367386609790*z3-40.147661210465904882*z4+342.08159265890118616*z5-1924.6585503289824747*z6+6697.5462336543921967*z7-11380.206971579947342*z8+7281.7777777777777778*z9;
	 else if ( z <= .41317591116653482557) res[0] = 0.55495825136985377258e-1-1.9954045460924933281*z+31.842869439989531566*z2-295.53309602992800583*z3+1751.8523387895340956*z4-6825.9184073410988144*z5+17190.008116337684193*z6-26070.453766345607804*z7+21387.793028420052658*z8-7281.7777777777777778*z9;
	 else if (z <= .58682408883346517443) res[0]= -5.0552273567844730368+109.32887016157253927*z-1045.8992902455457948*z2+5790.8132711456916207*z3-20344.110431614285835*z4+46652.423802928302734*z5-69098.246296949255063*z6+63433.022751509310967*z7-32768.*z8+7281.7777777777777778*z9;
	 else if (z <= .75000000000000000000) res[0] = 115.12660322810743222-1733.8750743359833380*z+11518.029011754604590*z2-44165.951051435753480*z3+107351.98778554328413*z4-170952.98921104980870*z5+178114.29678565342937*z6-117112.10953898518655*z7+44148.206971579947345*z8-7281.7777777777777778*z9;
	 else if(z <= .88302222155948901760) res[0] = -978.37339677189257666+11388.124925664016721*z-58465.970988245395545*z2+173562.04894856424670*z3-328104.01221445671602*z4+409655.01078895019136*z5-337981.70321434657063*z6+177799.89046101481346*z7-54155.793028420052661*z8+7281.7777777777777778*z9;
	 else if(z <= .96984631039295419203) res[0] = 3775.1114393095821432-37060.671129748843081*z+161002.10069787071451*z2-406369.18580927442303*z3+657031.84382374066480*z4-705986.05914739567384*z5+504308.21331499667623*z6-231002.27100754651194*z7+61583.695595825291858*z8-7281.7777777777777778*z9;
	 else if(z <= 1.) res[0] = -7280.7777777777780210+65536.000000000000220*z-262143.99999999999959*z2+611669.33333333333273*z3-917503.99999999999968*z4+917503.99999999999990*z5-611669.33333333333333*z6+262144.*z7-65536.*z8+7281.7777777777777778*z9;
	 else res[0] = .99999999999996338747-0.38084892777430476847e-13*z+0.38100254603956599662e-12*z2-0.53798512751231230287e-12*z3+0.29695038035979669252e-12*z4-0.45208610065940934494e-13*z5-0.86211581928152639460e-14*z6-0.34984443524634146274e-16*z7+0.39321600000000000000e-15*z8;
	 break;
  case 10:
	 if      (z<= 0) res[0] = 0.0;
	 else if ( z <= 0.24471741852423213942e-1) res[0] = 26214.400000000000000*z10;
	 else if (z <= 0.95491502812526287949e-1) res[0] = -0.40384875164793945666e-11+0.16502656577670215063e-8*z-0.30346002768152967394e-6*z2+0.33067803075241608518e-4*z3-0.23647133796461901013e-2*z4+.11595643958194340222*z5-3.9486509338409739117*z6+92.203161336119367790*z7-1412.9025105591736969*z8+12830.240592323261991*z9-26214.400000000000000*z10;
	 else if(z <= .20610737385376343542) res[0] = 0.33053440291072819838e-5-0.34613889315882092304e-3*z+0.16311435276225912962e-1*z2-.45548358741953987547*z3+8.3455418873216695582*z4-104.78853899503713546*z5+911.52979780112876273*z6-5386.0909981065753686*z7+20100.637954127342092*z8-37234.808434250520462*z9+26214.400000000000000*z10;
	 else if(z <= .34549150281252628795) res[0] =-0.72494865223342578718e-2+.35154771772884343771*z-7.6666856608031469320*z2+98.948976685545888113*z3-835.66989640124444868*z4+4809.2448942041588794*z5-18956.890164236811193*z6+49698.705936454644500*z7-80122.844904456821299*z8+70824.814388791403557*z9-26214.400000000000000*z10;
	 else if(z <= .50000000000000000000) res[0] = 1.2631596940688140388-36.419521042260577435*z+471.27363001388107786*z2-3597.7386834608254167*z3+17888.968343881403582*z4-60227.277582968936504*z5+137912.67572410236643*z6-209757.02946215209234*z7+201493.05587002070501*z8-110312.23463778237889*z9+26214.400000000000000*z10;
	 else if(z <= .65450849718747371205) res[0]= -49.936840305931193005+987.58047895773946460*z-8744.7263699861190370*z2+45554.261316539174767*z3-154143.03165611859655*z4+352649.52241703106350*z5-550215.32427589763360*z6+576674.97053784790771*z7-388330.94412997929500*z8+151831.76536221762111*z9-26214.400000000000000*z10;
	 else if(z <= .79389262614623656458) res[0] = 706.41315638872070864-10568.419174901716989*z+70707.257318275144570*z2-278157.28316680560265*z3+711384.62043728070230*z4-1234240.7730875343176*z5+1470245.1972753673968*z6-1187318.7353027093976*z7+622347.51540533418924*z8-191319.18561120859644*z9+26214.400000000000000*z10; 
	 else if(z <= .90450849718747371205) res[0] = -4507.7858074779671021+55110.474823368299615*z-301578.63281137226899*z2+972342.26927039137513*z3-2045126.9589535473005*z4+2932335.1420031392324*z5-2903326.8470495769608*z6+1960693.9090018934247*z7-864635.36204587265791*z8+224909.19156574947953*z9-26214.400000000000000*z10;
	 else if(z <= .97552825814757678606) res[0]= 14709.693783341442707-157352.75410787171601*z+755442.42687295305590*z2-2143961.0430694622613*z3+3984148.4028061700548*z4-5066628.3355196876724*z5+4466203.5870620838613*z6-2695050.3555994998386*z7+1065588.7371796498158*z8-249313.75940767673801*z9+26214.400000000000000*z10;
	 else if(z <= 1.) res[0]=-26213.399999999999546+262143.99999999999818*z-1179647.9999999999949*z2+3145727.9999999999930*z3-5505023.9999999999962*z4+6606028.7999999999996*z5-5505024.*z6+3145728.*z7-1179648.*z8+262144.*z9-26214.400000000000000*z10;
	 else res[0]= 1.0000000000004344917-0.23771387168790120845e-11*z+0.60627974858379033630e-11*z2-0.70198079820713437333e-11*z3+0.37343903197696794288e-11*z4-0.77687433196916569664e-12*z5+0.45868652506776221160e-13*z6-0.21801870076643113956e-13*z7+0.54625363142958494196e-16*z8+0.15728640000000000000e-14*z9;
	 break;
	 
  }
  
  if(flags & EVAL_FD) {
	 
	 switch(deg){
		
	 case 2:
		if(z<=0) res[1] = 0;
		else if(z <= .50000000000000000000) res[1] = 4.*z;
		else if(z <= 1.) res[1] = -4.*z+4.;
		else res[1] = 0;
		break;
	 case 3:
		if(z<=0) res[1] = 0;
		else if(z < .25000000000000000000) res[1] =16.*z2;
		else if(z < .75000000000000000000) res[1] =-16.*z2+16.*z-2.; 
		else if(z < 1.) res[1]=16.*z2-32.*z+16.;
		else res[1] = 0;
		break;
	 case 4:
		if(z<=0) res[1] = 0;
		else if (z < .14644660940672623781) res[1] = 64.*z3;
		else if (z < .50000000000000000000) res[1] = -64.*z3+56.235498012182875317*z2-8.2354980121828753150*z+.40202025355333863352;
		else if(z < .85355339059327376225) res[1] = 64.*z3-135.76450198781712468*z2+87.764501987817124684*z-15.597979746446661366;
		else if(z < 1.) res[1] = -64.*z3+192.*z2-192.*z+64;
		else res[1] = -0.53573911150189318881e-19*z+0.53573911150189318881e-19;
		break;
	 case 5:
		if(z<=0) res[1] = 0;
		else if(z < 0.95491502812526287948e-1) res[1] = 256.*z4;
		else if(z < .34549150281252628792) res[1] = -256.*z4+195.56659776005383771*z3-28.012422480060567434*z2+1.7832988800269188596*z-0.42572472504416375407e-1;
		else if(z < .65450849718747371208) res[1] = 256.*z4-512.*z3+338.67494832004037829*z2-82.674948320040378282*z+7.2523292150113563936;
		else if(z < .90450849718747371204) res[1] = -256.*z4+828.43340223994616228*z3-977.31262919989905429*z2+491.54175279993270286*z-86.705098312484227238;
		else if(z < 1.) res[1] = 256.*z4-1024.*z3+1536.*z2-1024.*z+256.;
		else res[1] = -0.20480000000000000000e-17*z3+0.61440000000000000000e-17*z2-0.18189084971874737120e-16*z+0.50335829998180960800e-17;
		break;
	 case 6:
		if(z<=0) res[1] = 0;
		else if(z < 0.66987298107780676620e-1) res[1] = 1024.*z5;
		else if(z < .25000000000000000000) res[1] = -1024.*z5+685.94993262367412860*z4-91.899865247348257152*z3+6.1561236693889910688*z2-.20619104571486249852*z+0.27624362092913055274e-2;
		else if(z < .50000000000000000000) res[1] = 1024.*z5-1874.0500673763258715*z4+1188.1001347526517428*z3-313.84387633061100891*z2+39.793808954285137496*z-1.9972375637907086945;
		else if(z < .75000000000000000000) res[1] = -1024.*z5+3245.9499326236741285*z4-3931.8998652473482572*z3+2246.1561236693889911*z2-600.20619104571486252*z+62.002762436209291333;
		else if(z < .93301270189221932341) res[1] = 1024.*z5-4434.0500673763258708*z4+7588.1001347526517424*z3-6393.8438763306110085*z2+2639.7938089542851374*z-423.99723756379070875;
		else if(z < 1.) res[1]=-1024.*z5+5120.*z4-10240.*z3+10240.*z2-5120.*z+1024.;
		else res[1] = -0.20480000000000000000e-16*z4+0.25102628121036749150e-15*z3-0.53790065639086644528e-15*z2+0.36489928722954095856e-15*z-0.16720452122920761250e-15;
		break;
	 case 7:
		if(z<=0) res[1] = 0;
		else if(z < 0.49515566048790436879e-1) res[1] = 4096.*z6;
		else if(z < .18825509907063323475) res[1] = -4096.*z6+2433.7891024301475535*z5-301.27611262551591290*z4+19.890476338175412836*z3-.73866614614862100315*z2+0.14630188939051013232e-1*z-0.12073701445297728894e-3;
		else if(z < .38873953302184279784) res[1] = 4096.*z6-6819.3255270896172000*z5+4053.5889156049067218*z4-1073.2102534334968694*z3+153.59767348686211900*z2-11.607210974185483160*z+.36452440624340911689;
		else if(z < .61126046697815720216) res[1] = -4096.*z6+12288.*z5-14515.843091137978500*z4+8551.6861822759570000*z3-2652.5856358641070972*z2+424.74254472612859720*z-27.906542304612471398;
		else if(z < .81174490092936676525) res[1] = 4096.*z6-17756.674472910382800*z5+31396.961280156820718*z4-28867.890138089958013*z3+14502.245135919639841*z2-3769.6854026705890334*z+399.40812200071269745; 
		else if(z < .95048443395120956303) res[1] = -4096.*z6+22142.210897569852444*z5-49572.330600474778142*z4+58767.322949862412704*z3-38850.832888583242326*z2+13553.950211440157839*z-1944.3206905514169707;
		else if(z < 1.) res[1] = 4096.*z6-24576.*z5+61440.*z4-81920.*z3+61440.*z2-24575.999999999999992*z+4095.9999999999999974;
		else res[1] = 0.98304000000000000000e-16*z5+0.18311550115532740512e-15*z4-0.16948860706548015553e-14*z3-0.79553371588893437385e-16*z2+0.12749926930926765871e-14*z-0.36175677281249441306e-15;
		break;
	 case 8:
		if(z<=0) res[1] = 0;
		else if(z < 0.38060233744356621938e-1) res[1] = 16384.*z7;
		else if(z < .14644660940672623781) res[1] = -16384.*z7+8730.1041753455445132*z6-996.80941657870539846*z5+63.231332322602133910*z4-2.4065992881653293110*z3+0.54957438861944815575e-1*z2-0.69723098969226935242e-3*z+0.37909677773566962996e-5; 
		else if(z < .30865828381745511415) res[1] = 16384.*z7-24861.233303931693010*z6+13761.203021253007166*z5-3538.8701395153120225*z4+525.10894800147537628*z3-46.296760547078916843*z2+2.2619867487495773094*z-0.47333694316756025672e-1;
		else if(z < .50000000000000000000) res[1] = -16384.*z7+45937.569204980891252*z6-51796.707614942673757*z5+30186.116839861764392*z4-9884.3876548190725192*z3+1881.4896535510527853*z2-196.08042863196215206*z+8.6983848113431433435;
		else if(z < .69134171618254488585) res[1] = 16384.*z7-68750.430795019108748*z6+120235.29238505732623*z5-113173.88316013823562*z4+61795.612345180927480*z3-19622.510346448947213*z2+3387.9195713680378490*z-247.30161518865685704;
		else if(z < .85355339059327376225) res[1] = -16384.*z7+89826.766696068306996*z6-208657.80319766284911*z5+265788.64540777432867*z4-200196.99252383643877*z3+89053.339899920151606*z2-21656.130034272647554*z+2226.1264183148314156;
		else if(z < .96193976625564337808) res[1] = 16384.*z7-105957.89582465445548*z6+292680.18436450543832*z5-447409.25312038775716*z4+408555.49192887392052*z3-222708.30838209390760*z2+67045.607281747613232*z-8589.8262441998840621;
		else if(z < 1.) res[1] = -16384.*z7+114688.*z6-344064.*z5+573440.*z4-573440.00000000000008*z3+344064.00000000000027*z2-114688.00000000000026*z+16384.000000000000064;
		else res[1] = -0.91750400000000000000e-15*z6+0.82562721698166651444e-14*z5-0.25920243729745365142e-13*z4+0.42452290749014990856e-13*z3-0.33750951349936464951e-13*z2+0.78448992052763070588e-14*z-0.54820642047180343371e-14;
		break;
	 case 9:
		if(z<=0) res[1] = 0;
		else if(z < 0.30153689607045807973e-1) res[1] = 65536.*z8;
		else if(z < .11697777844051098241) res[1] = -65536.000000000000001*z8+31618.435233397665141*z7-3336.9436866092395990*z6+201.24232832441230355*z5-7.5852483801169143815*z4+.18297858019711396128*z3-0.27587396560008565757e-2*z2+0.23767479798199473350e-4*z-0.89584651072079832878e-7;
		else if(z < .25000000000000000000)res[1] = 65536.000000000000001*z8-91041.655772639578736*z7+46882.823635580745377*z6-11547.951301973894849*z5+1710.4079632945059308*z4-160.59064484186361953*z3+9.4007119102159829364*z2-.31426112002093686972*z+0.45954539075066661973e-2;
		else if(z < .41317591116653482558) res[1] = -65536.000000000000001*z8+171102.34422736042125*z7-182493.17636441925463*z6+103140.04869802610516*z5-34129.592036705494071*z4+7007.4093551581363820*z3-886.59928808978401749*z2+63.685738879979063132*z-1.9954045460924933281;
		else if(z < .58682408883346517442) res[1] = 65536.000000000000001*z8-262144.*z7+444031.15926056517677*z6-414589.47778169553038*z5+233262.11901464151367*z4-81376.441726457143344*z3+17372.439813437074863*z2-2091.7985804910915896*z+109.32887016157253927;
		else if(z < .75000000000000000000) res[1] = -65536.000000000000001*z8+353185.65577263957875*z7-819784.76677289630585*z6+1068685.7807139205762*z5-854764.94605524904355*z4+429407.95114217313656*z3-132497.85315430726044*z2+23036.058023509209180*z-1733.8750743359833380;
		else if(z < .88302222155948901759) res[1] = 65536.000000000000001*z8-433246.34422736042125*z7+1244599.2332271036942*z6-2027890.2192860794238*z5+2048275.0539447509570*z4-1312416.0488578268641*z3+520686.14684569274013*z2-116931.94197649079109*z+11388.124925664016721;
		else if(z < .96984631039295419207) res[1] = -65536.000000000000001*z8+492669.56476660233489*z7-1617015.8970528255834*z6+3025849.2798899800574*z5-3529930.2957369783692*z4+2628127.3752949626591*z3-1219107.5574278232691*z2+322004.20139574142904*z-37060.671129748843081;
		else if(z < 1.) res[1] = 65536.000000000000001*z8-524288.*z7+1835008.*z6-3670016.*z5+4587519.9999999999995*z4-3670015.9999999999987*z3+1835007.9999999999982*z2-524287.99999999999918*z+65536.000000000000220;
		else res[1] = 0.31457280000000000000e-14*z7-0.24489110467243902392e-15*z6-0.51726949156891583675e-13*z5-0.22604305032970467247e-12*z4+0.11878015214391867700e-11*z3-0.16139553825369369086e-11*z2+0.76200509207913199320e-12*z-0.38084892777430476847e-13;
		break;
	 case 10:
		if(z<=0) res[1] = 0;
		else if(z < 0.24471741852423213942e-1) res[1] = 262144.*z9;
		else if(z < 0.95491502812526287948e-1) res[1] = -262144.*z9+115472.16533090935792*z8-11303.220084473389574*z7+645.42212935283557452*z6-23.691905603045843470*z5+.57978219790971701110*z4-0.94588535185847604044e-2*z3+0.99203409225724825554e-4*z2-0.60692005536305934786e-6*z+0.16502656577670215063e-8;
		else if(z < .20610737385376343542) res[1] = 262144.*z9-335113.27590825468416*z8+160805.10363301873673*z7-37702.636986746027579*z6+5469.1787868067725761*z5-523.94269497518567730*z4+33.382167549286678230*z3-1.3664507622586196264*z2+0.32622870552451825926e-1*z-0.34613889315882092304e-3;
		else if(z < .34549150281252628792) res[1] = -262144.*z9+637423.32949912263201*z8-640982.75923565457040*z7+347890.94155518251150*z6-113741.34098542086716*z5+24046.224471020794398*z4-3342.6795856049777948*z3+296.84693005663766433*z2-15.333371321606293864*z+.35154771772884343770;
		else if(z < .50000000000000000000) res[1] =262144.*z9-992810.11174004141001*z8+1611944.4469601656401*z7-1468299.2062350646463*z6+827476.05434461419858*z5-301136.38791484468250*z4+71555.873375525614328*z3-10793.216050382476250*z2+942.54726002776215574*z-36.419521042260577435;
		else if(z < .65450849718747371208) res[1] = -262144.*z9+1366485.8882599585900*z8-3106647.5530398343600*z7+4036724.7937649353540*z6-3301291.9456553858019*z5+1763247.6120851553176*z4-616572.12662447438620*z3+136662.78394961752431*z2-17489.452739972238074*z+987.58047895773946459;
		else if(z < .79389262614623656452) res[1] = 262144.*z9-1721872.6705008773680*z8+4978780.1232426735137*z7-8311231.1471189657839*z6+8821471.1836522043814*z5-6171203.8654376715880*z4+2845538.4817491228090*z3-834471.84950041680795*z2+141414.51463655028914*z-10568.419174901716989;
		else if(z < .90450849718747371204) res[1] = -262144.*z9+2024182.7240917453158*z8-6917082.8963669812632*z7+13724857.363013253973*z6-17419961.082297461765*z5+14661675.710015696162*z4-8180507.8358141892028*z3+2917026.8078111741257*z2-603157.26562274453794*z+55110.474823368299617;
		else if(z < .97552825814757678607) res[1] = 262144.*z9-2243823.8346690906421*z8+8524709.8974371985264*z7-18865352.489196498872*z6+26797221.522372503165*z5-25333141.677598438360*z4+15936593.611224680219*z3-6431883.1292083867842*z2+1510884.8537459061118*z-157352.75410787171601;
		else if(z < 1.) res[1] = -262144.*z9+2359296.*z8-9437184.*z7+22020096.*z6-33030144.*z5+33030144.*z4-22020095.999999999985*z3+9437183.9999999999790*z2-2359295.9999999999898*z+262143.99999999999818;
		else res[1] = 0.14155776000000000000e-13*z8+0.43700290514366795354e-15*z7-0.15261309053650179771e-12*z6+0.27521191504065732696e-12*z5-0.38843716598458284833e-11*z4+0.14937561279078717714e-10*z3-0.21059423946214031200e-10*z2+0.12125594971675806726e-10*z-0.23771387168790120845e-11;
		break;
      
	 }
    
    
	 //res[1] = res[1]/scale;
	 
    
    
  }
  
  if(flags & EVAL_SD) {
	 assert(0);     
  }
  ebiFunctionReturn(0);
}	

// ----------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SplinelendFunc1D"
int splineBlendFunc1D(int flags,int deg, double t, double scale, double * res){
  enum {	 EVAL_VL=1,	 EVAL_FD=2,	 EVAL_SD=4  };
  //spline blending function of arbitrary degree
  if (deg==3){
    double x= 3.0*t;
    int tx= (int)ceil(x); 
    t = x-tx+1;
    double t2 = t*t; double t3 = t2*t;
    switch(tx){
    case 1:{ 
      res[0] = (1.0/6.0)*((-t3+3.0*t2-3.0*t+1.0)+(4.0-6.0*t2+3.0*t3)+(1.0+3.0*t+3.0*t2-3.0*t3));
      res[1] = -3.0*t2/2.;
    } break;
    case 2:{
      res[0] = (1.0/6.0)*((-t3+3.0*t2-3.0*t+1.0)+(4.0-6.0*t2+3.0*t3));
      res[1] = (1.0-2.0*t2+2.0*t)*(-3.0)/2.0;
    } break;
    case 3:{
      res[0] = (1.0/6.0)*(-t3+3.0*t2-3.0*t+1.0);
      res[1] = (1.0+t2-2.0*t)*(-3.0)/2.0;
    }break;
    default: assert(0); break;
    }
    
    return(0);
  }
  else{
	 
    double *poucp = new double[2*deg];
    for (int i=0; i<deg; i++){ 
      poucp[i] = 1.0;
      poucp[i+deg] = 0.0;
    }
	 
    double x = deg*t; 
    int tx = (int)ceil(x);
    int e = (int)floor((deg-1)/2.0);
    double ex = ceil((deg-1)/2.0);
    double *Nj = new double[2*deg];
	 
    if(flags & EVAL_VL){
      
      res[0] = 0.0;
      SplineBasis(x, tx,deg,2*deg, Nj);
      for(int i=tx-(e+1); i<=tx+e; i++)
		  res[0]+=Nj[i+e]*poucp[i+e];      
      res+=1;
    }
    
    if(flags &EVAL_FD){
      res[0] = 0.0;
      SplineBasis(x, tx, deg-1,2*deg,Nj);
      for(int i=tx-(e+1); i<tx+ex; i++){      
		  res[0] += Nj[i+e]*(poucp[i+1+e] - poucp[i+e]);
      }
      res[0] *= (double(2*deg)/2.);
      res+=1;
    }
    if(flags & EVAL_SD){
      assert(0);
    }
    
    delete[] Nj;
    delete[] poucp;
  }
  return 0;
}


// ----------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SplineBasis"
inline int SplineBasis(double u, int k, int p, int n, double *N){ 
  ebiFunctionBegin;
  
  assert(k>0);
 
  double x = u-k+1;
  double x2 = x*x;
  double x3 = x2*x;

  for(int i=0; i<n; i++) N[i] = 0.0;
  
  if (p==2){
    N[k-1] = (1.0+x2-2.0*x)/2.0;
    N[k]   = (1.0-2.0*x2+2.0*x)/2.0;
    N[k+1] = x2/(2.0);
    return(0);
  }
  else if (p==3){
    N[k-1] = (-x3+3.0*x2-3.0*x+1.0)*(1.0/6.0);
    N[k]   = (4.0-6.0*x2+3.0*x3)*(1.0/6.0);
    N[k+1] = (1.0+3.0*x+3.0*x2-3.0*x3)*(1/6.0);
    N[k+2] = x3*(1/6.0);
    return(0);
  }
  else if (p==4){
    double x4 = x2*x2;
    N[k-1] = (x4-4.0*x3+6.0*x2-4.0*x+1.0)/24.0;
    N[k] =   (11.0-12.0*x-6.0*x2+12.0*x3-4.0*x4)/24.0;
    N[k+1] = (11.0+12.0*x-6.0*x2-12.0*x3+6.0*x4)/24.0;
    N[k+2] = (1.0+4.0*x+6.0*x2+4.0*x3-4.0*x4)/24.0;
    N[k+3] = x4/24.0;
    return(0);
  }
  else {
    k--;
    N[k+p] = 1.0;
    
    for (int d=1; d<=p; d++){
        N[k-d+p] = ((k+1-u)/d)*N[k-d+1+p];
        for (int i=k-d+1; i<=k-1; i++)
            N[i+p] = ((u-i)/d)*N[i+p] + ((i+d+1-u)/d)*N[i+1+p];
        N[k+p] = ((u-k)/d)*N[k+p];
    }
  }
  ebiFunctionReturn(0);
}


END_EBI_NAMESPACE

