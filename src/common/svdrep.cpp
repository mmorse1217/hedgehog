/*! \file */
#include "mobo_blas.h"
#include "mobo_lapack.h"
#include "svdrep.hpp"

BEGIN_EBI_NAMESPACE

using std::min;
using std::max;
using std::abs;

//int    SVDRep::_wssize = 4194304;
//double SVDRep::_wsbuf[4194304];

/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "SVDRep::construct"
ebiEC SVDRep::construct(double epsilon, const DblNumMat& K, const int type)
{
  ebiFunctionBegin;
  
  PetscBLASInt m = K.m();
  PetscBLASInt n = K.n();
  PetscBLASInt k = min(m, n);
  
  DblNumMat tU(m, k);
  DblNumVec tS(k);
  DblNumMat tVT(k, n);
  PetscBLASInt wssize = 4194304;

  //SVD
  if(1) {	
	 PetscBLASInt INFO;
	 char JOBU  = 'S';
	 char JOBVT = 'S';
	 ebiAssert( wssize >= max(3*min(m,n)+max(m,n), 5*min(m,n)));
	 //double* wsbuf = _wsbuf;
	 double* wsbuf = new double[wssize];
	 DGESVD(&JOBU, &JOBVT, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, &INFO);
	 ebiAssert(INFO==0);
	 delete[] wsbuf;
	 
  }
  if(0) {
	 PetscBLASInt INFO;
	 char JOBZ = 'S';

#if defined(EBI_linux_mpi) || defined(EBI_linux)
	  PetscBLASInt wssize = 4*min(m,n)*min(m,n) + max(m,n) + 9*min(m,n);
#endif
#if defined(EBI_solaris)
	  PetscBLASInt wssize = 3*min(m,n)*min(m,n) + max(max(m,n),4*min(m,n)*min(m,n)+4*min(m,n));
#endif
#if defined(EBI_alpha)
	 PetscBLASInt wssize = 4*min(m,n)*min(m,n) + max(m,n) + 9*min(m,n);
	 ebiAssert(0);
#endif
	 PetscBLASInt wisize = 8*min(m,n);
	 	 
	 double* wsbuf = new double[wssize];
	 PetscBLASInt*    wibuf = new PetscBLASInt[   wisize];
	 DGESDD(&JOBZ, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, wibuf, &INFO);
	 delete[] wsbuf;
	 delete[] wibuf;
	 
	 ebiAssert(INFO==0);
  }
  
  //cutoff
  double cutoff = 0.0;
  if (type == 0){ 
	 cutoff = epsilon*tS(0);
  }
  else {
	 cutoff = 0.0;
  }
  int cnt=0;
  while(cnt< k)
    if(abs(tS(cnt)) >= cutoff)
      cnt++;
    else
      break;
  
  _matU.resize(m, cnt);
  _matS.resize(cnt);	
  _matVT.resize(cnt,n);
  
  for(int i=0; i<m; i++)
    for(int j=0; j<cnt; j++)
      _matU(i,j) = tU(i,j);
  for(int i=0; i<cnt; i++)
    _matS(i) = tS(i);
  for(int i=0; i<cnt; i++)
    for(int j=0; j<n; j++)
      _matVT(i,j) = tVT(i,j);
  
  ebiFunctionReturn(0);
}

/* ********************************************************************** 
#undef __FUNCT__
#define __FUNCT__ "SVDRep::dgemv"
ebiEC SVDRep::dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol)
{
  ebiFunctionBegin;

  ebiAssert(Y.m() == _matU.m());
  ebiAssert(X.m() == _matVT.n());
  iC( dgemv(alpha, X.data(), beta, Y.data(), tol) );
  
  ebiFunctionReturn(0);
}

 ********************************************************************** 

#undef __FUNCT__
#define __FUNCT__ "SVDRep::dgemv"
ebiEC SVDRep::dgemv(double alpha, double* X, double beta, double* Y, double tol)
{
  ebiFunctionBegin;
  iA(0);
  int K = 1; //prevent matrix of zero size
  while(K<_matS.m() && _matS(K)>=tol) 
	 K++;
  //buf = VT(1:K,:) * X
  double* buf;// = _wsbuf;
  {
	 char TRANS = 'N';
	 int M = _matVT.m();
	 int N = _matVT.n();
	 double ALPHA = 1.0;
	 double BETA = 0.0;
	 int INC = 1;
	 DGEMV(&TRANS, &K, &N, &ALPHA, _matVT.data(), &M, X, &INC, &BETA, buf, &INC); //first K rows
  }
  // buf = S(1:K) .* buf;
  for(int i=0; i<K; i++)
	 buf[i] = buf[i] * _matS(i);
  // y = U(:,1:K) * buf
  {
	 char TRANS = 'N';
	 int M = _matU.m(); //int N = _matU.n();
	 double ALPHA = alpha;
	 double BETA = beta;
	 int INC = 1;
	 DGEMV(&TRANS, &M, &K, &ALPHA, _matU.data(), &M, buf, &INC, &BETA, Y, &INC);	
  }
  
  PetscLogFlops( 2 * (_matVT.m()*_matVT.n()+_matU.m()*_matU.n()) );
  ebiFunctionReturn(0);
}
*/

END_EBI_NAMESPACE
