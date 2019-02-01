#include "mobo_lapack.h"
#include "Dcuhre.hpp"
#include "memAlloc.hpp"

double glb_xtarg, glb_ytarg, glb_ztarg;
int glb_kt;
double glb_lambda;
int glb_nk;
int glb_k;
int glb_sdof;
int glb_tdof;

BEGIN_EBI_NAMESPACE

#undef __FUNCT__
#define __FUNCT__ "Dcuhre::setup"
int Dcuhre::setup(){
  ebiFunctionBegin;
  std::cout << "DCUHRE SETUP" << std::endl;

  int64_t wksze, maxpts, minpts, key;
  PetscBool flg = PETSC_FALSE;
  iC( PetscOptionsGetInt(NULL,  prefix().c_str(),  "-minpts", &minpts, &flg) );
  if (flg != true) { minpts = 5000; }
  
  iC( PetscOptionsGetInt(NULL,  prefix().c_str(),  "-maxpts", &maxpts, &flg) );
  if (flg != true) { maxpts = 50000000; }

  iC( PetscOptionsGetInt(NULL,  prefix().c_str(),  "-worksize", &wksze, &flg) );
  if (flg != true) { wksze = 999999999; }

  iC( PetscOptionsGetInt(NULL,  prefix().c_str(),  "-key", &key, &flg) ); 
  if (flg != true) { key = 0; }
  double epsabs, epsrel;
  double eps;

  iC( PetscOptionsGetReal(NULL, prefix().c_str(),  "-eps", &eps, &flg) );
  if (flg == true){
	 epsabs = eps;
	 epsrel = eps;
  }
  else {
	 iC( PetscOptionsGetReal(NULL, prefix().c_str(),  "-epsabs", &epsabs, &flg) );  
	 iC( PetscOptionsGetReal(NULL, prefix().c_str(),  "-epsrel", &epsrel, &flg) );
  }

  if (flg != true){
	 epsabs = 10e-12;
	 epsrel = 10e-12;
  } 

  /* This will need to change if kval != 4 and Nk != 20 at some point */
  /* Should probably only allocate if needed */

  //double eps = 1e-12; /* Slows things down but seems necessary */
  //_ndim = 3;
  _ndim = 1; // changed mjm 2/9
  _numfun = 0;
  _minpts = minpts;
  _maxpts = maxpts;
  _epsabs = epsabs;
  _epsrel = epsrel;
  _Alim=vecalloc((short int) _ndim); /*allocate memory for ndim double vector*/
  _Blim=vecalloc((short int) _ndim);
  //_Alim[0] = -1.0; _Alim[1] = -1.0; _Alim[2] = -1.0;// changed mjm 2/9
  //_Blim[0] = 1.0; _Blim[1] = 1.0; _Blim[2] = 1.0;// changed mjm 2/9
  _Alim[0] = 0.;
  _Blim[0] = 1.;
  _key = key;
  _restar = 0;
  _worksize = wksze;
  _isalloc = false;
  cerr << "TBLS EPS = " << _epsabs << " " << _epsrel << endl;
  ebiFunctionReturn(0);
}

Dcuhre::~Dcuhre(){
  if (_isalloc){
	 free(_Alim); free(_Blim);
	 free(_result);
	 free(_abserr);
	 free(_work);
  }
}

void Dcuhre::setLims(double a){
  setLims(a,-a);
}

void Dcuhre::setLims(double a, double b){
  setLims(a,a,a,b,b,b);
}

void Dcuhre::setLims(double ax, double ay, double az, double bx, double by, double bz){
  _Alim[0] = ax; _Alim[1] = ay; _Alim[2] = az;
  _Blim[0] = bx; _Blim[1] = by; _Blim[2] = bz;
}

void Dcuhre::printLims(){
  for (int i = 0; i < _ndim; i++){
	 std::cout << "A(" << i << ") = " << _Alim[i] << "  B(" << i << ") = " << _Blim[i] << std::endl;
  }
}

void Dcuhre::alloc(){
  //iA(_numfun != (long)0);
  //iA(_numfun == _nk*_sdof*_tdof);
  //iA(_nk == _k*(_k+1)*(_k+2)/6);
  _result=vecalloc((int) _numfun);  /* resultant approximation to the integral(s) */
  _abserr=vecalloc((int) _numfun);  /* An array of estimates for absolute errors */
  _work=vecalloc((int) _worksize);
  _isalloc = true;
} 

#undef __FUNCT__
#define __FUNCT__ "Dcuhre::eval"
int Dcuhre::eval(double x, double y, double z){
  ebiFunctionBegin;
  if (!_isalloc) alloc();
  glb_xtarg = x; glb_ytarg = y; glb_ztarg = z;
  glb_kt = _kt;
  glb_lambda = _lambda;
  glb_nk = _nk;
  glb_k = _k;
  glb_sdof = _sdof;
  glb_tdof = _tdof;
  
  DCUHRE(&_ndim,
          &_numfun,
          _Alim,
          _Blim,
          &_minpts,
          &_maxpts,
          (void *)func,
          &_epsabs,
          &_epsrel,
          &_key,
          &_worksize,
          &_restar,
          _result,
          _abserr,
          &_neval,
          &_ifail,
          _work);

  /*
     for (int i = 0; i < 3; i++){ 
     std::cout << _Alim[i] << " " << _Blim[i] << " ";
  }
  std::cout << std::endl;
  std::cout << _ifail << std::endl;
  */
  ebiFunctionReturn((int)_ifail);
}

END_EBI_NAMESPACE
