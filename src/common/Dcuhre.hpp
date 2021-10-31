#ifndef _VFMM_DCUHRE_HPP_
#define _VFMM_DCUHRE_HPP_

#include "ebi.hpp"
#include "ebiobject.hpp"

extern "C" {
  void func(long *ndim, double *pnt, long *Xnumfun, double *funvls);
}

BEGIN_EBI_NAMESPACE

class Dcuhre: public EbiObject
{
private:
  /* All of the following are for Gaussian Quadrature */
  /* All of this is set in allocateGaussQuad() if needed */
  long _ndim;      /* Number of variables in integrand(s) */
  long _numfun;    /* Number of components in the vector integrand function */
  double * _Alim;       /* Lower limit of integration */
  double * _Blim;       /* Upper limit of integration */
  long _minpts;    /* Minimum number of funsub calls */
  long _maxpts;    /* Maximum number of funsub calls */
  double _epsabs;  /* Requested absolute error */
  double _epsrel;  /* Requested relative error */
  long _key;       /* Selects integration rule */
  long _worksize;        /* Length of work array */
  long _restar;    /* restar = 0, this is the first attempt
							= 1, continuation of a previous attempt */
  double *_result;  /* An array of the approximation to the integrals */
  double *_abserr;  /* An array of estimates for absolute errors */
  long _neval;     /* Number of function evaluations used */
  long _ifail;     /* Exit condition indicator */
  double *_work;    /* Work array */

  int _kt;
  int _nk;
  double _lambda;
  int _k;
  int _sdof;
  int _tdof;

  bool _isalloc;

  void alloc();

public:
  Dcuhre(const string& n, const string& p): EbiObject(n,p) {;}
  ~Dcuhre();

  int setup();

  long& ndim() { return _ndim; }
  long& numfun() { return _numfun; }
  double*& A() { return _Alim; }
  double*& B() { return _Blim; }
  long& minPts() { return _minpts; }
  long& maxPts() { return _maxpts; }
  double& epsAbs() { return _epsabs; }
  double& epsRel() { return _epsrel; }
  long& key() { return _key; }
  long& worksize() { return _worksize; }
  long& restar() { return _restar; }
  double*& result() { return _result; }
  double*& absErr() { return _abserr; }
  long& numEval() { return _neval; }
  long& iFail() { return _ifail; }
  double*& work() { return _work; }

  int& kt() { return _kt; }
  int& nk() { return _nk; }
  int& k() { return _k; }
  int& sdof() { return _sdof; }
  int& tdof() { return _tdof; }
  
  double& lambda() { return _lambda; }

  bool& isalloc() { return _isalloc; }

  void setLims(double a);
  void setLims(double a, double b);
  void setLims(double ax, double ay, double az,double bx, double by, double bz);
  void printLims();

  int eval(double x, double y, double z);

  double res(int n) { iA( n >= 0 && n <= _numfun); return _result[n]; }

};

END_EBI_NAMESPACE

#endif
