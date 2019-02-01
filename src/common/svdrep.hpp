/*! \file */
#ifndef _SVDREP_HPP_
#define _SVDREP_HPP_

//#include "ebi_namespace.hpp"
#include "nummat.hpp"

BEGIN_EBI_NAMESPACE

using std::endl;

class SVDRep
{
public:
  SVDRep()   {}
  SVDRep(const SVDRep& C): _matU(C._matU), _matS(C._matS), _matVT(C._matVT)  {}
  ~SVDRep()  {}

  SVDRep& operator=(const SVDRep& c)  { _matU = c._matU; _matS = c._matS; _matVT = c._matVT; return *this; }
  //access
  DblNumMat& U() { return _matU; }
  DblNumVec& S() { return _matS; }
  DblNumMat& VT(){ return _matVT; }
  //ops
  int construct(double epsilon, const DblNumMat& M, const int type=0);
  
  //int dgemv(double alpha, const DblNumVec& X, double beta, DblNumVec& Y, double tol=0.0); // y <- a Mx + b y
  //int dgemv(double alpha, double* X, double beta, double* Y, double tol=0.0);
  
  int m() const { return _matU.m(); }
  int k() const { return _matS.m(); } //length
  int n() const { return _matVT.n(); }
  
protected:
  DblNumMat _matU;
  DblNumVec _matS;
  DblNumMat _matVT;
  //static int _wssize;
  //static double _wsbuf[];
  //static int _wisize;
  //static int _wibuf[];
};	

inline ostream& operator<<( ostream& os, SVDRep& svdrep)
{
  os<<svdrep.U().m()<<" "<<svdrep.S().m()<<" "<<svdrep.VT().n()<<endl;
  os<<svdrep.U()<<svdrep.S()<<svdrep.VT()<<endl;
  return os;
}

//int matvec(double a, const DblNumVec& X, double b, DblNumVec& Y); // y <- a Mx + b y
//int matvec(double a, double* X, double b, double* Y);


END_EBI_NAMESPACE
#endif
