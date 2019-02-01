/*! \file */
#ifndef _NUMMAT_HPP_
#define _NUMMAT_HPP_

extern "C" {
#include <string.h>
# include <unistd.h>
}
#include "numvec.hpp"
//#include <cstring>

BEGIN_EBI_NAMESPACE

template <class F>
class NumMat
{
public:
  int _m, _n;
  bool _owndata;
  bool __teardown;
  F* _data;
  Vec _v;
public:
  NumMat(){
	 _m = 0; _n = 0; _owndata = true; _data = NULL; _v = NULL;
  }
  NumMat(int m, int n): _m(m), _n(n), _owndata(true), __teardown(false), _v(NULL) {
	 if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
  }
  
  NumMat(int m, int n, bool owndata, F* data): _m(m), _n(n), _owndata(owndata),  __teardown(false), _v(NULL) {
	 if(_owndata) {
		if(_m>0 && _n>0) {
            _data = new F[_m*_n];
            assert( _data!=NULL );
            memset(_data, 0, _m*_n*sizeof(F));
        } else { 
            _data=NULL;
        }
		if(_m>0 && _n>0) memcpy( _data, data, _m*_n*sizeof(F) );
	 } else {
		_data = data;
	 }
  }
  NumMat(int m, Vec v): _m(m),  _owndata(false),  __teardown(true), _v(v) {
      int64_t size; 
      VecGetLocalSize(v, &size);
      _n = size/m;
      VecGetArray(v, &_data);
  }
  NumMat(int m, int n, Vec v): _m(m), _n(n), _owndata(false), __teardown(true), _v(v) {
      VecGetArray(v, &_data);
  }

  
  NumMat(const NumMat& C): _m(C._m), _n(C._n), _owndata(C._owndata), _v(C._v) {
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
		if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
	 } else {
		_data = C._data;
	 }
  }

  ~NumMat() { 
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 } 
  }
  NumMat& operator=(const NumMat& C) {
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
	 _m = C._m; 
     _n=C._n; 
     _owndata=C._owndata;
     _v = C._v;
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
		if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
	 } else {
		_data = C._data;
	 }
	 return *this;
  }

  void appendRows(const NumMat& C){
	 assert(_owndata);
	 int m = _m;
	 int n = _n;
	 if (m != 0) assert (m == C._m);
	 else resize(C._m, C._n);
	 if (n < (n + C._n) && n != 0){
		F* tempdata = new F[m*n];
		if (n > 0){
		  assert(tempdata != NULL);
		  memset(tempdata, 0, m*n*sizeof(F));
		}
		else tempdata=NULL;
		if (n > 0) memcpy(tempdata, _data, m*n*sizeof(F));
		resize(m, n + C._n);
		if (n > 0) memcpy(_data, tempdata, m*n*sizeof(F));
		memcpy(_data + (m*n), C._data, (C._m)*(C._n)*sizeof(F));
		delete[] tempdata; tempdata = NULL;
	 }
	 else {
		memcpy(_data + (m*n), C._data, (C._m)*(C._n)*sizeof(F));
	 }
  }
  
  void resize(int m, int n)  {
	 assert( _owndata==true );
     assert(_v == NULL);
	 if(_m!=m || _n!=n) {
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
		_m = m; _n = n;
		if(_m>0 && _n>0) {
		  _data = new F[_m*_n];
		  assert( _data!=NULL );
		  memset(_data, 0, _m*_n*sizeof(F));
		} else _data=NULL;
	 }
  }
  const F& operator()(int i, int j) const  { 
	 assert( i>=0 && i<_m && j>=0 && j<_n );
	 return _data[i+j*_m];
  }
  F& operator()(int i, int j)  { 
	 assert( i>=0 && i<_m && j>=0 && j<_n );
	 return _data[i+j*_m];
  }
  
  F* data() const { return _data; }
  F* clmdata(int j) { return &(_data[j*_m]); }
  int m() const { return _m; }
  int n() const { return _n; }
  void restore_local_vector();
};

template <class F> inline ostream& operator<<( ostream& os, const NumMat<F>& mat)
{
  os<<mat.m()<<" "<<mat.n()<<endl;
  os<<"mat=[";
  os.precision(15);
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=0; i<mat.m(); i++) {
	 for(int j=0; j<mat.n(); j++)
		os<<" "<<mat(i,j);
	 os<<";"<<endl;
  }
   os<<"];";

  return os;
}
template <class F> inline void setvalue(NumMat<F>& M, F val)
{
  for(int i=0; i<M.m(); i++)
	 for(int j=0; j<M.n(); j++)
		M(i,j) = val;
}
template <class F> inline void clear(NumMat<F>& M)
{
  memset(M.data(), 0, M.m()*M.n()*sizeof(F));
}

typedef NumMat<bool>   BolNumMat;
typedef NumMat<int>    IntNumMat;
typedef NumMat<double> DblNumMat; 
typedef NumMat<imaginary> ComplexNumMat;

NumMat<double> get_local_vector(int m, int n, Vec v);
template<>
void NumMat<double>::restore_local_vector();
template<>
NumMat<double>::~NumMat<double>();


/*
  void getColumn(int j, Vector<F>& vec)  {
  assert( j>=0 && j<n() );
  vec.resize(m());
  for(int i=0; i<m(); i++)
  vec(i) = (*this)(i,j);
  }
  void getRow(int i, Vector<F>& vec)  {
  assert( i>=0 && i<m() );
  vec.resize(n());
  for(int j=0; j<n(); j++)
  vec(j) = (*this)(i,j);
  }
  void setColumn(int j, Vector<F>& vec)  {
  assert( j>=0 && j<n() );
  assert( vec.length() == m() );
  for(int i=0; i<m(); i++)
  (*this)(i,j) = vec(i);
  }
  void setRow(int i, Vector<F>& vec)  {
  assert( i>=0 && i<m());
  assert( vec.length() == n());
  for(int j=0; j<n(); j++)
  (*this)(i,j) = vec(j);
  }
*/


END_EBI_NAMESPACE

#endif




