/*! \file */
#ifndef _OFFMAT_HPP_
#define _OFFMAT_HPP_

#include "offvec.hpp"

BEGIN_EBI_NAMESPACE

template <class F>
class OffMat
{
public:
  int _m, _n;
  int _s, _t;
  bool _owndata;
  F* _data;
public:
  OffMat(int m=0, int n=0, int s=0, int t=0): _m(m), _n(n), _s(s), _t(t), _owndata(true) {
	 if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
  }
  OffMat(int m, int n, int s, int t, bool owndata, F* data): _m(m), _n(n), _s(s), _t(t), _owndata(owndata) {
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
		if(_m>0 && _n>0) memcpy( _data, data, _m*_n*sizeof(F) );
	 } else {
		_data = data;
	 }
  }
  OffMat(const OffMat& C): _m(C._m), _n(C._n), _s(C._s), _t(C._t), _owndata(C._owndata) {
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
		if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
	 } else {
		_data = C._data;
	 }
  }
  ~OffMat() { 
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
  }
  OffMat& operator=(const OffMat& C) {
	 if(_owndata) { 
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
	 }
	 _m = C._m; _n=C._n; _s = C._s; _t=C._t; _owndata=C._owndata;
	 if(_owndata) {
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
		if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
	 } else {
		_data = C._data;
	 }
	 return *this;
  }
  void resize(int m, int n, int s, int t)  {
	 assert( _owndata==true );
	 if(_m!=m || _n!=n) {
		if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
		_m = m; _n = n; _s = s; _t = t;
		if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); memset(_data, 0, _m*_n*sizeof(F)); } else _data=NULL;
	 }
  }
  const F& operator()(int i, int j) const  {
 	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	 return _data[(i-_s) + (j-_t)*_m];
  }
  F& operator()(int i, int j)  {
	 assert( i>=_s && i<_m+_s && j>=_t && j<_n+_t );
	 return _data[(i-_s) + (j-_t)*_m];
  }
  
  F* data() const { return _data; }
  F* clmdata(int j) { return &(_data[j*_m]); }
  int m() const { return _m; }
  int n() const { return _n; }
  int s() const { return _s; }
  int t() const { return _t; }
};

template <class F> inline ostream& operator<<( ostream& os, const OffMat<F>& mat)
{
  os<<mat.m()<<" "<<mat.n()<<" "<<mat.s()<<" "<<mat.t()<<endl;
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=mat.s(); i<mat.s()+mat.m(); i++) {
	 for(int j=mat.t(); j<mat.t()+mat.n(); j++)
		os<<" "<<mat(i,j);
	 os<<endl;
  }
  return os;
}
template <class F> inline void setvalue(OffMat<F>& mat, F val)
{
  for(int i=mat.s(); i<mat.s()+mat.m(); i++)
	 for(int j=mat.t(); j<mat.t()+mat.n(); j++)
		mat(i,j) = val;
}
template <class F> inline void clear(OffMat<F>& M)
{
  memset(M.data(), 0, M.m()*M.n()*sizeof(F));
}

typedef OffMat<bool>   BolOffMat;
typedef OffMat<int>    IntOffMat;
typedef OffMat<double> DblOffMat;

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




