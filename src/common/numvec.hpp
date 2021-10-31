/*! \file */
#ifndef _NUMVEC_HPP_
#define _NUMVEC_HPP_

#include "ebi.hpp"
BEGIN_EBI_NAMESPACE

using std::ostream;
using std::ios_base;
using std::endl;

#include <complex>

typedef std::complex<double> imaginary ;

template <class F>
class NumVec
{
public:
  int  _m;
  bool _owndata;
  F* _data;
  Vec _v;
public:
  NumVec(int m=0): _m(m), _owndata(true), _v(NULL)  {
	 if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
  }
  NumVec(int m, bool owndata, F* data): _m(m), _owndata(owndata), _v(NULL) {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, data, _m*sizeof(F) );
	 } else {
		_data = data;
	 }
  }
  NumVec(int m, Vec v): _m(m), _owndata(false), _v(v) {
      VecGetArray(v, &_data);
  }
  NumVec(Vec v):  _owndata(false), _v(v) {
      int64_t size; 
      VecGetLocalSize(v, &size);
      _m = size;
      VecGetArray(v, &_data);
  }
  NumVec(const NumVec& C): _m(C._m), _owndata(C._owndata), _v(NULL) {
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data = C._data;
	 }
  }
  ~NumVec() {
	 if(_owndata) {
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
  }
  NumVec& operator=(const NumVec& C)  {
	 if(_owndata) { 
		if(_m>0) { delete[] _data; _data = NULL; }
	 }
	 _m = C._m; _owndata=C._owndata; _v = C._v;
	 if(_owndata) {
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL); memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
		if(_m>0) memcpy( _data, C._data, _m*sizeof(F) );
	 } else {
		_data =C._data;
	 }
	 return *this;
  }

  void append(const NumVec& C){
	 assert(_owndata);
	 int m = _m;
	 if (m < (_m + C._m)){
		F* tempdata = new F[m];
		if (_m > 0){
		  assert(tempdata != NULL);
		  memset(tempdata, 0, m*sizeof(F));
		}
		else tempdata=NULL;
		if (m > 0) memcpy(tempdata, _data, m*sizeof(F));
		resize(m + C._m);
		if (m > 0) memcpy(_data, tempdata, m*sizeof(F));
		memcpy(_data + (m), C._data, (C._m)*sizeof(F));
		delete[] tempdata; tempdata = NULL;
	 }
	 else {
		memcpy(_data + (m), C._data, (C._m)*sizeof(F));
	 }
  } 
  
  void resize(int m)  {
	 assert(_owndata==true);
     assert(_v == NULL);
	 if(m !=_m) {
		if(_m>0) { delete[] _data; _data = NULL; }
		_m = m;
		if(_m>0) { _data = new F[_m]; assert(_data!=NULL);  memset(_data, 0, _m*sizeof(F)); } else _data=NULL;
	 }
  }
  const F& operator()(int i) const  {	 assert(i>=0 && i<_m);
	 return _data[i]; 
  }
  F& operator()(int i)  {	 assert(i>=0 && i<_m);
	 return _data[i]; 
  }
  
  F* data() const { return _data; }
  int m () const { return _m; }
  F linfty( void ) const  { F cur=F(0); for(int i=0; i<_m; i++) cur=max(cur,abs(_data[i])); return cur; }
  F l2( void ) const  {
      F sum_of_squares=F(0);
      for(int i=0; i<_m; i++)
          sum_of_squares += pow(_data[i], 2);
      return sqrt(sum_of_squares);
  }
void restore_local_vector();
};

template <class F>
double dot(NumVec<F> a, NumVec<F> b){
    assert(a.m() == b.m());
    F dot_product = 0.;
    for(int i = 0; i < a.m(); i++)
        dot_product += a(i)*b(i);
    return dot_product;
}


template <class F> inline ostream& operator<<( ostream& os, const NumVec<F>& vec)
{
  os<<vec.m()<<endl;
  os<<"vec= [";
  os.precision(15);
  os.setf(ios_base::scientific, ios_base::floatfield);
  for(int i=0; i<vec.m(); i++)	 os<<" "<<vec(i);
  os<<"];";
  os<<endl;
  return os;
}
template <class F> inline void setvalue(NumVec<F>& V, F val)
{
  for(int i=0; i<V.m(); i++)
	 V(i) = val;
}
template <class F> inline void clear(NumVec<F>& V)
{
  memset(V.data(), 0, V.m()*sizeof(F));
}

typedef NumVec<bool>   BolNumVec;
typedef NumVec<int>    IntNumVec;
typedef NumVec<char>   ChrNumVec;
typedef NumVec<double> DblNumVec; 
typedef NumVec<imaginary> ComplexNumVec; 

NumVec<double> get_local_vector(int m, Vec v);
template<>
void NumVec<double>::restore_local_vector();

END_EBI_NAMESPACE
#endif


