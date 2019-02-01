/*! \file */
#ifndef _ROTATION_HPP_EBI_
#define _ROTATION_HPP_EBI_

#include "quaternion.hpp"

BEGIN_EBI_NAMESPACE

/* ********************************************************************** */
/// implements rotations need to represent the orientation 
///
class Rotation {
public:
  Rotation(double ain) : dof_(1), dim_(DIM2), Q_(Quaternion::ZERO),a_(ain) {} /// 2d constructro
  Rotation(Quaternion Qin) : dof_(4), dim_(DIM3), Q_(Qin),a_(0) {}  /// 3d constructor
  Rotation() : dof_(4), dim_(DIM3), Q_(Quaternion::ZERO),a_(0) {}  ///< DEFAULT IS 3D
  Rotation(DimType d) : dim_(d), Q_(Quaternion::ZERO),a_(0) { if(d==DIM2) dof_=1; else dof_=4;}

  void toDouble(double *arr) const;
  void toDoubleAdd(double *arr) const;
  void setToZero() {Q_=Quaternion::ZERO;a_=0;}

  /// assignment
  Rotation & operator=( const Rotation &R);
  Rotation & operator-=( Rotation &R);
  Rotation & operator+=( Rotation &R);  
  
  Quaternion & Q() {return Q_;};
  double & a() {return a_;}
  const Quaternion & Q() const {return Q_;};
  const double & a() const {return a_;}  
  const DimType dim() const {return dim_;}
  const int dof() const {return dof_;}                                  ///< degrees of freedom
  CPoint rotate( const CPoint &cen, const CPoint &x) const;             ///< rotation about cen
  CPoint operator*( const CPoint &v) const; ///< rotate vector v

  void rotateInertia( const DblNumMat &In, DblNumMat &Out) const; ///< Updates Inertia (tensor in 3D, nothing in 2D)

  void toRotationMatrix( DblNumMat &Rot) const;
  
private:
  int dof_;
  DimType dim_;
  Quaternion Q_;
  double a_;
  
  CPoint rotate2d( const CPoint &cen, const CPoint &p) const;
  CPoint rotate3d( const CPoint &cen, const CPoint &p) const;
};



/* ********************************************************************** */
/* ********************************************************************** */

/* ********************************************************************** */
inline Rotation &Rotation::operator-=( Rotation &R){
  Q_ -= R.Q();  a_ -= R.a(); return *this;
}

/* ********************************************************************** */
inline Rotation &Rotation::operator+=( Rotation &R){
  Q_ += R.Q();  a_ -= R.a(); return *this;
}

/* ********************************************************************** */
inline  CPoint Rotation::operator*( const CPoint &v) const {
  switch(dim_){
  case DIM2 : return rotate2d( CPoint(DIM2, 0.0), v); 
  case DIM3 : return Q_*v;
  }
}
/* ********************************************************************** */
inline Rotation& Rotation::operator=( const Rotation &R){
  dim_ = R.dim();
  dof_ = R.dof();
  Q_   = R.Q();
  a_   = R.a();
  return *this;
}

/* ********************************************************************** */
inline void Rotation::toDouble(double *arr) const {
  if(dim_==DIM2) { arr[0] = this->a_; }
  if(dim_==DIM3) { arr[0] = this->Q_.w; arr[1] = this->Q_.x; arr[2] = this->Q_.y; arr[3] = this->Q_.z;}
}
/* ********************************************************************** */
inline void Rotation::toDoubleAdd(double *arr) const {
  if(dim_==DIM2) { arr[0] += this->a_; }
  if(dim_==DIM3) { arr[0] += this->Q_.w; arr[1] += this->Q_.x; arr[2] += this->Q_.y; arr[3] += this->Q_.z;}
}

/* ********************************************************************** */
inline CPoint Rotation::rotate( const CPoint &cen, const CPoint &p) const {
  return dim_==DIM2 ? rotate2d(cen,p) : rotate3d(cen,p);
}

/* ********************************************************************** */
inline CPoint Rotation::rotate2d( const CPoint &cen, const CPoint &p) const {
  double tmp[2];
  tmp[0] =  cos(a())*(p[0] - cen[0]) - sin(a())*(p[1] - cen[1]);
  tmp[1] =  sin(a())*(p[0] - cen[0]) + cos(a())*(p[1] - cen[1]);
  CPoint ret(DIM2,tmp);
  return ret;
}

/* ********************************************************************** */
inline CPoint Rotation::rotate3d( const CPoint &cen, const CPoint &p) const { return  Q()*(p-cen);  }  


/* ********************************************************************** */
inline void Rotation::toRotationMatrix( DblNumMat &Rot) const
{
  if( this->dim() == DIM3 ) this->Q().toRotationMatrix( Rot);
  if( this->dim() == DIM2 ){
	 assert( Rot.m()==DIM2 &&  Rot.n() == DIM2 );
	 Rot(0,0) =  cos(a_);    Rot(0,1) = -sin(a_);
	 Rot(1,0) =-Rot(0,1);    Rot(1,1) = Rot(0,0);
  }
}

/* ********************************************************************** */
// Auxiliary inlines
/* ********************************************************************** */

inline Rotation doubleToRotation(double *arr, DimType d){
  if(d==DIM2) { Rotation  R(DIM2); R.a() = arr[0]; return R; }
  if(d==DIM3) { Rotation  R(DIM3); Quaternion tmp( arr[0],arr[1],arr[2],arr[3]); R.Q()=tmp; return R;}
}


/* ********************************************************************** */
/// Given rotation, its center, and R
inline void RigidBodyTransformation( const CPoint &center, const CPoint &translation, const Rotation &R,  vector<CPoint> & points)
{
  for(vector<CPoint>::iterator pi=points.begin(); pi!=points.end(); pi++){
	 *pi = center + translation + R.rotate( center, *pi);
  }
}
inline void RigidBodyTransformation( const CPoint &center, const CPoint &translation, const Rotation &R, CPoint & point)
{
  point = center + translation + R.rotate( center, point);
}

/// Apply DblNumMat to a vector
inline CPoint operator*( const DblNumMat &M, const CPoint p)
{
  assert( M.n() == p.dim() );
  assert( M.m() == p.dim() );

  CPoint tmp = p;
  for( int i=0; i<M.m(); i++){
	 tmp(i) = 0;
	 for( int j=0; j<M.n(); j++) tmp(i) += M(i,j) * p(j);
  }
  return tmp;
}

/* ********************************************************************** */
/// Convert a cross product with CPoint to a DblNumMat operator
inline void CPointCrossToDblNumMat( const CPoint &v, DblNumMat &V){
  assert(v.dim()!=DIM2);
  assert(V.m() == V.n());
  assert(V.m() == DIM3);
  V(0,0)=V(1,1)=V(2,2)=0;
  V(0,1)=-v[2];
  V(0,2)= v[1];
  V(1,0)= v[2];
  V(1,2)=-v[0];
  V(2,0)=-v[1];
  V(2,1)= v[0];
}  


/* ********************************************************************** */
inline void Rotation::rotateInertia( const DblNumMat &In, DblNumMat &Out) const
{
  if( dim_==DIM2){	 Out = In; }
  if( dim_==DIM3){
	 DblNumMat R(DIM3,DIM3),Wrk(DIM3,DIM3);
	 toRotationMatrix(R);
	 DblNumMatMult( R, In, Wrk);
	 DblNumMatTranspose(R);
	 DblNumMatMult(Wrk,R,Out);
  }
}


  
  


END_EBI_NAMESPACE

#endif // _ROTATION_HPP_EBI_
