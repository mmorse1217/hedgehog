/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#include "vecmatop.hpp"
#include "kernel3d.hpp"

BEGIN_EBI_NAMESPACE

//double Kernel3d::_mindif = 1e-8;

// ---------------------------------------------------------------------- 
uint Kernel3d::srcDOF() const
{
  int dof = 0;
  switch(_kernelType) {
	 //laplace kernels
  case KNL_LAP_S_U: dof = 1; break;
  case KNL_LAP_D_U: dof = 1; break;
  case KNL_LAP_E_U: dof = 1; break;
  case KNL_LAP_I  : dof = 1; break;
	 //Modified Laplacian
  case KNL_MODHEL_S_U: dof = 1; break;
  case KNL_MODHEL_D_U: dof = 1; break;
	 //stokes kernels
  case KNL_STK_F_U: dof = 4; break;
  case KNL_STK_S_U: dof = 3; break;
  case KNL_STK_S_P: dof = 3; break;
  case KNL_STK_D_U: dof = 3; break;
  case KNL_STK_D_P: dof = 3; break;
  case KNL_STK_R_U: dof = 3; break;
  case KNL_STK_R_P: dof = 3; break;
  case KNL_STK_I  : dof = 3; break;
  case KNL_STK_E  : dof = 3; break;
	 //Unsteady Stokes kernels
  case KNL_UNSTK_F_U: dof = 4; break;
  case KNL_UNSTK_S_U: dof = 3; break;
  case KNL_UNSTK_S_P: dof = 3; break;
  case KNL_UNSTK_D_U: dof = 3; break;
  case KNL_UNSTK_D_P: dof = 3; break;
	 //navier kernels:	 //case KNL_NAV_F_U: dof = 3; break; //used for fmm
  case KNL_NAV_S_U: dof = 3; break;
  case KNL_NAV_D_U: dof = 3; break;
  case KNL_NAV_R_U: dof = 3; break;
  case KNL_NAV_I  : dof = 3; break;
  case KNL_NAV_E  : dof = 3; break;
	 //others
  case KNL_SQRTLAP: dof = 1; break;
  case KNL_EXP:     dof = 1; break;
	 //error
  case KNL_ERR:     dof = 0; break;
  }
  return dof;
}

// ---------------------------------------------------------------------- 
uint Kernel3d::trgDOF() const
{
  int dof = 0;
  switch(_kernelType) {
	 //laplace kernels
  case KNL_LAP_S_U: dof = 1; break;
  case KNL_LAP_D_U: dof = 1; break;
  case KNL_LAP_E_U: dof = 3; break;
  case KNL_LAP_I  : dof = 1; break;
	 //Modified Laplacian
  case KNL_MODHEL_S_U: dof = 1; break;
  case KNL_MODHEL_D_U: dof = 1; break;
	 //stokes kernels
  case KNL_STK_F_U: dof = 3; break;
  case KNL_STK_S_U: dof = 3; break;
  case KNL_STK_S_P: dof = 1; break;
  case KNL_STK_D_U: dof = 3; break;
  case KNL_STK_D_P: dof = 1; break;
  case KNL_STK_R_U: dof = 3; break;
  case KNL_STK_R_P: dof = 1; break;
  case KNL_STK_I  : dof = 3; break;
  case KNL_STK_E  : dof = 3; break;
	 //Unsteady Stokes
  case KNL_UNSTK_F_U: dof = 3; break;	
  case KNL_UNSTK_S_U: dof = 3; break;
  case KNL_UNSTK_S_P: dof = 1; break;
  case KNL_UNSTK_D_U: dof = 3; break;
  case KNL_UNSTK_D_P: dof = 1; break;
	 //navier kernels:	 //  case KNL_NAV_F_U: dof = 3; break; //used for fmm
  case KNL_NAV_S_U: dof = 3; break;
  case KNL_NAV_D_U: dof = 3; break;
  case KNL_NAV_R_U: dof = 3; break;
  case KNL_NAV_I  : dof = 3; break;
  case KNL_NAV_E  : dof = 3; break;
	 //others
  case KNL_SQRTLAP: dof = 1; break;
  case KNL_EXP    : dof = 1; break;
	 //error
  case KNL_ERR:     dof = 0; break;
  }
  return dof;
}

// ---------------------------------------------------------------------- 
bool Kernel3d::homogeneous() const
{
  bool ret = false;
  switch(_kernelType) {
	 //laplace kernels
  case KNL_LAP_S_U: ret = true; break;
  case KNL_LAP_E_U: ret = true; break;
	 //Modified Laplace
  case KNL_MODHEL_S_U: ret = false; break;
	 //stokes kernels
  case KNL_STK_F_U: ret = true; break;
  case KNL_STK_S_U: ret = true; break;
  case KNL_UNSTK_F_U: ret = false; break;
	 //navier kernels
  case KNL_NAV_S_U: ret = true; break;
	 //others:
  case KNL_SQRTLAP: ret = true; break;
  case KNL_EXP    : ret = false; break;
  default: assert(0);
  }
  return ret;
}

// ---------------------------------------------------------------------- 
void Kernel3d::homogeneousDeg(vector<double>& degreeVec) const
{
  switch(_kernelType) {
	 //laplace
  case KNL_LAP_S_U: degreeVec.resize(1); degreeVec[0]=1; break;
  case KNL_LAP_E_U: degreeVec.resize(3); degreeVec[0]=1; degreeVec[1]=1; degreeVec[2]=1;  break;
	 //Mod Lap
  case KNL_MODHEL_S_U: degreeVec.resize(0); break;
	 //stokes
  case KNL_STK_F_U: degreeVec.resize(4); degreeVec[0]=1; degreeVec[1]=1; degreeVec[2]=1; degreeVec[3]=2; break;
  case KNL_UNSTK_F_U: degreeVec.resize(0); break;
	 //navier
  case KNL_NAV_S_U: degreeVec.resize(3); degreeVec[0]=1; degreeVec[1]=1; degreeVec[2]=1; break;
	 //others
  case KNL_SQRTLAP: degreeVec.resize(1); degreeVec[0]=0.5; break;
  case KNL_EXP    : degreeVec.resize(0); break;
  default: assert(0);
  }
  return;
}

double Kernel3d::density_unscaling_value(double position_scale){
    switch(_kernelType){
        case KNL_LAP_S_U:
            return pow(position_scale, -1);

        case KNL_LAP_D_U:
            return pow(position_scale, -2);

        case KNL_STK_S_U:
            return pow(position_scale, -1);

        case KNL_STK_D_U:
            return pow(position_scale, -2);
        
        case KNL_STK_S_P:
            return pow(position_scale, -2);

        case KNL_STK_D_P:
            return pow(position_scale, -3);

        case KNL_MODHEL_S_U:
            return pow(position_scale, -1);

        case KNL_MODHEL_D_U:
            return pow(position_scale, -2);

        case KNL_NAV_S_U:
            return pow(position_scale, -1);

        case KNL_NAV_D_U:
            return pow(position_scale, -2);

        default: 
            assert(0);
    }

}

// ---------------------------------------------------------------------- 
DblNumMat Kernel3d::get_identity_matrix(int n){

    DblNumMat id_matrix(n, n); // matrix is initialized to all 0.0
    for(int i = 0; i < n; i++){
        id_matrix(i,i) = 1.0;
    }

    return id_matrix;
}

// ---------------------------------------------------------------------- 
// MJM TODO rename source to target and target to source here
DblNumMat Kernel3d::get_cross_product_matrix(DblNumVec source, DblNumVec target){
    // need to be the same size
    ebiAssert(source.m() == target.m());

    // need to be in R^3 or this doesn't make sense
    ebiAssert( source.m() == 3); //hence target.m() == 3
    DblNumMat cross_product_matrix(3, 3);

    // compute target - source  
    double x = target(0) - source(0);
    double y = target(1) - source(1);
    double z = target(2) - source(2);

    // See KNL_STK_E in Kernel3d.cpp kernel(...) for original
    // Compute this in a loop, but this is never used; might as well do it
    // explicitly
    cross_product_matrix(0, 0) =  0; cross_product_matrix(0, 1) =  z; cross_product_matrix(0, 2) = -y;
    cross_product_matrix(1, 0) = -z; cross_product_matrix(1, 1) =  0; cross_product_matrix(1, 2) =  x;
    cross_product_matrix(2, 0) =  y; cross_product_matrix(2, 1) = -x; cross_product_matrix(2, 2) =  0;

    return cross_product_matrix;
}




// kernel function should be split up
// ---------------------------------------------------------------------- 
int Kernel3d::kernel(const DblNumMat& srcPos, const DblNumMat& srcNor, const DblNumMat& trgPos, DblNumMat& inter)
{
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  ebiAssert(srcPos.m()==dim() && srcNor.m()==dim() && trgPos.m()==dim());
  ebiAssert(srcPos.n()==srcNor.n());
  ebiAssert(srcPos.n()*srcDOF == inter.n());
  if(trgPos.n()*trgDOF != inter.m()){
      cout << "Kernel size mismatch:  trgPos.n(): " << trgPos.n() << ", trgDOF: " << trgDOF << ", inter.m(): " << inter.m() << endl;
  }
  ebiAssert(trgPos.n()*trgDOF == inter.m());
  
  if(       _kernelType==KNL_LAP_S_U) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 /*! One over four pi */
	 double OOFP = 1.0/(4.0*M_PI);
	 for(int i=0; i<trgPos.n(); i++) {
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if( r <_mindif) {
			 inter(i,j) = 0.0;
		  } else {
			 inter(i,j) = OOFP / r;
		  }
		}	
	 }
  } else if(_kernelType==KNL_LAP_D_U) {
	 //---------------------------------
	 /*! One over four pi */
	 double OOFP = 1.0/(4.0*M_PI);
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double nx = srcNor(0,j);
			 double ny = srcNor(1,j);
			 double nz = srcNor(2,j);
			 double rn = x*nx + y*ny + z*nz;
			 double r3 = r2*r;
			 inter(i,j) = - OOFP / r3 * rn;
		  }
		}
  } else if(_kernelType==KNL_LAP_E_U) {
	 //---------------------------------
	 /*! One over four pi */
	 double OOFP = 1.0/(4.0*M_PI);
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 int is = i*3; int js = j;
			 double r3 = r2*r;
             // Changed from (x,x,x) to (x,y,z) in the following lines:
			 inter(is  ,js  ) = - (OOFP * x) / r3;
			 inter(is+1,js  ) = - (OOFP * y) / r3;
			 inter(is+2,js  ) = - (OOFP * z) / r3;
		  }
		}
  } else if(_kernelType==KNL_LAP_I) {
	 //---------------------------------
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++)
		  inter(i,j) = 1.0;
  } else if(_kernelType==KNL_MODHEL_S_U) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 /*! One over four pi */
	 ebiAssert(_coefs.size()>=1);
	 double lambda = _coefs[1]; //std::cerr << lambda << std::endl;
	 double OOFP = 1.0/(4.0*M_PI);
	
	 for(int i=0; i<trgPos.n(); i++) {
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
	
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if( r <_mindif ) {
			 inter(i,j) = 0.0;
		  } else {
			 inter(i,j) = (OOFP / r) * exp(-lambda * r);
		  }
		}	
	 }
  } else if(_kernelType==KNL_MODHEL_D_U) {
	 //---------------------------------
	 double lambda = _coefs[1];
	 double OOFP = 1.0/(4.0*M_PI);
	 for(int i=0; i<trgPos.n(); i++) {
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
          double y = trgPos(1,i) - srcPos(1,j); 
          double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
          double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++) {
				for(int s=0;s<srcDOF;s++) {
				  inter(i*trgDOF+t, j*srcDOF+s) = 0.0;
				}
			 }
		  } else {
			 double nx = srcNor(0,j);
			 double ny = srcNor(1,j);
			 double nz = srcNor(2,j);
			 double rn = x*nx + y*ny + z*nz;
			 double r3 = r2*r;
			 inter(i,j) = - OOFP * (((1.0/r3 + lambda/r2)*exp(-lambda * r)) * rn);
		  }
		}
	 }  
  }  else if(_kernelType==KNL_STK_F_U) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 double OOFP = 1.0/(4.0*M_PI);
	 double OOEP = 1.0/(8.0*M_PI);
	 double oomu = 1.0/mu;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 double G = oomu * OOEP / r;
			 double H = oomu * OOEP / r3;
			 double A = OOFP / r3;
			 int is = i*3;			 int js = j*4;
			 inter(is,   js)   = G + H*x*x; inter(is,   js+1) =     H*x*y; inter(is,   js+2) =     H*x*z; inter(is,   js+3) = A*x;
			 inter(is+1, js)   =     H*y*x; inter(is+1, js+1) = G + H*y*y; inter(is+1, js+2) =     H*y*z; inter(is+1, js+3) = A*y;
			 inter(is+2, js)   =     H*z*x; inter(is+2, js+1) =     H*z*y; inter(is+2, js+2) = G + H*z*z; inter(is+2, js+3) = A*z;
		  }
		}
  } else if(_kernelType==KNL_STK_S_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 /*! One over eight pi */
	 double OOFP = 1.0/(4.0*M_PI);
	 double OOEP = OOFP/2.0;;
	 double oomu = 1.0/mu;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 double G = OOEP / r;
			 double H = OOEP / r3;
			 inter(i*3,   j*3)   = oomu*(G + H*x*x); inter(i*3,   j*3+1) = oomu*(    H*x*y); inter(i*3,   j*3+2) = oomu*(    H*x*z);
			 inter(i*3+1, j*3)   = oomu*(    H*y*x); inter(i*3+1, j*3+1) = oomu*(G + H*y*y); inter(i*3+1, j*3+2) = oomu*(    H*y*z);
			 inter(i*3+2, j*3)   = oomu*(    H*z*x); inter(i*3+2, j*3+1) = oomu*(    H*z*y); inter(i*3+2, j*3+2) = oomu*(G + H*z*z);

			 //inter(i*3,   j*3)   = 1.0; inter(i*3,   j*3+1) = 0.0     ; inter(i*3,   j*3+2) = 0.0;
			 //inter(i*3+1, j*3)   = 0.0     ; inter(i*3+1, j*3+1) = 1.0; inter(i*3+1, j*3+2) = 0.0;
			 //inter(i*3+2, j*3)   = 0.0     ; inter(i*3+2, j*3+1) = 0.0     ; inter(i*3+2, j*3+2) = 1.0;
		  }
		}
  } else if(_kernelType==KNL_STK_S_P) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 /*! One over four pi */
	 double OOFP = 1.0/(4.0*M_PI);
	 
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 inter(i  ,j*3  ) = OOFP*x/r3;			 inter(i  ,j*3+1) = OOFP*y/r3;			 inter(i  ,j*3+2) = OOFP*z/r3;
		  }
		}
  } else if(_kernelType==KNL_STK_D_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 /*! Six over eight pi */
	 double SOEP = 6.0/(8.0*M_PI);
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double nx = srcNor(0,j);
			 double ny = srcNor(1,j);
			 double nz = srcNor(2,j);
			 double rn = x*nx + y*ny + z*nz;
			 double r5 = r2*r2*r;
			 double C = - SOEP / r5;
			 inter(i*3,   j*3)   = C*rn*x*x;
             inter(i*3,   j*3+1) = C*rn*x*y;
             inter(i*3,   j*3+2) = C*rn*x*z;
			 inter(i*3+1, j*3)   = C*rn*y*x;				  inter(i*3+1, j*3+1) = C*rn*y*y;				  inter(i*3+1, j*3+2) = C*rn*y*z;
			 inter(i*3+2, j*3)   = C*rn*z*x;				  inter(i*3+2, j*3+1) = C*rn*z*y;				  inter(i*3+2, j*3+2) = C*rn*z*z;
		  }
		}
  } else if(_kernelType==KNL_STK_D_P) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 /*! One over two pi */
	 double OOTP = 1.0/(2.0*M_PI);
	 double coef = mu*OOTP;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double nx = srcNor(0,j);
			 double ny = srcNor(1,j);
			 double nz = srcNor(2,j);
			 double rn = x*nx + y*ny + z*nz;
			 double r3 = r2*r;
			 double r5 = r3*r2;
			 int is = i;			 int js = j*3;
			 inter(is  ,js  ) = coef*(nx/r3 - 3*rn*x/r5);
			 inter(is  ,js+1) = coef*(ny/r3 - 3*rn*y/r5);
			 inter(is  ,js+2) = coef*(nz/r3 - 3*rn*z/r5);
		  }
		}
  } else if(_kernelType==KNL_STK_R_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 double coef = 1.0/(8.0*M_PI*mu)/r3;
			 inter(i*3,   j*3)   = coef*0;			 inter(i*3,   j*3+1) = coef*z;			 inter(i*3,   j*3+2) = coef*(-y);
			 inter(i*3+1, j*3)   = coef*(-z);		 inter(i*3+1, j*3+1) = coef*0;			 inter(i*3+1, j*3+2) = coef*x;
			 inter(i*3+2, j*3)   = coef*y;			 inter(i*3+2, j*3+1) = coef*(-x);		 inter(i*3+2, j*3+2) = coef*0;
		  }
		}
  } else if(_kernelType==KNL_STK_R_P) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  inter(i,j*3  ) = 0;
		  inter(i,j*3+1) = 0;
		  inter(i,j*3+2) = 0;
		}
  } else if(_kernelType==KNL_STK_I) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  inter(i*3,   j*3  ) = 1;		 inter(i*3,   j*3+1) = 0;		 inter(i*3,   j*3+2) = 0;
		  inter(i*3+1, j*3  ) = 0;		 inter(i*3+1, j*3+1) = 1;		 inter(i*3+1, j*3+2) = 0;
		  inter(i*3+2, j*3  ) = 0;		 inter(i*3+2, j*3+1) = 0;		 inter(i*3+2, j*3+2) = 1;
		}
  } else if(_kernelType==KNL_STK_E) { //levi-civita tensor
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  inter(i*3,   j*3  ) = 0;		 inter(i*3,   j*3+1) = z;		 inter(i*3,   j*3+2) = -y;
		  inter(i*3+1, j*3  ) = -z;  	 inter(i*3+1, j*3+1) = 0;		 inter(i*3+1, j*3+2) = x;
		  inter(i*3+2, j*3  ) = y;		 inter(i*3+2, j*3+1) = -x;		 inter(i*3+2, j*3+2) = 0;
		}
  }
  else if(_kernelType==KNL_UNSTK_F_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 double lambda = _coefs[1];
	 double lambda2 = lambda*lambda;
	 /*! One over eight pi */
	 
	 double OOFP = 1.0/(4.0*M_PI*lambda2);
	 double oomu = 1.0/mu;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r*r2;
			 double explr = exp(-lambda*r);
			 double f1 = -1.0/r2 + (1.0/r2 + lambda/r) * explr;
			 double f2 = 2.0/r3 - (2.0/r3 + 2.0*lambda/r2 + lambda2/r) * explr;

			 double G = OOFP * (-f2 - f1/r);
			 double H = OOFP * (f2/r2 - f1/r3);
			 double A = OOFP / r3;
			 int is = i*3;			 int js = j*4;
			 inter(is,   js)   = G + H*x*x; inter(is,   js+1) =     H*x*y; inter(is,   js+2) =     H*x*z; inter(is,   js+3) = A*x;
			 inter(is+1, js)   =     H*y*x; inter(is+1, js+1) = G + H*y*y; inter(is+1, js+2) =     H*y*z; inter(is+1, js+3) = A*y;
			 inter(is+2, js)   =     H*z*x; inter(is+2, js+1) =     H*z*y; inter(is+2, js+2) = G + H*z*z; inter(is+2, js+3) = A*z;

			 //inter(i*3,   j*3)   = 1.0; inter(i*3,   j*3+1) = 0.0     ; inter(i*3,   j*3+2) = 0.0;
			 //inter(i*3+1, j*3)   = 0.0     ; inter(i*3+1, j*3+1) = 1.0; inter(i*3+1, j*3+2) = 0.0;
			 //inter(i*3+2, j*3)   = 0.0     ; inter(i*3+2, j*3+1) = 0.0     ; inter(i*3+2, j*3+2) = 1.0;
		  }
		}
  }
  else if(_kernelType==KNL_UNSTK_S_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=1);
	 double mu = _coefs[0];
	 double lambda = _coefs[1];
	 double lambda2 = lambda*lambda;
	 /*! One over eight pi */
	 
	 double OOFP = 1.0/(4.0*M_PI*lambda2);
	 double oomu = 1.0/mu;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r*r2;
			 double explr = exp(-lambda*r);
			 double f1 = -1.0/r2 + (1.0/r2 + lambda/r) * explr;
			 double f2 = 2.0/r3 - (2.0/r3 + 2.0*lambda/r2 + lambda2/r) * explr;

			 double G = OOFP * (-f2 - f1/r);
			 double H = OOFP * (f2/r2 - f1/r3);
			 inter(i*3,   j*3)   = oomu*(G + H*x*x); inter(i*3,   j*3+1) = oomu*(    H*x*y); inter(i*3,   j*3+2) = oomu*(    H*x*z);
			 inter(i*3+1, j*3)   = oomu*(    H*y*x); inter(i*3+1, j*3+1) = oomu*(G + H*y*y); inter(i*3+1, j*3+2) = oomu*(    H*y*z);
			 inter(i*3+2, j*3)   = oomu*(    H*z*x); inter(i*3+2, j*3+1) = oomu*(    H*z*y); inter(i*3+2, j*3+2) = oomu*(G + H*z*z);

			 //inter(i*3,   j*3)   = 1.0; inter(i*3,   j*3+1) = 0.0     ; inter(i*3,   j*3+2) = 0.0;
			 //inter(i*3+1, j*3)   = 0.0     ; inter(i*3+1, j*3+1) = 1.0; inter(i*3+1, j*3+2) = 0.0;
			 //inter(i*3+2, j*3)   = 0.0     ; inter(i*3+2, j*3+1) = 0.0     ; inter(i*3+2, j*3+2) = 1.0;
		  }
		}
  }


  else if(_kernelType==KNL_NAV_S_U) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 ebiAssert(_coefs.size()>=2);
	 double mu = _coefs[0];
	 double ve = _coefs[1];
	 double sc1 = (3.0-4.0*ve)/(16.0*M_PI*(1.0-ve));
	 double sc2 = 1.0/(16.0*M_PI*(1.0-ve));
	 double oomu = 1.0/mu;
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 double G = sc1 / r;
			 double H = sc2 / r3;
			 inter(i*3,   j*3)   = oomu*(G + H*x*x);  inter(i*3,   j*3+1) = oomu*(    H*x*y);  inter(i*3,   j*3+2) = oomu*(    H*x*z);
			 inter(i*3+1, j*3)   = oomu*(    H*y*x);  inter(i*3+1, j*3+1) = oomu*(G + H*y*y);  inter(i*3+1, j*3+2) = oomu*(    H*y*z);
			 inter(i*3+2, j*3)   = oomu*(    H*z*x);  inter(i*3+2, j*3+1) = oomu*(    H*z*y);  inter(i*3+2, j*3+2) = oomu*(G + H*z*z);
		  }
		}
  } else if(_kernelType==KNL_NAV_D_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=2);
	 double mu = _coefs[0];
	 double ve = _coefs[1];
	 double dc1 = -(1-2.0*ve)/(8.0*M_PI*(1.0-ve));
	 double dc2 =  (1-2.0*ve)/(8.0*M_PI*(1.0-ve));
	 double dc3 = -3.0/(8.0*M_PI*(1.0-ve));
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double nx = srcNor(0,j);
			 double ny = srcNor(1,j);
			 double nz = srcNor(2,j);
			 double rn = x*nx + y*ny + z*nz;
			 double r3 = r2*r;
			 double r5 = r3*r2;
			 double A = dc1 / r3;
             double B = dc2 / r3;
             double C = dc3 / r5;
			 double rx = x;
             double ry = y;
             double rz = z;
			 //&&&&&&&
			 inter(i*3,   j*3)   = A*(rn+nx*rx) + B*(rx*nx) + C*rn*rx*rx;
			 inter(i*3,   j*3+1) = A*(   nx*ry) + B*(rx*ny) + C*rn*rx*ry;
			 inter(i*3,   j*3+2) = A*(   nx*rz) + B*(rx*nz) + C*rn*rx*rz;
			 inter(i*3+1, j*3)   = A*(   ny*rx) + B*(ry*nx) + C*rn*ry*rx;
			 inter(i*3+1, j*3+1) = A*(rn+ny*ry) + B*(ry*ny) + C*rn*ry*ry;
			 inter(i*3+1, j*3+2) = A*(   ny*rz) + B*(ry*nz) + C*rn*ry*rz;
			 inter(i*3+2, j*3)   = A*(   nz*rx) + B*(rz*nx) + C*rn*rz*rx;
			 inter(i*3+2, j*3+1) = A*(   nz*ry) + B*(rz*ny) + C*rn*rz*ry;
			 inter(i*3+2, j*3+2) = A*(rn+nz*rz) + B*(rz*nz) + C*rn*rz*rz;
		  }
		}
  } else if(_kernelType==KNL_NAV_R_U) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=2);
	 double mu = _coefs[0];
	 double ve = _coefs[1];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 for(int t=0;t<trgDOF;t++)
				for(int s=0;s<srcDOF;s++)
				  { inter(i*trgDOF+t, j*srcDOF+s) = 0.0; }
		  } else {
			 double r3 = r2*r;
			 double coef = 1.0/(8.0*M_PI*mu)/r3;
			 inter(i*3,   j*3)   = coef*0;			 inter(i*3,   j*3+1) = coef*z;			 inter(i*3,   j*3+2) = coef*(-y);
			 inter(i*3+1, j*3)   = coef*(-z);		 inter(i*3+1, j*3+1) = coef*0;			 inter(i*3+1, j*3+2) = coef*x;
			 inter(i*3+2, j*3)   = coef*y;			 inter(i*3+2, j*3+1) = coef*(-x);		 inter(i*3+2, j*3+2) = coef*0;
		  }
		}
  } else if(_kernelType==KNL_NAV_I) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=2);
	 double mu = _coefs[0];
	 double ve = _coefs[1];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  inter(i*3,   j*3  ) = 1;		 inter(i*3,   j*3+1) = 0;		 inter(i*3,   j*3+2) = 0;
		  inter(i*3+1, j*3  ) = 0;		 inter(i*3+1, j*3+1) = 1;		 inter(i*3+1, j*3+2) = 0;
		  inter(i*3+2, j*3  ) = 0;		 inter(i*3+2, j*3+1) = 0;		 inter(i*3+2, j*3+2) = 1;
		}
  } else if(_kernelType==KNL_NAV_E) {
	 //---------------------------------
	 ebiAssert(_coefs.size()>=2);
	 double mu = _coefs[0];
	 double ve = _coefs[1];
	 for(int i=0; i<trgPos.n(); i++)
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  inter(i*3,   j*3  ) = 0;		 inter(i*3,   j*3+1) = z;		 inter(i*3,   j*3+2) = -y;
		  inter(i*3+1, j*3  ) = -z;  	 inter(i*3+1, j*3+1) = 0;		 inter(i*3+1, j*3+2) = x;
		  inter(i*3+2, j*3  ) = y;		 inter(i*3+2, j*3+1) = -x;		 inter(i*3+2, j*3+2) = 0;
		}
  } else if(_kernelType==KNL_SQRTLAP) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 for(int i=0; i<trgPos.n(); i++) {
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 inter(i,j) = 0;
		  } else {
			 inter(i,j) = 1 / sqrt(r);
		  }
		}
	 }
  } else if(_kernelType==KNL_EXP) {
	 //---------------------------------
	 for(int i=0; i<trgPos.n(); i++) {
		for(int j=0; j<srcPos.n(); j++) {
		  double x = trgPos(0,i) - srcPos(0,j);
		  double y = trgPos(1,i) - srcPos(1,j);
		  double z = trgPos(2,i) - srcPos(2,j);
		  double r2 = x*x + y*y + z*z;
		  double r = sqrt(r2);
		  if(r<_mindif) {
			 inter(i,j) = 0;
		  } else {
			 inter(i,j) = exp(-r2);
		  }
		}
	 }
  } else if(_kernelType==KNL_ERR) {
	 //----------------------------------------------------------------------------------
	 //---------------------------------
	 ebiAssert(0);
  }
  return 0;
}

END_EBI_NAMESPACE
