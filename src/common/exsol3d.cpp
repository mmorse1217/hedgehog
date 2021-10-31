#include "vecmatop.hpp"
#include "vec3t.hpp"
#include "exsol3d.hpp"

BEGIN_EBI_NAMESPACE

#undef __FUNCT__
#define __FUNCT__ "Exsol3d::tdof"
int Exsol3d::tdof(int qt)
{
  int dof = 0;
  if(       _et==KNL_LAP_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 1; break;
	 case QNT_RHS: dof = 1; break;
	 case QNT_MAX_U: dof = 1; break;
	 case QNT_MAX_RHS: dof = 1; break;
	 }
  } else if(       _et==KNL_MODHEL_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 1; break;
	 case QNT_RHS: dof = 1; break;
	 case QNT_MAX_U: dof = 1; break;
	 case QNT_MAX_RHS: dof = 1; break;
	 }
  } else if(_et==KNL_STK_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 3; break;
	 case QNT_P: dof = 1; break;
	 case QNT_RHS: dof = 3; break;
	 case QNT_MAX_U: dof = 3; break;
	 case QNT_MAX_RHS: dof = 3; break;
	 }
  } else if(_et==KNL_NAV_S_U) {
	 switch(qt) {
	 case QNT_U: dof = 3; break;
	 case QNT_RHS: dof = 3; break;
	 }
  } else {
	 cerr<<"error"<<endl;
  }
  return dof;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Exsol3d::quantity"
int Exsol3d::quantity(int qt, const DblNumMat& trgpos, DblNumVec& trgval)
{
  ebiFunctionBegin;
  if (_ct != CHS_EMPTY) { ebiAssert(trgpos.n()*tdof(qt)==trgval.m()); }
  double L = 250.0;
  if(       _et==KNL_LAP_S_U) {
	 if(       _ct==CHS_LAP_ZRO) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		}
	 } else if(       _ct==CHS_LAP_CST) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 1;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		}
	 } else if(_ct==CHS_LAP_XYZ) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = 1 + x*y*z;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		}
	 } else if(_ct==CHS_LAP_X2) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);
			 trgval(i) = x*x;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = -2.0;
		  }
		}
	 } else if(_ct==CHS_LAP_ESQ) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 //trgval(i) = exp(-sqrt(2.0)*x*y*z)*sin(1.0*M_PI*(x*x+y*y+z*z))+ (1.0/396.0)*pow((x+y+z),12.0);
			 trgval(i) = exp(sqrt(2.0)*PI*x)*sin(PI*(y+z)) + 1.0/6.0*(x*x*x + y*y*y + z*z*z);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);

			 //double T = sqrt(2.0);
			 //double P = M_PI;
			 //double I = 1.0;
			 //double r2 = x*x + y*y + z*z;
			 //double xyz = x*y*z;
			 //double tpr2 = I*P*r2;
			 //double etp = exp(-T*xyz);

			 //double val = T*T*y*y*z*z*etp*sin(tpr2)-12.0*T*y*z*etp*cos(tpr2)*I*P*x-4.0*I*I*etp*sin(tpr2)*P*P*x*x+6.0*I*etp*cos(tpr2)*P+T*T*x*x*z*z*etp*sin(tpr2)-4.0*I*I*etp*sin(tpr2)*P*P*y*y+T*T*x*x*y*y*etp*sin(tpr2)-4.0*I*I*etp*sin(tpr2)*P*P*z*z;
			 // val = val + pow((x+y+z),10.0);
			 //trgval(i) = -val;
			 trgval(i) = -(x+y+z);
			 //trgval(i) = -x;
		  }
		}
	 } else if(_ct==CHS_LAP_X2Y2Z2) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);
			 double y = trgpos(1,i);
			 double z = trgpos(2,i);
			 trgval(i) = x*x + y*y + z*z;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = -6.0;
		  }
		}
	 } else if(_ct==CHS_LAP_ANA){
		//------------------
		if(       qt==QNT_U) {
		  double RV = 0.9;
		  double OOFP = 0.25*(1.0/M_PI);
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i) - RV;
			 double y = trgpos(1,i);
			 double z = trgpos(2,i);
			 double r = sqrt(x*x + y*y + z*z);
			 trgval(i) = 0.0;
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 x = trgpos(0,i) + RV;
			 y = trgpos(1,i);
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 x = trgpos(0,i);
			 y = trgpos(1,i) - RV;
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 x = trgpos(0,i);
			 y = trgpos(1,i) + RV;
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 x = trgpos(0,i);
			 y = trgpos(1,i);
			 z = trgpos(2,i) - RV;
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 x = trgpos(0,i);
			 y = trgpos(1,i);
			 z = trgpos(2,i) + RV;
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 /*
			 x = trgpos(0,i) - RV/sqrt(2.0);
			 y = trgpos(1,i) - RV/sqrt(2.0);
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i) - RV/sqrt(2.0);
			 y = trgpos(1,i) + RV/sqrt(2.0);
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i) + RV/sqrt(2.0);
			 y = trgpos(1,i) - RV/sqrt(2.0);
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i) + RV/sqrt(2.0);
			 y = trgpos(1,i) + RV/sqrt(2.0);
			 z = trgpos(2,i);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 */
		  

			 /*
			 x = trgpos(0,i) - RV/sqrt(2.0);
			 y = trgpos(1,i);
			 z = trgpos(2,i) - RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i) - RV/sqrt(2.0);
			 y = trgpos(1,i);
			 z = trgpos(2,i) + RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 
			 x = trgpos(0,i) + RV/sqrt(2.0);
			 y = trgpos(1,i);
			 z = trgpos(2,i) - RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i) + RV/sqrt(2.0);
			 y = trgpos(1,i);
			 z = trgpos(2,i) + RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }



			 x = trgpos(0,i);
			 y = trgpos(1,i) - RV/sqrt(2.0);
			 z = trgpos(2,i) - RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i);
			 y = trgpos(1,i) - RV/sqrt(2.0);
			 z = trgpos(2,i) + RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 
			 x = trgpos(0,i);
			 y = trgpos(1,i) + RV/sqrt(2.0);
			 z = trgpos(2,i) - RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }

			 x = trgpos(0,i);
			 y = trgpos(1,i) + RV/sqrt(2.0);
			 z = trgpos(2,i) + RV/sqrt(2.0);
			 r = sqrt(x*x + y*y + z*z);
			 if (r > 1e-12){
				trgval(i) += -OOFP*(1.0/r);
			 }
			 */

			 
		 			 
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0.0;
		  }
		}
	 } else if(_ct==CHS_HARPER) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(i) = (sin(PI*x)*sin(PI*y)*sin(PI*z));
		  }
		} else if(qt==QNT_RHS) {
		  double cnst = PI;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(i) = 3.0*cnst*cnst*(sin(cnst*x)*sin(cnst*y)*sin(cnst*z));
		  }
		}
	 } else if (_ct== CHS_LAP_ATAN1 || _ct== CHS_LAP_ATAN2 || _ct== CHS_LAP_ATAN5 || _ct== CHS_LAP_ATAN10 || _ct== CHS_LAP_ATAN20 || _ct == CHS_LAP_ATAN30 || _ct == CHS_LAP_ATAN40 || _ct == CHS_LAP_ATAN80|| _ct == CHS_LAP_TOR_ATAN10){
		double scale = 0.0;
		if (_ct==CHS_LAP_ATAN1) scale = 1.0;
		else if (_ct==CHS_LAP_ATAN2) scale = 2.0;
		else if (_ct==CHS_LAP_ATAN5) scale = 5.0;
		else if (_ct==CHS_LAP_ATAN10 || _ct==CHS_LAP_TOR_ATAN10) scale = 10.0;
		else if (_ct==CHS_LAP_ATAN20) scale = 20.0;
		else if (_ct==CHS_LAP_ATAN30) scale = 30.0;
		else if (_ct==CHS_LAP_ATAN40) scale = 40.0;
		else if (_ct==CHS_LAP_ATAN80) scale = 80.0;
		else { ebiAssert(0); }

		double R = 0.8;
		if (_ct == CHS_LAP_TOR_ATAN10) { R = 0.48; }
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -(1.0/M_PI)*(atan(scale*(R*R - r2)) + M_PI/2.0);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 double r2mR2 = r2 - R*R;
			 double fac = (scale*scale)*(r2mR2*r2mR2) + 1.0;
			 trgval(i) =  (((8.0/M_PI)*(scale*scale*scale)*(r2)*(r2mR2))/(fac*fac) - (6.0/M_PI)*(scale/fac));
		  }
		}
	 }	else if(_ct==CHS_LAP_EXTHAR) {
		double ll = 2.5;
		double c = 0.0;
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 z = z - c;
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = exp(-ll*r2);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 z = z - c;
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -exp(-ll*r2)*((4.0*ll*ll*r2) - 6.0*ll);
			 //std::cout << x << " " << y << " " << z << " " << trgval(i) << endl;
		  }
		}
	 } else if(       _ct==CHS_LAP_FREE_SPACE) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = exp(-L*r2);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += exp(-L*r2);
				  }
				}
			 }
			 //cerr << i << " " << trgval(i) << endl;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
				  }
				}
			 }	
			 //trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z;
			 //trgval(i) = 1.0 + x + y + z;
			 //if (abs(trgval(i)) < 10e-64) trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = exp(0.0);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += exp(-L*r2);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		}else if(qt== QNT_MAX_RHS){
		  trgval(0) = -exp(0)*((0.0) - 6.0*L);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		}
	 } else if(       _ct==CHS_LAP_SPH) {
		double RSPH = coefs()[0];
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 if (sqrt(r2) < RSPH) trgval(i) = 1.0/(6.0) * (RSPH*RSPH - r2) + RSPH*RSPH/(3.0);
			 else trgval(i) = RSPH*RSPH*RSPH/(3.0*sqrt(r2));
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 if (sqrt(r2) < RSPH) trgval(i) = 1.0;
			 else trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = abs(1.0/(6.0) * (RSPH*RSPH) + RSPH*RSPH/(3.0));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_COL) {
		Point3 C1(-0.3125, -0.0625, 0.3125);
		Point3 C2(-0.0625, 0.3125, -0.3125);
		Point3 C3(0.3125, -0.3125, -0.0625);
		double RVAL = coefs()[0];
		double MVAL = coefs()[1];
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double al = 4.0*M_PI*MVAL;

			 double x1 = x-C1(0); double y1 = y-C1(1); double z1 = z-C1(2);
			 double r_1= sqrt(x1*x1 + y1*y1 + z1*z1);	 
			 x1 = x-C2(0); y1 = y-C2(1); z1 = z-C2(2);
			 double r_2= sqrt(x1*x1 + y1*y1 + z1*z1);
			 x1 = x-C3(0); y1 = y-C3(1); z1 = z-C3(2);
			 double r_3= sqrt(x1*x1 + y1*y1 + z1*z1);

			 double al2 = al*al;
			 double al3 = al*al2;
			 double al4 = al*al3;
			 double al5 = al*al4;
			 double al6 = al*al5;
			 double al7 = al*al6;
			 double tmp = 0.0;
			 double RS1 = r_1/RVAL;
			 double RS2 = r_2/RVAL;
			 double RS3 = r_3/RVAL;
			 double val = 0.0;
			 
			 if (RS1 < 1.0) {
				double RS = RS1;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);
				
				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.0/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));
		
				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }
			 else if (RS2 < 1.0) {
				double RS = RS2;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);

				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.00/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));

				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }	 
			 else if (RS3 < 1.0) {
				double RS = RS3;
				tmp = pow(RS,6.0)/84.0 - pow(RS,5.0)/30.0 + pow(RS,4.0)/40.0;
				tmp = tmp + 60.0/(al6) - 9.0/(al4) - 1.0/120.0 + 120.0/(al6*RS);
				
				tmp = tmp + (-120.0/(al6*RS) - 9.0/(al4) + 300.00/(al6))*cos(al*RS);
				tmp = tmp + (36.0*RS/(al4) + pow(RS,2.0)/(2.0*al2))*cos(al*RS);
				tmp = tmp + (-30.0*pow(RS,2.0)/(al4) - pow(RS,3.0)/(al2))*cos(al*RS);
				tmp = tmp + (pow(RS,4.0)/(2.0*al2))*(cos(al*RS));

				tmp = tmp + (12.0/(al5*RS) - 360.0/(al7*RS) - 96.0/(al5) + 120.0*RS/(al5))*(sin(al*RS));
				tmp = tmp + (-3.0*RS/(al3) + 8.0*pow(RS,2.0)/al3 - 5.0*pow(RS,3.0)/al3)*(sin(al*RS));

				//Add in RS2
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
		
				val = -tmp/RVAL;
			 }
			 else{
				tmp = (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS1);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS2);
				tmp = tmp + (-1.0/210.0 - 12.0/(al4) + 360.0/(al6))/(RS3);
				val = -tmp/RVAL;
			 }
			 trgval(i) = val;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double al = 4.0*M_PI*MVAL;

			 double x1 = x-C1(0); double y1 = y-C1(1); double z1 = z-C1(2);
			 double r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_1= sqrt(r2)/RVAL;
			 	 
			 x1 = x-C2(0); y1 = y-C2(1); z1 = z-C2(2);
			 r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_2= sqrt(r2)/RVAL;
			 
			 x1 = x-C3(0); y1 = y-C3(1); z1 = z-C3(2);
			 r2 = (x1*x1 + y1*y1 + z1*z1);
			 double r_3= sqrt(r2)/RVAL;

			 double val = 0.0;
			 if (r_1 < 1.0) {
				double r_12 = r_1*r_1;
				val = ((r_1 - r_12)*sin(al*r_1/2.0));
				val = val*val;
			 }
			 else if (r_2 < 1.0) {
				double r_22 = r_2*r_2;
				val = ((r_2 - r_22)*sin(al*r_2/2.0));
				val = val*val;
			 }
			 else if (r_3 < 1.0) {
				double r_32 = r_3*r_3;
				val = ((r_3 - r_32)*sin(al*r_3/2.0));
				val = val*val;
			 }
			 else val = 0.0;
			 val = val/(RVAL*RVAL*RVAL);
			 trgval(i) = val;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = abs((-1.0/120.0 - (6.0/(pow(4.0*M_PI*MVAL,4.0))))/RVAL + (-1.0/105.0 - 24.0/(pow(4.0*M_PI*MVAL, 4.0)) + 720.0/(pow(4.0*M_PI*MVAL, 6.0)))/(sqrt(0.59375)));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_PER) {
		//------------------
		double CON = coefs()[0];
		double PK = coefs()[1];
		double adjval  = 0.0;
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*sin(M_PI*x*CON)*sin(M_PI*y*CON)*sin(M_PI*z*CON);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -PK*(3.0*(CON*CON)*M_PI*M_PI + adjval)*sin(M_PI*x*CON)*sin(M_PI*y*CON)*sin(M_PI*z*CON);
		  }
		  //#warning need to change max vals
		}  else if(qt== QNT_MAX_U){
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
	 } else if(       _ct==CHS_LAP_DIR) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -sin(M_PI*(1.0+x)/2.0)*sin(M_PI*(1.0+y)/2.0)*sin(M_PI*(1.0+z)/2.0);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = -(3.0*(0.5*0.5)*M_PI*M_PI)*sin(M_PI*(1.0+x)/2.0)*sin(M_PI*(1.0+y)/2.0)*sin(M_PI*(1.0+z)/2.0);
		  }
		  //#warning need to change max vals
		}  else if(qt== QNT_MAX_U){
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = 1.0;
		}
   } else if(       _ct==CHS_LAP_POLY_TEST) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0.0;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z + x*x*x*x;
		  }
		}
	 } else if (_ct== CHS_LAP_ATAN1 || _ct== CHS_LAP_ATAN2 || _ct== CHS_LAP_ATAN5 || _ct== CHS_LAP_ATAN10 || _ct== CHS_LAP_ATAN20 || _ct == CHS_LAP_ATAN30 || _ct == CHS_LAP_ATAN40 || _ct == CHS_LAP_ATAN80|| _ct == CHS_LAP_TOR_ATAN10){
		double scale = 0.0;
		if (_ct==CHS_LAP_ATAN1) scale = 1.0;
		else if (_ct==CHS_LAP_ATAN2) scale = 2.0;
		else if (_ct==CHS_LAP_ATAN5) scale = 5.0;
		else if (_ct==CHS_LAP_ATAN10 || _ct==CHS_LAP_TOR_ATAN10) scale = 10.0;
		else if (_ct==CHS_LAP_ATAN20) scale = 20.0;
		else if (_ct==CHS_LAP_ATAN30) scale = 30.0;
		else if (_ct==CHS_LAP_ATAN40) scale = 40.0;
		else if (_ct==CHS_LAP_ATAN80) scale = 80.0;
		else { ebiAssert(0); }

		double R = 0.8;
		if (_ct == CHS_LAP_TOR_ATAN10) { R = 0.48; }
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -(1.0/M_PI)*(atan(scale*(R*R - r2)) + M_PI/2.0);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 double r2mR2 = r2 - R*R;
			 double fac = (scale*scale)*(r2mR2*r2mR2) + 1.0;
			 trgval(i) =  (((8.0/M_PI)*(scale*scale*scale)*(r2)*(r2mR2))/(fac*fac) - (6.0/M_PI)*(scale/fac));
		  }
		}
	 }	 else if (_ct == CHS_LAP_HAR) {
		double p = M_PI;
		double q = -10.0;
		if (qt == QNT_U){
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0.0;
			 for (int c = -1; c <= 1; c = c+2){
				for (int dd = 0; dd < 3; dd++){
				  double zs = (dd == 0 ? 0.2 : (dd == 1 ? 0.6 : 1.0));
				  double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i) + c*zs;			 //double c0 = sqrt(2.0);
				  double r2 = x*x + y*y + z*z;
				  //double r = sqrt(r2);	
				  //trgval(i) += exp(q*r2)* + sin(p*x)*sin(p*y)*sin(p*z);
				  trgval(i) += sin(p*x)*sin(p*y)*sin(p*z) + exp(q*r2);
				}
			 }
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) += sin(p*x)*sin(p*y)*sin(p*z) + exp(q*r2);
		  }
		} else if (qt == QNT_RHS){
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0.0;
			 for (int c = -1; c <= 1; c = c+2){
				for (int dd = 0; dd < 3; dd++){
				  double zs = (dd == 0 ? 0.2 : (dd == 1 ? 0.6 : 1.0));
				  double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i) + c*zs;			 //double c0 = sqrt(2.0);
				  //double r2 = x*x + y*y + z*z;
				  //double r = sqrt(r2);	
				  //trgval(i) += -((4.0*q*q*r2) + 6.0*q)*exp(q*r2) + (3.0*p*p)*sin(p*x)*sin(p*y)*sin(p*z);
				  trgval(i) += 3*p*p*sin(p*x)*sin(p*y)*sin(p*z) - 4*q*q*x*x*exp(q*(x*x + y*y + z*z)) - 4*q*q*y*y*exp(q*(x*x + y*y + z*z)) - 4*q*q*z*z*exp(q*(x*x + y*y + z*z)) - 6*q*exp(q*(x*x + y*y + z*z));
				}
			 }
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 //double r2 = x*x + y*y + z*z;
			 trgval(i) += 3*p*p*sin(p*x)*sin(p*y)*sin(p*z) - 4*q*q*x*x*exp(q*(x*x + y*y + z*z)) - 4*q*q*y*y*exp(q*(x*x + y*y + z*z)) - 4*q*q*z*z*exp(q*(x*x + y*y + z*z)) - 6*q*exp(q*(x*x + y*y + z*z));
		  }
		}
	 } else if ( _ct == CHS_EMPTY){
		//Do nothing
	 }
	 else {
		cerr << _ct << endl;
		//------------------
		ebiAssert(0);
	 }
  } else if(       _et==KNL_MODHEL_S_U) {
	 double lambda = coefs()[1];
	 double alpha = lambda*lambda;
	 if(       _ct==CHS_MODHEL_CST) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		
			 //double y = trgpos(1,i);		
			 //double z = trgpos(2,i);
			 trgval(i) = exp(lambda*x);
			 //double r2 = x*x + y*y + z*z;
			 //double r = sqrt(r2);
			 //if (r > 1e-12)
			 //trgval(i) = exp(-lambda * r)/r;
			 //else
			 //trgval(i) = 0.0;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		}
	 } else if(_ct==CHS_MODHEL_XYZ) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);      double y = trgpos(1,i);    double z = trgpos(2,i);
			 trgval(i) = 1 + x*y*z;
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);      double y = trgpos(1,i);    double z = trgpos(2,i);
			 trgval(i) = alpha*(1 + x*y*z);
		  }
		}
	 }	 else if (_ct == CHS_MODHEL_HAR) {
		double p = M_PI;
		double q = -10.0;
		if (qt == QNT_U){
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) += sin(p*x)*sin(p*y)*sin(p*z) + exp(q*r2);
		  }
		} else if (qt == QNT_RHS){
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) += 3*p*p*sin(p*x)*sin(p*y)*sin(p*z) - 4*q*q*x*x*exp(q*(x*x + y*y + z*z)) - 4*q*q*y*y*exp(q*(x*x + y*y + z*z)) - 4*q*q*z*z*exp(q*(x*x + y*y + z*z)) - 6*q*exp(q*(x*x + y*y + z*z));
			 trgval(i) += alpha*(sin(p*x)*sin(p*y)*sin(p*z) + exp(q*r2));
		  }
		}
  } else if(_ct==CHS_MODHEL_ESQ) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);      double y = trgpos(1,i);    double z = trgpos(2,i);        //double c0 = sqrt(2.0);
			 trgval(i) = exp(sqrt(2.0)*PI*(x))*sin(PI*(y+z)) + 1.0/6.0*(x*x*x+y*y*y+z*z*z);
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);      double y = trgpos(1,i);    double z = trgpos(2,i);        //double c0 = sqrt(2.0);
			 trgval(i) = -(x+y+z) + alpha*(exp(sqrt(2.0)*PI*(x))*sin(PI*(y+z)) + 1.0/6.0*(x*x*x+y*y*y+z*z*z));
		  }
		}
  } else if(       _ct==CHS_MODHEL_FREE_SPACE) {
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = exp(-L*r2);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += exp(-L*r2);
				  }
				}
			 }
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 trgval(i) = -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 trgval(i) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
				  }
				}
			 }
			 //trgval(i) = 1 + x + y + z + x*x + x*y + x*z + y*y + y*z + z*z + x*x*x + x*x*y + x*x*z + x*y*y + x*y*z + x*z*z + y*y*y + y*y*z + y*z*z + z*z*z + x*x*x*x + x*x*x*x*x + x*x*x*x*x*x + x*x*x*x*x*x*x + x*y*z*x*y*z*x*y*z + sin(7.0*x);
			 //trgval(i) = 1.0;
			 //if (abs(trgval(i)) < 10e-64) trgval(i) = 0.0;
		  }
		} else if(qt== QNT_MAX_U){
		  trgval(0) = exp(0.0);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += exp(-L*r2);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		} else if(qt== QNT_MAX_RHS){
		  trgval(0) = -exp(0)*((0.0) - 6.0*L - alpha);
		  for (int a = -1; a < 2; a++) {
			 for (int b = -1; b < 2; b++) {
				for (int c = -1; c < 2; c++) {
				  double x1 = a*(0.075);	
				  double y1 = b*(0.075);
				  double z1 = c*(0.075);
				  double r2 = x1*x1 + y1*y1 + z1*z1;
				  trgval(0) += -exp(-r2*L)*((4.0*L*L*r2) - 6.0*L - alpha);
				}
			 }
		  }
		  trgval(0) = abs(trgval(0));
		}
	 }
	 else {
		iA(0);
	 }
  } else if(_et==KNL_STK_S_U) {
	 //----------------------------------------------
	 if(       _ct==CHS_STK_FREE_SPACE) {
		double l = L/2.0;
		double l2 = l*l;
		//------------------
		if(       qt==QNT_U) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 int trgdof = tdof(qt);
			 for (int d = 0; d < trgdof; d++){
				if (d == 0){      trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(z - y); }
				else if (d == 1){ trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(x - z); }
				else if (d == 2){ trgval(i*trgdof + d) = 2.0*l*exp(-l*r2)*(y - x); }
			 }
			 
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 for (int d = 0; d < trgdof; d++){
						if (d == 0){      trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(z1 - y1); }
						else if (d == 1){ trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(x1 - z1); }
						else if (d == 2){ trgval(i*trgdof + d) += 2.0*l*exp(-l*r2)*(y1 - x1); }
					 }
				  }
				}
			 }
			 
		  }
		} else if(qt==QNT_RHS) {
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);
			 double r2 = x*x + y*y + z*z;
			 int trgdof = tdof(qt);
			 for (int d = 0; d < trgdof; d++){
				if (d == 0){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(z - y)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(z - y);
				}
				else if (d == 1){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(x - z)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(x - z);
				}
				else if (d == 2){
				  trgval(i*trgdof + d) = -4.0*l2*exp(-l*r2)*(y - x)*(2.0*l*r2 - 5.0);
				  //trgval(i*trgdof + d) = -(y - x);
				}
			 }
			 for (int a = -1; a < 2; a++) {
				for (int b = -1; b < 2; b++) {
				  for (int c = -1; c < 2; c++) {
					 double x1 = x +  a*(0.075);	
					 double y1 = y +  b*(0.075);
					 double z1 = z +  c*(0.075);
					 r2 = x1*x1 + y1*y1 + z1*z1;
					 for (int d = 0; d < trgdof; d++){
						if (d == 0){		
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(z1 - y1)*(2.0*l*r2 - 5.0);
						}
						else if (d == 1){
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(x1 - z1)*(2.0*l*r2 - 5.0);
						}
						else if (d == 2){
						  trgval(i*trgdof + d) += -4.0*l2*exp(-l*r2)*(y1 - x1)*(2.0*l*r2 - 5.0);
						}
					 }
				  }
				}
			 }
		  }
		} else if(qt== QNT_MAX_U){
		  //#warning need to fix
		  trgval(0) = 1.0;
		} else if(qt== QNT_MAX_RHS){
		  //#warning need to fix
		  trgval(0) = 1.0;
		} else {
		  //------------------
		  ebiAssert(0);
		}
	 } else if(       _ct==CHS_STK_CST) {
		//------------------
		if(       qt==QNT_U) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 1; //towards z direction
		  }
		} else if(qt==QNT_P) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 } else if(_ct==CHS_STK_STK) {
		//------------------
		if(       qt==QNT_U) {
		  DblNumMat srcpos(3,1); //one point at origin
		  Kernel3d knl(KNL_STK_S_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_P) {
		  DblNumMat srcpos(3,1); //one point 
		  Kernel3d knl(KNL_STK_S_P, coefs());
		  DblNumMat inter(1*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 } else if(_ct==CHS_STK_ROT) {
		//------------------
		if(       qt==QNT_U) {
		  DblNumMat srcpos(3,1); //one point 
		  Kernel3d knl(KNL_STK_R_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_P) {
		  DblNumMat srcpos(3,1); //one point 
		  Kernel3d knl(KNL_STK_R_P, coefs());
		  DblNumMat inter(1*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
 		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 } else if(_ct==CHS_STK_PB2) {
		//------------------
		if(       qt==QNT_U) {
		  int numsrc = 8;
		  DblNumMat srcpos(3, numsrc);
		  int cnt = 0;
		  for(int i=0; i<2; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) {
			 srcpos(0,cnt) = -0.4 + i;			 srcpos(1,cnt) = -0.4 + j;			 srcpos(2,cnt) = -0.4 + k;
			 cnt++;
		  }
		  Kernel3d knl(KNL_STK_S_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3*numsrc);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3*numsrc);		  setvalue(srcden, 1.0); //all 1 everywhere
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_P) {
		  int numsrc = 8;
		  DblNumMat srcpos(3, numsrc);
		  int cnt = 0;
		  for(int i=0; i<2; i++) for(int j=0; j<2; j++) for(int k=0; k<2; k++) {
			 srcpos(0,cnt) = -0.4 + i;			 srcpos(1,cnt) = -0.4 + j;			 srcpos(2,cnt) = -0.4 + k;
			 cnt++;
		  }
		  Kernel3d knl(KNL_STK_S_P, coefs());
		  DblNumMat inter(1*trgpos.n(), 3*numsrc);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3*numsrc);
		  setvalue(srcden, 1.0);
          //all 1 everywhere
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;
             trgval(c++) = 0;
             trgval(c++) = 0;

		  }
		}
	 } else if(_ct==CHS_STK_RBR) {
		//------------------
		if(       qt==QNT_U) {
		  int c=0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);
             double y = trgpos(1,i);
             double z = trgpos(2,i);
			 //double c0 = sqrt(2.0);

			 trgval(c++) = z-y;
			 trgval(c++) = x-z;
			 trgval(c++) = y-x;

		  }
		} else if(qt==QNT_P) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;
             trgval(c++) = 0;
             trgval(c++) = 0;

		  }
		}
	 }
	 else if(_ct==CHS_STK_HAR) {
		L = 1.0; //resetting L here
		double p = M_PI/4.0;
		double q = -10.0;
		//double L2 = L*L;
		//------------------
		if(       qt==QNT_U) {
		  int cnt=0;
		  for(int i=0; i<trgpos.n(); i++) {
			 for (int a = -1; a <= 1; a = a+2){
				for (int b = -1; b <= 1; b = b+2){
				  for (int c = -1; c <= 1; c = c+2){
					 double x = trgpos(0,i) + a*0.25;		double y = trgpos(1,i) + b*0.25;		double z = trgpos(2,i) + c*0.25;			 //double c0 = sqrt(2.0);
					 double r2 = x*x + y*y + z*z;
					 //double r = sqrt(r2);	
					 trgval(3*cnt + 0) += -exp(q*r2)*(y - z) - sin(p*(y - z));
					 trgval(3*cnt + 1) +=  exp(q*r2)*(x - z) + sin(p*(x - z));
					 trgval(3*cnt + 2) += -exp(q*r2)*(x - y) - sin(p*(x - y));
				  }
				}
			 }
			 cnt++;
		  }
		} else if(qt==QNT_P) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		} else if(qt==QNT_RHS) {  
		  int cnt = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 for (int a = -1; a <= 1; a = a+2){
				for (int b = -1; b <= 1; b = b+2){
				  for (int c = -1; c <= 1; c = c+2){
					 double x = trgpos(0,i) + a*0.25;		double y = trgpos(1,i) + b*0.25;		double z = trgpos(2,i) + c*0.25;			 //double c0 = sqrt(2.0);
					 double r2 = x*x + y*y + z*z;
					 //double r = sqrt(r2);
					 trgval(3*cnt + 0) += 4*q*y*exp(q*(r2)) - 2*(p*p)*sin(p*(y - z)) - 4*q*z*exp(q*(r2)) + 6*q*exp(q*(r2))*(y - z) + 4*(q*q)*(x*x)*exp(q*(r2))*(y - z) + 4*(q*q)*(y*y)*exp(q*(r2))*(y - z) + 4*(q*q)*(z*z)*exp(q*(r2))*(y - z);
					 trgval(3*cnt + 1) += 2*(p*p)*sin(p*(x - z)) - 4*q*x*exp(q*(r2)) + 4*q*z*exp(q*(r2)) - 6*q*exp(q*(r2))*(x - z) - 4*(q*q)*(x*x)*exp(q*(r2))*(x - z) - 4*(q*q)*(y*y)*exp(q*(r2))*(x - z) - 4*(q*q)*(z*z)*exp(q*(r2))*(x - z);
					 trgval(3*cnt + 2) += 4*q*x*exp(q*(r2)) - 2*(p*p)*sin(p*(x - y)) - 4*q*y*exp(q*(r2)) + 6*q*exp(q*(r2))*(x - y) + 4*(q*q)*(x*x)*exp(q*(r2))*(x - y) + 4*(q*q)*(y*y)*exp(q*(r2))*(x - y) + 4*(q*q)*(z*z)*exp(q*(r2))*(x - y);
					 
					 //trgval(c++) =  2.0*q*exp(q*r2)*(y - z)*(2.0*r2 + 5.0) - sin(p*(z))*(p*p)*(z-y) - 2.0*p*cos(p*z);
					 //trgval(c++) = -2.0*q*exp(q*r2)*(x - z)*(2.0*r2 + 5.0) + sin(p*(x))*(p*p)*(x-z) - 2.0*p*cos(p*x);
					 //trgval(c++) =  2.0*q*exp(q*r2)*(x - y)*(2.0*r2 + 5.0) - sin(p*(y))*(p*p)*(y-x) - 2.0*p*cos(p*y);
				  }
				}
			 }
			 cnt++;
		  }
		}
	 }
	 else if(_ct==CHS_STK_ESQ) {
		double p = 1.0;
		//------------------
		if(       qt==QNT_U) {
		  int c=0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(c++) = -sin(p*(z))*(y-z);
			 trgval(c++) =  sin(p*(x))*(x-z);
			 trgval(c++) = -sin(p*(y))*(x-y);
		  }
		} else if(qt==QNT_P) {
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(i) = 0;
		  }
		} else if(qt==QNT_RHS) {  
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(c++) =  -sin(p*(z))*(p*p)*(y-z) - 2.0*p*cos(p*z);
			 trgval(c++) =   sin(p*(x))*(p*p)*(x-z) - 2.0*p*cos(p*x);
			 trgval(c++) =  -sin(p*(y))*(p*p)*(y-x) - 2.0*p*cos(p*y);
		  }
		}
	 } else if(_ct==CHS_STK_SIN) {
		double p = 1.0;
		iA(0); // NOT DIV FREE
		//------------------
		if(       qt==QNT_U) {
		  int c=0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(c++) = -p*sin(p*y - p*z)*sin(p*x);
			 trgval(c++) =  p*sin(p*x - p*z)*sin(p*y);
			 trgval(c++) = -p*sin(p*x - p*y)*sin(p*z);
			 /*
			 if (fabs(y-z) < 1e-12) { trgval(c++) = 1.0; }
			 else { trgval(c++) = -sin(M_PI*(p*y - p*z))/(M_PI*(p*y-p*z)); }

			 if (fabs(x-z) < 1e-12) { trgval(c++) = 1.0; }
			 else { trgval(c++) =  sin(M_PI*(p*x - p*z))/(M_PI*(p*x-p*z)); }

			 if (fabs(x-y) < 1e-12) { trgval(c++) = 1.0; }
			 else { trgval(c++) = -sin(M_PI*(p*x - p*y))/(M_PI*(p*x-p*y)); }
			 */
		  }
		} else if(qt==QNT_RHS) {
		  int c=0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 trgval(c++) =  -3.0*(p*p*p)*sin(p*y - p*z)*sin(p*x);
			 trgval(c++) =   3.0*(p*p*p)*sin(p*x - p*z)*sin(p*y);
			 trgval(c++) =  -3.0*(p*p*p)*sin(p*x - p*y)*sin(p*z);
			 /*
			 if (fabs(y-z) < 1e-12) { trgval(c++) = 0.0; }
			 else { trgval(c++) = (4*(p*p)*sin(M_PI*(p*y - p*z)))/(M_PI*(p*y - p*z)*(p*y - p*z)*(p*y - p*z)) - (4*(p*p)*cos(M_PI*(p*y - p*z)))/((p*y - p*z)*(p*y - p*z)) - (2*M_PI*(p*p)*sin(M_PI*(p*y - p*z)))/(p*y - p*z); }

			 if (fabs(x-z) < 1e-12) { trgval(c++) = 0.0; }
			 else { trgval(c++) = (4*(p*p)*cos(M_PI*(p*x - p*z)))/((p*x - p*z)*(p*x - p*z)) - (4*(p*p)*sin(M_PI*(p*x - p*z)))/(M_PI*(p*x - p*z)*(p*x - p*z)*(p*x - p*z)) + (2*M_PI*(p*p)*sin(M_PI*(p*x - p*z)))/(p*x - p*z); }

			 if (fabs(x-y) < 1e-12) { trgval(c++) = 0.0; }
			 else { trgval(c++) = (4*(p*p)*sin(M_PI*(p*x - p*y)))/(M_PI*(p*x - p*y)*(p*x - p*y)*(p*x - p*y)) - (4*(p*p)*cos(M_PI*(p*x - p*y)))/((p*x - p*y)*(p*x - p*y)) - (2*M_PI*(p*p)*sin(M_PI*(p*x - p*y)))/(p*x - p*y); }
			 */
		  }	
		}
	 } else if(_ct==CHS_STK_SPH) {
		//------------------
		if(       qt==QNT_U) {
		  int c=0;
		  double R = 0.8;
		  double U0 = 1.0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 double r = sqrt(r2);              double iRr = R/r;
			 double Fr = 1.0 - 1.5*iRr + 0.5*iRr*iRr*iRr;
			 double Gr = 1.0 - 0.75*iRr - 0.25*iRr*iRr*iRr;
			 
			 double theta = acos(z/r);
			 
			 if (fabs(x) < 1e-16) { x = 1e-16; }
			 double phi = atan(y/x);
			 double Vr = U0*Fr*cos(theta);     double Vth = -U0*Gr*sin(theta);
			 //cout<<"Vr = "<<Vr<<"Vth = "<<Vth<<endl;
			 
			 trgval(c++) = cos(phi)*(Vr*sin(theta) + Vth*cos(theta));
			 trgval(c++) = sin(phi)*(Vr*sin(theta) + Vth*cos(theta));
			 trgval(c++) = Vr*cos(theta) - Vth*sin(theta) - U0;
		  }
		} else if(qt==QNT_P) {
		  
		  double R = 0.8;
		  double mu = 1.0;
		  double U0 = 1.0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 double r = sqrt(r2);              double iRr = R/r;             double theta = acos(z/r);
			 double Pr = -(1.5*mu*U0/R)*iRr*iRr;
          
			 trgval(i) =  Pr*cos(theta);
			 
		  }
		}     else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 }  else if (_ct== CHS_STK_ATAN1 || _ct== CHS_STK_ATAN2 || _ct== CHS_STK_ATAN5 || _ct== CHS_STK_ATAN10 || _ct== CHS_STK_ATAN20 || _ct == CHS_STK_ATAN30 || _ct == CHS_STK_ATAN40 || _ct == CHS_STK_ATAN80){
		double scale = 0.0;
		if (_ct==CHS_STK_ATAN1) scale = 1.0;
		else if (_ct==CHS_STK_ATAN2) scale = 2.0;
		else if (_ct==CHS_STK_ATAN5) scale = 5.0;
		else if (_ct==CHS_STK_ATAN10) scale = 10.0;
		else if (_ct==CHS_STK_ATAN20) scale = 20.0;
		else if (_ct==CHS_STK_ATAN30) scale = 30.0;
		else if (_ct==CHS_STK_ATAN40) scale = 40.0;
		else if (_ct==CHS_STK_ATAN80) scale = 80.0;
		else { ebiAssert(0); }

		double R = 0.8;
		//------------------
		if(       qt==QNT_U) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 trgval(c++) = -(1.0/M_PI)*(atan(scale*(R*R - r2)) + M_PI/2.0)*(z - y);
			 trgval(c++) = -(1.0/M_PI)*(atan(scale*(R*R - r2)) + M_PI/2.0)*(x - z);
			 trgval(c++) = -(1.0/M_PI)*(atan(scale*(R*R - r2)) + M_PI/2.0)*(y - x);
		  }
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 double x = trgpos(0,i);		double y = trgpos(1,i);		double z = trgpos(2,i);			 //double c0 = sqrt(2.0);
			 double r2 = x*x + y*y + z*z;
			 double r2mR2 = r2 - R*R;
			 double fac = (scale*scale)*(r2mR2*r2mR2) + 1.0;
			 double s3 = scale*scale*scale;
			 trgval(c++) = -((10.0*(z - y))/M_PI)*(scale/fac) + (8*(s3)*(r2)*(z - y)*(r2mR2))/(M_PI*(fac)*(fac));
			 trgval(c++) = -((10.0*(x - z))/M_PI)*(scale/fac) + (8*(s3)*(r2)*(x - z)*(r2mR2))/(M_PI*(fac)*(fac));
			 trgval(c++) = -((10.0*(y - x))/M_PI)*(scale/fac) + (8*(s3)*(r2)*(y - x)*(r2mR2))/(M_PI*(fac)*(fac));
		  }
		}
	 } else {
		//------------------
		ebiAssert(0);
	 }
	 
  } else if(_et==KNL_NAV_S_U) {
	 //----------------------------------------------
	 //------------------

	 if(       _ct==CHS_NAV_CST) {
		//------------------
		if(       qt==QNT_U) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 1; //towards z direction
		  }
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 } else if(_ct==CHS_NAV_NAV) {
		//------------------
		if(       qt==QNT_U) {
		  DblNumMat srcpos(3,1); //one point at origin
		  Kernel3d knl(KNL_NAV_S_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
	 } else if(_ct==CHS_NAV_ROT) {
		//------------------
		if(       qt==QNT_U) {
		  DblNumMat srcpos(3,1); //one point 
		  Kernel3d knl(KNL_NAV_R_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3); srcden(0) = 1; //strength (1,0,0)
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
 		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;		trgval(c++) = 0;		trgval(c++) = 0;
		  }
		}
		
     } else if(_ct==CHS_NAV_PB3) {
		//------------------
		if(       qt==QNT_U) {
		  int numsrc = 8;
		  DblNumMat srcpos(3, numsrc);
		  int cnt = 0;
		  for(int i=0; i<2; i++)
              for(int j=0; j<2; j++)
                  for(int k=0; k<2; k++) {
			        srcpos(0,cnt) = -0.4 + i;
                    srcpos(1,cnt) = -0.4 + j;
                    srcpos(2,cnt) = -0.4 + k;
			 cnt++;
		  }
		  Kernel3d knl(KNL_NAV_S_U, coefs());
		  DblNumMat inter(3*trgpos.n(), 3*numsrc);
		  iC( knl.kernel(srcpos, srcpos, trgpos, inter) );
		  DblNumVec srcden(3*numsrc);
          setvalue(srcden, 1.0); //all 1 everywhere
		  iC( dgemv(1.0, inter, srcden, 0.0, trgval) );
		} else if(qt==QNT_RHS) {
		  int c = 0;
		  for(int i=0; i<trgpos.n(); i++) {
			 trgval(c++) = 0;
             trgval(c++) = 0;
             trgval(c++) = 0;

		  }
		}
     else{
         cerr << "KernelNotYetImplementedError: " << _et << endl;
        ebiAssert(0);
     }
  }
  }
  ebiFunctionReturn(0);
}

END_EBI_NAMESPACE
