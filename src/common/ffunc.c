#include <math.h>
#include <assert.h>
#include "ffunc.h"

void func(long *ndim, double *pnt, long *Xnumfun, double *funvls)
{
  /*printf("KT = %d %d\n", glb_kt, glb_k);*/
  double OOFP = 1.0/(4.0*M_PI);
  double OOEP = 1.0/(8.0*M_PI);
  
  double x = pnt[0];
  double y = pnt[1];
  double z = pnt[2];

  double rx = (x - glb_xtarg);
  double ry = (y - glb_ytarg);
  double rz = (z - glb_ztarg);
  double r2 = rx*rx + ry*ry + rz*rz;
  double r = sqrt(r2);
  double oor = 1.0/r;
  /*double oor3 = 1.0/(r*r*r);*/
    
  if (r > 1e-15){
	 int ii,jj;
	 double tmpvals[glb_nk];
	 
	 for (ii=0; ii < glb_nk; ii++){
		tmpvals[ii] = 0.0;
	 }

	 if (glb_k <= 6) {
		tmpvals[0] = 1.0;
		int cnt = 1;
		for (ii=1; ii <= glb_k - 1; ii++){
		  for (jj = (1 + ii*(ii+1)*(ii+2)/6); jj <= ((ii+1)*(ii+2)*(ii+3)/6) - (ii+1); jj++){
			 tmpvals[jj-1] = x*tmpvals[jj-(ii*(ii+1)/2)-1];
			 cnt++;
		  }
		  for (jj = ((ii+1)*(ii+2)*(ii+3)/6 - ii); jj <= ((ii+1)*(ii+2)*(ii+3)/6 - 1); jj++){
			 tmpvals[jj-1] = y*tmpvals[jj-(ii*(ii+3)/2)-1];
			 cnt++;
		  }
		  jj = (ii+1)*(ii+2)*(ii+3)/6;
		  tmpvals[jj-1] = z*tmpvals[(ii*(ii+1)*(ii+2)/6)-1];
		  cnt++;
		}
		assert(cnt == glb_nk);
	 }
	 else { /* Use Cheb Polys */
		double T[3*glb_k];
		T[0] = 1.0; T[glb_k] = 1.0;   T[2*glb_k] = 1.0;
		T[1] = x;   T[glb_k + 1] = y; T[2*glb_k + 1] = z;
		for (ii = 2; ii < glb_k; ii++){
		  for (jj = 0; jj < 3; jj++){
			 int tmp_int = jj*glb_k + ii;
			 T[tmp_int] = 2.0*T[jj*glb_k + 1]*T[tmp_int-1] - T[tmp_int-2];
		  }
		}
		int cnt = 0;
		int JJ;
		for (JJ=0; JJ < glb_k; JJ++){
		  for (ii=JJ; ii >= 0; ii--){
			 for (jj = (JJ-ii); jj >= 0; jj--){
				int kk = (JJ-(ii+jj));
				tmpvals[cnt] = T[ii]*T[glb_k+jj]*T[2*glb_k + kk];
				cnt++;
			 }
		  }
		}
		assert(cnt == glb_nk);
	 }
	 /* Make basic function values with proper dof */
	 int s, t;
	 int std = glb_sdof*glb_tdof;
	 for (ii = 0; ii < glb_nk; ii++){
		for (s = 0; s < glb_sdof; s++){
		  for (t = 0; t < glb_tdof; t++){
			 funvls[ii*std + s*glb_tdof + t] = tmpvals[ii];
		  }
		}
	 }

	 /* Everything below here is kernel specific */
	 
	 if (glb_kt == 111){ /* Laplacian */
		assert(std == 1);
		double stdvals[std];
		stdvals[0] = oor * OOFP;

		for (ii = 0; ii < glb_nk; ii++){
		  for (jj = 0; jj < std; jj++){
			 funvls[ii*std + jj] *= stdvals[jj];
		  }
		}
	 }
	 else if (glb_kt == 211) { /* Helmholtz */
		assert(std == 1);
		double stdvals[std];
		stdvals[0] = oor * OOFP * exp(-(double)(glb_lambda)*r);
		for (ii = 0; ii < glb_nk; ii++){
		  for (jj = 0; jj < std; jj++){
			 funvls[ii*std + jj] *= stdvals[jj];
		  }
		}			
	 }
	 else if (glb_kt == 311){
		double r3 = r2*r;
		double G = OOEP / r;
		double H = OOEP / r3;
		double rx2 = rx*rx;  double rxry = rx*ry; double rxrz = rx*rz;
		                     double ry2 = ry*ry;  double ryrz = ry*rz;
									                     double rz2 = rz*rz;

		double Hrxry = H*rxry;     double Hrxrz = H*rxrz;     double Hryrz = H*ryrz;
		double GHrx2 = G + H*rx2;  double GHry2 = G + H*ry2;  double GHrz2 = G + H*rz2;

		assert(std == 9);
		double stdvals[std];
		stdvals[0] = GHrx2;   stdvals[1] = Hrxry;   stdvals[2] = Hrxrz;
		stdvals[3] = Hrxry;   stdvals[4] = GHry2;   stdvals[5] = Hryrz;
		stdvals[6] = Hrxrz;   stdvals[7] = Hryrz;   stdvals[8] = GHrz2;
  
		for (ii = 0; ii < glb_nk; ii++){
		  for (jj = 0; jj < std; jj++){
			 funvls[ii*std + jj] *= stdvals[jj];
		  }
		}
	 }
	 else {
		/*printf("Either you selected wrong kernel or not yet implemented\n");*/
		assert(0);
	 }
  }
  else {
	 int ii;
	 for (ii = 0; ii < glb_nk*glb_sdof*glb_tdof; ii++){
		funvls[ii] = 0.0;
	 }
  }
}
	
