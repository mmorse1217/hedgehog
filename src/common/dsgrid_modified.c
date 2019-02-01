/************************************************************************
*                                                                       *
*                Copyright (C)  2000                                    *
*        University Corporation for Atmospheric Research                *
*                All Rights Reserved                                    *
*                                                                       *
*    The use of this Software is governed by a License Agreement.       *
*                                                                       *
************************************************************************/

#include "ngmath/include/ngmath.h"
#include "ngmath/include/dstypes.h"
#include "ngmath/include/dsproto_mod.h"
#include "ngmath/include/dsuhead_mod.h"
#include "ngmath/include/dsgvars_mod.h"

#include <stdio.h>

#define PIH 1.5707963

/*
 *  Interpolate randomly-spaced 3D input data to a regularly spaced grid.
 */
void c_dsgrid3d_mod(int n, double x[], double y[], double z[], double u[],
                  int nx, int ny, int nz, double xo[], double yo[], 
						 double zo[], double ds_output[], int *ier)
/*
 *    Arguments
 *    ---------
 *        n - The number of input data points.
 *        x - An array of X coordinate values of the input data points.
 *        y - An array of Y coordinate values of the input data points.
 *        z - An array of Z coordinate values of the input data points.
 *        u - The functional value at coordinate (x,y,z).
 *       nx - The dimension of the array xo containing the X coordinate 
 *            values for the output grid.
 *       ny - The dimension of the array yo containing the Y coordinate 
 *            values for the output grid.
 *       nz - The dimension of the array zo containing the z coordinate 
 *            values for the output grid.
 *       xo - The array containing the X coordinate values for the output 
 *            grid (must be monotone increasing, but need not be equally 
 *            spaced.
 *       yo - The array containing the Y coordinate values for the output 
 *            grid (must be monotone increasing, but need not be equally 
 *            spaced.
 *       zo - The array containing the Z coordinate values for the output 
 *            grid (must be monotone increasing, but need not be equally 
 *            spaced.
 *     *ier - An error return value ( = 0 is no error).
 *
 *   Return value
 *   ------------
 *      A pointer to the first element of a linear array that is
 *      laid out as a 3D array (i.e. the last subscript varies fastest)
 *      of dimension: nx by ny by nz.
 */
{
  int       i, j, k;
  static    double perror = 1.;
  double    xc, yc, zc;

  int       iw, it;
  double    normalization_factor, interpolated_value, weight_sum;


  DSpointd3 q;

/*
 *  Get memory and do initialization.
 */
  *ier = 0;

  DSpointd3 ds_input_points[n];
  double ds_weights[n];
  int ds_permutation_vector[n];
  double ds_distances[n];
  
  /* dsinit */
  double xmn, ymn, zmn, xmx, ymx, zmx, tlm;
  int ds_shadowing = 0, ds_error_status = 0, ds_set_max_dist = 1;
  double ds_scale = 1.0, ds_max_dist = 0, ds_epsilon_test;
  if (n < 3) {
	 *ier = ds_error_status;
    return;
  }

  if ((nx <= 0) || (ny <=0) || (nz <= 0)) {
	 *ier = ds_error_status;
    return;
  }

  xmn = xo[0];
  ymn = yo[0];
  zmn = zo[0];
  xmx = xo[0];
  ymx = yo[0];
  zmx = zo[0];
  for (i = 0; i < nx; i++) {
    xmn = MIN(xo[i], xmn);
    xmx = MAX(xo[i], xmx);
  }
  for (i = 0; i < ny; i++) {
    ymn = MIN(yo[i], ymn);
    ymx = MAX(yo[i], ymx);
  }
  for (i = 0; i < nz; i++) {
    zmn = MIN(zo[i], zmn);
    zmx = MAX(zo[i], zmx);
  }
  for (i = 0; i < n; i++) {
    xmn = MIN(x[i], xmn);
    xmx = MAX(x[i], xmx);
    ymn = MIN(y[i], ymn);
    ymx = MAX(y[i], ymx);
    zmn = MIN(z[i], zmn);
    zmx = MAX(z[i], zmx);
  }
  
/*
 *  Find the maximum span of the three coordinates.
 */
  tlm = MAX(MAX((xmx - xmn), ymx - ymn), zmx - zmn);

/*
 *  Set the maximum distance for point inclusion to include all
 *  points, if this has not been specifically set by the user.
 */
  if (ds_set_max_dist == 1){
	 ds_max_dist = 2.*((xmx - xmn) + (ymx - ymn) + (zmx - zmn));
  }
  
/*
 *  Scale and store the input values.
 */
  ds_scale = 1./tlm;
  for (i = 0; i < n; i++) {
    ds_input_points[i].x = x[i] * ds_scale; 
    ds_input_points[i].y = y[i] * ds_scale; 
    ds_input_points[i].z = z[i] * ds_scale; 
  }

  ds_epsilon_test = 1.E-5;

  for (i = 0; i < nx*ny*nz; i++) {
    ds_output[i] = ds_mod_missing_value;
  }

  if (*ier != 0) return(&perror);
    
/*
 *  Loop over the output grid
 */
  for (i = 0; i < nx; i++) {
    xc = xo[i]; 
    for (j = 0; j < ny; j++) {
      yc = yo[j]; 
      for (k = 0; k < nz; k++) {
        zc = zo[k]; 

		  /*
			*  Compute the distances from (xc, yc, zc) to all the input points.
			*/
        q.x = xc * ds_scale;
        q.y = yc * ds_scale;
        q.z = zc * ds_scale;
		  dsdist_mod(n, ds_input_points, q, ds_distances, ds_scale, ds_max_dist);
  
		  /*
			*  Calculate the interpolated value.
			*/
        if (ds_shadowing) {

			 /*
			  *  Calculate interpolated value when shadowing option is on.
			  */
			 /*
			 double svalue(int num_points, double *values, 
			        double xc, double yc, double zc)
			 */
			 int num_points = n;
			 
			 int       ia;

			 for (iw = 0; iw < num_points; iw++) {
				if ( (ds_distances[iw] < ds_epsilon_test) && (ds_distances[iw] >= 0.) ) {
				  for (it = 0; it < num_points; it++) {
					 ds_weights[it] = 0.;
				  }
				  ds_weights[iw] = 1.;
				  goto label_svalue;
				}
			 }

			 for (it = 0; it < num_points; it++) {
				ds_permutation_vector[it] = it;
			 }
			 /*
			  *  Sort a linear array ar in place in ascending order and
			  *  sort a companion integer array in the same order.
			  */
			 double v;
			 int iii, jjj, hhh, iv;
			 n = num_points;
			 for (hhh = 1; hhh <= n/9; hhh = 3*hhh+1);

			 for ( ; hhh > 0; hhh /= 3) {
				for (iii = hhh; iii < n; iii++) {
				  v = ds_distances[iii];
				  iv = ds_permutation_vector[iii];
				  jjj = iii;
				  while (jjj > hhh-1 && ds_distances[jjj-hhh] > v) {
					 ds_distances[jjj] = ds_distances[jjj-hhh];
					 ds_permutation_vector[jjj] = ds_permutation_vector[jjj-hhh];
					 jjj -= hhh;
				  }
				  ds_distances[jjj] = v;
				  ds_permutation_vector[jjj] = iv;
				}
			 }
		  
			 
			 /*
			  * Claculate weights when shadowing is on.
			  */
			 /*int n = num_points;*/
			 int       lpw, lpa;
			 double    minimum_angle, angle, shadow_scale;
			 DSpointd3 p1, p2, p3;
  
			 lpw = ds_permutation_vector[0];
			 if (ds_distances[lpw] > 0.) {
				ds_weights[lpw] = 1./dist_pow_mod(ds_distances[0]);
			 }
			 else {
				ds_weights[lpw] = 0.;
			 }	
  
			 for (iw = 1; iw < n; iw++) {
				minimum_angle = PIH;
				lpw = ds_permutation_vector[n-iw];
				p1.x = ds_input_points[lpw].x;
				p1.y = ds_input_points[lpw].y;
				p1.z = ds_input_points[lpw].z;
				p2.x = q.x;
				p2.y = q.y;
				p2.z = q.z;
				for (ia = 0; ia < n - iw; ia++) {
				  lpa = ds_permutation_vector[ia];
				  p3.x = ds_input_points[lpa].x;
				  p3.y = ds_input_points[lpa].y;
				  p3.z = ds_input_points[lpa].z;
				  angle = dsangd_mod(p1,p2,p3);
				  if (angle < minimum_angle) minimum_angle = angle;
				}
				shadow_scale = tan(0.5*minimum_angle);
	 
				if (ds_distances[n-iw] > 0.) {
				  ds_weights[lpw] = shadow_scale / dist_pow_mod(ds_distances[n-iw]);
				}
				else {
				  ds_weights[lpw] = 0.;
				}
			 }
		  
  
			 for (weight_sum = 0., iw = 0; iw < num_points; iw++) {
				weight_sum += ds_weights[iw];
			 }

			 if ((int)(weight_sum) == 0) {
				printf("svalue weight_sum = 0 error\n");
				exit(-1);
			 }
			 normalization_factor = 1./weight_sum;
			 for (iw = 0; iw < num_points; iw++) {
				ds_weights[iw] *= normalization_factor;
			 }

		  label_svalue:
			 for (interpolated_value = 0., iw = 0; iw < num_points; iw++) {
				interpolated_value += u[iw] * ds_weights[iw];
			 }
			 
			 
          ds_output[nz*ny*i + nz*j + k] = interpolated_value;
        }
        else {
			 /*ivalue */
			 /*
			  *  Calculate interpolated value when the shadowing option is
			  *  turned off.
			  */
			 int num_points = n;
			 for (iw = 0; iw < num_points; iw++) {
				if ( (ds_distances[iw] < ds_epsilon_test) && (ds_distances[iw] >= 0.) ) {
				  for (it = 0; it < num_points; it++) {
					 ds_weights[it] = 0.;
				  }
				  ds_weights[iw] = 1.;
				  goto label_ivalue;
				}
			 }

			 n = num_points;
				for (it = 0; it < n; it++) {
				  if (ds_distances[it] > 0.) {
					 ds_weights[it] = 1./ dist_pow_mod(ds_distances[it]);
				  }
				  else {
					 ds_weights[it] = 0.;
				  }
				}
			 

			 
			 for (weight_sum = 0., iw = 0; iw < num_points; iw++) {
				weight_sum += ds_weights[iw];
			 }

			 if ((int)(weight_sum) != 0) {
				normalization_factor = 1./weight_sum;
				for (iw = 0; iw < num_points; iw++) {
				  ds_weights[iw] *= normalization_factor;
				}
			 }
			 else {

			 }

		  label_ivalue:
			 for (interpolated_value = 0., iw = 0; iw < num_points; iw++) {
				interpolated_value += u[iw] * ds_weights[iw];
			 }

			 if ((int)(weight_sum) == 0) { interpolated_value = 0.0; }
			 
			 ds_output[nz*ny*i + nz*j + k] = interpolated_value;
		  }
      }
    }
  }
  
}


/*
 *  Calculate interpolated value when the shadowing option is
 *  turned off.
 */
double ivalue_mod(int num_points, double *values)
{
  int       iw, it;
  double    normalization_factor, interpolated_value, weight_sum;

  for (iw = 0; iw < num_points; iw++) {
    if ( (ds_distances[iw] < ds_epsilon_test) && (ds_distances[iw] >= 0.) ) {
      for (it = 0; it < num_points; it++) {
        ds_weights[it] = 0.;
      }
      ds_weights[iw] = 1.;
      goto label_1;
    }
  }

  dweights(num_points);
  
  for (weight_sum = 0., iw = 0; iw < num_points; iw++) {
    weight_sum += ds_weights[iw];
  }

  if ((int)(weight_sum) == 0) {
    DSErrorHnd(14, "ivalue", stderr, "\n");    
    return(ds_mod_missing_value);
  }
  normalization_factor = 1./weight_sum;
  for (iw = 0; iw < num_points; iw++) {
    ds_weights[iw] *= normalization_factor;
  }
 
label_1:
  for (interpolated_value = 0., iw = 0; iw < num_points; iw++) {
    interpolated_value += values[iw] * ds_weights[iw];
  }

  return(interpolated_value);
  
}

/*
 *  Calculate interpolated value when shadowing option is on.
 */
double svalue_mod(int num_points, double *values, 
                   double xc, double yc, double zc)
{
  int       iw, it;
  double    normalization_factor, interpolated_value, weight_sum;

  for (iw = 0; iw < num_points; iw++) {
    if ( (ds_distances[iw] < ds_epsilon_test) && (ds_distances[iw] >= 0.) ) {
      for (it = 0; it < num_points; it++) {
        ds_weights[it] = 0.;
      }
      ds_weights[iw] = 1.;
      goto label_1;
    }
  }

  for (it = 0; it < num_points; it++) {
    ds_permutation_vector[it] = it;
  }
  dssortd(num_points, ds_distances, ds_permutation_vector);
  sweights(num_points, xc, yc, zc);
  
  for (weight_sum = 0., iw = 0; iw < num_points; iw++) {
    weight_sum += ds_weights[iw];
  }

  if ((int)(weight_sum) == 0) {
    DSErrorHnd(14, "svalue", stderr, "\n");
    return(ds_mod_missing_value);
  }
  normalization_factor = 1./weight_sum;
  for (iw = 0; iw < num_points; iw++) {
    ds_weights[iw] *= normalization_factor;
  }
 
label_1:
  for (interpolated_value = 0., iw = 0; iw < num_points; iw++) {
    interpolated_value += values[iw] * ds_weights[iw];
  }
  return(interpolated_value);
  
}

/*
 *  3D interpolation in point mode.
 */
/*
void c_dspnt3d_mod(int n, double xi[], double yi[], double zi[], double ui[],
               int m, double xo[], double yo[], double zo[], double uo[],
						 int *ier);

void c_dspnt3d_mod(int n, double xi[], double yi[], double zi[], double ui[],
               int m, double xo[], double yo[], double zo[], double uo[],
               int *ier)
{
  printf("NEED TO FIX IF WANT TO USE\n"); exit(0);
  
  int    i;
  double xt[1], yt[1], zt[1];

  for (i = 0; i < m; i++) {
    xt[0] = xo[i];
    yt[0] = yo[i];
    zt[0] = zo[i];
    uo[i] = *c_dsgrid3d(n, xi, yi, zi, ui, 1, 1, 1, xt, yt, zt, ier);
  }
}
*/

/*
 *  Calculate powers of the distances.
 */
double dist_pow_mod(double dist)
{
  double dtmp;

  if ((int)(ds_mod_exponent) == 3) {
    dtmp = dist*dist*dist;
  }
  else if ((int)(ds_mod_exponent) == 1) {
    dtmp = dist;
  }
  else if ((int)(ds_mod_exponent*2.0) == 1) {
    dtmp = sqrt(dist);
  }
  else if ((int)(ds_mod_exponent) == 2) {
    dtmp = dist*dist;
  }
  else if ((int)(ds_mod_exponent) == 4) {
    dtmp = dist*dist;
    dtmp = dtmp*dtmp;
  }
  else if ((int)(ds_mod_exponent) == 5) {
    dtmp = dist*dist*dist*dist*dist;
  }
  else if ((int)(ds_mod_exponent) == 6) {
    dtmp = dist*dist*dist;
    dtmp = dtmp*dtmp;
  }
  else if ((int)(ds_mod_exponent) == 7) {
    dtmp = dist*dist*dist*dist*dist*dist*dist;
  }
  else if ((int)(ds_mod_exponent) == 8) {
    dtmp = dist*dist;
    dtmp = dtmp*dtmp;
    dtmp = dtmp*dtmp;
  }
  else if ((int)(ds_mod_exponent) == 9) {
    dtmp = dist*dist*dist;
    dtmp = dtmp*dtmp*dtmp;
  }
  else if ((int)(ds_mod_exponent) == 10) {
    dtmp = dist*dist*dist*dist*dist;
    dtmp = dtmp*dtmp;
  }
  else {
    dtmp = pow(dist, ds_mod_exponent);
  }
  if (dtmp < (double) 1.E-30) {
    return ( (double) 1.E-30);
  }
  else {
    return (dtmp);
  }
}

/*
 *  Given three points in three space a, b, c, find the angle
 *  between the vector from b to a and the vector from b to c.
 */
double dsangd_mod (DSpointd3 a, DSpointd3 b, DSpointd3 c)
{
  DSpointd3 vector_1,vector_2;
  double    cosd;

  vector_1.x = a.x - b.x;
  vector_1.y = a.y - b.y;
  vector_1.z = a.z - b.z;
  
  vector_2.x = c.x - b.x;
  vector_2.y = c.y - b.y;
  vector_2.z = c.z - b.z;

  cosd = dotd_mod(vector_1,vector_2)/(magd_mod(vector_1)*magd_mod(vector_2));
  if (cosd >  1.) cosd =  1.; 
  if (cosd < -1.) cosd = -1.; 

  return(acos(cosd));
}

/*
 *  Find the magnitude of a vector.
 */
double magd_mod(DSpointd3 p)
{
  return(sqrt(p.x*p.x + p.y*p.y + p.z*p.z));
}

/*
 *  Find the dot product of two vectors.
 */
double dotd_mod(DSpointd3 p, DSpointd3 q)
{
  return(p.x*q.x + p.y*q.y + p.z*q.z);
}


/*
 *  Computes distances between a set of input points and
 *  a given point.  The distances are returned in a sorted array.
 */
void dsdist_mod(int n, DSpointd3 p[], DSpointd3 q, double *sdist, double ds_scale, double ds_max_dist)
/*
 *        n - The number of points in the p array.
 *        p - An array of n points.
 *        q - An individual point.
 *    sdist - A pointer to an array of n distances to point q.
 */
{
  DSpointd3 dp;
  int       i;
  double    dtmp;

  for (i = 0; i < n; i++) {
    dp.x = p[i].x - q.x; 
    dp.y = p[i].y - q.y; 
    dp.z = p[i].z - q.z; 
    dtmp = magd_mod(dp);
    if (dtmp <= ds_max_dist * ds_scale) {
      sdist[i] = dtmp;
    }
    else {
      sdist[i] = -1.;
    }
  }
}

