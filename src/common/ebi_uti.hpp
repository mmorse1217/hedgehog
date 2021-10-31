/*! \file */
#ifndef _EBI_UTI_HPP_
#define _EBI_UTI_HPP_

#include "ebi_namespace.hpp"
#include "ebi_include.hpp"

BEGIN_EBI_NAMESPACE

using std::string;

/***********************************************************************************************/
/***********************************************************************************************/
enum DimType {  //DIM0=0,
  DIM2=2,
  DIM3=3
};

#define SCL_EPS DBL_EPSILON
#define SCL_MAX DBL_MAX
#define SCL_MIN DBL_MIN

#define MAXDIM 3
#define PI M_PI

//complex number
typedef std::complex<double> cpx;

int radmult_spacing(const double spac, double& rad);
int bis3dspacing2bdsurfspacing(const double bisspac, double& bdspac);
double error_estimate(int n, double target_accuracy);
/*
 * q = number of quadrature points 1d
 * delta = distance from point to boundary
 * phi = density value on-surface at closest point
 * max_jacobian = max jacobian value over closest patch
 */
double error_estimate2(int q, double delta, double phi, double max_jacobian);
double error_estimate3(int q, double delta, double phi, double max_jacobian);
double error_estimate_fit(double q, double delta, double phi, double max_jacobian);
double near_zone_approx_size(int q,  double phi, double max_jacobian,
        double target_accuracy);
double near_zone_approx_size_fit(int q,  double phi, double max_jacobian,
        double target_accuracy);
bool is_quad_accurate_at_target(int q, double delta, double phi, 
        double max_jacobian, double target_accuracy);
bool is_quad_accurate_at_target_fit(int q, double delta, double phi, 
        double max_jacobian, double target_accuracy);

bool AE2C(double A, double B, int maxUlps=1);

// ---------------------------------------------------------------------- 
inline int pow2(int l) { if(l >=0) { return (1<<l); } else { return (int)(pow(2.0,(double)(l))); } }
// ---------------------------------------------------------------------- 
int file2string(MPI_Comm comm, const char* infile, string& str);
int string2file(MPI_Comm comm, string& str, const char* outfile);


END_EBI_NAMESPACE

#endif



