#include "common/ebi.hpp"
#include "utils.hpp"
BEGIN_EBI_NAMESPACE

using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::endl;

#define MAXSTRLEN 10000
double error_estimate(int n, double target_accuracy){
    double temp= 1.;
    return -1./(2*n)*log(target_accuracy)*temp;
}
double working_error_estimate(int q, double delta, double phi, double max_jacobian){
    //double phi = 1.;
    //double max_jacobian = 2.;
    //double r= dist_to_patch; // TODO fix
    //double delta = max_jacobian*r; // TODO fix
    //double delta = dist_to_patch; // TODO fix
    double r = 1./(max_jacobian)*delta;
    double temp = 1./log(q);
    if(max_jacobian > 1.)
        temp *= pow(max_jacobian,3);
    else
        temp *= 1./pow(max_jacobian,3);
    //temp *= 3;
    //return .5*phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-2*q*r);
    return phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-4*q*r)*temp; // works for flat patch and both sides of curved patch
    //return phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-4*q*r)*temp*1./(pow(max_jacobian,3)); 
}
double error_constant(int q, double phi, double max_jacobian){
    //double temp = 1./log(q)*(4*M_PI*M_PI);
    double temp = 1./log(q);
    double fudge_factor= 2*100.;
    temp *= fudge_factor;
    //if(max_jacobian > 1.)
        //temp *= pow(max_jacobian,2);
    //else
        //temp *= 1./pow(max_jacobian,2);
    return phi*sqrt(M_PI)/tgamma(4)*2*q*temp*(max_jacobian);

}
double error_estimate2(int q, double delta, double phi, double max_jacobian){
    //double phi = 1.;
    //double max_jacobian = 2.;
    //double r= dist_to_patch; // TODO fix
    //double delta = max_jacobian*r; // TODO fix
    //double delta = dist_to_patch; // TODO fix
    //temp *= 3;
    //return .5*phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-2*q*r);
    
    double r = delta/(max_jacobian);
    double constant = error_constant(q,phi,max_jacobian);
    return constant*delta*exp(-4*q*r); // works for flat patch and both sides of curved patch
    //return phi*max_jacobian*sqrt(M_PI)/tgamma(4)*sqrt(2*q*2*q)*delta*exp(-4*q*r)*temp*1./(pow(max_jacobian,3)); 
}
double error_estimate3(int q, double delta, double phi, double max_jacobian){
    delta = fabs(delta);
    cout << "jac:" <<  max_jacobian << endl;
    //double r = 2*delta/max_jacobian;
    //return 16./(4*M_PI*tgamma(3.5)*pow(r, 4))*r*sqrt((M_PI*r))*exp(-r*(q-1));//*1./pow(max_jacobian,3);
    // vaguely working
    double r = 4*delta/max_jacobian;
    return 16./(4*M_PI*tgamma(3.5)*pow(r, 4))*r*sqrt((M_PI*r))*exp(-r*(q-1));//*1./pow(max_jacobian,3);
    
    
    //return 8./(tgamma(3.5)*pow(2*2*delta, 3.5))*delta*1./(2*q)*exp(-2*r*q)*max_jacobian;//*1./pow(max_jacobian,3);
    //return 8./(tgamma(3.5)*pow(2*2*delta, 3.5))*(delta/(2.*q))*exp(-2*r*q)*max_jacobian;//*1./pow(max_jacobian,3);
}
double error_estimate_fit(double q, double delta, double phi, double max_jacobian){
    
    double H = max_jacobian;
    //vector<double> c = {  -0.7194,  -0.6945,  0.0071};
    //vector<double> c = {-0.82207252,0.21929627,-0.51105122,0.08497819,0.00208908,0.08778897  };
     //double s =  -0.81417011 ;
    /*vector<double> c = { 0.07562064, -0.04479058, 0.00443401, 0.02537857, 
        0.88204874, -0.04963572, 0.45872447, -0.02664137};
    double s = 0.117006644786;
    vector<double> values = { 
        H, H*H*H, q*delta, log10(q), log10(q)*H*H, log10(delta), log10(delta)*H*H
    };*/
   //vector<double> c = { -0.06206106,-0.03847751,-0.82671431,0.55786321  };
    ////double s =  0.64899458 ;
    //vector<double> c = { 0.43882073,0.89574204,-0.07052429,0.01043488  };
     //double s =  0.30365267 ;
     //vector<double> c = { -0.43954992,-0.89789702,0.02052794,0.01246691  };
      //double s =  -0.30196938 ;
      //vector<double> c = { 0.43669588,0.89421950,-0.02620255,0.09271910,0.01956943,-0.00135144  };
       //double s =  0.37761694 ;
       //vector<double> c = { 0.45774468,0.88653077,0.02107063,0.04271104,0.04756030,-0.00167631  };
        //double s =  0.15809456 ;
        vector<double> c = { 0.45805428,0.88738207,0.02135794,0.04775302,-0.00168032  };
        double s =  0.13446854 ;
    vector<double> values = { 
        //sqrt(q*delta), q, log(delta), log10(q), sqrt(q*delta)*H*H
        //q*delta, log10(q), log10(delta), delta*H,q*H
        q*delta, log10(q), delta*H,q*H
    };
    double error_exponent = 0.;
    for (int i = 1; i < c.size(); i++) {
        double error_term =c[i]*values[i-1]; 
        error_exponent += error_term;
    }
    error_exponent += s; 
    error_exponent *= -1.;
    error_exponent /= c[0];
 
    return pow(10., error_exponent );
}

double near_zone_approx_size_fit(int q, double phi, double max_jacobian,
        double target_accuracy){
    vector<double> c = { 0.45805428,0.88738207,0.02135794,0.04775302,-0.00168032  };
    double s =  0.13446854 ;
    double alpha = c[0];
    double beta = c[1];
    double lambda = c[2];
    double eta = c[3];
    double kappa= c[4];

    double H = max_jacobian;
    //cout << 1./(beta*q + eta*H) <<", " << alpha*log10(target_accuracy)  << ", " <<  lambda*log10(q) << ", " <<  kappa*q*H << endl;
    return -1./(beta*q + eta*H)*(alpha*log10(target_accuracy) - lambda*log10(q) - kappa*q*H);

}
bool is_quad_accurate_at_target_fit(int q, double delta, double phi, 
        double max_jacobian, double target_accuracy){
    //double constant = error_constant(q, phi, max_jacobian);
    return error_estimate_fit(q,delta, phi, max_jacobian) < target_accuracy;
    //return delta > near_zone_approx_size(q, phi, max_jacobian, target_accuracy);
}

double near_zone_approx_size(int q, double phi, double max_jacobian,
        double target_accuracy){
    double constant = error_constant(q, phi, max_jacobian);
    //return -max_jacobian*q/4.*log(constant*target_accuracy*1e4);
    return max_jacobian/(q*4.)*log(constant/target_accuracy);
}
bool is_quad_accurate_at_target(int q, double delta, double phi, 
        double max_jacobian, double target_accuracy){
    //double constant = error_constant(q, phi, max_jacobian);
    return error_estimate2(q,delta, phi, max_jacobian) < target_accuracy;
    //return delta > near_zone_approx_size(q, phi, max_jacobian, target_accuracy);
}

//Tale surface discretization and return radmult - this is optimum setup
#undef __FUNCT__
#define __FUNCT__ "radmult_spacing"
int radmult_spacing(const double spc, double& rad) 
{	
  ebiFunctionBegin;

  // roughly recovers the values below
  //rad = 2./sqrt(spc);
  double constant = Options::get_double_from_petsc_opts("-pou_radius_constant");
  //recovers the values Lexing used in the 2006 paper
  rad = constant/sqrt(spc);
  //rad = 1.3/sqrt(spc);
  //rad = 1.6/sqrt(spc);
    // MJM TODO REMOVE!!!!!!!!!!!
  /*
  if (spc <= 0.2) { rad = 4; }
  if (spc <= 0.1) { rad = 7; }
  if (spc <= 0.05) { rad = 9; }
  if (spc <= 0.025) { rad = 13; }
  if (spc <= 0.0125) { rad = 17; }
  if (spc <= 0.00625) { rad = 25; }
  if (spc <= 0.003125) { rad = 35; }
  if (spc <= 0.0015625) { cerr << "radmult not designated for such a highly-refined surface yet" << endl; iA(0); }
  */
  /*
  if (spc <= 0.2) { rad = 10; }
  if (spc <= 0.1) { rad = 13; }
  if (spc <= 0.05) { rad = 13; }
  if (spc <= 0.025) { rad = 15; }
  if (spc <= 0.0125) { rad = 19; }
  if (spc <= 0.00625) { rad = 27; }
  if (spc <= 0.003125) { rad = 37; }
  if (spc <= 0.0015625) { cerr << "radmult not designated for such a highly-refined surface yet" << endl; iA(0); }
  PetscBool flg = PETSC_FALSE;
  int64_t opt_override; 
  PetscOptionsGetInt(NULL, "", "-bis3d_re_radmult", &opt_override, &flg);
  if(flg){
    rad = opt_override;
  }
  */
  ebiFunctionReturn(0);
}				

//Take surface discretization and return spacing var for bdsurf
#undef __FUNCT__
#define __FUNCT__ "bis3dspacing2bdsurfspacing"
int bis3dspacing2bdsurfspacing(const double bisspac, double& bdspac){
  ebiFunctionBegin;
  if (bisspac <= 0.2) { bdspac = 0.032; }
  if (bisspac <= 0.1) { bdspac = 0.016; }
  if (bisspac <= 0.05) { bdspac = 0.008; }
  if (bisspac <= 0.025) { bdspac = 0.004; }
  if (bisspac <= 0.0125) { bdspac = 0.002; }
  if (bisspac <= 0.00625) { bdspac = 0.001; }
  if (bisspac <= 0.003125) { bdspac = 0.0005; }
  if (bisspac <= 0.0015625) { cerr << "bd spacingt not designated for such a highly-refined surface yet" << endl; iA(0); }
  ebiFunctionReturn(0);
}

// Usable AlmostEqual function
//From: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "AE2C"
bool AE2C(double AA, double BB, int maxUlps)
{
  /* Cast these double to floats.  Really need to rewrite this for doubles,
	* but the values we are dealing with generally involve low-digit accuracy (radii of nodes
	* that are never too deep in octree), so this should be fine
	*/
  float A = (float)(AA); float B = (float)(BB);
  int MAXU = maxUlps;
  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  assert(MAXU > 0 && MAXU < 4 * 1024 * 1024);
  int aInt = *(int*)&A;
  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0) { aInt = 0x80000000 - aInt; }
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt = *(int*)&B;
  if (bInt < 0) { bInt = 0x80000000 - bInt; }
  int intDiff = abs(aInt - bInt);
  bool ret = false;
  if (intDiff <= MAXU) { ret =  true; }

  return ret;
}


/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "file2string"
int file2string(MPI_Comm comm, const char* infile, string& str)
{
  ebiFunctionBegin;

  int mpiRank;
  iC( MPI_Comm_rank(comm, &mpiRank) );
  if(mpiRank == 0) {
	 cerr << infile << endl;
	 //-----------------
	 ifstream fin(infile);
	 ostringstream sout;
	 
	 ebiAssert(fin.good());
	 
	 char tt[MAXSTRLEN];
	 //int i=0;
	 while(fin.eof()==false) {
		fin.getline(tt,MAXSTRLEN);
		sout<<tt<<endl;
	 }
	 //copy
	 char* buf = new char[sout.str().length()+1];
	 sout.str().copy(buf, sout.str().length());
	 buf[sout.str().length()] = 0;
	 //bcast length
	 int length = sout.str().length() + 1; // cstyle
	 iC( MPI_Bcast(&length, 1, MPI_INT, 0, comm) );
	 //bcast data
	 iC( MPI_Bcast(buf, length, MPI_CHAR, 0, comm) );
	 str = buf;
	 delete[] buf;
	 fin.close();
  } else {
	 //-----------------
	 int length;
	 iC( MPI_Bcast(&length, 1, MPI_INT, 0, comm) );
	 char* buf = new char[length];
	 iC( MPI_Bcast(buf, length, MPI_CHAR, 0, comm) );
	 str = buf;
	 delete[] buf;
  }
  
  ebiFunctionReturn(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "string2file"
int string2file(MPI_Comm comm, string& str, const char* outfile)
{
  ebiFunctionBegin;

  int mpiRank;
  iC( MPI_Comm_rank(comm, &mpiRank) );
  if(mpiRank == 0) {
	 //-----------------
	 ofstream fout(outfile);
	 istringstream sin(str);
	 ebiAssert(sin.good());
	 
	 string tt; //char tt;
	 while(sin.eof()==false) {
		sin>>tt;
		fout<<tt;
	 }
	 fout.close();
  }
  
  ebiFunctionReturn(0);
}

END_EBI_NAMESPACE

