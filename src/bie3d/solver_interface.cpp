#include "solver_interface.hpp"
#include "vec2t.hpp"
#include <omp.h>

BEGIN_EBI_NAMESPACE

using std::pair;
using std::min;
using std::max;
using std::abs;
using std::cerr;
using std::ofstream;

SolverInterface::SolverInterface(const string& n, const string& p):
  EbiObject(n,p), _patch_samples(NULL), _refined_patch_samples(NULL) {;}

SolverInterface::~SolverInterface()
{
  if(_patch_samples!=NULL) delete _patch_samples;
  if(_refined_patch_samples!=NULL) delete _refined_patch_samples;
}
/*void SolverInterface::setFromOptions()
{
  ebiFunctionBegin;
  
  ebiFunctionReturn(0);
}*/



void SolutionDensity::localSize(int& la, int& lb, int& lm)
{
  int num_poles;
  if (_domain_boundedness == DOM_UNBND){
      num_poles = _number_of_poles;

  } else { 
      num_poles = _number_of_poles - 1;
  }
  
  la = _local_num_points * _kernel.get_sdof();
  lb = (this->mpiRank()==0) ? num_poles*_kernel.get_pdof() : 0;
  lm = la + lb;

}

// ---------------------------------------------------------------------- 
void SolutionDensity::globalSize(int& ga, int& gb, int& gm)
{

  int num_poles;
  if (_domain_boundedness == DOM_UNBND){
      num_poles = _number_of_poles;
  } else { 
      num_poles = _number_of_poles - 1;
  }

  ga = _global_num_points * _kernel.get_sdof();
  gb = num_poles*_kernel.get_pdof();
  gm = ga + gb;
}

END_EBI_NAMESPACE
