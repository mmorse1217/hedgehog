#ifndef _BE3DOV_HPP_
#define _BE3DOV_HPP_

#include "solver_interface.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "solver_utils.hpp"

BEGIN_EBI_NAMESPACE

using std::pair;

//-------------------------------------------------
class Evaluator: public EbiObject
{
public:
  //-----------------------------

protected:
  Kernel3d _knl;

  PatchSamples* _patch_samples;
  
public:
  Evaluator(const string& n, const string& p): EbiObject(n,p)  {;}
  
  virtual ~Evaluator() {;}
  virtual int setup()=0;
  virtual int eval(Vec, Vec)=0;

  Kernel3d& knl() { return _knl; }

  int source_dof() { return _knl.srcDOF(); }  
  int target_dof() { return _knl.trgDOF(); }  
  
  void set_surface_discretization(PatchSamples* patch_samples){
      _patch_samples = patch_samples;
  }

};

END_EBI_NAMESPACE

#endif
