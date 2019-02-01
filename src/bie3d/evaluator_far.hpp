#ifndef _BE3DOVFAR_HPP_
#define _BE3DOVFAR_HPP_

#include "evaluator.hpp"

BEGIN_EBI_NAMESPACE

// ---------------------------------------------------------------------- 
class EvaluatorFar: public Evaluator
{
protected:
  //PARAMS(REQ)
  Vec _target_3d_position;

public:

  unique_ptr<FMM> fmm;

  EvaluatorFar(const string& n, const string& p): Evaluator(n,p), _target_3d_position(NULL) {;} 
  
  ~EvaluatorFar() {;}
  //accessors
  Vec& target_3d_position() { return _target_3d_position; } 
  //int setFromOptions(); //nothing to do
  int setup();
  int eval(Vec den, Vec val);
};

END_EBI_NAMESPACE

#endif
