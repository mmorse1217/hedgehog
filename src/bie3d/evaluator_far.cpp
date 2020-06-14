#include "evaluator_far.hpp"
#include "solver_utils.hpp"
#include "common/stats.hpp"

BEGIN_EBI_NAMESPACE

int EvaluatorFar::setup()
{
  ebiFunctionBegin;
  ebiAssert(_target_3d_position!=NULL);

  PatchSamples* patch_samples = this->_patch_samples;
  fmm = unique_ptr<FMM>(new PvFMM(patch_samples->sample_point_3d_position(),
                 patch_samples->sample_point_normal(),
                 _target_3d_position, this->knl()));
  

  ebiFunctionReturn(0);
}

int EvaluatorFar::eval(Vec density, Vec val)
{
  ebiFunctionBegin;

  PatchSamples* patch_samples = this->_patch_samples;
  double zero = 0.0;
   VecSet( val, zero);
  // duplicate input density and scale by the weights
  Vec scaled_density;
   VecDuplicate(density, &scaled_density);
   denscale(this->source_dof(), patch_samples->sample_point_combined_weight(), density, scaled_density);

  // put DblNumVec interfaces on the local part of scaled density array and output
  // array for FMM
  
  // Evaluate
  fmm->evaluate(scaled_density, val);
  //((PvFMM*)fmm.get())->evaluate_direct(scaled_density, val);


  VecDestroy(&scaled_density);
  ebiFunctionReturn(0);
}


END_EBI_NAMESPACE
