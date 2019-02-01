#ifndef _BE3DOVRON_HPP_
#define _BE3DOVRON_HPP_

#include "evaluator.hpp"
#include "collocation_patch_samples.hpp"

BEGIN_EBI_NAMESPACE

// ---------------------------------------------------------------------- 
class EvaluatorOnSurface: public Evaluator
{
//public:
  //----------------------------------

  class UnitVectorDensities { //cstmdfy data
  protected:
    // DZ  rename _data, consider adding operator[] instead of accessor
	 vector<Vec> _untvec;
  public:
	 UnitVectorDensities() {;}
	 ~UnitVectorDensities() { clear(); }
	 vector<Vec>& untvec() { return _untvec; }
	 void clear() {
		for(uint i=0; i<_untvec.size(); i++) {
		  if(_untvec[i]!=NULL) { VecDestroy(&_untvec[i]); _untvec[i]=NULL; }
		}
	 }
  };
protected:
  //PARAMS(REQ)  
  Vec _target_3d_position; 
  Vec _target_as_face_point; 
  //Vec _density_at_targets;
  
  //COMPONENTS
  CollocationPatchSamples* _collocation_data;  
  UnitVectorDensities _singularity_cancellation_data; 
  
  //TEMPORARY VARIABLE
  Vec _dat;
  vector<DblNumMat> _refined_datvec;
  vector<DblNumMat> _refined_positions;
  vector<DblNumMat> _refined_normals;
  int64_t _surface_interpolation_num_samples;
  int64_t _refinement_factor;
  
public:
  FMM* fmm;

  EvaluatorOnSurface(const string& n, const string& p):
	 Evaluator(n,p),  _dat(NULL), fmm(NULL)  {
         _collocation_data = new CollocationPatchSamples();
     } 
  EvaluatorOnSurface(Kernel3d kernel, Vec target_3d_position, Vec target_as_face_point, 
          PatchSamples* patch_samples);
  ~EvaluatorOnSurface() {
	 if(fmm!=NULL) delete fmm;
	 if(_collocation_data!=NULL) delete _collocation_data;
  }
  int setFromOptions();
  int setup();
  int setup_no_fmm();
  int eval(Vec den, Vec val);
  
  //accessor
  Vec& target_3d_position() { return _target_3d_position; } 
  Vec& target_as_face_point() { return _target_as_face_point; }
  //eval
  int singular_evaluation(Vec den, Vec val); //called by eval() and setup_unit_vector_densities()  //int singular_evaluationbmdfy(Vec den, Vec val); //bruno type modification

  int apply_singularity_cancellation(Vec den, Vec val); //constant modification
  
protected:
  //setup
  //int distribute_collocation_points(Vec face_point, Vec pos, CollocationPointData&);
  bool use_singularity_cancellation(Kernel3d& knl); //use constant modification or not

  int setup_unit_vector_densities();

  
  int subtract_inaccurate_part(int pi, int sdof, double* xy, int floating_POU_radius_scale,
          DblNumMat& pos, DblNumMat& nor, DblNumVec& jaw, DblNumMat& dat);
  int add_singular_quadrature_part(int pi, int sdof, double* xy, int floating_POU_radius_scale,
          DblNumMat& pos, DblNumMat& nor, DblNumVec& jaw, DblNumMat& dat);
  //int RADMULT();  //int RAD(int pi);
  double eta(double dist);
};

END_EBI_NAMESPACE

#endif
