#ifndef _BIS3DOVGENERIC_HPP_
#define _BIS3DOVGENERIC_HPP_

#include "common/kernel3d.hpp"
#include "solver_interface.hpp"
#include "evaluator_on_surface.hpp"
#include "evaluator_qbkix.hpp"
#include "bdry3d/p4est_interface.hpp"

BEGIN_EBI_NAMESPACE

//------------------------------------------------------
//------------------------------------------------------
enum EvaluationType {
    SINGULAR_EVAL = 0,
    INTERIOR_EXTRAPOLATION = 1,
    EXTRAPOLATION_AVERAGE = 2,
    INTERPOLATION_ACROSS_SURFACE = 3,
    SMOOTH_QUADRATURE = 4,
    NONE = 5

};

class SolverGMRESDoubleLayer: public SolverInterface
{
protected:
    p4est_t* _p4est; 
  KSP _ksp;
  Mat _mat;
  Mat _pole_matrix;
  Mat _constraint_matrix; 
  Mat precond_schur_complement_inverse;

  DblNumMat _poles;
  Evaluator* _on_surface_evaluator;

  Equation_type _equation_type;
  EvaluationType _evaluation_type;

  bool _is_incompressible;

  int _nits;

  Kernel3d _problem_kernel;
  // Evaluator specific variables
  Vec _target_in_out;
  Vec _interpolation_directions;

public:
  bool _compute_refined_surface;
  SolverGMRESDoubleLayer(const string& n, const string& p);
  SolverGMRESDoubleLayer(PatchSurf* surface);
  ~SolverGMRESDoubleLayer();
  
  int setup();
  int solve(Vec b0, Vec x);
  int setFromOptions();
  int mmult(SolutionDensity& solution_in, SolutionDensity& solution_out);
  int pcmult(SolutionDensity& solution_in, SolutionDensity& solution_out);
  static int mmultWrapper(Mat, Vec, Vec);
  static int pcmultWrapper(PC, Vec, Vec);

  int core_evaluation(SolutionDensity& solution_in, Vec out, Kernel_variable qt, Vec tp,
          Evaluator *evaluator);

  int fareval(Vec tp, int qt, Vec den, Vec val);
  int roneval(Vec tp, Vec tb, int qt, Vec den, Vec val);
  //int neaeval(Vec tp, Vec ti, Vec bp, Vec bb, int qt, Vec den, Vec val, Vec reg);

  int neaeval(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        int quantity, 
        Vec density, 
        Vec potential); 


  int neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        int quantity, 
        Vec density, 
        Vec potential); 

  int neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential); 
  int neaeval_extrapolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec dist_from_boundary,  //closest on surface sample for interpolation
        Vec interp_point_spacing,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential); 

  int neaeval_extrapolate_average(
        Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec dist_from_boundary,  //closest on surface sample for interpolation
        Vec interp_point_spacing,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential);
  int neaeval_extrapolate_average(
        Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential);
  int neaeval_interpolate(Vec near_targets,  // targets in Omega_2
        Vec intermediate_targets, //targets in Omega_1
        Vec target_in_out_near, // whether the Omega_2 targets are in/out
        Vec closest_sample_3d_position,  //closest on surface sample for interpolation
        Vec closest_sample_as_face_point,
        Vec interpolation_directions,
        int quantity, 
        Vec density, 
        Vec potential); 

  void evaluate(Vec target_points, 
          Vec density,
          Vec &potential,
          NumVec<OnSurfacePoint>& on_surface_points, 
          Kernel_variable kernel_variable,
          bool precomputed_closest_points=false);
  void evaluate_subvector(
          // global arrays for the whole problem
          Vec target_points, 
          Vec density,
          Vec &potential,
          Vec closest_on_surface_points_3d_position,
          Vec closest_on_surface_points_as_face_point,
          Vec interpolation_directions,
          Vec target_far_field,
          Vec target_interpolant_spacing,
          Vec target_in_out,
          NumVec<OnSurfacePoint>& on_surface_points,
          // local indices for the current subvector.
          vector<int64_t>sub_vector_indices,
          Kernel_variable kernel_variable,
          EvaluationType evaluation_type);
  
  int& nits() { return _nits; }  
  
  void set_evaluation_type(EvaluationType e){
      _evaluation_type =  e;
  }
  Equation_type equation_type(){
      return _equation_type;
  }

  const Kernel3d problem_kernel(){
    return _problem_kernel;
  }

Vec singular_correction(Vec sources, 
        Vec targets, Vec targets_as_face_points, Vec potential, Kernel3d kernel);

void populate_qbx_data(const NumVec<OnSurfacePoint> on_surface_points,
        Vec& target_in_out, Vec& closest_sample_3d_position,
    Vec& closest_on_surface_points_as_face_point,
    Vec& interpolation_directions,
    Vec& target_far_field,
    Vec& target_interpolant_spacing
    );

};
Vec greens_identity(
        MPI_Comm comm,
        Kernel3d problem_kernel, 
        Vec singularity_location,
        Vec singularity_strength,
        Vec target_points,
        SolverGMRESDoubleLayer* solver);

Vec greens_identity(
        MPI_Comm comm,
        Vec dirichlet_data,
        Vec neumann_data,
        Vec target_points,
        SolverGMRESDoubleLayer* solver);
Vec greens_identity(
        MPI_Comm comm,
        Vec dirichlet_data,
        Vec neumann_data,
        NumVec<OnSurfacePoint> on_surface_points_targets,
        Vec target_points,
        SolverGMRESDoubleLayer* solver);
END_EBI_NAMESPACE

#endif
