#ifndef _EVAL_NEAR_EXTRAPOLATION_HPP_
#define _EVAL_NEAR_EXTRAPOLATION_HPP_

#include "evaluator.hpp"
#include "evaluator_near.hpp"
#include "collocation_patch_samples.hpp"

BEGIN_EBI_NAMESPACE

class EvaluatorQBKIX: public Evaluator {
    protected:
        PatchSamples* _refined_patch_samples;
        unique_ptr<FMM> _fmm;
        CollocationPatchSamples* _collocation_data;
        CollocationPatchSamples* _refined_collocation_data;
        
        Vec _target_3d_position;
        Vec _target_as_face_point;
        Vec _closest_sample_3d_position;
        Vec _closest_sample_as_face_point;
        Vec _target_in_out;
        Vec _aux_interpolation_points;
        Vec _interpolation_directions;
        Vec _target_far_field;
        Vec _target_interpolant_spacing;
        
        ExpansionType _expansion_type;

    public:

    EvaluatorQBKIX(Kernel3d kernel,
            PatchSamples* patch_samples,
            PatchSamples* refined_patch_samples,
            Vec& target_3d_position,
            Vec& closest_sample_3d_position,
            Vec& closest_sample_as_face_point,
            Vec& target_in_out,
            Vec interpolation_directions = NULL,
            Vec target_far_field = NULL,
            Vec target_interpolant_spacing = NULL
            ):

        Evaluator("", ""),
        _refined_patch_samples(refined_patch_samples),
        _target_3d_position(target_3d_position),
        _closest_sample_3d_position(closest_sample_3d_position),
        _closest_sample_as_face_point(closest_sample_as_face_point),
        _target_in_out(target_in_out),
        _expansion_type(EXTRAPOLATE_ONE_SIDE),
        _interpolation_directions(interpolation_directions),
        _target_far_field(target_far_field),
        _target_interpolant_spacing(target_interpolant_spacing)
    
    {
        this->_patch_samples = patch_samples;
        this->_knl = kernel;
        _collocation_data = new CollocationPatchSamples();
        _refined_collocation_data = new CollocationPatchSamples();
    }

    ~EvaluatorQBKIX(){
        delete _collocation_data;
        delete _refined_collocation_data;
        if(_aux_interpolation_points){
            VecDestroy(&_aux_interpolation_points);

        }
        //delete _fmm;
    }

    int setup();
    int eval(Vec density, Vec potential);

    Vec compute_interpolation_target_potential(Vec density);
    Vec compute_refined_density(Vec density);

    void lagrange_extrapolation(DblNumMat target_3d_position, 
            DblNumMat closest_sample_3d_position, 
            vector<double> interpolation_nodes,
            DblNumMat interpolation_point_potential_local,
            double h,
            double (*eval_point)(double, double),
            DblNumMat& final_potential); 
    void lagrange_extrapolation(
            vector<double> interpolation_nodes,
            DblNumMat interpolation_point_potential_local,
            double h,
            DblNumVec eval_points,
            DblNumMat& final_potential); 

};

double extrapolation_eval_point_blendsurf_legacy(double distance_from_target_to_closest_sample, 
        double node_spacing);
double extrapolation_eval_point_blendsurf(double distance_from_target_to_closest_sample, 
        double node_spacing, double exp_distance_to_boundary);
double extrapolation_eval_point_qbkix(double distance_from_target_to_closest_sample, 
        double node_spacing, double exp_distance_to_boundary);

void evaluate_lagrange_interpolant(
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        double h,
        int target_dof,
        double eval_points,
        DblNumVec& final_potential); 


Vec compute_refined_density(Vec density, 
        PatchSamples* patch_samples, 
        PatchSamples* refined_patch_samples, 
        CollocationPatchSamples* refined_collocation_data,
        Kernel3d kernel);

void lagrange_extrapolation_bary(DblNumMat target_3d_position, 
        DblNumMat closest_sample_3d_position, 
        vector<double> interpolation_nodes,
        DblNumMat interpolation_point_potential_local,
        DblNumMat node_spacing,
        DblNumMat expansion_distance_to_boundary,
        int target_dof,
        double (*eval_point)(double, double, double),
        DblNumMat& final_potential); 

END_EBI_NAMESPACE

#endif
