#ifndef _EVAL_NEAR_EXTRAPOLATION_AVERAGE_HPP_
#define _EVAL_NEAR_EXTRAPOLATION_AVERAGE_HPP_

#include "evaluator_qbkix.hpp"
BEGIN_EBI_NAMESPACE

class EvaluatorQBKIXAverage: public EvaluatorQBKIX{
    private:
        Vec _interior_interpolation_points;
        Vec _exterior_interpolation_points;
    public:
    EvaluatorQBKIXAverage(Kernel3d kernel,
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
    EvaluatorQBKIX(kernel,
            patch_samples,
            refined_patch_samples,
            target_3d_position,
            closest_sample_3d_position,
            closest_sample_as_face_point,
            target_in_out,
            interpolation_directions,
            target_far_field,
            target_interpolant_spacing){;}
    /*~EvaluatorQBKIXAverage(){
        if(_interior_interpolation_points)
            delete _interior_interpolation_points;
        if(_exterior_interpolation_points)
            delete _exterior_interpolation_points;
    }*/

    int setup();
    int eval(Vec density, Vec potential);
    int eval_2(Vec density, Vec potential);
    Vec construct_check_points();

};

END_EBI_NAMESPACE
#endif
