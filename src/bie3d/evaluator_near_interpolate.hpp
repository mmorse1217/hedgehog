#ifndef _EVAL_NEAR_INTERPOLATION_AVERAGE_HPP_
#define _EVAL_NEAR_INTERPOLATION_AVERAGE_HPP_

#include "evaluator_qbkix.hpp"

// MJM TODO redesign EvaluatorNear constructor and member variables so that this
// can be directly extended from EvaluatorNear
BEGIN_EBI_NAMESPACE
class EvaluatorNearInterpolate: public EvaluatorQBKIX{

    public:
    EvaluatorNearInterpolate(Kernel3d kernel,
            PatchSamples* patch_samples,
            PatchSamples* refined_patch_samples,
            Vec& target_3d_position,
            Vec& closest_sample_3d_position,
            Vec& closest_sample_as_face_point,
            Vec& target_in_out,
            Vec interpolation_directions = NULL):

        EvaluatorQBKIX(kernel, 
                patch_samples,
                refined_patch_samples,
                target_3d_position,
                closest_sample_3d_position,
                closest_sample_as_face_point,
                target_in_out,
                interpolation_directions) {
            _expansion_type = INTERPOLATE_ACROSS_SURFACE;
        }

    int setup();
    int eval(Vec density, Vec potential);
    static vector<double> get_interpolation_nodes();
};
Vec evaluate_smooth_quadrature(MPI_Comm comm, Vec refined_samples, Vec refined_normals, 
        Vec targets, Kernel3d kernel, Vec refined_density);
END_EBI_NAMESPACE
#endif
