#ifndef __ERROR_ESTIMATE_HPP__
#define __ERROR_ESTIMATE_HPP__

#include "common/numtns.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"
#include <omp.h>
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "bdry3d/face_map_subpatch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "fmm3d/pvfmm_bis_interface.hpp"
#include "common/ebiobject.hpp"
#include "common/stats.hpp"
#include "common/numvec.hpp"


BEGIN_EBI_NAMESPACE

namespace ErrorEstimate {
    class CircleArc {
        double _curvature;
        double _r;
        imaginary theta(imaginary t){
            double a = .5;
            return a*t/_r;
        }

        public: 
        CircleArc(double curvature): _curvature(curvature) {
            // to prevent division by zero. Use a curvature that will produce a
            // panel that will function as flat
            if(fabs(curvature) <=1e-10){
                double sign = (_curvature > 0) - (_curvature < 0);
                if(sign == 0.)
                    sign = 1.;
                _curvature = sign*1e-10;
            }
            _r = 1./_curvature;
        }
        imaginary operator()(imaginary t){
            imaginary I(0., 1.);
            imaginary theta_t = theta(t);
            return I*_r*(cos(theta_t) + I*sin(theta_t)-1.);
        }

        imaginary derivative(imaginary t){

            double a = .5;
            imaginary I(0., 1.);
            imaginary theta_t = a*t/_r;
            return I*a*(-sin(theta_t) + I*cos(theta_t));
        }

        imaginary second_derivative(imaginary t){
            double a = .5;
            imaginary I(0., 1.);
            imaginary theta_t = a*t/_r;
            return -a*a/_r*I*( cos(theta_t) + I*sin(theta_t));

        }
    };


    enum CurvatureDirection {
        FIRST = 0,
        SECOND = 1
    };

    ComplexNumVec compute_pullback(
            DblNumMat target_points,
            NumVec<OnSurfacePoint> closest_points,
            CurvatureDirection curvature_direction,
            PatchSurfFaceMap* face_map);
    
    imaginary compute_pullback_under_patch(
            Point3 target,
            OnSurfacePoint closest_point,
            CurvatureDirection curvature_direction,
            FaceMapSubPatch* patch);
    
    DblNumMat push_forward(
            ComplexNumVec pullbacks,
            NumVec<OnSurfacePoint> closest_points,
            CurvatureDirection curvature_direction,
            PatchSurfFaceMap* face_map);

    /*
     * computes error estimate at a complex target due to density_values and the
     * curvature of arc. Note that if density_values is a vector-valued, it
     * computes the maximum error estimate over each component
     *
     * @param CircleArc     arc     1d approx of patch centered at target point
     *                              curvature needs to be rescaled by patch
     *                              length for accurate estimate
     * @param imaginary     target  target point position relative to circle arc
     *                              rescaled by patch length
     * @param int           quadrature_order   
     * @param DblNumMat     density_values      values of the density along
     *                                          [-1.1] of the arc. should be
     *                                          quadrature_order number of
     *                                          values
     * @return double               error estimate at target point 
     */
double evaluate_error_estimate_on_patch( FaceMapSubPatch* patch,
        CurvatureDirection curvature_direction, Point3 target, 
        OnSurfacePoint closest_point,
        int quadrature_order, DblNumMat density_values);
    /*
     * Compute the error estimate at target due to patch and density_values.
     * @param FaceMapSubPatch* patch            patch to integrate
     * @param Point3        target              target point at which to
     *                                          estimate the quad error
     * @param OnSurfacePoint closest_point      closest_point to target on patch
     *
     * @param int           quadrature_order   
     * @param DblNumMat     density_values      values of the density at
     *                                          quadrature points. should be q^2
     *                                          points.
     * @param DblNumMat     uv_values           uv coordinates of the density
     *                                          values on patch
     * @return double               error estimate at target point 
     */
    double evaluate_error_estimate( FaceMapSubPatch* patch, Point3 target, 
            OnSurfacePoint closest_point, int quadrature_order, 
            DblNumMat uv_values, DblNumMat density_values);
}

END_EBI_NAMESPACE
#endif
