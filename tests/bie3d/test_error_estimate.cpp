#include "../catch.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/p4est_interface.hpp"
#include "bie3d/markgrid.hpp"
#include "common/stats.hpp"
#include "bdry3d/on_surface_point.hpp"
#include "common/utils.hpp"
#include "../utils/evaluation_utils.hpp"
#include "common/vtk_writer.hpp"
#include "bie3d/error_estimate.hpp"

using namespace Ebi;
using ErrorEstimate::CurvatureDirection;
using Markgrid::NearFieldMap;

void test_pullback_and_pushforward(string filename, bool inside_domain, CurvatureDirection curvature_direction){
    // compute the pullback of a target point; push pullback forward through
    // panel, and check error \|target - P(P^{-1}(target))\|; should be less
    // than \eps_opt
    unique_ptr<PatchSurfFaceMap> face_map( new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;

    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/"+filename+".wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/"+filename+".wrl");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-bis3d_spacing", ".2");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "12");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", "1.");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    // sample surface
    unique_ptr<PatchSamples> samples(new PatchSamples("", ""));
    samples->bdry() = face_map.get();
    vector<int> patches(face_map->patches().size(), 0);
    samples->patch_partition() = patches;
    samples->setup(); 

    Vec targets;
    VecDuplicate(samples->sample_point_3d_position(), &targets);
    VecCopy(samples->sample_point_3d_position(), targets);

    // move target points off surface along the normal
    double d = .01;
    if(inside_domain){
        d *= -1;
    }
    VecAXPY(targets, d, samples->sample_point_normal()); 

    // Compute pullback of targets
    DblNumMat targets_local(3, targets);
    NumVec<OnSurfacePoint> closest_points = samples->sample_point_as_on_surface_point();
    for (int i = 0; i < closest_points.m(); i++) {
        closest_points(i).distance_from_target = fabs(d);
        if(inside_domain)
            closest_points(i).inside_domain = INSIDE;
        else
            closest_points(i).inside_domain = OUTSIDE;
    }

    ComplexNumVec pullbacks = compute_pullback(targets_local, closest_points, 
            curvature_direction, face_map.get());
    
    // Push forward the computed pullback t.
    DblNumMat push_forward = ErrorEstimate::push_forward(pullbacks, closest_points, 
            curvature_direction, face_map.get());
    //write_qbkix_points_to_vtk(targets_local, closest_points, 0, "output/pullback_true_targets_");
    //write_qbkix_points_to_vtk(push_forward, closest_points, 0, "output/pullback_approx_targets_");

    //write_face_map_mesh_to_vtk( face_map.get(), 0,"output/pullback_surface_", 12);



    // Check |targets - P(P^{-1}(t))| < \eps = 1e-7
    assert(push_forward.n() == targets_local.n());
    for (int i = 0; i < push_forward.n(); i++) {
        Point3 estimated_target(push_forward.clmdata(i));
        Point3 true_target(targets_local.clmdata(i));
        CHECK((estimated_target - true_target).l2() <= 1e-7);
    }

}


TEST_CASE("Test pullback computation", "[error-estimate][pullback][newton]"){
    SECTION("compute pullback on flat panel"){
        // pullback of a target on a flat panel should be equal to the target
        
        unique_ptr<PatchSurfFaceMap> face_map( new PatchSurfFaceMap("BD3D_", "bd3d_"));
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;

        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/flat_patch.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");

        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();


        assert(face_map->patches().size() == 1); // only a single flat patch
        auto patch = face_map->subpatch(0);
        int n = 11;
        double step_size = 2./(n-1);
        double curvature = 0.;


        for (int k = -1; k < 2; k+=2) {
            DblNumMat targets(3, n*n);
            double d = -.1*k;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int index = i*n + j;
                    Point3 target(-1. + i*step_size, -1. + j*step_size, d);
                    for (int d = 0; d < DIM; d++) {
                       targets(d,index) = target(d); 
                    }

                }
            }
            NumVec<OnSurfacePoint> closest_points(n*n);
            Markgrid::mark_near_field(targets, face_map.get(), closest_points);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int index = i*n + j;
                    // Compute pullback of target point
                    Point3 target(targets.clmdata(index));
                    OnSurfacePoint closest_point = closest_points(index);

                    imaginary pullback = ErrorEstimate::compute_pullback_under_patch(
                            target, closest_point, CurvatureDirection::FIRST, patch);

                    // curvature = 0 =>imaginary component should be equal to
                    // the physical distance from the patch
                    CHECK( fabs(d - pullback.imag()) <= 1e-6);
                    pullback = ErrorEstimate::compute_pullback_under_patch(
                            target, closest_point, CurvatureDirection::SECOND, patch);

                    // curvature = 0 =>imaginary component should be equal to
                    // the physical distance from the patch
                    CHECK( fabs(d - pullback.imag()) <= 1e-6);
                }
            }


        }

    } SECTION("compute pullback on curved panel"){
        // check that curve(pullback) = target, inside/outside domain, for each
        // curvature direction, and for several domains
        //vector<string> domains = {"cube", "pipe", "ppp"};
        vector<string> domains = {"cube"};
        vector<int> inside_outisde = {0, 1};
        vector<CurvatureDirection> curvatures= {CurvatureDirection::FIRST, CurvatureDirection::SECOND};
        
        for(auto const& domain: domains){
            for(auto const& curvature : curvatures){
                for(auto const& in_out : inside_outisde){
                    test_pullback_and_pushforward(domain, in_out, curvature);
                }
            }
        }
    }
}


TEST_CASE("Test error estimate", "[error-estimate][cases]"){
    SECTION(""){
        PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/flat_patch.wrl");
        PetscOptionsSetValue(NULL, "-poly_coeffs_file", "wrl_files/poly/flat_patch.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
        Options::set_value_petsc_opts("-bis3d_spacing", ".06666");


        unique_ptr<PatchSurfFaceMap> face_map( new PatchSurfFaceMap("BD3D_", "bd3d_"));
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();

        assert(face_map->patches().size() == 1); // only a single flat patch
        auto patch = face_map->subpatch(0);
        int n = 11;
        double step_size = 2./(n-1);
        double curvature = 0.;

        unique_ptr<PatchSamples> samples(new PatchSamples("", ""));
        samples->bdry() = face_map.get();
        vector<int> patches(face_map->patches().size(), 0);
        samples->patch_partition() = patches;
        samples->setup(); 

        Vec targets;
        VecDuplicate(samples->sample_point_3d_position(), &targets);
        VecCopy(samples->sample_point_3d_position(), targets);

        // move target points off surface along the normal
        double d = .0001;
        bool inside_domain = true;
        if(inside_domain){
            d *= -1;
        }
        VecAXPY(targets, d, samples->sample_point_normal()); 
        DblNumMat targets_local(DIM, targets);
        
        NumVec<OnSurfacePoint> closest_points = samples->sample_point_as_on_surface_point();
        DblNumMat uv_values = samples->sample_point_parametric_preimage(0);
        
        DblNumMat density_values(1, uv_values.n()); 
        setvalue(density_values, 1.);
        int q = int(1./Options::get_double_from_petsc_opts("-bis3d_spacing")) +1;
        for (int i = 0; i < targets_local.n(); i++) {
            Point3 target(targets_local.clmdata(i));
            OnSurfacePoint closest_point = closest_points(i);
            closest_point.distance_from_target = d;
            double error = ErrorEstimate::evaluate_error_estimate(patch, target, 
                    closest_point, q, uv_values, density_values);
            cout << error << endl;
        }

    } SECTION(""){

    }
}
