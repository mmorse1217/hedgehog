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
#include "bie3d/evaluator_far.hpp"

using namespace hedgehog;
using ErrorEstimate::CurvatureDirection;
using Markgrid::NearFieldMap;

void perturb_flat_patch2(unique_ptr<PatchSurfFaceMap>& face_map, double tx, double ty){
    assert(face_map->num_patches() == 1);
    auto coefficients= face_map->face_map().control_points_write(0);
    // 4x4 grid of coefficients
    Point3 dx(0.,0.,tx);
    Point3 dy(0.,0.,ty);
    for (int i = 0; i < 4; i+=3) {
        for (int j = 1; j < 3; j++) {
            (*coefficients)(i*4 + j) += dx/2.; 
        }
    }
    for (int i = 1; i < 3; i++) {
        for (int j = 0; j < 4; j+=3) {
            (*coefficients)(i*4 + j) += dy/2.; 
        }
    }
    for (int i = 0; i < 4; i+=3) {
        for (int j = 0; j < 4; j+=3) {
            (*coefficients)(i*4 + j) += (dy+dx)*.5; 
            
        }
    }
        
    
}

double compute_integral2(PatchSurfFaceMap* face_map, Point3 target_point){
    //setup sampling and quadrature evaluator with current spacing
    unique_ptr<PatchSamples> samples (new PatchSamples("", ""));
    int num_patches = face_map->num_patches();
    vector<int> patch_partition(num_patches, 0);
    samples->bdry() = face_map;
    samples->patch_partition() = patch_partition;
    samples->setup();
    
    Vec target;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1*DIM, target);
    PetscInt idx[3] = {0,1,2};

    VecSetValues(target, 3, idx, target_point.array(), INSERT_VALUES);
    VecAssemblyBegin(target);
    VecAssemblyEnd(target);
    
    Kernel3d kernel(121, vector<double>(2,1.));
    unique_ptr<EvaluatorFar> evaluator(new EvaluatorFar("", ""));
    evaluator->knl()                  = kernel;
    evaluator->target_3d_position()   = target;
    evaluator->set_surface_discretization(samples.get());
    evaluator->setFromOptions();
    evaluator->setup();
    
    Vec density;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, samples->local_num_sample_points(), density);
    Point3 charge_location(0., 5., 0.);
    {
        DblNumMat density_local(1,density);
        DblNumMat sample_points_local(DIM, samples->sample_point_3d_position());
        for(int i = 0; i < density_local.n(); i++){
            Point3 y(sample_points_local.clmdata(i));
            density_local(0,i) = 1./(y-charge_location).length();
        }
    }
    Vec potential;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, 1, potential);
    VecSet(potential, 0.);
    // compute the integral
    evaluator->eval(density,  potential);

    double integral = 0;
    {
        DblNumMat p_local(1, potential);
        integral = p_local(0,0);
    }
    VecDestroy(&target);
    VecDestroy(&density);
    VecDestroy(&potential);
    return integral;
}

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
        Options::set_value_petsc_opts("-target_accuracy", "1e-9");
        Options::set_value_petsc_opts("-bis3d_spacing", ".06666");


        unique_ptr<PatchSurfFaceMap> face_map( new PatchSurfFaceMap("BD3D_", "bd3d_"));
        face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        face_map->setFromOptions();
        face_map->setup();
        double shift = 1.;
        perturb_flat_patch2(face_map, shift, shift);
        face_map->refine_test();
       
        unique_ptr<PatchSurfFaceMap> refined_face_map( new PatchSurfFaceMap("BD3D_", "bd3d_"));
        refined_face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        refined_face_map->setFromOptions();
        refined_face_map->setup();
        perturb_flat_patch2(refined_face_map, shift, shift);
        refined_face_map->refine_uniform(4);
        
        
        cout << "num_patches: "<< face_map->patches().size() << endl;
        assert(face_map->patches().size() == 1); // only a single flat patch
        auto patch = face_map->subpatch(0);

        unique_ptr<PatchSamples> samples(new PatchSamples("", ""));
        samples->bdry() = face_map.get();
        vector<int> patches(face_map->patches().size(), 0);
        samples->patch_partition() = patches;
        samples->setup(); 

        Vec targets;
        VecDuplicate(samples->sample_point_3d_position(), &targets);
        VecCopy(samples->sample_point_3d_position(), targets);

        // move target points off surface along the normal
        //double d = .05;
        double d = .3;
        bool inside_domain = true;
        VecAXPY(targets, d, samples->sample_point_normal()); 
        DblNumMat targets_local(DIM, targets);
        
        NumVec<OnSurfacePoint> closest_points = samples->sample_point_as_on_surface_point();
        DblNumMat uv_values = samples->sample_point_parametric_preimage(0);
        
        DblNumMat density_values(1, uv_values.n()); 
        setvalue(density_values, 1.);
        int q = int(1./Options::get_double_from_petsc_opts("-bis3d_spacing")) +1;
        for (int i = q*q/2; i < q*q/2+ 1; i++){//targets_local.n(); i++) {
        //for (int i = 0; i < targets_local.n(); i++) {
            Point3 target(targets_local.clmdata(i));
            OnSurfacePoint closest_point = closest_points(i);
            closest_point.distance_from_target = d;
            double error = ErrorEstimate::evaluate_error_estimate(patch, target, 
                    closest_point, q, uv_values, density_values);
            double computed_integral = compute_integral2(face_map.get(), target);
            double true_integral= compute_integral2(refined_face_map.get(), target);
            cout << "true_integral: " << true_integral << endl;
            cout << "computed_integral: " << computed_integral << endl;
            cout << "estimated error: " << error << endl;
            double true_error = fabs(true_integral-computed_integral)/fabs(true_integral);
            cout << "true error: " << true_error << endl;
            CHECK(fabs(true_error - error)/fabs(true_error) <=10); 
        }

    } SECTION(""){

    }
}

