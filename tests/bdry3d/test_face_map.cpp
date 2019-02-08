
#include "../catch.hpp"
#include "common/nummat.hpp"
#include "common/utils.hpp"
#include "bie3d/solver_utils.hpp"
#include "bie3d/markgrid.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bdry3d/patch_samples.hpp"

using namespace Ebi;

// Initialize polynomial density \phi(y), y = (x,y,z) \in \partial\Omega
// phi(y) = (1,1,1)*(\sum_k=0^poly_degree x^k*y^k*z^k)
void compute_polynomial_density(int poly_degree,
        DblNumMat sample_point_local, DblNumMat& density_local){
    int sdof = density_local.m();
    int num_local_points = sample_point_local.n();
    assert(sample_point_local.n() == density_local.n());

    for(int i =0; i < num_local_points; i++){
        for(int d = 0; d < sdof; d++){
            density_local(d,i) = 1. + sample_point_local(0,i)*sample_point_local(1,i)*sample_point_local(2,i);
            /*
            for(int kk = 0; kk <= poly_degree; kk++){
                double kkth_term = 1.;
                for(int j = 0; j < DIM; j++)
                    kkth_term *= pow(sample_point_local(j,i), kk);

                density_local(d,i) += kkth_term;
            }*/
        }
    }
}

void test_gauss_bonnet(string domain, int euler_characteristic){
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/"+domain+".wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/explicit_torus_patches.poly");
    Options::set_value_petsc_opts("-bis3d_spacing", ".03333");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    
    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    //face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();
    
    unique_ptr<PatchSamples> patch_samples(new PatchSamples("",""));
    vector<int> patch_partition(face_map->patches().size(), 0 );
    patch_samples->bdry() = face_map.get();
    patch_samples->patch_partition() = patch_partition;
    patch_samples->setup();

    int num_samples = patch_samples->local_num_sample_points();
    int num_samples_per_patch = patch_samples->num_sample_points()[0]; // valid for face-map only, all patches have same number of points
    int num_patches = face_map->num_patches(); 

    DblNumVec gaussian_curvature(num_samples);
    DblNumVec quadrature_weights(patch_samples->sample_point_combined_weight());
    DblNumMat uv_coords(2, patch_samples->sample_point_parametric_preimage());
    int iter = 0;
    
    // Evaluate Gaussian curvature at samples
    for (int pi = 0; pi < num_patches; pi++) {
        for (int i = 0; i < num_samples_per_patch*num_samples_per_patch; i++) {
            Point2 uv(uv_coords.clmdata(i));
            auto patch = face_map->subpatch(pi);
            gaussian_curvature(iter++) = patch->gaussian_curvature(uv);
        }
    }

    // compute the integral
    double integral = dot(gaussian_curvature, quadrature_weights);
    cout << "integral: " << integral << endl;
    //int euler_characteristic = 2;
    CHECK(fabs(integral - 2.*M_PI*euler_characteristic)+ 1e-3*fabs(2.*M_PI*euler_characteristic) <= 1e-3);
}

TEST_CASE("Test face-map density refinement","[face-map][geom]"){
    // 0: Initialize surfaces and patch samplings 
    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    PatchSamples* patch_samples = new PatchSamples("","");
    PatchSamples* refined_patch_samples = new PatchSamples("","");

    // Setup it up
    PetscOptionsSetValue(NULL, "-bd3d_filename", "wrl_files/cube.wrl");
    PetscOptionsSetValue(NULL, "-bd3d_meshfile", "wrl_files/cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-3");
    //Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    face_map->setFromOptions();
    face_map->setup();
    face_map->refine_test();

    vector<int> patch_partition(face_map->patches().size(), 0 );

    patch_samples->bdry() = face_map;
    patch_samples->patch_partition() = patch_partition;
    patch_samples->setup();

    refined_patch_samples->bdry() = face_map;
    refined_patch_samples->patch_partition() = patch_partition;
    refined_patch_samples->setup(true);

    /*SECTION("constant density"){

        // Laplace kernel first
        Kernel3d kernel(121, vector<double>(1,1.0));

        int num_local_points = patch_samples->local_num_sample_points();
        int num_local_refined_points = refined_patch_samples->local_num_sample_points();
        
        Vec density;
        VecCreateMPI(PETSC_COMM_WORLD,
            num_local_points*kernel.get_sdof(),
            PETSC_DETERMINE,
            &density);
        
        // constant density
        const double one = 1.0;
        VecSet(density, one);
        
        // Interpolate to refined surface
        Vec refined_density = patch_samples->refine_density(
                    kernel.get_sdof(),
                    0, // TODO remove variable unused
                    density,
                    refined_patch_samples);
        
        // Check that we've interpolated 1 exactly
        // (we're using high order interpolation so we should be able to get the
        // constant function exactly on all patches, i.e. every point is 1)
        DblNumMat refined_density_local = get_local_vector(kernel.get_sdof(), 
                num_local_refined_points,
                refined_density);
        for(int i = 0; i < refined_density_local.n(); i++){
            for(int d = 0; d < kernel.get_sdof(); d++)
                CHECK((fabs(refined_density_local(d,i) -1.) <= 1e-14));
        }
        VecDestroy(&density);
        VecDestroy(&refined_density);
    }
    SECTION("polynomial density"){

        // Laplace kernel first
        Kernel3d kernel(121, vector<double>(1,1.0));

        int num_local_points = patch_samples->local_num_sample_points();
        int num_local_refined_points = refined_patch_samples->local_num_sample_points();
        
        Vec density;
e1       VecCreateMPI(PETSC_COMM_WORLD,
            num_local_points*kernel.get_sdof(),
            PETSC_DETERMINE,
            &density);
        
        int max_degree = 2;
        for(int k = 0; k < max_degree; k++){

            DblNumMat density_local = get_local_vector(kernel.get_sdof(), 
                    num_local_points,
                    density);
            DblNumMat sample_point_local = get_local_vector(DIM,
                    num_local_points,
                    patch_samples->sample_point_3d_position());
            
            // compute the density
            compute_polynomial_density(k, sample_point_local, density_local);
            
            density_local.restore_local_vector();
            sample_point_local.restore_local_vector();



            DblNumMat refined_sample_point_local= get_local_vector(
                    DIM,
                    num_local_refined_points,
                    refined_patch_samples->sample_point_3d_position());

            DblNumMat refined_sample_point_preimage= get_local_vector(
                    2,
                    num_local_refined_points,
                    refined_patch_samples->sample_point_parametric_preimage());



            // Interpolate to refined surface
            Vec refined_density = patch_samples->refine_density(
                    kernel.get_sdof(),
                    0, // TODO remove variable unused
                    density,
                    refined_patch_samples);

            // Check that we've interpolated 1 exactly
            // (we're using high order interpolation so we should be able to get
            // polynomials of degree n exactly on all patches if face-map uses
            // polynomial patches of bidegree n)
            Vec true_refined_density;
            VecCreateMPI(PETSC_COMM_WORLD,
                num_local_refined_points*kernel.get_sdof(),
                PETSC_DETERMINE,
                &true_refined_density);
        
            DblNumMat true_refined_density_local = get_local_vector(
                    kernel.get_sdof(), 
                    num_local_refined_points,
                    true_refined_density);


            
            // Compute true refined_density (same as above)
            compute_polynomial_density(k, refined_sample_point_local, true_refined_density_local);
            
            DblNumMat refined_density_local = get_local_vector(
                    kernel.get_sdof(), 
                    num_local_refined_points,
                    refined_density);
            for(int i = 0; i < refined_density_local.n(); i++){
                for(int d = 0; d < kernel.get_sdof(); d++){
                    CHECK(fabs(refined_density_local(d,i) - true_refined_density_local(d,i)) <= 1e-2);
                }
            }
            refined_sample_point_local.restore_local_vector();
            refined_sample_point_preimage.restore_local_vector();
            refined_density_local.restore_local_vector();
            true_refined_density_local.restore_local_vector();
            VecDestroy(&true_refined_density);
            VecDestroy(&refined_density);
        }
        VecDestroy(&density);

    }
    SECTION("Test qbkix point generation per patch"){

        // sample the surface; these are the near targets
        Vec near_targets;
        VecDuplicate(patch_samples->sample_point_3d_position(), &near_targets);
        VecCopy(patch_samples->sample_point_3d_position(), near_targets);
        
        // for on-surface eval, the target is the closest on-surface point
        Vec closest_on_surface_points;
        VecDuplicate(near_targets, &closest_on_surface_points);
        VecCopy(near_targets, closest_on_surface_points);
        
        // the interpolation direction is the interior pointing normals
        Vec interpolation_directions;
        VecDuplicate(patch_samples->sample_point_normal(), &interpolation_directions);
        VecCopy(patch_samples->sample_point_normal(), interpolation_directions);
        VecScale(interpolation_directions, -1.);
        
        // load options file value into a big vector :(
        int num_local_points = patch_samples->local_num_sample_points();
        Vec far_field;
        Vec interpolant_spacing;
        VecCreateMPI(
                MPI_COMM_WORLD,
                num_local_points,
                PETSC_DETERMINE,
                &far_field);
        VecCreateMPI(
                MPI_COMM_WORLD,
                num_local_points,
                PETSC_DETERMINE,
                &interpolant_spacing);
        VecSet(far_field, Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
        VecSet(interpolant_spacing, Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));




        VecDestroy(&near_targets);
        VecDestroy(&closest_on_surface_points);
    }*/
    SECTION("Test FaceMapSubPatch bounding box code"){
        for(int pi = 0; pi < face_map->patches().size(); pi++){
            FaceMapSubPatch* patch = face_map->subpatch(pi);
            
            // find bounding box on patch
            Point3 min, max;
            patch->bounding_box(min, max);

            // sample the patch
            DblNumMat samples_on_pith_patch = 
                patch_samples->sample_point_3d_position(pi);

            // each on-surface sample should be inside the patch's bounding box
            for(int ti =0; ti < samples_on_pith_patch.n(); ti++){
                Point3 target(samples_on_pith_patch.clmdata(ti));
                CHECK(min < target);
                CHECK(max > target);
            }
        }
        
    }
    /*
    SECTION("Test single qbkix point contained in bounding box"){
    }*/
}

TEST_CASE("Test patch size after refinement", "[face-map][geom][test]"){
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
    // Load patch + polynomial from file
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/flat_patch.wrl");
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/flat_patch.wrl");
    Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/flat_patch.poly");
    PatchSurfFaceMap* face_map= new PatchSurfFaceMap("BD3D_", "bd3d_");
    face_map->_surface_type = PatchSurfFaceMap::POLYNOMIAL;


    SECTION("Test FaceMapSubPatch characteristic length"){
        // Single patch defined on [-1,1] x [-1, 1] x 0
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
        face_map->setFromOptions();
        face_map->setup();
        face_map->refine_test();
        REQUIRE(face_map->num_patches() == 1);
        cout << "im confused " << face_map->num_patches() << endl;
        auto patch = face_map->subpatch(0);
        CHECK(fabs(patch->characteristic_length() - 2.) < 1e-9);


    }
    SECTION("Test FaceMapSubPatch refined characteristic length"){
        for(int i =0; i <5; i++){
            //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", to_string(i));
            face_map->setFromOptions();
            face_map->setup();
            //face_map->refine_test();
            face_map->refine_uniform(i);
            int num_patches = pow(4,i);
            REQUIRE(face_map->num_patches() == num_patches);
            cout << "im confused " << face_map->num_patches() << endl;
            for(int pi =0; pi < num_patches; pi++){
            auto patch = face_map->subpatch(pi);
            cout << patch->characteristic_length() << "," << sqrt(4./double(num_patches)) << endl;
            CHECK(fabs(patch->characteristic_length() - sqrt(4./double(num_patches))) < 1e-6);
            }
        }
    }
}

TEST_CASE("Test Gauss-Bonnet Theorem", "[diff-geom][face-map][geom]"){
 test_gauss_bonnet("cube",2);
 test_gauss_bonnet("pipe",2);
 test_gauss_bonnet("newtorus",0);
}
