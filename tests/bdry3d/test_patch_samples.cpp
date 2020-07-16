#include "../catch.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "common/utils.hpp"
#include "common/stats.hpp"
#include "../utils/evaluation_utils.hpp"
#include "../utils/regression_utils.hpp"
using namespace Ebi;


void setup_face_map(unique_ptr<PatchSurfFaceMap>& surface, 
        unique_ptr<PatchSamples>& samples){
    surface = std::move(unique_ptr<PatchSurfFaceMap> (new PatchSurfFaceMap("BD3D_", "bd3d_")));
    surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
    surface->_coarse = true;
    surface->setFromOptions();
    surface->setup();
    surface->refine_test();

    vector<int> partition(surface->patches().size(), 0);
    samples = std::move(unique_ptr<PatchSamples>(new PatchSamples("", "")));
    samples->bdry() = surface.get();
    samples->patch_partition() = partition;
    samples->setup();
}
void setup_blended(unique_ptr<PatchSurfBlended>& surface, 
        unique_ptr<PatchSamples>& samples){
    surface = std::move(unique_ptr<PatchSurfBlended> (new PatchSurfBlended("BD3D_", "bd3d_")));
    surface->setFromOptions();
    surface->setup();

    vector<int> partition(surface->patches().size(), 0);
    samples = std::move(unique_ptr<PatchSamples>(new PatchSamples("", "")));
    samples->bdry() = surface.get();
    samples->patch_partition() = partition;
    samples->setup();
}

string build_prefix(unique_ptr<PatchSamples>& samples){
    string prefix = "PatchSamples/";
    prefix += Test::get_domain() + "/";
    if(dynamic_cast<PatchSurfFaceMap*>(samples->bdry())){
        prefix += "face_map/";
        if(Options::get_double_from_petsc_opts("-bd3d_facemap_adaptive")){
            assert(fabs(Options:: get_double_from_petsc_opts("-bd3d_facemap_fit_accuracy") - 1e-4) <=1e-16);
            prefix += "adaptive/";
        } else {
            prefix += "no_ref/";

        }
    } else if(dynamic_cast<PatchSurfBlended*>(samples->bdry())){
        prefix += "blended/";
    } else if(dynamic_cast<PatchSurfAnalytic*>(samples->bdry())){
        prefix += "analytic/";
    } else { 
        assert(0);//????
    }
    return prefix;
}

void dump_regression_data(unique_ptr<PatchSamples>& samples){
    string prefix = build_prefix(samples);

    Regression::dump(samples->sample_point_3d_position(), prefix+"sample_point_3d_position.reg");
    Regression::dump(samples->sample_point_normal(), prefix+"sample_point_normal.reg");
    Regression::dump(samples->sample_point_jacobian(), prefix+"sample_point_jacobian.reg");
    Regression::dump(samples->sample_point_quad_weight(), prefix+"sample_point_quad_weight.reg");
    Regression::dump(samples->sample_point_far_field(), prefix+"sample_point_far_field.reg");
    Regression::dump(samples->sample_point_blend_func_value(), prefix+"sample_point_blend_func_value.reg");
    Regression::dump(samples->sample_point_interpolant_spacing(), prefix+"sample_point_interpolant_spacing.reg");
    Regression::dump(samples->sample_point_combined_weight(), prefix+"sample_point_combined_weight.reg");
}

void test_vec(Vec vec, string file){
    Vec vec_computed;
    VecDuplicate(vec, &vec_computed);
    Regression::load(vec_computed, file);
    Regression::compare(vec_computed, vec);
}

void compare_to_regression_data(unique_ptr<PatchSamples>& samples){
    string prefix = build_prefix(samples);
    test_vec(samples->sample_point_3d_position(), prefix+"sample_point_3d_position.reg");
    test_vec(samples->sample_point_normal(), prefix+"sample_point_normal.reg");
    test_vec(samples->sample_point_jacobian(), prefix+"sample_point_jacobian.reg");
    test_vec(samples->sample_point_quad_weight(), prefix+"sample_point_quad_weight.reg");
    test_vec(samples->sample_point_far_field(), prefix+"sample_point_far_field.reg");
    test_vec(samples->sample_point_blend_func_value(), prefix+"sample_point_blend_func_value.reg");
    test_vec(samples->sample_point_interpolant_spacing(), prefix+"sample_point_interpolant_spacing.reg");
    test_vec(samples->sample_point_combined_weight(), prefix+"sample_point_combined_weight.reg");

}


TEST_CASE("Regression test initialization for PatchSamples", "[patch-samples][regression-init]"){

    SECTION("explicit torus, non-adaptive, no patch refinement"){
        cout << "Are you sure you want to regenerate regression data? Old state will be lost. [y to confirm]" << endl;
        char a;
        cin >> a;
        if(a != 'y'){
            exit(0);
        }
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/newtorus.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/newtorus.wrl");
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/explicit_torus_patches.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
        Options::set_value_petsc_opts("-bdtype", "2");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");

        unique_ptr<PatchSurfFaceMap> surface;
        unique_ptr<PatchSamples> samples;
        setup_face_map(surface, samples);
        REQUIRE(surface);
        REQUIRE(samples);
        dump_regression_data(samples);


    } SECTION("blended pipe"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".125");
        Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
        Options::set_value_petsc_opts("-bdtype", "1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

        unique_ptr<PatchSurfBlended> surface;
        unique_ptr<PatchSamples> samples;
        setup_blended(surface, samples);
        REQUIRE(surface);
        REQUIRE(samples);
        dump_regression_data(samples);


    }SECTION("face-map fitting blended cube"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
        Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
        Options::set_value_petsc_opts("-bdtype", "1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        
        

        unique_ptr<PatchSurfFaceMap> surface;
        unique_ptr<PatchSamples> samples;
        setup_face_map(surface, samples);
        dump_regression_data(samples);


    }
}

TEST_CASE("Regression test for PatchSamples", "[patch-samples][regression]"){
   
    SECTION("explicit torus, non-adaptive, no patch refinement"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/newtorus.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/newtorus.wrl");
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_meshes/poly/explicit_torus_patches.poly");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "3");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".1");
        Options::set_value_petsc_opts("-bdtype", "2");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        
    unique_ptr<PatchSamples> samples;
        unique_ptr<PatchSurfFaceMap> surface;
        setup_face_map(surface, samples);
        compare_to_regression_data(samples);
       
    } SECTION("blended pipe"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/pipe.wrl");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".125");
        Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
        Options::set_value_petsc_opts("-bdtype", "1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

    unique_ptr<PatchSamples> samples;
        unique_ptr<PatchSurfBlended> surface;
        setup_blended(surface, samples);
        compare_to_regression_data(samples);


    } 
     SECTION("face-map fitting blended cube"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
        Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
        Options::set_value_petsc_opts("-bdtype", "1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        
        

        unique_ptr<PatchSurfFaceMap> surface;
        unique_ptr<PatchSamples> samples;
        setup_face_map(surface, samples);
        compare_to_regression_data(samples);


    }
}
/*TEST_CASE("Regression test for PatchSamples", "[patch-samples][regression]"){
}*/
