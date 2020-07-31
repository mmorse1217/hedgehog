#include "../catch.hpp"
#include "../utils/evaluation_utils.hpp"
#include "../utils/regression_utils.hpp"
#include "bdry3d/patch_samples.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "common/stats.hpp"
#include "common/utils.hpp"
#include <vector>
using namespace hedgehog;
/****************************************************************************/
// Define configuration options for regression test cases 
/****************************************************************************/
void set_reg_options_patch_samples_blended_pipe() {
  Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/pipe.wrl");
  Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/pipe.wrl");
  Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
  Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
  Options::set_value_petsc_opts("-bis3d_spacing", ".125");
  Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
  Options::set_value_petsc_opts("-bdtype", "1");
  Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
  Options::set_value_petsc_opts("-bdsurf_interpolate", "0");
}

void set_reg_options_patch_samples_poly_torus() {
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
}
void set_reg_options_patch_samples_poly_cube() {

  Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl");
  Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
  Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
  Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
  Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
  Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".0125");
  Options::set_value_petsc_opts("-bdtype", "2");
  Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
}

/****************************************************************************/
// Explicit lists of vectors to store and compare against for regression tests
/****************************************************************************/
vector<Vec> get_regression_vectors(PatchSamples* samples){
  vector<Vec> computed_data = {
      samples->sample_point_3d_position(),
      samples->sample_point_normal(),
      samples->sample_point_jacobian(),
      samples->sample_point_quad_weight(),
      samples->sample_point_far_field(),
      samples->sample_point_blend_func_value(),
      samples->sample_point_interpolant_spacing(),
      samples->sample_point_combined_weight(),
  };
  return computed_data;
}


vector<string> get_regression_file_names(PatchSamples* samples){
  vector<string> file_names = {

      "sample_point_3d_position.reg",
      "sample_point_normal.reg",
      "sample_point_jacobian.reg",
      "sample_point_quad_weight.reg",
      "sample_point_far_field.reg",
      "sample_point_blend_func_value.reg",
      "sample_point_interpolant_spacing.reg",
      "sample_point_combined_weight.reg"

  };
  return file_names;
}



TEST_CASE("Regression test initialization for PatchSamples",
          "[patch-samples][regression-init]") {

  unique_ptr<PatchSurf> surface;
  unique_ptr<PatchSamples> samples;
  SECTION("explicit torus, non-adaptive, no patch refinement") {
    set_reg_options_patch_samples_poly_torus();

    Regression::setup_face_map(surface, samples);
  }
  SECTION("blended pipe") {
    set_reg_options_patch_samples_blended_pipe();
    Regression::setup_blended(surface, samples);
  }
  SECTION("face-map fitting blended cube") {
    set_reg_options_patch_samples_poly_cube();
    Regression::setup_face_map(surface, samples);
  }
  REQUIRE(surface);
  REQUIRE(samples);
  
  vector<Vec> computed_data = get_regression_vectors(samples.get());
  vector<string> file_names = get_regression_file_names(samples.get());

  Regression::dump_regression_data(computed_data, file_names, "PatchSamples");
}

TEST_CASE("Regression test for PatchSamples", "[patch-samples][regression]") {

  unique_ptr<PatchSamples> samples;
  unique_ptr<PatchSurf> surface;

  SECTION("explicit torus, non-adaptive, no patch refinement") {

    set_reg_options_patch_samples_poly_torus();
    Regression::setup_face_map(surface, samples);
  }
  SECTION("blended pipe") {
    set_reg_options_patch_samples_blended_pipe();
    Regression::setup_blended(surface, samples);
  }
  SECTION("face-map fitting blended cube") {
    set_reg_options_patch_samples_poly_cube();
    Regression::setup_face_map(surface, samples);
  }
  REQUIRE(surface);
  REQUIRE(samples);
  
  vector<Vec> computed_data = get_regression_vectors(samples.get());
  vector<string> file_names = get_regression_file_names(samples.get());

  Regression::compare_to_regression_data(computed_data, file_names, "PatchSamples");
}
