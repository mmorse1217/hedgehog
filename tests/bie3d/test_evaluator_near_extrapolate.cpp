
#include "common/kernel3d.hpp"
#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_surf_analytic.hpp"
#include "bie3d/solver_utils.hpp"
#include <sampling.hpp>
#include "common/nummat.hpp"
#include "common/vtk_writer.hpp"
#include "../utils/evaluation_utils.hpp"
using namespace Ebi;

TEST_CASE("Test qbkix extrapolation of known functions", "[critical][eval][near][qbkix]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    SECTION("Test single flat patch, extrapolate analytic function"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/small_flat_patch.wrl"); // single small_flat_patch
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/small_flat_patch.wrl");
        Options::set_value_petsc_opts("-poly_coeffs_file", "wrl_files/poly/small_flat_patch.poly");

        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::POLYNOMIAL;
        surface->setFromOptions();
        surface->setup();
        surface->_coarse = true;
        surface->refine_test();

        test_qbkix_extrapolation_no_quad(surface.get(), eval_x, true);
        test_qbkix_extrapolation_no_quad(surface.get(), eval_exp_sin_sin, true);
    }
}

TEST_CASE("Test qbkix extrapolation of known functions slow", "[eval][near][qbkix]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    Options::set_value_petsc_opts("-bis3d_spacing", ".25");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".25");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "1");
    //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");

    SECTION("Test single curved flat patch with adaptive refinement, extrapolate analytic functions"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/gaussian.wrl"); // single small_flat_patch
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/gaussian.wrl");
        Options::set_value_petsc_opts("-analytic_function", "4");
        Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "6");
        Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-6");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-dump_qbkix_points", "1");
        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::ANALYTIC;
        surface->setFromOptions();
        surface->setup();
        surface->_coarse = true;
        surface->refine();

        test_qbkix_extrapolation_no_quad(surface.get(), eval_x, true);
        test_qbkix_extrapolation_no_quad(surface.get(), eval_exp_sin_sin, true);

    }

    SECTION("Test cube, extrapolate analytic function"){
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // single cube
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");

        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::BLENDED;
        surface->setFromOptions();
        surface->setup();
        surface->_coarse = true;
        surface->refine_test();

        Options::set_value_petsc_opts("-bis3d_spacing", ".25");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
        test_qbkix_extrapolation_no_quad(surface.get(), eval_x, false);
    }
}
