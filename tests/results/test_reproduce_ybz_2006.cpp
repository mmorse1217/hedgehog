
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
#include "common/stats.hpp"
using namespace hedgehog;

TEST_CASE("", "[ybz]"){
}
TEST_CASE("Reproduce table 1 in Ying, Biros, Zorin 2006 paper", "[ybz][table-1]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".5");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "16");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", ".032");

    
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "12");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");

    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    string output_folder = "output/test_reproduce_ybx_2006/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement = 5; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 


            /*check_potential(surface, EvaluationType::EXTRAPOLATION_AVERAGE,  
                    "qbx_", output_folder, i);
            */
            check_potential(surface, EvaluationType::SINGULAR_EVAL,  
                    "singular_eval_", output_folder, i);
            // half the surface spacing
            double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            double h_interp = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h_interp/2.)); 
            // and run it again

        }
    }
}
