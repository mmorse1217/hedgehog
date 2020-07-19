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
using namespace Ebi;

void check_density(unique_ptr<PatchSurfAnalytic>& surface, 
        EvaluationType solver_matvec_type,  string eval_type, 
        string output_folder, int i=0){
    int num_coarse_patches = surface->patches().size();
    stats.add_result("num coarse patches", num_coarse_patches);

    TestConfig test;
    // data generated from Y^1_1 i.e. (1,1)-spherical harmonic
    test.bc_type = BoundaryDataType::SPHERICAL_HARMONIC;

    // Evaluate density at the collocation points
    test.target_type   = TargetType::COLLOCATION_POINTS;

    // Solve with GMRES using solver_matvec_type
    test.solver_matvec_type = solver_matvec_type;
    test.solution_scheme = SolutionScheme::GMRES_SOLVE;

    // compare density against true spherical harmonic solution
    test.evaluation_scheme = EvaluationScheme::CHECK_DENSITY;

    string filename = eval_type + string("_ref_lvl_")+to_string(i);
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;

    test.time_coarse_fmm = true;
    test.dump_values = true;
    run_test(surface.get(),test);

    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.dump_key_values("", stats._file_prefix);
    stats.print_results();
    stats.clear();

}

BEGIN_EBI_NAMESPACE
void check_potential(PatchSurf* surface, 
        EvaluationType solver_matvec_type,  SolutionScheme solution_scheme, string eval_type, 
        string output_folder, int i, BoundaryDataType boundary_condition =BoundaryDataType::CONSTANT_DENSITY){
    int num_coarse_patches = surface->patches().size();
    stats.add_result("num coarse patches", num_coarse_patches);

    TestConfig test;
    test.bc_type = boundary_condition;

    if(boundary_condition == BoundaryDataType::HARMONIC){
        test.singularity_type= SingularityType::SPHERE;
        //test.sphere_radius_bc = .89;
        test.sphere_radius_bc = 1.;
    } else if (boundary_condition != BoundaryDataType::CONSTANT_DENSITY){
        assert(0);
    }
    
    // Solve with GMRES using solver_matvec_type
    test.solver_matvec_type = solver_matvec_type;
    
    test.solution_scheme = solution_scheme;

    if(solution_scheme == SolutionScheme::GMRES_SOLVE){
        // solve, then evaluate far field error with smooth quadrature
        test.target_type   = TargetType::SPHERE_FAR;
        test.sphere_radius_targets = .2;
        test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
        /*test.target_type   = TargetType::COLLOCATION_POINTS;
        test.evaluation_scheme = EvaluationScheme::ON_QBKIX;*/

    } else if(solution_scheme == SolutionScheme::EXPLICIT_DENSITY){
        // we know the density exactly, test singular quadrature error
        test.target_type   = TargetType::COLLOCATION_POINTS;

        if(solver_matvec_type == EvaluationType::EXTRAPOLATION_AVERAGE){
            test.evaluation_scheme = EvaluationScheme::ON_QBKIX;
        } else if(solver_matvec_type == EvaluationType::SINGULAR_EVAL){
            test.evaluation_scheme = EvaluationScheme::ON_SINGULAR;
        } else {
            assert(0);
        }

    } else {
        assert(0);
    }


    string filename = eval_type + string("_ref_lvl_")+to_string(i);
    string full_path= output_folder + filename + string(".vtp");


    stats._file_prefix = output_folder + filename;
    test.error_filename = full_path;

    test.dump_values = true;
    test.time_coarse_fmm = true;
    run_test(surface,test);

    stats.add_result("dist to boundary", Options::get_double_from_petsc_opts("-boundary_distance_ratio"));
    stats.add_result("check pt spacing", Options::get_double_from_petsc_opts("-interpolation_spacing_ratio"));
    stats.add_result("coarse spacing", Options::get_double_from_petsc_opts("-bis3d_spacing"));
    stats.add_result("upsampled spacing", Options::get_double_from_petsc_opts("-bis3d_rfdspacing"));
    stats.add_result("num check points", Options::get_double_from_petsc_opts("-near_interpolation_num_samples"));
    stats.add_result("fft refinement factor", Options::get_double_from_petsc_opts("-dnref"));
    stats.dump_key_values("", stats._file_prefix);
    if(solution_scheme == SolutionScheme::GMRES_SOLVE){
    cout << "T_"<< eval_type <<  "/T_FMM: " << stats._results["total matvec time"]/stats._results["FMM benchmark time"]/stats._results["number of GMRES iterations"] << endl;
    } else {
    cout << "T_"<< eval_type <<  "/T_FMM: " << stats._results["total matvec time"]/stats._results["FMM benchmark time"] << endl;
    }
    stats.print_results(); 
    stats.clear();

}


END_EBI_NAMESPACE
TEST_CASE("Test singular quad vs. qbkix solver on analytic surfaces", 
        "[results][qbx-vs-singular-quad][density]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "0"); // analytic surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/sphere.wrl"); // single sphere
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/sphere.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03125");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".15");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".5");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".08333");
    Options::set_value_petsc_opts("-qbkix_convergence_type", "classic");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");
    
    // These give 1e-7 abs error...
    //Options::set_value_petsc_opts("-bis3d_spacing", ".03125");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0078125");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    //double alpha = .42857099;
    //double alpha = .259; // best value for 6x upsampling
    //double alpha = 0.421052;
    double alpha= 0.210526902;
    //double alpha = .6;
    double check_pt_h = alpha*sqrt(h);///1.5;
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h));

    unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    string output_folder = "output/test_qbkix_vs_singular_quad/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =5;
        for(int i = 0; i < num_levels_refinement; i++){
            /*
             * QBKIX cost: Let N_c be the number of coarse points and N_u be
             * the number upsampled points. L is the order of the interpolant on
             * one side of the boundary
             *
             * QBKIX cost is:
             * density interpolation (2d barycentric) from N_c points evaluated at N_u points + 
             * PvFMM on N_u sources and 2*L*N_c targets + 
             * N_c * L + L^2 for barycentric extrapolation
             *
             * Singular quadrature cost: Let N_c be the number of coarse points,
             * and N_u be the number of upsampled points. Let k be the FFT
             * upsampling amount along 1 dimension
             *
             * Singular quadrature cost is:
             * density interpolation from N_c points to k^2*N_c points via 2D FFT +
             * compute bspline interpolant from  k^2*N_c points evaluated at N_c\sqrt(N_c) points  + 
             * evaluating the \sqrt(N_c) singular quadrature at N_c points +
             * PvFMM on N_c sources and N_c targets + 
             * 2*N_c to compute the correction
             *
             * We need dump:
             * N_c  - num coarse points
             * N_s  - num upsampled points
             * k - FFT refinement factor (should be dnref)
             * L - qbkix interpolant order
             *
             * The rest of the cost can be compared in python; it's the sum of
             * the terms listed above.
             */

            unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 


            check_density(surface, EvaluationType::EXTRAPOLATION_AVERAGE,  
                    "qbx_", output_folder, i);
            check_density(surface, EvaluationType::SINGULAR_EVAL,  
                    "singular_eval_", output_folder, i);
            
            // half the surface spacing
            double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            
            double check_pt_h = alpha*sqrt(h/2.);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h));

            // and run it again

        }
    }
}


TEST_CASE("Test singular quad vs. qbkix eval on blendsurf laplace", 
        "[results][qbx-vs-singular-quad][laplace]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".075");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".6");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".15");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_pouctrl","1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");


    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    //double alpha = .42857099; // 36 x upsampling
    double alpha = 1.12105;
//1.12105, 0.186842
    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    //double check_pt_h = .928*sqrt(h);
    //double check_pt_h = 1./4.*(h); // initial spacing is too coarse for \sqrt(h)
    double check_pt_h = alpha*sqrt(h);
    double tt = 6.;
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));
    Options::set_value_petsc_opts("-pou_radius_constant", "1.");
    
    /*
    double rad; 
    radmult_spacing(h, rad) ;
    cout << "h: " << h << ", rad: " << rad << ", " << h*rad << ", "  << log2(h*rad) << endl;
    cout << "max level: " << -int(log2(h*rad)) << endl;
    Options::set_value_petsc_opts("-bis3d_maxlevel", to_string(-int(log2(h*rad)))); 
    */

    string output_folder = "output/test_qbkix_vs_singular_quad/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =5; 
        //int num_levels_refinement =3; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 

            check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::EXPLICIT_DENSITY,
              "qbx_", output_folder, i);
            check_potential(surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::EXPLICIT_DENSITY,
                    "singular_eval_", output_folder, i);

            // half the surface spacing
            h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/2./10));
            /* 
            radmult_spacing(h, rad) ;
            cout << "h: " << h << ", rad: " << rad << ", " << h*rad << ", "  << log2(h*rad) << endl;
            cout << "max level: " << -int(log2(h*rad)) << endl;
            Options::set_value_petsc_opts("-bis3d_maxlevel", to_string(-int(log2(h*rad)))); 
            */
            double check_pt_h = alpha*sqrt(h/2.);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));

            // and run it again

        }
    }
}

TEST_CASE("Test singular quad vs. qbkix eval on blendsurf navier constant density", 
        "[results][qbx-vs-singular-quad][navier-const-den]"){
    Options::set_value_petsc_opts("-kt", "511"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".075");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_pouctrl","1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");


    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    //double alpha = .42857099; // 36 x upsampling
    //double alpha = 0.2691271301;
    //double alpha = 1.0157/5.;
    //double alpha = 1.2;
    double alpha = .8;
    //alpha = 1.0157/6.;
    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    //double check_pt_h = .928*sqrt(h);
    //double check_pt_h = 1./4.*(h); // initial spacing is too coarse for \sqrt(h)
    double check_pt_h = alpha*sqrt(h);
    double tt = 6.;
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));
    Options::set_value_petsc_opts("-pou_radius_constant", "1.0");
    
    /*
    double rad; 
    radmult_spacing(h, rad) ;
    cout << "h: " << h << ", rad: " << rad << ", " << h*rad << ", "  << log2(h*rad) << endl;
    cout << "max level: " << -int(log2(h*rad)) << endl;
    Options::set_value_petsc_opts("-bis3d_maxlevel", to_string(-int(log2(h*rad)))); 
    */

    string output_folder = "output/test_qbkix_vs_singular_quad/const_density/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =5; 
        //int num_levels_refinement =3; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 

            check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::EXPLICIT_DENSITY,
              "qbx_", output_folder, i);
            check_potential(surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::EXPLICIT_DENSITY,
                    "singular_eval_", output_folder, i);

            // half the surface spacing
            h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/2./10));
            /* 
            radmult_spacing(h, rad) ;
            cout << "h: " << h << ", rad: " << rad << ", " << h*rad << ", "  << log2(h*rad) << endl;
            cout << "max level: " << -int(log2(h*rad)) << endl;
            Options::set_value_petsc_opts("-bis3d_maxlevel", to_string(-int(log2(h*rad)))); 
            */
            double check_pt_h = alpha*sqrt(h/2.);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));

            // and run it again

        }
    }
}


TEST_CASE("Test singular quad vs. qbkix solver on blendsurf navier", 
        "[results][qbx-vs-singular-quad][navier]"){




    Options::set_value_petsc_opts("-kt", "511"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl"); // blob
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".075");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".2");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03333");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".02");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
    //Options::set_value_petsc_opts("-bis3d_ksp_max_it", "100");

    Options::set_value_petsc_opts("-pou_radius_constant", "1.0");

    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
    //Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    //double alpha = .42857099;
    //double alpha = .183333; // upsampling 36x
    //double alpha = 0.2263153093; // upsampling 4x
    double alpha = 1.0157/4.;
    double tt = 6.;
    double check_pt_h = alpha*sqrt(h);
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));

    string output_folder = "output/test_qbkix_vs_singular_quad/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =5; 
        //int num_levels_refinement =4; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 

           check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::GMRES_SOLVE,
                    "qbx_", output_folder, i, BoundaryDataType::HARMONIC);
            if(i < 5){
                check_potential(surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::GMRES_SOLVE,
                        "singular_eval_", output_folder, i, BoundaryDataType::HARMONIC);
            }

            // half the surface spacing
            h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/2./10));

            double check_pt_h = alpha*sqrt(h/2.);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));

            // and run it again

        }
    }
}

TEST_CASE("Test singular quad vs. qbkix solver on blendsurf laplace solve", 
        "[results][qbx-vs-singular-quad][laplace-solve]"){




    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl"); // blob
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".075");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".2");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03333");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".02");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
    //Options::set_value_petsc_opts("-bis3d_ksp_max_it", "100");

    Options::set_value_petsc_opts("-pou_radius_constant", "1.0");

    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
    //Options::set_value_petsc_opts("-qbkix_convergence_type","adaptive");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    //double alpha = .42857099;
    //double alpha = .183333; // upsampling 36x
    //double alpha = 0.2263153093; // upsampling 4x
    double alpha = 1.12105/5.;
    alpha = 1.12105/7.; // 1e-10
    alpha = 1.12105/6.5; 
    //1.12105, 0.186842gt
    double tt = 6.;
    double check_pt_h = alpha*sqrt(h);
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));

    string output_folder = "output/test_qbkix_vs_singular_quad/solve/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";
    SECTION("test solver convergence using QBKIX for singular evaluation"){
        int num_levels_refinement =5; 
        //int num_levels_refinement =4; 
        for(int i = 0; i < num_levels_refinement; i++){

            unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 

           check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::GMRES_SOLVE,
                    "qbx_", output_folder, i, BoundaryDataType::HARMONIC);
            if(i < 5){
                check_potential(surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::GMRES_SOLVE,
                        "singular_eval_", output_folder, i, BoundaryDataType::HARMONIC);
            }

            // half the surface spacing
            h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            double h_up = Options::get_double_from_petsc_opts("-bis3d_rfdspacing");
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(h/2.));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(h_up/2.)); 
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/2./10));

            double check_pt_h = alpha*sqrt(h/2.);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h/tt));

            // and run it again

        }
    }
}


TEST_CASE("Test singular quad on blendsurf vs. qbkix on face-map", 
        "[results][qbx-vs-singular-quad][mixed-bdry]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    //Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "26");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive","0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");

    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-pou_radius_constant", "1.1");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));

    double blended_spacing = .35;
    double face_map_spacing= .090909;
    string output_folder = "output/test_qbkix_vs_singular_quad/face_map_vs_blended/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";

    SECTION("test solver convergence using QBKIX for singular evaluation"){
            Options::set_value_petsc_opts("-bdtype", "2"); // face-map surface
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
            Options::set_value_petsc_opts("-bdsurf_interpolate", "0");
            
            unique_ptr<PatchSurfFaceMap> face_map_surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
            face_map_surface->_surface_type = PatchSurfFaceMap::BLENDED;
            face_map_surface->_coarse = true;
            face_map_surface->setFromOptions();
            face_map_surface->setup();
        int num_levels_refinement =5; 
        //int num_levels_refinement =1; 
        for(int i = 0; i < num_levels_refinement; i++){
            Options::set_value_petsc_opts("-bdtype", "1"); // blended surface

            Options::set_value_petsc_opts("-bis3d_spacing", to_string(blended_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(blended_spacing/2.));
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(blended_spacing/10.));
            Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");

            unique_ptr<PatchSurfBlended> blended_surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            blended_surface->setFromOptions();
            blended_surface->setup();

            check_potential(blended_surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::GMRES_SOLVE,
                    "singular_eval_", output_folder, i, BoundaryDataType::HARMONIC);
            
            if(i > 0){
            face_map_surface->refine_uniform(1);
            } else { face_map_surface->refine_test();}

            Options::set_value_petsc_opts("-bdtype", "2"); // face-map surface
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
            Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

            double check_pt_spacing = .075;
            //double check_pt_spacing = .135;
            double tt = 6.;
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_spacing));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_spacing/tt));

            check_potential(face_map_surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::GMRES_SOLVE,
                    "qbx_", output_folder, i, BoundaryDataType::HARMONIC);

            // half the blended surface spacing
            blended_spacing /= 2.;
            
            // and run it again

        }
    }
}

TEST_CASE("Test singular quad on blendsurf vs. qbkix on face-map, navier", 
        "[results][qbx-vs-singular-quad][mixed-bdry-navier]"){
    Options::set_value_petsc_opts("-kt", "511"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "16");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive","0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");

    Options::set_value_petsc_opts("-bis3d_ptsmax", "250");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-pou_radius_constant", "1.1");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/10));

    double blended_spacing = .35;
    //double face_map_spacing= .090909;
    double face_map_spacing= .07142857;

    string output_folder = "output/test_qbkix_vs_singular_quad/face_map_vs_blended/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/";

    SECTION("test solver convergence using QBKIX for singular evaluation"){
            Options::set_value_petsc_opts("-bdtype", "2"); // face-map surface
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
            Options::set_value_petsc_opts("-bdsurf_interpolate", "0");
            
            unique_ptr<PatchSurfFaceMap> face_map_surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
            face_map_surface->_surface_type = PatchSurfFaceMap::BLENDED;
            face_map_surface->_coarse = true;
            face_map_surface->setFromOptions();
            face_map_surface->setup();
        int num_levels_refinement =5; 
        //int num_levels_refinement =1; 
        for(int i = 0; i < num_levels_refinement; i++){
            Options::set_value_petsc_opts("-bdtype", "1"); // blended surface

            Options::set_value_petsc_opts("-bis3d_spacing", to_string(blended_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(blended_spacing/2.));
            Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(blended_spacing/10.));
            Options::set_value_petsc_opts("-bdsurf_interpolate", "1");
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");

            unique_ptr<PatchSurfBlended> blended_surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            blended_surface->setFromOptions();
            blended_surface->setup();

            check_potential(blended_surface.get(), EvaluationType::SINGULAR_EVAL,  SolutionScheme::GMRES_SOLVE,
                    "singular_eval_", output_folder, i, BoundaryDataType::HARMONIC);
            
            if(i > 0){
            face_map_surface->refine_uniform(1);
            } else { face_map_surface->refine_test();}

            Options::set_value_petsc_opts("-bdtype", "2"); // face-map surface
            Options::set_value_petsc_opts("-bis3d_spacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bis3d_rfdspacing", to_string(face_map_spacing));
            Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
            Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

            double check_pt_spacing = .08;
            double tt = 6.;
            //double check_pt_spacing = .075;
            //double check_pt_spacing = .135;
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_spacing));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_spacing/tt));

            check_potential(face_map_surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::GMRES_SOLVE,
                    "qbx_", output_folder, i, BoundaryDataType::HARMONIC);

            // half the blended surface spacing
            blended_spacing /= 2.;
            
            // and run it again

        }
    }
}



TEST_CASE("Parameter qbkix error on sphere", 
        "[param][analytic]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "0"); // analytic surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/sphere.wrl"); // single sphere
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/sphere.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".5");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0833");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    Options::set_value_petsc_opts("-qbkix_convergence_type", "classic");
    Options::set_value_petsc_opts("-bis3d_ksp_max_it" ,"50");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-pou_radius_constant", "1.0");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    double alpha = 1.;
    double check_pt_h = alpha*sqrt(h);///1.5;
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h));

    unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    string output_folder = "output/test_qbkix_vs_singular_quad/";
    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/param-sweep/";
    SECTION("test param sweep QBKIX on sphere"){
        int iter=20;
        double step = 1./double(iter-1);
        double alpha = 1.;
        for(int i = 0; i < iter-1; i++){
            

            unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 


            check_density(surface, EvaluationType::EXTRAPOLATION_AVERAGE,  
                    "qbx_", output_folder, i);
            
            // half the surface spacing
            double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
            alpha -= step; 
            //double check_pt_h=alpha;//*sqrt(h);
            
            double check_pt_h = alpha*sqrt(h);
            Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(check_pt_h));
            Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(check_pt_h));

            // and run it again

        }
    }
}

TEST_CASE("Parameter qbkix error on blendsurf cube", 
        "[qbkix][param][solver][blended]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    //Options::set_value_petsc_opts("-kt", "511"); // not a Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "3");

    Options::set_value_petsc_opts("-bd3d_facemap_adaptive","0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor","1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".015625");
    Options::set_value_petsc_opts("-bis3d_spacing", ".3");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".075");
    //Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    Options::set_value_petsc_opts("-bd3d_facemap_patch_order", "20");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-bis3d_ksp_max_it", "25");


    // Best check point spacing for blendsurf: .928*\sqrt(h)
//  .910526 .227632 .019
//1.01579 .169298 .0182
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-bis3d_np", "10");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    double alpha = 1.2;
    double beta = alpha/6.;
    double check_pt_h = alpha*sqrt(h);
    Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(alpha*sqrt(h)));
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(beta*sqrt(h)));

    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    string output_folder = "output/test_qbkix_vs_singular_quad/";

    output_folder += Test::get_domain() + "/" + Test::get_kernel() + "/param-sweep/";
    int iter=20;
    double step = .2/double(iter-1);
    for(int i = 0; i < iter-1; i++){
        unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
        surface->setFromOptions();
        surface->setup();
        /*
           stats.print_results(); 
           unique_ptr<PatchSurfFaceMap> face_map_surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
           face_map_surface->_surface_type = PatchSurfFaceMap::BLENDED;
           face_map_surface->_coarse = true;
           face_map_surface->setFromOptions();
           face_map_surface->setup();
           */


        check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::GMRES_SOLVE,
                "qbx_", output_folder, i, BoundaryDataType::HARMONIC);
            //check_potential(surface.get(), EvaluationType::EXTRAPOLATION_AVERAGE,  SolutionScheme::EXPLICIT_DENSITY,
              //"qbx_", output_folder, i);

        // half the surface spacing
        double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
        alpha -= step; 
        double check_pt_h=alpha*sqrt(h);
        double beta = alpha/6.;
            cout << "iteration: " << i << "sqrt(h), alpha, beta: " << sqrt(h) << ", " << alpha << ", " << beta << endl;
        Options::set_value_petsc_opts("-boundary_distance_ratio", to_string(alpha*sqrt(h)));
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", to_string(beta*sqrt(h)));

        // and run it again

    }
}



void dump(Vec targets, Vec density, Vec potential){
        TestConfig test = {"output/error.vtp"};
    tear_down(targets, density, potential, test, 1, 1e-6, 1e-6, true);

}
TEST_CASE("Test singular quad with constant density", 
        "[debug]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03125");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
            
    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
            surface->setFromOptions();
            surface->setup();
            stats.print_results(); 

    vector<int> patch_partition(surface->patches().size(), 0);  //All in one processor
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = surface.get();
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom");
    solver->set_evaluation_type(EvaluationType::EXTRAPOLATION_AVERAGE);
    solver->_compute_refined_surface = true;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .5;
    
    solver->setFromOptions();
    solver->setup(); // memroy leak <<==========!!!!!

    Vec targets;
    VecDuplicate(solver->patch_samples()->sample_point_3d_position(),
            &targets);
    VecCopy(solver->patch_samples()->sample_point_3d_position(),
            targets);
    Vec targets_face_point;
    VecDuplicate(solver->patch_samples()->sample_as_face_point(),
            &targets_face_point);
    VecCopy(solver->patch_samples()->sample_as_face_point(),
            targets_face_point);
    Vec density;
    Vec potential;
    Petsc::create_mpi_vec(solver->mpiComm(), 
            Petsc::get_vec_size(targets)/DIM,
            density);
    VecDuplicate(density, &potential);
    VecSet(density, 1.);
        solver->roneval(targets,
                     targets_face_point,
                     VAR_U, density, potential);
        Vec t;
        VecDuplicate(density, &t);
        VecSet(t, .5);
    Vec error;
    VecDuplicate(t, &error);
    VecCopy(t, error);
    int minus_one = -1.;
    VecAXPY( error, minus_one,  potential);
    VecAbs(error);
    DblNumMat error_local = get_local_vector(1, Petsc::get_vec_size(error), error);
    for(int i = 0; i < Petsc::get_vec_size(error); i++){
        error_local(0,i) = error_local(0,i) < 1e-16 ? 1e-16 : error_local(0,i);
    }
    error_local.restore_local_vector(); 
    VecLog(error); // WARNING CHANGE TO BASE 10
    VecScale(error, 1./log(10.));
    double m;
    VecMax(error, NULL,&m);
    cout << "max error: " << m << endl;
    //dump(targets, t, potential);
}

TEST_CASE("Test blendsurf timing", "[debug][blendsurf][timing]"){
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // blob
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
    Options::set_value_petsc_opts("-LL", "4");
    Options::set_value_petsc_opts("-dnref", "10");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "1");

    double h = Options::get_double_from_petsc_opts("-bis3d_spacing");
    Options::set_value_petsc_opts("-bdsurf_interpolant_spacing", to_string(h/4));
    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    stats.start_timer("blended setup");
    surface->setup();
    stats.stop_timer("blended setup");
   
    int num_iter = 20;
    stats.start_timer("blended eval");
    for(auto const& patch: surface->patches()){
        for (int i = 0; i < num_iter; i++) {
            for (int j = 0; j < num_iter; j++) {
                //Point2 xy(i/double(num_iter),j/double(num_iter));
                double jac;
                Point3 ret[3];
                patch->estimate_jacobian(&jac);
                //patch->xy_to_patch_coords(xy.array(), PatchSamples::EVAL_FD|PatchSamples::EVAL_VL,  (double*)ret);
                //jac += 1.;
            }
        }
    }
    stats.stop_timer("blended eval");
    stats.add_result("num_evaluations", num_iter);

    stats.print_results(); 
    stats.dump_key_values("", "output/blended_interp_timing.vtp");
    //stats.dump_key_values("", "output/blended_eval_timing.vtp");

}

TEST_CASE("test pou", "[pou]"){
    int POU_type = 3;
    int POU_degree= 5;
    int num_iter = 500;
    double eps = 1e-3;
    DblNumVec t(num_iter);
    DblNumVec POU(num_iter);

    for(int p = 1; p < 11; p++){
        eps = pow(10,-p);
        // test that for all all t \in [0, \eps], pou(t) = 1
        for(int i =0; i < num_iter; i++){
            t(i) = eps*double(i)/double(num_iter-1);
            double a;  
            pou1d(PatchSamples::EVAL_VL, t(i), eps, 1.0-eps, &a, POU_type, POU_degree);
            CHECK(a == 1.); // not floating point, we want bitwise exact value of 1
        }

        // test that for all all t \in [1-\eps, 1], pou(t) = 0 
        for(int i =0; i < num_iter; i++){
            t(i) = eps*double(i)/double(num_iter-1) + 1.;
            double a;  
            pou1d(PatchSamples::EVAL_VL, t(i), eps, 1.0-eps, &a, POU_type, POU_degree);
            CHECK(a == 0.); // not floating point, we want bitwise exact value of 0
        }
    }


    for(int i =0; i < num_iter; i++){
        t(i) = double(i)/double(num_iter-1);
    }
    for(int i =0; i < num_iter; i++){
        double a;  
        double eps = 1e-3;//pow(10,-i);
        pou1d(PatchSamples::EVAL_VL, t(i), eps, 1.0-eps, &a, POU_type, POU_degree);
        POU(i) =a;
   }
    ofstream f;
    f.open("pou_func_values.txt");
    for(int i =0; i < num_iter; i++){
    f << t(i) ;
        f << " ";
    }
    for(int i =0; i < num_iter; i++){
    f << POU(i) << " ";
    }
    f.close();


}

