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


Vec run_test_get_potential(PatchSurf* surface, 
        TestConfig test_type){
    Vec true_potential;
    Vec computed_potential;
    Vec targets;
    Vec boundary_data;
    Vec solved_density;

    auto solver = setup_solver(surface, test_type.solver_matvec_type);

    NumVec<OnSurfacePoint> closest_points;
    setup_target_data(solver, test_type, targets, closest_points, 
            computed_potential, true_potential);

    setup_problem_data(solver, test_type, boundary_data, solved_density);
    
    solve_and_evaluate(solver,
            boundary_data,
            targets,
            closest_points, solved_density, 
            computed_potential,
            test_type);


    VecDestroy(&true_potential);
    VecDestroy(&solved_density);
    VecDestroy(&boundary_data);
    VecDestroy(&targets);

    return computed_potential;
}

string build_prefix(unique_ptr<SolverGMRESDoubleLayer>& solver){
    string prefix = "SolverGMRESDoubleLayer/";
    if(dynamic_cast<PatchSurfFaceMap*>(solver->bdry())){
        prefix += "face_map/";
        if(Options::get_double_from_petsc_opts("-bd3d_facemap_adaptive")){
            assert(fabs(Options:: get_double_from_petsc_opts("-bd3d_facemap_fit_accuracy") - 1e-4) <=1e-16);
            prefix += "adaptive/";
        } else {
            prefix += "no_ref/";

        }
    } else if(dynamic_cast<PatchSurfBlended*>(solver->bdry())){
        prefix += "blended/";
    } else if(dynamic_cast<PatchSurfAnalytic*>(solver->bdry())){
        prefix += "analytic/";
    } else { 
        assert(0);//????
    }
    return prefix;
}
void dump_regression_data(vector<Vec> computed_data, vector<string> file_names){
    assert(computed_data.size() == file_names.size());
}




TEST_CASE("Test singular quad solver on analytic surfaces", "[solver][analytic]"){
    SECTION("Test single sphere Laplace interior"){
        Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
        Options::set_value_petsc_opts("-dom", "0"); // interior problem
        Options::set_value_petsc_opts("-bdtype", "0"); // analytic surface
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/sphere.wrl"); // single sphere
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/sphere.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".05");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".0625");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".03125");
        /*Options::set_value_petsc_opts("-bis3d_spacing", ".125");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".0625");
        //Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
        //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03125");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".0625");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".03125");
        */

        unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
        surface->setFromOptions();
        surface->setup();
        
        //test_constant_boundary_data(surface.get(), true);
        //test_constant_boundary_data(surface.get());

        TestConfig test;
        test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
        test.target_type   = TargetType::GRID;
        test.evaluation_scheme = EvaluationScheme::AUTOEVAL_QBKIX;
        test.solver_matvec_type = EvaluationType::EXTRAPOLATION_AVERAGE;
        //test.solution_scheme   = SolutionScheme::EXPLICIT_DENSITY;
        test.solution_scheme   = SolutionScheme::GMRES_SOLVE;

        test.dump_values = false;
        run_test(surface.get(),test);

    }
    /*SECTION("Test single sphere Stokes exterior"){
        Options::set_value_petsc_opts("-kt", "311"); // Stokes problem
        Options::set_value_petsc_opts("-dom", "1"); // exterior problem
        Options::set_value_petsc_opts("-bdtype", "0"); // analytic surface
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/sphereflip.wrl"); // single sphereflip
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/sphereflip.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", ".0625");
        Options::set_value_petsc_opts("-bis3d_rfdspacing", ".03125");
        
        unique_ptr<PatchSurfAnalytic> surface(new PatchSurfAnalytic("BD3D_", "bd3d_"));
        surface->setFromOptions();
        surface->setup();
        
        test_solver_constant_boundary_data(surface.get());

    }*/

}
/* The Rosetta Stone to get Blended solver convergence: NOTE deprecated
 rfdhs = {0.192:0.096, 0.096:0.048, 0.048:0.024, 0.024:0.008, 0.012:0.004, 0.006:0.002, 0.003:0.001}
 hsradmult = {0.192:4, 0.096:7, 0.048:9, 0.024:13, 0.012:17, 0.006:25, 0.003:35} 
 hsbdspacing = {0.192:0.032, 0.096:0.016, 0.048:0.008, 0.024:0.004, 0.012:0.002, 0.006:0.001, 0.003:0.0005}
*/
TEST_CASE("Test solver on blended surfaces", "[solver][blended]"){
    SECTION("Test cube Laplace interior"){
        Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
        Options::set_value_petsc_opts("-dom", "0"); // interior problem
        Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_meshes/wrl/cube.wrl"); // single cube
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_meshes/wrl/cube.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", ".192");

        // NOTE BUG figure out why blendsurf is so damn slow now
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        
        unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
        surface->setFromOptions();
        surface->setup();
        
        //test_constant_boundary_data(surface.get());

        TestConfig test;
        test.bc_type = BoundaryDataType::CONSTANT_DENSITY;
        test.target_type   = TargetType::GRID;
        test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
        //test.solution_scheme   = SolutionScheme::EXPLICIT_DENSITY;
        test.solution_scheme   = SolutionScheme::GMRES_SOLVE;

        test.dump_values = true;
        run_test(surface.get(),test);

    }

    SECTION("Test pipe"){

    }
    SECTION("Test propeller exterior only"){

    }

}
TEST_CASE("Test solver on face-map surfaces", "[solver][face-map]"){
    SECTION("Test cube Laplace interior"){
        Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
        Options::set_value_petsc_opts("-dom", "0"); // interior problem
        Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // single cube
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", ".2");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");
        Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        //Options::set_value_petsc_opts("-bd3d_facemap_fit_accuracy", "1e-2");
        //Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "2");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    
        
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".08");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".08");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    Options::set_value_petsc_opts("-qbkix_convergence_type","classic");
    //Options::set_value_petsc_opts("-bis3d_spacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", "0.14285");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".125");
    //Options::set_value_petsc_opts("-bis3d_rfdspacing", ".125");
    Options::set_value_petsc_opts("-bis3d_spacing", ".071428");
    Options::set_value_petsc_opts("-bis3d_rfdspacing", ".071428");
    Options::set_value_petsc_opts("-dump_qbkix_points", "1");
    Options::set_value_petsc_opts("-dom", "0");
    Options::set_value_petsc_opts("-bdtype", "2");
    Options::set_value_petsc_opts("-bis3d_np", "16");
    Options::set_value_petsc_opts("-bis3d_ptsmax", "5000");
    Options::set_value_petsc_opts("-bdsurf_interpolate", "0");

        
        unique_ptr<PatchSurfFaceMap> surface(new PatchSurfFaceMap("BD3D_", "bd3d_"));
        surface->_surface_type = PatchSurfFaceMap::BLENDED;
        surface->setFromOptions();
        surface->setup();
        surface->_coarse = true;
        surface->refine_test();
        
        test_constant_boundary_data(surface.get(),true);

    }
}

TEST_CASE("test evaluate()", "[solver][debug]"){
        Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
        Options::set_value_petsc_opts("-dom", "0"); // interior problem
        Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // single cube
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", "1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".5");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".5");

    unique_ptr<PatchSurfFaceMap> face_map(new PatchSurfFaceMap("BD3D_", "bd3d_"));
    face_map->_surface_type = PatchSurfFaceMap::BLENDED;
    face_map->_coarse = true;
    face_map->setFromOptions();
    face_map->setup();
    //face_map->refine();
    //face_map->refine_uniform(1);
    face_map->refine_test();

    
    unique_ptr<SolverGMRESDoubleLayer> solver(new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = face_map.get();
    //solver->_compute_refined_surface = false;
    int num_patches = face_map->patches().size();
    vector<int> patch_partition(num_patches, 0);  //All in one processor
    
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom");
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->eqcoefs() = vector<double>(2,1.0); 
    solver->_compute_refined_surface = true;
    
    solver->setFromOptions();
    solver->setup();
    auto k = solver->problem_kernel(); 
    int num_sources = Petsc::get_vec_size(solver->patch_samples()->sample_point_3d_position())/DIM;
    int num_targets = 2;
    Vec targets;
    Vec density;
    Vec potential;
    Petsc::create_mpi_vec(solver->mpiComm(), num_targets*DIM, targets);
    Petsc::create_mpi_vec(solver->mpiComm(), num_sources*k.get_sdof(), density);
    Petsc::create_mpi_vec(solver->mpiComm(), num_targets*k.get_tdof(), potential);
    VecSet(density,1.);
    VecSet(targets,0.);
    DblNumMat targets_local = get_local_vector(DIM, num_targets, targets);
    DblNumMat sample_position_local= 
        get_local_vector(DIM, num_sources, solver->patch_samples()->sample_point_3d_position());
    targets_local(0,1) = sample_position_local(0,20);
    targets_local(1,1) = sample_position_local(1,20);
    targets_local(2,1) = sample_position_local(2,20);
    targets_local.restore_local_vector();
    sample_position_local.restore_local_vector();
    //VecScale(targets, .95);
    NumVec<OnSurfacePoint> on_surface_points(num_targets);
    solver->evaluate(targets,
            density,
            potential,
            on_surface_points,
            VAR_U, 
            false);
}

TEST_CASE("Debug solver memory leak", "[solver][debug][memory]"){
        Options::set_value_petsc_opts("-kt", "111");
        Options::set_value_petsc_opts("-dom", "0"); 
        Options::set_value_petsc_opts("-bdtype", "2"); 
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); 
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
        Options::set_value_petsc_opts("-bis3d_spacing", ".1");
        Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
        Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "3");
        Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
        Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    unique_ptr<PatchSurfFaceMap> surface = setup_face_map(PatchSurfFaceMap::BLENDED, 0, "uniform");
    unique_ptr<SolverGMRESDoubleLayer> solver = setup_solver(surface.get(),  EXTRAPOLATION_AVERAGE, true);
    auto k = solver->problem_kernel(); 
    
    Vec boundary_data;
    Vec density;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, solver->local_total_dof(), density);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, solver->local_total_dof(), boundary_data);
    
    PetscRandom seed;
    PetscRandomCreate(MPI_COMM_WORLD, &seed);
    for (int i = 0; i < 100; i++) {
        VecSetRandom(density, seed);
        VecSetRandom(boundary_data, seed);
        solver->solve(boundary_data, density);
        double norm;
        VecNorm(density, NORM_2, &norm);
        cout <<  "norm: " << norm<< endl;
    }
    
}

TEST_CASE("Debug evaluation memory leak", "[eval][debug][memory]"){
    Options::set_value_petsc_opts("-kt", "111");
    Options::set_value_petsc_opts("-dom", "0"); 
    Options::set_value_petsc_opts("-bdtype", "2"); 
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); 
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".05");
    Options::set_value_petsc_opts("-upsampling_type", "uniform");
    Options::set_value_petsc_opts("-uniform_upsampling_num_levels", "1");
    
    unique_ptr<PatchSurfFaceMap> surface = setup_face_map(PatchSurfFaceMap::BLENDED, 0, "uniform");
    unique_ptr<SolverGMRESDoubleLayer> solver = setup_solver(surface.get(),  EXTRAPOLATION_AVERAGE, true);
    auto k = solver->problem_kernel(); 

    Vec potential;
    Vec density;
    Vec targets;
    int num_targets = 1000;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, solver->local_total_dof(), density);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, solver->local_total_dof(), potential);
    //Petsc::create_mpi_vec(MPI_COMM_WORLD, num_targets*DIM, targets);
    targets = solver->patch_samples()->sample_point_3d_position();
    NumVec<OnSurfacePoint> closest_points = solver->patch_samples()->sample_point_as_on_surface_point();


    PetscRandom seed;
    PetscRandomCreate(MPI_COMM_WORLD, &seed);
    for (int i = 0; i < 100; i++) {
        VecSetRandom(density, seed);
        VecSetRandom(potential, seed);
        VecSetRandom(targets, seed);
        // no closest point check, could be a leak there
        solver->evaluate(targets, density, potential, closest_points, VAR_U);
        double norm;
        VecNorm(potential, NORM_2, &norm);
        cout <<  "norm: " << norm<< endl;
    }

}
 
TEST_CASE("stampede test", "[solver][stampede]"){
    PatchSurfFaceMap* surface;
    SolverGMRESDoubleLayer* solver;
    MPI_Comm comm;
    Vec solved_density;
    Vec boundary_data;

    comm = MPI_COMM_WORLD;
    //PetscOptionsInsertFile(comm, NULL, "opt/morse_cases_vessel_debug.opt", PETSC_TRUE);
    Options::set_value_petsc_opts("-kt","311");
    Options::set_value_petsc_opts("-dom", "0"); // 1 exterior problem, 0 interior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/vessel_section_scaling_debug.wrl"); // small vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/vessel_section_scaling_debug.wrl");
        Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/two_cube.wrl"); // single cube
        Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/two_cube.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/pipe.wrl"); // small vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/pipe.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".2"); // for test
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".01666");
    //Options::set_value_petsc_opts("-bis3d_ptsmax", "1000000000");
    
    surface = new PatchSurfFaceMap("BD3D_", "bd3d_");
    surface->_surface_type = PatchSurfFaceMap::BLENDED;
    surface->setFromOptions();
    surface->setup();
    surface->_coarse = true;
    surface->refine_test();
    //surface->refine_uniform(2);
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    
    solver = new SolverGMRESDoubleLayer(surface);
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->_compute_refined_surface = true;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    solver->setFromOptions();
    solver->setup();
    //std::cout<<"solver setup\n";

    // boundary data
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            boundary_data);
    VecSet(boundary_data, 1.);
    
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            solved_density);
    VecSet(solved_density, 0.);

    solver->solve(boundary_data, solved_density);

    // clean
    if(surface)
        delete surface;
    if(solver)
        delete solver;
    Petsc::destroy_vec(boundary_data);
    Petsc::destroy_vec(solved_density);
}
