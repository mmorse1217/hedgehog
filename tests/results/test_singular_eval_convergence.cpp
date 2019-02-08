#include "../catch.hpp"
#include "common/utils.hpp"
#include "../utils/evaluation_utils.hpp"
#include "bdry3d/patch_surf_blended.hpp"
#include "bdry3d/patch_surf_face_map.hpp"
#include "bdry3d/patch_surf_analytic.hpp"

using namespace Ebi;
TEST_CASE("Test singular eval convergence", "[results][singular-eval]"){
    Options::set_value_petsc_opts("-kt", "111"); // Laplace problem
    Options::set_value_petsc_opts("-dom", "0"); // interior problem
    Options::set_value_petsc_opts("-bdtype", "1"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // single cube
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".096");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "1");

    unique_ptr<PatchSurfBlended> surface(new PatchSurfBlended("BD3D_", "bd3d_"));
    surface->setFromOptions();
    surface->setup();
    
    auto solver = setup_solver(surface.get(), SINGULAR_EVAL);
    /*
    test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
    test.solve = false;
    test.target_type = TargetType::GRID;
    test.problem_type = ProblemType::CONSTANT_DENSITY;

    run_test(surface.get(),test);
    
    //TestType test;
    test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
    test.solve = false;
    test.target_type = TargetType::COLLOCATION_POINTS;
    test.problem_type = ProblemType::CONSTANT_DENSITY;
    run_test(surface.get(),test,false);*/

    TestConfig test;
    test.bc_type = BoundaryDataType::HARMONIC;
    test.target_type   = TargetType::GRID;
    test.evaluation_scheme = EvaluationScheme::SMOOTH_QUAD;
    test.solution_scheme   = SolutionScheme::GMRES_SOLVE;
    test.solver_matvec_type= SINGULAR_EVAL;
    run_test(surface.get(),test);


}
