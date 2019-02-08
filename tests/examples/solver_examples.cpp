
#include "../catch.hpp"
#include "bie3d/solver_gmres_double_layer.hpp"
using namespace Ebi;
TEST_CASE("Solver examples: blended surface", "[examples][solver][blended]"){
    // default options are listed in opt/test_cases.opt. Any options explcitly
    // set here have priority and override defaults
   /* 
    // Set up the solver
    vector<int> patch_partition(surface->patches().size(), 0);  //All in one processor
    unique_ptr<SolverGMRESDoubleLayer> solver( new SolverGMRESDoubleLayer("BIS3D_", "bis3d_"));
    solver->bdry() = surface;
    solver->patch_partition() = patch_partition;
    solver->dom() = Options::get_int_from_petsc_opts("-dom"); // TODO remove

    // specify which matvec to use inside GMRES
    solver->set_evaluation_type(eval_type); 
    solver->_compute_refined_surface = compute_refined_surface;

    // Pick the 
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    
    solver->setFromOptions();
    solver->setup(); // memory leak <<==========!!!!!
*/
}
TEST_CASE("Solver examples: face-map surface", "[examples][solver][face-map]"){

}
