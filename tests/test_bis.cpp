//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "common/ebi.hpp"
#include "common/ebi_petsc.hpp"
#include "common/utils.hpp"
#include "deps/petsc/petsc.hpp"
#include <geogram/basic/common.h>
#include "bdry3d/geogram_interface.hpp"
using namespace std;
using namespace Ebi;

int main( int argc, char**  argv ) {

    Petsc::PetscContext petsc_context(argc, argv, "opt/morse_cases.opt", "Unit test");
    // Needs to be called once to initialize geogram
    GEO::initialize();
    
    srand48( (long)time(NULL) );
    MPI_Comm comm = petsc_context.comm();
    
    // everything will explode if we use mpirun
    int mpi_size;
    MPI_Comm_size(comm, &mpi_size);  
    //assert(mpi_size==1);

    cerr << "starting tests... " << endl;
    int result = Catch::Session().run( argc, argv );
    cerr << "finished tests! " << endl;

    return result;
}
