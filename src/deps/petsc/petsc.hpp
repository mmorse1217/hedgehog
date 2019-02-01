#pragma once

#include <memory>
#include <iostream>
#include <array>
#include <mutex>

#include <petsc.h>

// forward definitions:
namespace Petsc
{
class Vector;
}
#include <deps/petsc/vector.hpp>


namespace Petsc
{

class PetscContext
{
  public:
    PetscContext( int argc,
                  char** argv,
                  std::string filename,
                  std::string help )
    {
        /*int ac = argc;
        char** av = new char* [static_cast<size_t>( argc ) + 1];
        for ( int i = 0; i < argc; i++ ) {
            av[i] = new char[strlen( argv[i] ) + 1];
            std::copy( argv[i], argv[i] + strlen( argv[i] ) + 1, av[i] );
        }
        av[ac] = NULL;*/
        PetscInitialize( &argc, &argv, filename.c_str(), help.c_str() );
        comm_ = PETSC_COMM_WORLD;
    }

    PetscContext( const int argc, const char** argv )
    {
        int ac = argc;
        char** av = new char* [static_cast<size_t>( argc ) + 1];
        for ( int i = 0; i < argc; i++ ) {
            av[i] = new char[strlen( argv[i] ) + 1];
            std::copy( argv[i], argv[i] + strlen( argv[i] ) + 1, av[i] );
        }
        av[ac] = NULL;
        PetscInitialize( &ac, &av, PETSC_NULL, PETSC_NULL );
        comm_ = PETSC_COMM_WORLD;
    }

    MPI_Comm comm() const { return comm_; }

    int rank() const
    {
        int r;
        MPI_Comm_rank( this->comm_, &r );
        return r;
    }

    PetscContext( const PetscContext& ) = delete;

    ~PetscContext()
    {
        PetscFinalize();
    }

    // private:
    MPI_Comm comm_;
};
}
