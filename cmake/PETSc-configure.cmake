cmake_minimum_required(VERSION 3.1)
project(petsc-download)
include(ExternalProject)
    # Compile PETSc 
    ExternalProject_Add(
        PETSc
        URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.14.5.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/_deps/PETSc
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/PETSc
        CONFIGURE_COMMAND cd <SOURCE_DIR> && ./configure --with-pic=0 --with-shared-libraries=0 --prefix=<INSTALL_DIR> --with-debugging=0 --with-64-bit-indices=1 --with-single-library=1 
        BUILD_COMMAND make -C <SOURCE_DIR>
        INSTALL_COMMAND make -C <SOURCE_DIR> install
    )
    
