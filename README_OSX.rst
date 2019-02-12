* The following instruction works on Mac OSX Yosemite and they should work with some variation for other Macs with a later OS version. Hopefully.
* Depending on the version of software already installed on your machine, you can skip some of the steps below.

Pre build
-----------
::

  $ hg clone https://subversive.cims.nyu.edu/geonum/kifmm/mobo-temp
  $ cd mobo-temp
  $ export MOBO_DIR=`pwd`
  $ mkdir path/to/libs/

Install necessary dependecies
-----------------------------
:: 
  $ brew install gcc47

Download and Install MPICH 1.2.1
--------------------------------
:: 
  $ cd path/to/libs/
  $ wget http://www.mpich.org/static/downloads/1.2.1/mpich2-1.2.1.tar.gz
  $ tar xvf mpich2-1.2.1.tar.gz -C .
  $ cd mpich2-1.2.1
  $ ./configure CC=gcc-4.7 CXX=g++-4.7 F90=gfortran-4.7 --enable-sharedlibs=gcc -prefix=path/to/libs/mpich2-1.2.1
  $ make 
  $ make install 
  $ export PATH=path/to/libs/mpich2-1.2.1/bin:$PATH
  $ rm ../mpich2-1.2.1.tar.gz
  $ touch /.mpd.conf && chmod 600 .mpd.conf

Download and Build PETSC
-----------
::
  $ cd path/to/libs/
  $ mkdir petsc-3.7.2
  $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.2.tar.gz
  $ tar xvf petsc-3.7.2.tar.gz -C .
  $ export MPI_HOME=path/to/libs/mpich2-1.2.1
  $ cd petsc-3.7.2/
  $ ./configure --with-cc=gcc-4.7 --with-cxx=g++-4.7 --with-fc=gfortran-4.7  PETSC_ARCH=mac-osx PETSC_DIR=path/to/libs/petsc-3.7.2 --with-debugging=1 --with-64-bit-indices=1 --with-single-library=1  --with-blas-lapack-dir=/usr/lib/ --with-mpi-dir=$MPI_HOME
  # NOTE previous step assumes you're using the OSX installed BLAS/LAPACK located in /usr/lib/
  $ make all test
  $ export PETSC_ARCH=mac-osx
  $ export PETSC_DIR=path/to/libs/petsc-3.7.2
  $ rm ../petsc-3.7.2.tar.gz

Download and Build FFTW
----------
::
  $ cd /path/to/libs
  $ wget http://www.fftw.org/fftw-3.3.4.tar.gz 
  $ mkdir fftw-3.3.4
  $ export FFTW_DIR=path/to/libs/fftw-3.3.4
  $ cd fftw-3.3.4 
  $ ./configure --prefix=`pwd`  --enable-openmp
  $ make
  $ make install
  $ rm ../fftw-3.3.4.tar.gz

Download and build PvFMM
------------------------
::
  $ cd /path/to/libs
  $ git clone https://github.com/dmalhotra/pvfmm
  $ cd pvfmm
  $ autoreconf --install
  $ ./configure --with-fftw=${FFTW_DIR} CC=gcc-4.7 MPICC=mpicc F77=gfortran-4.7
  # NOTE may or may not need to specify --prefix here.
  $ #if ./configure complains or doesn't complete successfully, check that
    dependencies like gcc, gfortran, FFTW, BLAS, MPI, etc. are properly linked
    in the script (see https://github.com/dmalhotra/pvfmm/blob/develop/INSTALL
    for proper flags to override defaults).
  $ make
  $ make install 
  $ export PVFMM_DIR=path/to/libs/pvfmm
  $ export PVFMM_INC= ${PVFMM_DIR}/include
  $ export PVFMM_LIB= ${PVFMM_DIR}/lib

Download and Build p4est 
----------
::
  $ cd /path/to/libs
  $ wget http://p4est.github.io/release/p4est-1.1.tar.gz
  $ tar xvf p4est-1.1.tar.gz
  $ cd p4est-1.1
  $ ./configure CC=gcc-4.7 FC=gfortran-4.7 F77=gfortran-4.7 --prefix=`pwd`
  $ make
  $ make install
  $ export P4EST_DIR=`pwd`


Note: blendsurf needs GLUT and apple's Accelerate framework and blendsurf and face-map need OpenGl 
Might need some path tweaking to link properly depending on version :( (see
makefile.in in each directory)
Build Blendsurf
---------------
:: 
  $ cd ${MOBO_DIR}/blendsurf3
  $ make
  $ export BLENDUSRF_DIR=${MOBO_DIR}/blendsurf3

Build Face-map 
---------------
:: 
  $ cd ${MOBO_DIR}/face_map
  $ make
  $ export FACEMAP_DIR=${MOBO_DIR}/face_map

Build Hedgehog 
----------
::
  $ export MACHINE_NAME=osx-yosemite
  $ cd ${MOBO_DIR}/pvfmm/pvfmm-utils
  $ make lib
  $ cd ${MOBO_DIR}
  $ source ${MOBO_DIR}/config/checkenv.rc
  $ # will warn you if you have an unset variable needed for compilation
  $ cd ${MOBO_DIR}/src/ebi/
  $ make BOPT=O

Build an example
----------------
:: 

  $ cd harper
  $ make  BOPT=O on_surface_eval_test_with_solve
  $ on_surface_eval_test_with_solve --kernel laplace --layer single  --bounded --domain cube --face-map

Later use
---------
You should add above defined variables to some rc file that you can source (e.g. your .bashrc or machine specific rc file in config dir). A modified version of following is a good starting point: ::

   export MOBO_DIR=<mobo_dir>
   source ${MOBO_DIR}/config/default.rc

 Note: if Ebi build complains about a missing environment variable, check the Dockerfile at the project root for a set of working build instructions on Ubuntu 14.04
 Another note: you may need to adjust your fortran compiler in machine_makefile/makefile.osx-yosemite.in

Compile Geogram
---------------
::
 cd geogram_1.6.2/
 mkdir build && cd build
 cmake -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8 -DVORPALINE_PLATFORM=Linux64-gcc -DGEOGRAM_LIB_ONLY=ON -DGEOGRAM_WITH_GRAPHICS=OFF ..
 #for intel repalce gcc/g++ with full path to icc and icpc; Linux64-gcc with Linux64-icc
 # possibly bugs in cmake/platforms/Linux-icc.cmake; may need to change the
 # value of CMAKE_C_COMPILER/CMAKE_CXX_COMPILER to full path to icc/icpc
 make 

