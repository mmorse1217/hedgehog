* Prerequisites packages are petsc2.2.1 and fftw3 both of which are included in the reposoitory.
* The following instruction works on energon machines in cims and they should work with some variation for other machines.
* Depending on the version of software already installed on your machine, you can skip some of the steps below.

Pre build
-----------
::

  $ hg clone https://arahimian@bitbucket.org/arahimian/legacy-mobo
  $ cd legacy-mobo
  $ export MOBO_DIR=`pwd`

Build PETSC
-----------
petsc config file gets the libraries from different locations if CC/CXX is not set explicitly::

  $ export PETSC_DIR=${MOBO_DIR}/petsc-2.2.1
  $ export PETSC_ARCH=linux-gnu
  $ cd ${PETSC_DIR}
  $ ./config/configure.py --download-fblaslapack --with-mpi=0 -CC=gcc47 -CXX=g++47 --with-fc=0 --with-x=0
  $ make BOPT=O

Build FFTW
----------
::

  $ export FFTW_DIR=${MOBO_DIR}/fftw-3.3.4
  $ cd ${FFTW_DIR}
  $ ./configure --prefix=`pwd`  --enable-openmp
  $ make
  $ make install

Build MOBO
----------
The following assumes that makefile and config files for the MACHINE_NAME exists.
::

  $ export MACHINE_NAME=`hostname -s|sed s/[[:digit:]]*//g`
  $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PETSC_DIR}/lib/libO/${PETSC_ARCH}
  $ source ${MOBO_DIR}/config/checkenv.rc
  $ # If checkenv complains about a variable, look at ${MOBO_DIR}/config/default.rc for the expected way to set it.
  $ make BOPT=O

Build Blendsurf
---------------
::

    $ export BLENDSURF_DIR=${MOBO_DIR}/blendsurf3
    $ cd ${BLENDSURF_DIR}
    $ make

Build PvFMM
-----------
::
    
    $ cs ${MOBO_DIR}
    $ git clone https://github.com/dmalhotra/pvfmm
    $ export PVFMM_DIR=${MOBO_DIR}/pvfmm
    $ cd ${PVFMM_DIR}
    $ ./configure
    $ #if ./configure complains or doesn't complete successfully, check that
    dependencies like gcc, gfortran, FFTW, BLAS, MPI, etc. are properly linked
    in the script (see https://github.com/dmalhotra/pvfmm/blob/develop/INSTALL
    for proper flags to override defaults).
    $ make 
    $ make install
    $ make all-examples
    $ #to compile PvFMM examples; test in examples/bin/example*

Build an example
----------------
:: 

  $ cd harper
  $ make  BOPT=O ebstt0_tol

Later use
---------
You should add above defined variables to some rc file that you can source (e.g. your .bashrc or machine specific rc file in config dir). A modified version of following is a good starting point: ::

 export MOBO_DIR=<mobo_dir>
 source ${MOBO_DIR}/config/default.rc

Docker
------
Clone the repo, in the top level directory, run 
:: 
    docker build -t hedgehog .

and wait for the container to be built. Some tests will be run at various points
of the build to ensure things are being set up appropriately. Once the container is built, run

::
    docker run -it -v`pwd`:/code/ hedgehog /bin/bash
    
