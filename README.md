# `hedgehog`

`hedgehog` is a library to solve constant coefficient elliptic partial differential equations on complex geometries.

## Using `hedgehog`

### Prerequisites
To install `hedgehog`, you will need at least the following dependences:
 - GCC v7.5.0 or higher (or another compiler supporting OpenMP)
 - CMake v3.17.3 or higher
 - MPI (Tested with OpenMPI 4.0.5 compiled with Intel icpc 19.1.2)
 - [Patchwork](https://github.com/mmorse1217/patchwork)
 - VTK 7.0.0 or higher
These must be installed before compiling `hedgehog`

`hedgehog` depends on the following libraries:
 - PETSc v3.7.2 or higher
 - Eigen v3.3.7 or higher
 - FFTW v3.3.4 or higher
 - [p4est](https://github.com/cburstedde/p4est) 1.1 (newer releases untested)
 - [PVFMM](https://github.com/dmalhotra/pvfmm)
 - [STKFMM](https://github.com/wenyan4work/STKFMM) (and implicitly [the PVFMM new_BC branch](https://github.com/wenyan4work/pvfmm))
 - [Geogram](https://github.com/alicevision/geogram)
 - [Blendsurf](https://github.com/mmorse1217/blendsurf)
 - [Nanospline](https://github.com/qnzhou/nanospline)

`hedgehog` attempts to download and install these dependencies for you.

### Downloading `hedgehog`
The best to get `hedgehog` at the moment is by cloning from Github or HTTPS
```bash
$ git clone https://github.com/mmorse1217/`hedgehog`.git
```
or with SSH keys, if you have them configured
```bash
$ git clone git@github.com:mmorse1217/`hedgehog`.git
```
### Using `hedgehog`
There are two main ways to use `hedgehog`: compiling directly with CMake or with a Docker image.
Docker is the preferred method of quickly compiling and running `hedgehog`, but CMake is often preferred when trying to incorporate `hedgehog` into a new project.

### Installing `hedgehog` with CMake
In an ideal world, running the following commands from within the `hedgehog`
directory should work:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
and optionally,
```bash
$ make install
```

The CMakeLists.txt first attempts to find dependencies already installed by the user, then
downloads and install missing dependencies.
This has been tested on Ubuntu 20.04 and RedHat 8.2. 
However, sometimes CMake has trouble with one stage or the other.
It first checks `/usr/local/include` and `/usr/local/lib` for headers and libraries, then checks certain environment variables.
If the find module is unable to locate your preinstalled version, or you would just like to install things yourself, setting the following environment variables for each dependency will ensure CMake finds the proper package
(unless otherwise specified, setting an environment variable `PACKAGE_DIR` will tell CMake that `package` has headers in `${PACKAGE_DIR}/include/` and libraries in `${PACKAGE_DIR}/lib`):

 - PETSc: `${PETSC_DIR}` and `${PETSC_ARCH}` (see [here](https://www.mcs.anl.gov/petsc/documentation/installation.html#envvars) for more details).
    if PETSc is an in-source build (i.e. run without `make install`), `PETSC_ARCH` is needed and must match the "arch" that PETSc determined during the compilation process
 - Eigen: `${EIGEN3_ROOT}` or `${EIGEN3_ROOT_DIR}`
 - FFTW: `${FFTW_DIR}` (note: it needs float and double version of FFTW)
 - p4est: `${P4EST_DIR}` or `${p4est_DIR}`. Also requires compilation of the `sc` library.
 - PVFMM: `${PVFMM_DIR}` or `${PVFMM_DIR}`. 
 - STKFMM: `${STKFMM_DIR}` or `${STKFMM_DIR}`.
 - Geogram: `${GEOGRAM_INSTALL_PREFIX}`
 - Blendsurf: `${BLENDSURF_DIR}` or `${Blendsurf_DIR}`
 - Nanospline: `${NANOSPLINE_DIR}` or `${Nanospline_DIR}`, with headers in `${NANOSPLINE_DIR}/include/nanospline` and libraries in `${NANOSPLINE_DIR}/lib` 

To force cmake to use a certain compiler over another, use 
```bash
$ CC=<your-c-compiler> CXX=<your-c++-compiler> cmake ..
```
instead of the standard line. 
There is a known issue with PETSc autodownloading. 
With Python 3.8.6, there seems to be a strange issue with the PETSc install script failing with an error from `shtuils`, looking something like:
```
Errno 11: Resource temporarily unavailable
```
Fortunately, PETSc is downloaded and compiled. Simply setting:
```
export PETSC_DIR=<path-to-hedgehog>/build/_deps/PETSc/
```
and correctly setting `PETSC_ARCH` will point CMake to the new installation (look in `<path-to-hedgehog>/build/_deps/PETSc/src/PETSc` for a folder with a architecture-related name or look at the PETSc build configure output).

### Building and running `hedgehog` with Docker
TODO: double check this section
To build a Docker image from scratch, run:
```bash
$ make image-build
```
To download the pre-built container from Docker Hub and skip building dependencies
```bash
$ docker pull mmorse1217/hedgehog
$ docker tag mmorse1217/hedgehog hedgehog
```
To start the container:
```bash
$ make container-create
$ make container-start
$ make container-exec
root@<container-hash> $
```
To compile hedgehog:
```bash
root@<container-hash> $ mkdir build
root@<container-hash> $ cd build
root@<container-hash> $ cmake ..
root@<container-hash> $ make 
```
To run the unit tests:
```
root@<container-hash> $ cd ..
root@<container-hash> $ build/test/test_hedgehog [regression]
```
