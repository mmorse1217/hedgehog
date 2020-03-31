cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

#sets compiler to match $HEDGEHOG_CC/CXX
execute_process(COMMAND which $ENV{HEDGEHOG_CC} OUTPUT_VARIABLE cc)
execute_process(COMMAND which $ENV{HEDGEHOG_CXX} OUTPUT_VARIABLE cpp)
message(STATUS "$ENV{CC}")
set(ENV{CC}  $ENV{HEDGEHOG_CC} )
set(ENV{CXX} $ENV{HEDGEHOG_CXX})
message(STATUS "$ENV{CC}")

project(Hedgehog)
set(Hedgehog_VERSION 0.1)
list(APPEND CMAKE_MODULE_PATH "${Hedgehog_SOURCE_DIR}/cmake")

#find_library(MPI_LIB
#    NAMES  libmpi.dylib libmpich.dylib
#    PATHS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH} $ENV{MPI_HOME} 
#    PATH_SUFFIXES lib
#    SHARED
#    )
#find_library(PETSC_LIB
#    NAMES  libpetsc.dylib 
#    PATHS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH} 
#    PATH_SUFFIXES lib
#    SHARED
#    )
#find_library(P4EST_LIB
#    NAMES  libp4est.dylib 
#    PATHS $ENV{P4EST_DIR}
#    PATH_SUFFIXES lib
#    SHARED
#    )
#find_library(SC_LIB
#    NAMES  libsc.dylib 
#    PATHS $ENV{P4EST_DIR}
#    PATH_SUFFIXES lib
#    SHARED
#    )
#find_library(BLENDSURF_LIB
#    NAMES  libblend.a 
#    PATHS $ENV{BLENDSURF_DIR}
#    STATIC
#    )
#find_library(FACEMAP_LIB
#    NAMES  libpatchwork.a 
#    PATHS $ENV{FACEMAP_DIR}
#    PATH_SUFFIXES cbuild/lib
#    STATIC
#    )
#find_library(PVFMM_LIB
#    NAMES  libpvfmm.a 
#    PATHS $ENV{PVFMM_DIR}
#    PATH_SUFFIXES lib/pvfmm
#    STATIC
#    )
#find_library(PVFMM_INTERFACE_LIB
#    NAMES  libpvfmmwrap.a 
#    PATHS $ENV{HEDGEHOG_DIR}
#    PATH_SUFFIXES utils/pvfmm-utils
#    STATIC
#    )
#find_library(GEOGRAM_LIB
#    NAMES  libgeogram.a 
#    PATHS $ENV{GEOGRAM_DIR}
#    PATH_SUFFIXES build/lib
#    STATIC
#    )
#find_library(FFTW_LIB
#    NAMES  libfftw3.a 
#    PATHS $ENV{FFTW_DIR}
#    PATH_SUFFIXES lib
#    STATIC
#    )

find_package(OpenMP REQUIRED)
find_package(Blendsurf REQUIRED)
#find_package(patchwork REQUIRED)
find_package(pvfmmutils REQUIRED)
find_package(pvfmm REQUIRED)
find_package(p4est REQUIRED)
find_package(VTK 7 REQUIRED)
find_package(PETSc REQUIRED)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#find_package(expat REQUIRED)

FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

add_subdirectory(src)
add_subdirectory(tests)

