cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

project(Hedgehog)
set(Hedgehog_VERSION 0.1)
list(APPEND CMAKE_MODULE_PATH "${Hedgehog_SOURCE_DIR}/cmake")
option(BUILD_SHARED_LIBS "Build shared libraries" OFF) 

include(FetchContent)
include(ExternalProject)

find_package(OpenMP REQUIRED)
find_package(Patchwork REQUIRED)
find_package(VTK 7 REQUIRED)


find_package(PETSc)
include(cmake/PETSc.cmake)
#
find_package(Eigen3)
include(cmake/Eigen3.cmake)
#
find_package(STKFMM)
include(cmake/STKFMM.cmake)
#
find_package(FFTW)
include(cmake/FFTW.cmake)
#
find_package(Nanospline)
include(cmake/Nanospline.cmake)
#
find_package(p4est)
include(cmake/p4est.cmake)
#
find_package(PVFMM)
include(cmake/PVFMM.cmake)
#
find_package(Blendsurf)
include(cmake/Blendsurf.cmake)
#
find_package(Geogram)
include(cmake/Geogram.cmake)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

add_subdirectory(src)
add_subdirectory(tests)

