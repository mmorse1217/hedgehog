
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(PVFMMUtils)
set(PVFMMUtils_VERSION 0.0.1)

set(pvfmmutils_SRC
    pvfmm_interface.cc
    )
set(pvfmmutils_INC
    pvfmm_interface.h
    pvfmm_laplace_kernel.cpp
    pvfmm_mod_helmholtz_kernel.cpp
    pvfmm_navier_kernel.cpp
    pvfmm_stks_kernel.tpp
    )

list(APPEND CMAKE_MODULE_PATH "${PVFMMUtils_SOURCE_DIR}/cmake")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PVFMMUtils_SOURCE_DIR}/lib")
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(PVFMM REQUIRED)

add_library(pvfmmutils STATIC ${pvfmmutils_SRC})
target_include_directories(pvfmmutils PUBLIC
   # including headers when building from source
   $<BUILD_INTERFACE:${pvfmmutils_SOURCE_DIR}/include> 
   # including headers when installing 
   $<INSTALL_INTERFACE:include> 
   )
target_include_directories(pvfmmutils PUBLIC 
    ${MPI_CXX_INCLUDE_DIRS}
    )
message(STATUS "MPI dirs: " ${MPI_CXX_INCLUDE_DIRS})
message(STATUS "PVFMM dirs: " ${PVFMM_INCLUDE_DIR}/pvfmm)
message(STATUS "PVFMM dirs: " ${PVFMM_LIB_DIR}/pvfmm/${PVFMM_STATIC_LIB})
target_link_libraries(pvfmmutils PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX
    PVFMM::PVFMM)


set(PVFMMUtils_INCLUDE_DEST "include/pvfmmutils-${PVFMMUtils_VERSION}")
set(PVFMMUtils_LIB_DEST "lib/pvfmmutils-${PVFMMUtils_VERSION}")
# Install pvfmmutils in CMAKE_INSTALL_PREFIX (defaults to /usr/local on linux)
install(
    # targets to install
    TARGETS pvfmmutils 
    # name of the CMake "export group" containing the targets we want to install
    EXPORT pvfmmutilsTargets
    # Dynamic, static library and include destination locations after running
    # "make install"
    LIBRARY DESTINATION ${PVFMMUtils_LIB_DEST}
    ARCHIVE DESTINATION ${PVFMMUtils_LIB_DEST} 
    INCLUDES DESTINATION ${PVFMMUtils_INCLUDE_DEST}
    )
install(
    # export group we want to save (matches name defined above containing the
    # install targets)
    EXPORT pvfmmutilsTargets
    # CMake file in which to store the export group info
    FILE  pvfmmutilsTargets.cmake
    # namespace prepends all targets in the export group
    NAMESPACE PVFMMUtils::
    # where to place the resulting file (here we're putting it with the library)
    DESTINATION ${PVFMMUtils_LIB_DEST}
    )

# Copies the files in ${PVFMMUtils_INC} to ${PVFMMUtils_INCLUDE_DEST}
install(FILES ${pvfmmutils_INC} DESTINATION ${PVFMMUtils_INCLUDE_DEST})

# Write a "version file" in case some wants to only load a particular version of
# pvfmmutils @BUG this fails at the moment.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    pvfmmutilsConfigVersion.cmake
    VERSION ${PVFMMUtils_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# Copies the resulting CMake config files to the installed library directory
install(FILES 
    "cmake/pvfmmutils-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/pvfmmutilsConfigVersion.cmake"
    DESTINATION ${PVFMMUtils_LIB_DEST})

