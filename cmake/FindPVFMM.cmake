# - Try to find the PVFMM library
# Once done this will define
#
#  PVFMM_FOUND - system has PVFMM
#  PVFMM_INCLUDE_DIR - PVFMM include directory
#  PVFMM_LIB - PVFMM library directory
#  PVFMM_LIBRARIES - PVFMM libraries to link

if(PVFMM_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/pvfmm*")
file(GLOB inc_glob "/usr/local/include/pvfmm*")
find_library(PVFMM_LIB 
    NAMES libpvfmm.a pvfmm
    HINTS
        ${PVFMM_DIR}
        ${PVFMM_DIR}
        $ENV{PVFMM_DIR}
        $ENV{PVFMM_DIR}
        ENV PVFMM_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(PVFMM_INCLUDE_DIR pvfmm.hpp
    HINTS
        ${PVFMM_DIR}
        ${PVFMM_DIR}
        $ENV{PVFMM_DIR}
        $ENV{PVFMM_DIR}
        ENV PVFMM_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PVFMM
    "\nPVFMM not found --- You can download it using:\n\tgit clone 
    https://github.com/dmalhotra/pvfmm\n and setting the PVFMM_DIR environment variable accordingly"
    PVFMM_LIB PVFMM_INCLUDE_DIR)
mark_as_advanced(PVFMM_INCLUDE_DIR PVFMM_LIB)

set(PVFMM_INCLUDE_DIRS ${PVFMM_INCLUDE_DIR})
set(PVFMM_LIBRARIES ${PVFMM_LIB})

if(PVFMM_FOUND AND NOT TARGET PVFMM::PVFMM)
    add_library(PVFMM::PVFMM UNKNOWN IMPORTED)
    # Interface include directory
    set_target_properties(PVFMM::PVFMM PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${PVFMM_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(PVFMM::PVFMM PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${PVFMM_LIBRARIES}"
        )
endif()
