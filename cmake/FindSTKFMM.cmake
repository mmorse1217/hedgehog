# - Try to find the STKFMM library # Once done this will define
#
#  STKFMM_FOUND - system has STKFMM
#  STKFMM_INCLUDE_DIR - STKFMM include directory
#  STKFMM_LIB - STKFMM library directory
#  STKFMM_LIBRARIES - STKFMM libraries to link

if(STKFMM_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/*")
file(GLOB inc_glob "/usr/local/include/*")
find_library(STKFMM_LIB 
    NAMES libSTKFMM_STATIC.a STKFMM
    HINTS
        ${STKFMM_DIR}
        ${STKFMM_DIR}
        $ENV{STKFMM_DIR}
        $ENV{STKFMM_DIR}
        ENV STKFMM_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(STKFMM_INCLUDE_DIR STKFMM.h
    HINTS
        ${STKFMM_DIR}
        ${STKFMM_DIR}
        $ENV{STKFMM_DIR}
        $ENV{STKFMM_DIR}
        ENV STKFMM_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        STKFMM
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(STKFMM
    "\nSTKFMM not found --- You can download it using:\n\tgit clone https://github.com/wenyan4work/stkfmm\n and setting the STKFMM_DIR environment variable accordingly"
    STKFMM_LIB STKFMM_INCLUDE_DIR)
mark_as_advanced(STKFMM_INCLUDE_DIR STKFMM_LIB)

set(STKFMM_INCLUDE_DIRS ${STKFMM_INCLUDE_DIR})
set(STKFMM_LIBRARIES ${STKFMM_LIB})

if(STKFMM_FOUND AND NOT TARGET STKFMM::STKFMM)
    add_library(STKFMM::STKFMM UNKNOWN IMPORTED)
    # Interface include directory
    set_target_properties(STKFMM::STKFMM PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${STKFMM_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(STKFMM::STKFMM PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${STKFMM_LIBRARIES}"
        )
endif()
