# - Try to find the p4est library
# Once done this will define
#
#  p4est_FOUND - system has p4est
#  p4est_INCLUDE_DIR - p4est include directory
#  p4est_LIB - p4est library directory
#  p4est_LIBRARIES - p4est libraries to link

if(p4est_FOUND)
    return()
endif()


file(GLOB lib_glob "/usr/local/lib/p4est*")
file(GLOB inc_glob "/usr/local/include/p4est*")
find_library(p4est_LIB 
    NAMES libp4est.a p4est 
    HINTS
        ${p4est_DIR}
        ${P4EST_DIR}
        $ENV{p4est_DIR}
        $ENV{P4EST_DIR}
        ENV P4EST_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_library(sc_LIB 
    NAMES libsc.a sc
    HINTS
        ${p4est_DIR}
        ${P4EST_DIR}
        $ENV{p4est_DIR}
        $ENV{P4EST_DIR}
        ENV P4EST_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(p4est_INCLUDE_DIR p4est.h
    HINTS
        ${p4est_DIR}
        ${P4EST_DIR}
        $ENV{p4est_DIR}
        $ENV{P4EST_DIR}
        ENV P4EST_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p4est
    "\np4est not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/p4est\n and setting the P4EST_DIR environment variable accordingly"
    p4est_LIB p4est_INCLUDE_DIR)
mark_as_advanced(p4est_INCLUDE_DIR p4est_LIB)



#set(p4est_LIB "${p4est_INCLUDE_DIR}/../lib/")
#set(p4est_LIBRARIES "${p4est_INCLUDE_DIR}/../lib/libp4est.a" "${p4est_INCLUDE_DIR}/../lib/libsc.a" )
#list(APPEND CMAKE_MODULE_PATH "${p4est_INCLUDE_DIR}/../cmake")
set(p4est_INCLUDE_DIRS ${p4est_INCLUDE_DIR})
set(p4est_LIBRARIES ${p4est_LIB})
set(sc_LIBRARIES ${sc_LIB})



if(p4est_FOUND AND NOT TARGET p4est::p4est)
    add_library(p4est::p4est UNKNOWN IMPORTED)
    add_library(p4est::sc UNKNOWN IMPORTED)
    # Interface include directory
    set_target_properties(p4est::p4est PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${p4est_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(p4est::p4est PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${p4est_LIBRARIES}"
        )
    set_target_properties(p4est::sc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${p4est_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(p4est::sc PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${sc_LIBRARIES}"
        )
endif()
