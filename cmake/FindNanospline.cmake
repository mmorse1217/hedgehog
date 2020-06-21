# - Try to find the Nanospline library
# Once done this will define
#
#  Nanospline_FOUND - system has Nanospline
#  Nanospline_INCLUDE_DIR - Nanospline include directory
#  Nanospline_LIB - Nanospline library directory
#  Nanospline_LIBRARIES - Nanospline libraries to link

if(Nanospline_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/nanospline*")
file(GLOB inc_glob "/usr/local/include/nanospline*")
find_library(Nanospline_LIB 
    NAMES nanospline
    HINTS
        ${Nanospline_DIR}
        ${NANOSPLINE_DIR}
        $ENV{Nanospline_DIR}
        $ENV{NANOSPLINE_DIR}
        ENV NANOSPLINE_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
        build
)
find_path(Nanospline_INCLUDE_DIR BSpline.h
    HINTS
        ${Nanospline_DIR}
        ${NANOSPLINE_DIR}
        $ENV{Nanospline_DIR}
        $ENV{NANOSPLINE_DIR}
        ENV NANOSPLINE_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include/nanospline
)

message(STATUS "NANOSPLINE_DIR: " $ENV{NANOSPLINE_DIR})
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Nanospline
    "\nNanospline not found --- You can download it using:\n\tgit clone https://github.com/qnzhou/nanospline\n and setting the NANOSPLINE_DIR environment variable accordingly"
    Nanospline_LIB Nanospline_INCLUDE_DIR)
mark_as_advanced(Nanospline_INCLUDE_DIR Nanospline_LIB)

set(Nanospline_INCLUDE_DIRS ${Nanospline_INCLUDE_DIR})
set(Nanospline_LIBRARIES ${Nanospline_LIB})

if(Nanospline_FOUND AND NOT TARGET Nanospline::Nanospline)
    add_library(Nanospline::Nanospline UNKNOWN IMPORTED)
    # Interface include directory
    set_target_properties(Nanospline::Nanospline PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Nanospline_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(Nanospline::Nanospline PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${Nanospline_LIBRARIES}"
        )
endif()
