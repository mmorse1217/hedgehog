# - Try to find the Blendsurf library
# Once done this will define
#
#  Blendsurf_FOUND - system has Blendsurf
#  Blendsurf_INCLUDE_DIR - Blendsurf include directory
#  Blendsurf_LIB - Blendsurf library directory
#  Blendsurf_LIBRARIES - Blendsurf libraries to link

if(Blendsurf_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/blendsurf*")
file(GLOB inc_glob "/usr/local/include/blendsurf*")
find_library(Blendsurf_LIB 
    NAMES blendsurf
    HINTS
        ${Blendsurf_DIR}
        ${BLENDSURF_DIR}
        $ENV{Blendsurf_DIR}
        $ENV{BLENDSURF_DIR}
        ENV BLENDSURF_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(Blendsurf_INCLUDE_DIR bdsurf.hpp
    HINTS
        ${Blendsurf_DIR}
        ${BLENDSURF_DIR}
        $ENV{Blendsurf_DIR}
        $ENV{BLENDSURF_DIR}
        ENV BLENDSURF_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Blendsurf
    "\nBlendsurf not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/blendsurf\n and setting the BLENDSURF_DIR environment variable accordingly"
    Blendsurf_LIB Blendsurf_INCLUDE_DIR)
mark_as_advanced(Blendsurf_INCLUDE_DIR Blendsurf_LIB)

set(Blendsurf_INCLUDE_DIRS ${Blendsurf_INCLUDE_DIR})
set(Blendsurf_LIBRARIES ${Blendsurf_LIB})
#set(Blendsurf_LIBRARIES "${Blendsurf_LIB}/../lib/libblendsurf.a"  )
#set(Blendsurf_render_LIBRARIES "${Blendsurf_LIB}/../lib/libblendsurf_render.a"  )
#set(Blendsurf_render_INCLUDE_DIR "${Blendsurf_LIB}/../vis/"  )

if(Blendsurf_FOUND AND NOT TARGET Blendsurf::Blendsurf)
    add_library(Blendsurf::Blendsurf UNKNOWN IMPORTED)
    # Interface include directory
    set_target_properties(Blendsurf::Blendsurf PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Blendsurf_INCLUDE_DIRS}"
        )

    # Link to library file
    set_target_properties(Blendsurf::Blendsurf PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${Blendsurf_LIBRARIES}"
        )
endif()
