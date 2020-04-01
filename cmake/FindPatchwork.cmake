# - Try to find the Patchwork library
# Once done this will define
#
#  Patchwork_FOUND - system has Patchwork
#  Patchwork_INCLUDE_DIR - Patchwork include directory
#  Patchwork_LIB - Patchwork library directory
#  Patchwork_LIBRARIES - Patchwork libraries to link

if(Patchwork_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/patchwork*")
file(GLOB inc_glob "/usr/local/include/patchwork*")
find_library(Patchwork_LIB 
    NAMES patchwork
    HINTS
        ${Patchwork_DIR}
        ${PATCHWORK_DIR}
        $ENV{Patchwork_DIR}
        $ENV{PATCHWORK_DIR}
        ENV PATCHWORK_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(Patchwork_INCLUDE_DIR face_map.hpp
    HINTS
        ${Patchwork_DIR}
        ${PATCHWORK_DIR}
        $ENV{Patchwork_DIR}
        $ENV{PATCHWORK_DIR}
        ENV PATCHWORK_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Patchwork
    "\nPatchwork not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/patchwork\n and setting the PATCHWORK_DIR environment variable accordingly"
    Patchwork_LIB Patchwork_INCLUDE_DIR)
mark_as_advanced(Patchwork_INCLUDE_DIR Patchwork_LIB)

set(Patchwork_INCLUDE_DIRS ${Patchwork_INCLUDE_DIR})
set(Patchwork_LIBRARIES ${Patchwork_LIB})
#set(Patchwork_LIBRARIES "${Patchwork_LIB}/../lib/libpatchwork.a"  )
#set(Patchwork_render_LIBRARIES "${Patchwork_LIB}/../lib/libpatchwork_render.a"  )
#set(Patchwork_render_INCLUDE_DIR "${Patchwork_LIB}/../vis/"  )

if(Patchwork_FOUND AND NOT TARGET Patchwork::Patchwork)
    add_library(Patchwork::Patchwork STATIC IMPORTED)
  set_target_properties(Patchwork::Patchwork PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES  "${Patchwork_INCLUDE_DIRS}"
    IMPORTED_LOCATION "${Patchwork_LIBRARIES}")
endif()
