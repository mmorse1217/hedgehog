# - Try to find the PVFMMUtils library
# Once done this will define
#
#  PVFMMUtils_FOUND - system has PVFMMUtils
#  PVFMMUtils_INCLUDE_DIR - PVFMMUtils include directory
#  PVFMMUtils_LIB - PVFMMUtils library directory
#  PVFMMUtils_LIBRARIES - PVFMMUtils libraries to link

if(PVFMMUtils_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/pvfmmutils*")
file(GLOB inc_glob "/usr/local/include/pvfmmutils*")
find_library(PVFMMUtils_LIB 
    NAMES pvfmmutils
    HINTS
        ${PVFMMUtils_DIR}
        ${PVFMMUTILS_DIR}
        $ENV{PVFMMUtils_DIR}
        $ENV{PVFMMUTILS_DIR}
        ENV PVFMMUTILS_DIR
    PATHS
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(PVFMMUtils_INCLUDE_DIR pvfmm_interface.h
    HINTS
        ${PVFMMUtils_DIR}
        ${PVFMMUTILS_DIR}
        $ENV{PVFMMUtils_DIR}
        $ENV{PVFMMUTILS_DIR}
        ENV PVFMMUTILS_DIR
    PATHS
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PVFMMUtils
    "\nPVFMMUtils not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/pvfmmutils\n and setting the PVFMMUTILS_DIR environment variable accordingly"
    PVFMMUtils_LIB PVFMMUtils_INCLUDE_DIR)
mark_as_advanced(PVFMMUtils_INCLUDE_DIR PVFMMUtils_LIB)

set(PVFMMUtils_INCLUDE_DIRS ${PVFMMUtils_INCLUDE_DIR})
set(PVFMMUtils_LIBRARIES ${PVFMMUtils_LIB})

if(PVFMMUtils_FOUND AND NOT TARGET PVFMMUtils::PVFMMUtils)
    add_library(PVFMMUtils::PVFMMUtils STATIC IMPORTED)
  set_target_properties(PVFMMUtils::PVFMMUtils PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES  "${PVFMMUtils_INCLUDE_DIRS}"
    IMPORTED_LOCATION "${PVFMMUtils_LIBRARIES}")
endif()
