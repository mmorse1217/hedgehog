# - Try to find the FFTW library
# Once done this will define
#
#  FFTW_FOUND - system has FFTW
#  FFTW_INCLUDE_DIR - FFTW include directory
#  FFTW_LIB - FFTW library directory
#  FFTW_LIBRARIES - FFTW libraries to link

if(FFTW_FOUND)
    return()
endif()

file(GLOB lib_glob "/usr/local/lib/fftw*")
file(GLOB inc_glob "/usr/local/include/fftw*")
find_library(FFTW_LIB 
    NAMES fftw3
    HINTS
        ${FFTW_DIR}
        ${FFTW_DIR}
        $ENV{FFTW_DIR}
        $ENV{FFTW_DIR}
        ENV FFTW_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    PATH_SUFFIXES 
        lib
)
find_path(FFTW_INCLUDE_DIR fftw3.h
    HINTS
        ${FFTW_DIR}
        ${FFTW_DIR}
        $ENV{FFTW_DIR}
        $ENV{FFTW_DIR}
        ENV FFTW_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW
    "\nFFTW not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/fftw\n and setting the FFTW_DIR environment variable accordingly"
    FFTW_LIB FFTW_INCLUDE_DIR)
mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIB)

set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(FFTW_LIBRARIES ${FFTW_LIB})
#set(FFTW_LIBRARIES "${FFTW_LIB}/../lib/libfftw.a"  )
#set(FFTW_render_LIBRARIES "${FFTW_LIB}/../lib/libfftw_render.a"  )
#set(FFTW_render_INCLUDE_DIR "${FFTW_LIB}/../vis/"  )

if(FFTW_FOUND AND NOT TARGET FFTW::FFTW)
    add_library(FFTW::FFTW STATIC IMPORTED)
  set_target_properties(FFTW::FFTW PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES  "${FFTW_INCLUDE_DIRS}"
    IMPORTED_LOCATION "${FFTW_LIBRARIES}")
endif()
