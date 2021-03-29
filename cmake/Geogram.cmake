include_guard()
if (NOT TARGET Geogram::geogram)
    if(BUILD_SHARED_LIBS)
        set(GEOGRAM_LIB_SUFFIX "so")
    else()
        set(GEOGRAM_LIB_SUFFIX "a")
    endif()
    message("Downloading Geogram...")
    # Download and compile Geogram 
    ExternalProject_Add(
        Geogram 
        GIT_REPOSITORY https://github.com/alicevision/geogram.git
        GIT_TAG        8b2ae6148c7ab1564fa2700673b4275296ce80d3
        GIT_SHALLOW TRUE
        PREFIX ${CMAKE_BINARY_DIR}/_deps/Geogram
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/Geogram
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DGEOGRAM_WITH_VORPALINE=OFF -DVORPALINE_PLATFORM=Linux64-gcc -DGEOGRAM_WITH_EXPLORAGRAM=OFF -DGEOGRAM_LIB_ONLY=ON -DGEOGRAM_WITH_GRAPHICS=OFF
    )
    
    # make an imported library target 
    add_library(Geogram::geogram INTERFACE IMPORTED)

    # Touch include directories for  Geogram to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(Geogram INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include/geogram1)
    set(GEOGRAM_INCLUDE_DIR ${INSTALL_DIR}/include/geogram1)


    target_include_directories(Geogram::geogram
        INTERFACE ${GEOGRAM_INCLUDE_DIR} 
        )
    set(GEOGRAM_LIB ${INSTALL_DIR}/lib/libgeogram.${GEOGRAM_LIB_SUFFIX})
    # set link library locations
    target_link_libraries(Geogram::geogram
        INTERFACE "${GEOGRAM_LIB}"
        )
    
    # Set explicit include/library variables (just in case)
    set(GEOGRAM_INCLUDE_DIRS 
        "${GEOGRAM_INCLUDE_DIR}"
        )
    set(GEOGRAM_LIBRARIES 
        "${GEOGRAM_LIB}"
         )
endif()
