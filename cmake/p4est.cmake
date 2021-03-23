include_guard()
if (NOT TARGET p4est::p4est)
    if(BUILD_SHARED_LIBS)
        set(p4est_LIB_SUFFIX "so")
    else()
        set(p4est_LIB_SUFFIX "a")
    endif()


    # Compile p4est
    ExternalProject_Add(
        p4est
        URL http://p4est.github.io/release/p4est-1.1.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/_deps/p4est
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/p4est
        CONFIGURE_COMMAND <SOURCE_DIR>/configure CC=mpicc CXX=mpic++ F77=mpif77 FC=mpif90 --enable-mpi --prefix=<INSTALL_DIR>
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )
    add_library(p4est::p4est UNKNOWN IMPORTED)
    add_library(p4est::sc UNKNOWN IMPORTED)
    
    # Touch include directory to pass generation phase without causing an error
    # because the directory doesn't exist
    ExternalProject_Get_Property(p4est INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
    
    # Set explicit include/library variables (just in case)
    set(p4est_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(sc_LIBRARIES ${INSTALL_DIR}/lib/libsc.${p4est_LIB_SUFFIX})
    set(p4est_LIBRARIES ${INSTALL_DIR}/lib/libp4est.${p4est_LIB_SUFFIX})

    # Interface include directory and link to library file for p4est...
    set_target_properties(p4est::p4est PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${p4est_INCLUDE_DIRS}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${p4est_LIBRARIES}"
        )
    # ... and sc
    set_target_properties(p4est::sc PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${p4est_INCLUDE_DIRS}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${sc_LIBRARIES}"
        )

endif() 
        
