include_guard()
if (NOT TARGET PETSc::PETSc)
    if(BUILD_SHARED_LIBS)
        set(PETSc_LIB_SUFFIX "so")
    else()
        set(PETSc_LIB_SUFFIX "a")
    endif()


    # Compile PETSc 
    ExternalProject_Add(
        PETSc
        URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.2.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/_deps/PETSc
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/PETSc
        CONFIGURE_COMMAND cd <SOURCE_DIR> && ./configure --with-pic=true --with-shared-libraries=0 --prefix=<INSTALL_DIR> --with-debugging=0 --with-64-bit-indices=1 --with-single-library=1 
        BUILD_COMMAND make -C <SOURCE_DIR>
        INSTALL_COMMAND make -C <SOURCE_DIR> install
    )
    
    # make an imported library target 
    #add_library(PETSc::PETSc INTERFACE IMPORTED)

    # Touch include directories for PETScfloat and PETScdouble to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(PETSc INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
    set(PETSc_INCLUDE_DIR ${INSTALL_DIR}/include)


#    target_include_directories(PETSc::PETSc
#        INTERFACE ${PETSc_FLOAT_INCLUDE_DIR} 
#        INTERFACE ${PETSc_DOUBLE_INCLUDE_DIR}
#        )
#
#    # set link library locations
#    target_link_libraries(PETSc::PETSc
#        INTERFACE ${CMAKE_BINARY_DIR}/_deps/PETScFloat/lib/libfftw3f.${PETSc_LIB_SUFFIX}
#        INTERFACE ${CMAKE_BINARY_DIR}/_deps/PETScFloat/lib/libfftw3f_threads.${PETSc_LIB_SUFFIX}
#        INTERFACE ${CMAKE_BINARY_DIR}/_deps/PETScDouble/lib/libfftw3.${PETSc_LIB_SUFFIX}
#        INTERFACE ${CMAKE_BINARY_DIR}/_deps/PETScDouble/lib/libfftw3_threads.${PETSc_LIB_SUFFIX}
#        )
#    
#    # Set explicit include/library variables (just in case)
#    set(PETSc_INCLUDE_DIRS 
#        "${CMAKE_BINARY_DIR}/_deps/PETScDouble/include" 
#        "${CMAKE_BINARY_DIR}/_deps/PETScFloat/include"
#        )
#    set(PETSc_LIBRARIES 
#         "${CMAKE_BINARY_DIR}/_deps/PETScFloat/lib/libfftw3f.${PETSc_LIB_SUFFIX}"
#         "${CMAKE_BINARY_DIR}/_deps/PETScFloat/lib/libfftw3f_threads.${PETSc_LIB_SUFFIX}"
#         "${CMAKE_BINARY_DIR}/_deps/PETScDouble/lib/libfftw3.${PETSc_LIB_SUFFIX}"
#         "${CMAKE_BINARY_DIR}/_deps/PETScDouble/lib/libfftw3_threads.${PETSc_LIB_SUFFIX}"
#         )
endif() 
