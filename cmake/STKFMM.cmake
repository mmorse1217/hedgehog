include_guard()

if (NOT TARGET STKFMM::STKFMM)
    if(BUILD_SHARED_LIBS)
        set(STKFMM_LIB_SUFFIX "so")
    else()
        set(STKFMM_LIB_SUFFIX "a")
    endif()
    
    find_package(Eigen3)
    include(cmake/Eigen3.cmake)

    message("Downloading STKFMM and PVFMM...")
    # Download and compile PVFMM 
    set(PVFMM_DIR ${CMAKE_BINARY_DIR}/_deps/PVFMM_newBC)
    set(STKFMM_DIR ${CMAKE_BINARY_DIR}/_deps/STKFMM)
    ExternalProject_Add(
        PVFMM_newBC
        GIT_REPOSITORY https://github.com/wenyan4work/pvfmm.git
        GIT_TAG        origin/new_BC
        GIT_SHALLOW TRUE
        PREFIX ${PVFMM_DIR} 
        INSTALL_DIR ${PVFMM_DIR}
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    )
   
    # Download and compile STKFMM 
    ExternalProject_Add(
        STKFMM
        GIT_REPOSITORY https://github.com/wenyan4work/STKFMM.git
        GIT_TAG        origin/master
        GIT_SHALLOW TRUE
        PREFIX ${STKFMM_DIR}
        INSTALL_DIR ${STKFMM_DIR}
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_PREFIX_PATH=${PVFMM_DIR} -Dpvfmm_DIR=${PVFMM_DIR} -DCMAKE_CXX_COMPILER=mpicxx 
        DEPENDS PVFMM_newBC 
        )
#cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/libs/pvfmm_newBC/install -Dpvfmm_DIR=/libs/pvfmm_newBC/install -DCMAKE_CXX_COMPILER=mpicxx .. && \
    
    # make an imported library target 
    add_library(PVFMM_newBC::PVFMM_newBC INTERFACE IMPORTED)

    # Touch include directories for PVFMMfloat and PVFMMdouble to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(PVFMM_newBC INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include/pvfmm)
    set(PVFMM_INCLUDE_DIR ${INSTALL_DIR}/include/pvfmm)


    target_include_directories(PVFMM_newBC::PVFMM_newBC
        INTERFACE ${PVFMM_INCLUDE_DIR} 
        )
    set(PVFMM_LIB ${INSTALL_DIR}/lib/libpvfmm.a)
    # set link library locations
    target_link_libraries(PVFMM_newBC::PVFMM_newBC
        INTERFACE "${PVFMM_LIB}"
        )
    
    # Set explicit include/library variables (just in case)
    set(PVFMM_INCLUDE_DIRS 
        "${PVFMM_INCLUDE_DIR}"
        )
    set(PVFMM_LIBRARIES 
        "${PVFMM_LIB}"
         )
    add_library(STKFMM::STKFMM INTERFACE IMPORTED)

    # Touch include directories for PVFMMfloat and PVFMMdouble to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(STKFMM INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include/STKFMM)
    set(STKFMM_INCLUDE_DIR ${INSTALL_DIR}/include/STKFMM)


    target_include_directories(STKFMM::STKFMM
        INTERFACE ${STKFMM_INCLUDE_DIR} 
        )
    set(STKFMM_LIB ${INSTALL_DIR}/lib/libSTKFMM_STATIC.a)
    # set link library locations
    target_link_libraries(STKFMM::STKFMM
        INTERFACE "${STKFMM_LIB}"
        )
    
    # Set explicit include/library variables (just in case)
    set(STKFMM_INCLUDE_DIRS 
        "${STKFMM_INCLUDE_DIR}"
        )
    set(STKFMM_LIBRARIES 
        "${STKFMM_LIB}"
         )
endif() 
