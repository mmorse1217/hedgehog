include_guard()


if (NOT TARGET PVFMM::PVFMM)
    if(BUILD_SHARED_LIBS)
        set(PVFMM_LIB_SUFFIX "so")
    else()
        set(PVFMM_LIB_SUFFIX "a")
    endif()


    message("Downloading PVFMM...")
    # Download and compile PVFMM 
    ExternalProject_Add(
        PVFMM
        GIT_REPOSITORY https://github.com/dmalhotra/pvfmm.git
        GIT_TAG        ffec8376dac7e2df134e56c1a37f22051ec483bb 
        GIT_SHALLOW TRUE
        PREFIX ${CMAKE_BINARY_DIR}/_deps/PVFMM
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/PVFMM
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    )
    
    # make an imported library target 
    add_library(PVFMM::PVFMM INTERFACE IMPORTED)

    # Touch include directories for PVFMMfloat and PVFMMdouble to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(PVFMM INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
    set(PVFMM_INCLUDE_DIR ${INSTALL_DIR}/include/pvfmm)


    target_include_directories(PVFMM::PVFMM
        INTERFACE ${PVFMM_INCLUDE_DIR} 
        )
    set(PVFMM_LIB ${INSTALL_DIR}/lib/libpvfmm.${PVFMM_LIB_SUFFIX})
    # set link library locations
    target_link_libraries(PVFMM::PVFMM
        INTERFACE "${PVFMM_LIB}"
        )
    
    # Set explicit include/library variables (just in case)
    set(PVFMM_INCLUDE_DIRS 
        "${PVFMM_INCLUDE_DIR}"
        )
    set(PVFMM_LIBRARIES 
        "${PVFMM_LIB}"
         )
endif() 
        
# TODO this code works if pvfmm_config.h is installed in
# ${pvfmm_SOURCE_DIR}/include instead of ${pvfmm_SOURCE_DIR}.
# when the issue is fixed in the core repo switch to configure-time setup of
# pvfmm 
#if (NOT TARGET PVFMM::PVFMM)
#    include(FetchContent)
#    FetchContent_Declare(
#        PVFMM
#        SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/PVFMM
#        BINARY_DIR ${CMAKE_BINARY_DIR}/_deps/PVFMM
#        GIT_REPOSITORY https://github.com/dmalhotra/pvfmm.git
#        GIT_TAG        ffec8376dac7e2df134e56c1a37f22051ec483bb 
#        GIT_SHALLOW TRUE
#    )
#    
#    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL  "Release or debug mode")
#
#    FetchContent_GetProperties(PVFMM)
#    if(NOT pvfmm_POPULATED)
#        FetchContent_Populate(PVFMM)
#
#        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#
#        add_subdirectory(${pvfmm_SOURCE_DIR} ${pvfmm_BINARY_DIR})
#        add_library(PVFMM::PVFMM INTERFACE IMPORTED)
#        target_include_directories(PVFMM::PVFMM SYSTEM INTERFACE
#            ${pvfmm_SOURCE_DIR}/include)
#        target_link_libraries(PVFMM::PVFMM INTERFACE
#            ${pvfmm_BINARY_DIR}/lib/libpvfmm.a)
#    endif()
#endif()
