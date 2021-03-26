include_guard()

if (NOT TARGET Blendsurf::Blendsurf)
    message("Downloading Blendsurf...")
    FetchContent_Declare(
        Blendsurf
        SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/Blendsurf
        BINARY_DIR ${CMAKE_BINARY_DIR}/_deps/Blendsurf
        GIT_REPOSITORY https://github.com/mmorse1217/blendsurf.git
        GIT_TAG        f13344873ec97650f9a1d43cffa5892a41570130
        GIT_SHALLOW TRUE
    )
    
    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL  "Release or debug mode")

    FetchContent_GetProperties(Blendsurf)
    if(NOT blendsurf_POPULATED)
        FetchContent_Populate(Blendsurf)
        add_subdirectory(${blendsurf_SOURCE_DIR} ${blendsurf_BINARY_DIR})
        add_library(Blendsurf::Blendsurf INTERFACE IMPORTED)
        target_include_directories(Blendsurf::Blendsurf SYSTEM INTERFACE
            ${blendsurf_SOURCE_DIR}/include)
        target_link_libraries(Blendsurf::Blendsurf INTERFACE
            ${blendsurf_BINARY_DIR}/lib/libblendsurf.a)
    endif()
endif()
