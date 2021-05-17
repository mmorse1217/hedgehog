
include_guard()
if (NOT TARGET LeanVTK::LeanVTK)
    message("Downloading LeanVTK...")
    FetchContent_Declare(
        LeanVTK
        SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/LeanVTK
        GIT_REPOSITORY https://github.com/mmorse1217/lean-vtk.git
    )

    FetchContent_GetProperties(LeanVTK)
    if(NOT LeanVTK_POPULATED)
        FetchContent_Populate(LeanVTK)
        # Build LeanVTK in LeanVTK_SOURCE_DIR and put build artifacts in
        # LeanVTK_BINARY_DIR
        message(${leanvtk_SOURCE_DIR} "   " ${leanvtk_BINARY_DIR})
        add_subdirectory(${leanvtk_SOURCE_DIR} ${leanvtk_BINARY_DIR})

        # Make link/include target
        add_library(LeanVTK::LeanVTK INTERFACE IMPORTED)
        target_include_directories(LeanVTK::LeanVTK SYSTEM INTERFACE
            ${leanvtk_SOURCE_DIR}/include)
        target_link_libraries(LeanVTK::LeanVTK INTERFACE
            ${leanvtk_SOURCE_DIR}/lib/libLeanVTK.a)
    endif()
endif()
