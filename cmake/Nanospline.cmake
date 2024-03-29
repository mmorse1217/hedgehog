include_guard()

if (NOT TARGET Nanospline::Nanospline)
    message("Downloading Nanospline...")
    FetchContent_Declare(
        Nanospline
        SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/Nanospline
        GIT_REPOSITORY https://github.com/qnzhou/nanospline.git
        GIT_TAG        origin/main
    )
    
    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL  "Release or debug mode")
    set(NANOSPLINE_HEADER_ONLY OFF CACHE INTERNAL  "Nanospline header only mode")
    set(NANOSPLINE_BUILD_TESTS OFF CACHE INTERNAL  "Build nanospline unit tests")

    FetchContent_GetProperties(Nanospline)
    if(NOT nanospline_POPULATED)
        FetchContent_Populate(Nanospline)
        add_subdirectory(${nanospline_SOURCE_DIR} ${nanospline_BINARY_DIR})
        add_library(Nanospline::Nanospline INTERFACE IMPORTED)
        target_include_directories(Nanospline::Nanospline SYSTEM INTERFACE
            ${nanospline_SOURCE_DIR})
        target_link_libraries(Nanospline::Nanospline INTERFACE
            ${nanospline_BINARY_DIR}/libnanospline.a nanospline_
            nanospline_forward_decl)
    endif()
endif()
