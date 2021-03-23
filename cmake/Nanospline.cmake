include_guard()

if (NOT TARGET Nanospline::Nanospline)
    FetchContent_Declare(
        Nanospline
        SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/Nanospline
        #BINARY_DIR ${CMAKE_BINARY_DIR}/_deps/Nanospline/build
        #INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/Nanospline/build
        GIT_REPOSITORY https://github.com/qnzhou/nanospline.git
        GIT_TAG        7213a53f9f79a6b98ebeef0bae8c9aed0fe7a721
        GIT_SHALLOW TRUE
    )
    
    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL  "Release or debug mode")
    set(NANOSPLINE_HEADER_ONLY OFF CACHE INTERNAL  "Nanospline header only mode")
    set(NANOSPLINE_BUILD_TESTS OFF CACHE INTERNAL  "Build nanospline unit tests")

    FetchContent_GetProperties(Nanospline)
    if(NOT nanospline_POPULATED)
        FetchContent_Populate(Nanospline)
        message(${nanospline_SOURCE_DIR} " " ${nanospline_BINARY_DIR})
        add_subdirectory(${nanospline_SOURCE_DIR} ${nanospline_BINARY_DIR})
        add_library(Nanospline::Nanospline INTERFACE IMPORTED)
        target_include_directories(Nanospline::Nanospline SYSTEM INTERFACE
            ${nanospline_SOURCE_DIR})
        target_link_libraries(Nanospline::Nanospline INTERFACE
            ${nanospline_BINARY_DIR}/libnanospline.a)
    endif()
endif()
