include_guard()
if (NOT TARGET FFTW::FFTW)
    message("Downloading FFTW...")
    if(BUILD_SHARED_LIBS)
        set(FFTW_LIB_SUFFIX "so")
    else()
        set(FFTW_LIB_SUFFIX "a")
    endif()


    # Compile FFTW for floats
    ExternalProject_Add(
        FFTWFloat
        URL http://www.fftw.org/fftw-3.3.9.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/_deps/FFTWFloat
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/FFTWFloat
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DENABLE_FLOAT=On -DENABLE_THREADS=On -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    )
    
    # Compile FFTW for doubles
    ExternalProject_Add(
        FFTWDouble
        URL http://www.fftw.org/fftw-3.3.9.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/_deps/FFTWDouble
        INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/FFTWDouble
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DENABLE_THREADS=On -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    )
    # make an imported library target 
    add_library(FFTW::FFTW INTERFACE IMPORTED)

    # Touch include directories for FFTWfloat and FFTWdouble to pass 
    # generation phase without causing an error because the directory doesn't exist
    ExternalProject_Get_Property(FFTWFloat INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
    set(FFTW_FLOAT_INCLUDE_DIR ${INSTALL_DIR}/include)

    ExternalProject_Get_Property(FFTWDouble INSTALL_DIR)
    file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
    set(FFTW_DOUBLE_INCLUDE_DIR ${INSTALL_DIR}/include)

    target_include_directories(FFTW::FFTW
        INTERFACE ${FFTW_FLOAT_INCLUDE_DIR} 
        INTERFACE ${FFTW_DOUBLE_INCLUDE_DIR}
        )

    # set link library locations
    target_link_libraries(FFTW::FFTW
        INTERFACE ${CMAKE_BINARY_DIR}/_deps/FFTWFloat/lib/libfftw3f.${FFTW_LIB_SUFFIX}
        INTERFACE ${CMAKE_BINARY_DIR}/_deps/FFTWFloat/lib/libfftw3f_threads.${FFTW_LIB_SUFFIX}
        INTERFACE ${CMAKE_BINARY_DIR}/_deps/FFTWDouble/lib/libfftw3.${FFTW_LIB_SUFFIX}
        INTERFACE ${CMAKE_BINARY_DIR}/_deps/FFTWDouble/lib/libfftw3_threads.${FFTW_LIB_SUFFIX}
        )
    
    # Set explicit include/library variables (just in case)
    set(FFTW_INCLUDE_DIRS 
        "${CMAKE_BINARY_DIR}/_deps/FFTWDouble/include" 
        "${CMAKE_BINARY_DIR}/_deps/FFTWFloat/include"
        )
    set(FFTW_LIBRARIES 
         "${CMAKE_BINARY_DIR}/_deps/FFTWFloat/lib/libfftw3f.${FFTW_LIB_SUFFIX}"
         "${CMAKE_BINARY_DIR}/_deps/FFTWFloat/lib/libfftw3f_threads.${FFTW_LIB_SUFFIX}"
         "${CMAKE_BINARY_DIR}/_deps/FFTWDouble/lib/libfftw3.${FFTW_LIB_SUFFIX}"
         "${CMAKE_BINARY_DIR}/_deps/FFTWDouble/lib/libfftw3_threads.${FFTW_LIB_SUFFIX}"
         )
endif() 
        
