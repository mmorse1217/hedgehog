include_guard()
if (NOT DEFINED PETSC_LIBRARIES)
    # Download and unpack PETSc at configure time
    message("Downloading PETSc...")
    set(CMAKE_PETSC_DIR ${CMAKE_BINARY_DIR}/_deps/PETSc)
    configure_file(cmake/PETSc-configure.cmake ${CMAKE_PETSC_DIR}/CMakeLists.txt)
    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY "${CMAKE_PETSC_DIR}"
        )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build .
        WORKING_DIRECTORY "${CMAKE_PETSC_DIR}"
        )

    # lines 16-38, taken from Jed Brown's FindPETSc.cmake in cmake/
    # Make a temporary makefile to probe the PETSc configuration
    set (petsc_config_makefile "${CMAKE_PETSC_DIR}/Makefile.petsc")
    set (petsc_conf_variables "${CMAKE_PETSC_DIR}/lib/petsc/conf/variables")
    file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPETSc.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")

    # macro to query values of stored petsc variables
    macro (PETSC_GET_VARIABLE name var)
        set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
        execute_process (COMMAND make --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
            OUTPUT_VARIABLE ${var}
            RESULT_VARIABLE petsc_return)
    endmacro (PETSC_GET_VARIABLE)
    petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
    petsc_get_variable (LIBNAME            libname)
    petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
    petsc_get_variable (PETSC_CC_INCLUDES            petsc_include)

    # Set explicit include/library variables to match those set in
    # FindPETSc.cmake
    # TODO Figure out how to expose a proper Petsc::Petsc target in
    # FindPETSc.cmake; seems broken at the moment.
    set(PETSC_LIBRARIES 
        "${libname}"
        "${petsc_libs_external}"
        )
    set(PETSC_INCLUDES
        "${petsc_include}"
        )
    set(PETSC_INCLUDE_DIRS 
        "${petsc_include}"
        )

    message("PETSC_LIBRARIES " ${PETSC_LIBRARIES})


endif() 
