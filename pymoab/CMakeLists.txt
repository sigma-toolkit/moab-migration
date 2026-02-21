find_program(CYTHON "cython")

#find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)

if (MPI_FOUND)
    message(STATUS "PyMOAB: MPI support detected in MOAB")
    
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import mpi4py; import os; print(os.path.join(os.path.dirname(mpi4py.__file__), 'include'))"
        OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE MPI4PY_RESULT
    )
    if(NOT MPI4PY_RESULT EQUAL 0)
        get_filename_component(PY_BIN_DIR ${Python_EXECUTABLE} DIRECTORY)
        get_filename_component(PY_PREFIX ${PY_BIN_DIR} DIRECTORY)
        
        execute_process(
            COMMAND ${Python_EXECUTABLE} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
            OUTPUT_VARIABLE PY_VER
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        unset(MPI4PY_INCLUDE_DIR)
        unset(MPI4PY_INCLUDE_DIR CACHE)

        find_path(MPI4PY_INCLUDE_DIR 
            NAMES mpi4py/mpi4py.h
            PATHS 
                "${PY_PREFIX}/lib/python${PY_VER}/site-packages/mpi4py/include"
                "${PY_PREFIX}/Lib/site-packages/mpi4py/include"
                "${Python_SITELIB}/mpi4py/include"
                "${Python_SITEARCH}/mpi4py/include"
            NO_DEFAULT_PATH
        )
        if(MPI4PY_INCLUDE_DIR)
            set(MPI4PY_RESULT 0)
            message(STATUS "PyMOAB: Found mpi4py include via fallback: ${MPI4PY_INCLUDE_DIR}")
        endif()
    endif()

    if(NOT MPI4PY_RESULT EQUAL 0)
        message(WARNING "PyMOAB: mpi4py not found but MPI is enabled in MOAB. parallelcomm module will not be built.")
        set(MPI4PY_FOUND "False")
    else()
        include_directories(${MPI4PY_INCLUDE_DIR})
        include_directories(
            ${CMAKE_SOURCE_DIR}/src/parallel/moab
            ${CMAKE_SOURCE_DIR}/src/parallel
            ${CMAKE_BINARY_DIR}/src/parallel
        )
        set(PyMOAB_MPI_INCLUDES "-I${CMAKE_SOURCE_DIR}/src/parallel/moab;-I${CMAKE_SOURCE_DIR}/src/parallel;-I${CMAKE_BINARY_DIR}/src/parallel")
        set(MPI4PY_FOUND "True")
        set(MOAB_MPI_BASIC "True")
        if(MOAB_HAVE_HDF5_PARALLEL)
            set(MOAB_MPI_IO "True")
        else()
            set(MOAB_MPI_IO "False")
        endif()
    endif()
else()
    message(STATUS "PyMOAB: MPI not found - parallelcomm module will not be built")
    set(MPI4PY_FOUND "False")
    set(MOAB_MPI_BASIC "False")
    set(MOAB_MPI_IO "False")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py.in ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)

if(MPI4PY_FOUND)
    message(STATUS "PyMOAB: Generating parallelcomm_defs.pxd with ParallelComm declarations")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/parallelcomm_defs.pxd.in ${CMAKE_CURRENT_BINARY_DIR}/parallelcomm_defs.pxd @ONLY)
else()
    message(STATUS "PyMOAB: Generating empty parallelcomm_defs.pxd (MPI disabled)")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/parallelcomm_defs.pxd "")
endif()

file(GLOB ALL_CYTHON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.pyx")

# Filter out parallelcomm if MPI not available
set(CYTHON_FILES "${ALL_CYTHON_FILES}")
if(NOT MPI4PY_FOUND)
    list(REMOVE_ITEM CYTHON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/parallelcomm.pyx")
    message(STATUS "PyMOAB: Excluding parallelcomm.pyx (MPI disabled or mpi4py not found)")
endif()

# Run a for loop to create a target for each Cython file
foreach(CYTHON_FILE ${CYTHON_FILES})
    # Extract the base filename (without extension)
    get_filename_component(TARGET_NAME ${CYTHON_FILE} NAME_WE)

    # Determine include directories for this module
    set(MODULE_INCLUDES "")
    set(MODULE_MPI_INCLUDES "")
    if(${TARGET_NAME} STREQUAL "parallelcomm")
        if(MPI4PY_FOUND)
            set(MODULE_INCLUDES "-I${MPI4PY_INCLUDE_DIR}")
            set(MODULE_MPI_INCLUDES ${PyMOAB_MPI_INCLUDES})
            # Add Python site-packages to Cython includes so it can find mpi4py.pxd in isolated builds
            if(Python_SITELIB)
                set(MODULE_MPI_INCLUDES ${MODULE_MPI_INCLUDES} -I ${Python_SITELIB})
            endif()
            if(Python_SITEARCH)
                set(MODULE_MPI_INCLUDES ${MODULE_MPI_INCLUDES} -I ${Python_SITEARCH})
            endif()
        else()
            continue()
        endif()
    endif()

    # Build C++ file
    add_custom_command(
        OUTPUT ${TARGET_NAME}.cxx
        DEPENDS ${CYTHON_FILE}
        VERBATIM
        COMMAND "${CYTHON}" --cplus
                -I ${CMAKE_SOURCE_DIR}/src
                -I ${CMAKE_BINARY_DIR}/src
                -I ${CMAKE_CURRENT_SOURCE_DIR}
                -I ${CMAKE_CURRENT_BINARY_DIR}
                -I ${Python_NumPy_INCLUDE_DIRS}
                ${MODULE_INCLUDES}
                ${MODULE_MPI_INCLUDES}
                "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pyx"
                --output-file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
    )

    # Creating a Python module
    python_add_library(
        ${TARGET_NAME} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
        WITH_SOABI
    )

    # Add the macro definition for the specific target
    target_compile_definitions(${TARGET_NAME} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

    # Add the MOAB library as a dependency
    add_dependencies(${TARGET_NAME} MOAB)

    # Link the target with the MOAB library
    target_link_libraries(${TARGET_NAME} PRIVATE MOAB)

    # Set the install RPATH
    set_target_properties(${TARGET_NAME} PROPERTIES INSTALL_RPATH "${PYMOAB_PYTHON_MODULE_RPATH}")

    # Install the target in the pymoab directory
    install(TARGETS ${TARGET_NAME} DESTINATION .)
endforeach()

# Install generated config.py
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.py DESTINATION . RENAME config.py)
