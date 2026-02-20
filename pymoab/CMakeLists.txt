find_program(CYTHON "cython")

# Find Python 
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Base include directories (always needed)
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)

# Check for MPI support in MOAB
if (MPI_FOUND)
    message(STATUS "PyMOAB: MPI support detected in MOAB")
    
    # Get mpi4py include directory (required when MPI is enabled)
    find_package(mpi4py QUIET)
    if(NOT mpi4py_FOUND)
        message(WARNING "PyMOAB: mpi4py not found but MPI is enabled in MOAB. parallelcomm module will not be built.")
        set(MPI4PY_FOUND "False")
    else()
        execute_process(
            COMMAND ${Python_EXECUTABLE} -c "import mpi4py; import os; print(os.path.join(os.path.dirname(mpi4py.__file__), 'include'))"
            OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        include_directories(${MPI4PY_INCLUDE_DIR})
        include_directories(
            ${CMAKE_SOURCE_DIR}/src/parallel/moab
            ${CMAKE_SOURCE_DIR}/src/parallel
            ${CMAKE_BINARY_DIR}/src/parallel
        )
        set(PyMOAB_MPI_INCLUDES "-I${CMAKE_SOURCE_DIR}/src/parallel/moab;-I${CMAKE_SOURCE_DIR}/src/parallel;-I${CMAKE_BINARY_DIR}/src/parallel")
        set(MPI4PY_FOUND "True")
    endif()
else()
    message(STATUS "PyMOAB: MPI not found - parallelcomm module will not be built")
    set(MPI4PY_FOUND "False")
endif()

# Enhanced MPI capability detection
if(MPI4PY_FOUND)
    message(STATUS "PyMOAB: Testing MPI capabilities...")
    
    # Test basic MPI functionality
    try_compile(HAS_BASIC_MPI ${CMAKE_BINARY_DIR}/test_mpi
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test_mpi.cpp
        CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${MPI4PY_INCLUDE_DIR}
        LINK_LIBRARIES=mpi4py
    )
    
    # Test MPI I/O functionality
    try_compile(HAS_MPI_IO ${CMAKE_BINARY_DIR}/test_mpi_io
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test_mpi_io.cpp
        CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${MPI4PY_INCLUDE_DIR}
        LINK_LIBRARIES=mpi4py
    )
    
    # Set capability flags
    if(HAS_BASIC_MPI)
        set(MOAB_MPI_BASIC "True" CACHE INTERNAL "Basic MPI support available")
    else()
        set(MOAB_MPI_BASIC "False" CACHE INTERNAL "Basic MPI support available")
    endif()
    
    if(HAS_MPI_IO)
        set(MOAB_MPI_IO "True" CACHE INTERNAL "MPI I/O support available")
    else()
        set(MOAB_MPI_IO "False" CACHE INTERNAL "MPI I/O support available")
    endif()
    
    message(STATUS "PyMOAB: MPI capabilities - Basic: ${MOAB_MPI_BASIC}, I/O: ${MOAB_MPI_IO}")
else()
    set(MOAB_MPI_BASIC "False" CACHE INTERNAL "Basic MPI support available")
    set(MOAB_MPI_IO "False" CACHE INTERNAL "MPI I/O support available")
endif()

# Export MPI4PY_FOUND to be used in Cython
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py.in ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)

# Generate moab_mpi.pxd based on MPI availability
if(MPI4PY_FOUND)
    message(STATUS "PyMOAB: Generating moab_mpi.pxd with MPI support")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/moab_mpi.pxd.in ${CMAKE_CURRENT_BINARY_DIR}/moab_mpi.pxd @ONLY)
else()
    message(STATUS "PyMOAB: Generating empty moab_mpi.pxd (MPI disabled)")
    # Create an empty moab_mpi.pxd file (must be truly empty for Cython)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/moab_mpi.pxd "")
endif()

# Generate parallelcomm_defs.pxd based on MPI availability and capabilities
if(MPI4PY_FOUND)
    message(STATUS "PyMOAB: Generating parallelcomm_defs.pxd with ParallelComm declarations")
    
    # Set capability flags for template
    set(MOAB_MPI_FULL "${MOAB_MPI_FULL}")
    set(MOAB_MPI_BASIC "${MOAB_MPI_BASIC}")
    set(MOAB_MPI_IO "${MOAB_MPI_IO}")
    
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/parallelcomm_defs.pxd.in ${CMAKE_CURRENT_BINARY_DIR}/parallelcomm_defs.pxd @ONLY)
else()
    message(STATUS "PyMOAB: Generating empty parallelcomm_defs.pxd (MPI disabled)")
    # Create an empty parallelcomm_defs.pxd file (must be truly empty for Cython)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/parallelcomm_defs.pxd "")
endif()

# Get list of Cython files (exclude parallelcomm if MPI not available)
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
            set(MODULE_INCLUDES ${MPI4PY_INCLUDE_DIR})
            set(MODULE_MPI_INCLUDES ${PyMOAB_MPI_INCLUDES})
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
