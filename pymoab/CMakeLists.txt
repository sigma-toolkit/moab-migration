# Find necessary packages
find_package(
    Python
    COMPONENTS Interpreter Development.Module
    REQUIRED
)
find_program(CYTHON "cython")

# Include the source and build directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Get a list of all Cython files in the pymoab directory
file(GLOB CYTHON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.pyx")

if (ENABLE_PYMOAB)
    find_package(Python COMPONENTS NumPy REQUIRED)

    # Get a list of all Python files in the pymoab directory
    file(GLOB PY_FILES "${CMAKE_SOURCE_DIR}/pymoab/*.py")
    file(GLOB PXD_FILES "${CMAKE_SOURCE_DIR}/pymoab/*.pxd")

    # Set the PyMOAB prefix
    if (NOT PYMOAB_PREFIX)
        set(PYMOAB_PREFIX ${CMAKE_INSTALL_PREFIX})
    endif()
    message(STATUS "PyMOAB prefix: ${PYMOAB_PREFIX} (DEPRECATED)")

    # Get the Python path directory name
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_paths()['purelib'].split('/')[-1])"
                OUTPUT_VARIABLE python_path
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    
    # Set the where to copy the files
    set(PYMOAB_LIB "/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/${python_path}/pymoab")
    
    # Copy the Python files to the Current binary directory
    file(COPY ${PY_FILES} ${PXD_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${PYMOAB_LIB})

    # Install the files to the PyMOAB prefix
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PYMOAB_LIB} DESTINATION ${PYMOAB_PREFIX}/${PYMOAB_LIB})
    
    message(WARNING "ENABLE_PYMOAB is deprecated and will be removed in the next release. Consider installing PyMOAB via pip instead.")
endif()

# Set the MOAB library name
set(MOAB_LIBS 
    ${CMAKE_BINARY_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}MOAB${CMAKE_SHARED_LIBRARY_SUFFIX}
)

# Run a for loop to create a target for each Cython file
foreach(CYTHON_FILE ${CYTHON_FILES})
    # Extract the base filename (without extension)
    get_filename_component(TARGET_NAME ${CYTHON_FILE} NAME_WE)

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
                "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pyx" 
                --output-file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
    )

    # Creating a Python module
    python_add_library(
        ${TARGET_NAME} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
        WITH_SOABI
    )

    if (SKBUILD)
        # Link the target with the MOAB library
        target_link_libraries(${TARGET_NAME} PRIVATE ${MOAB_LIBS})
        
        # Install the target in the pymoab directory
        install(TARGETS ${TARGET_NAME} DESTINATION pymoab)
    else()
        # Install the target in the deprecated pymoab directory
        install(TARGETS ${TARGET_NAME} DESTINATION ${PYMOAB_PREFIX}/${PYMOAB_LIB})
        
        # Create shortcut in current binary directory for testing
        add_custom_target(pymoab-local-install_${TARGET_NAME} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib/${TARGET_NAME} ${CMAKE_CURRENT_BINARY_DIR}/${PYMOAB_LIB}/${TARGET_NAME}
        DEPENDS ${TARGET_NAME})
    endif()
endforeach()
