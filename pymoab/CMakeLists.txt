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

    # Link the target with the MOAB library
    target_link_libraries(${TARGET_NAME} PRIVATE ${MOAB_LIBS})

    # Install the target in the pymoab directory
    install(TARGETS ${TARGET_NAME} DESTINATION pymoab)
endforeach()


if (ENABLE_PYMOAB)
    message(WARNING "ENABLE_PYMOAB is deprecated and will be removed in the next release.")
    message(WARNING "Consider installing PyMOAB via pip instead.")

    if (PYMOAB_PREFIX)
        message(WARNING "PYMOAB_PREFIX is deprecated and will be removed in the next release.")
        message(WARNING "Consider installing PyMOAB via pip instead.")
        # Install PyMOAB in the specified prefix
        execute_process(
            COMMAND 
            ${PYTHON_EXECUTABLE} -m pip install ../. -C wheel.cmake=false -C wheel.install-dir=${PYMOAB_PREFIX}
            RESULT_VARIABLE PYMOAB_INSTALL_RESULT
        )
    else()
        # Install PyMOAB in the default prefix
        execute_process(
            COMMAND 
            ${PYTHON_EXECUTABLE} -m pip install ../. -C wheel.cmake=false
            RESULT_VARIABLE PYMOAB_INSTALL_RESULT
        )
    endif()

    if(NOT PYMOAB_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "Installation of PyMOAB failed. Stopping CMake build.")
    endif()
endif()
