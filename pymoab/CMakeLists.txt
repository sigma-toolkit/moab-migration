# Find necessary packages
find_package(
    Python
    COMPONENTS Interpreter Development.Module NumPy
    REQUIRED
)
find_program(CYTHON "cython")

# Include the source and build directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Python_NumPy_INCLUDE_DIRS}
)

# Get a list of all Cython files in the pymoab directory
file(GLOB CYTHON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.pyx")

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
                -I ${Python_NumPy_INCLUDE_DIRS}
                "${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pyx" 
                --output-file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
    )

    # Creating a Python module
    python_add_library(
        ${TARGET_NAME} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.cxx"
        WITH_SOABI
    )

    # Add the MOAB library as a dependency
    add_dependencies(${TARGET_NAME} MOAB)

    # Link the target with the MOAB library
    target_link_libraries(${TARGET_NAME} PRIVATE MOAB)

    # Install the target in the pymoab directory
    install(TARGETS ${TARGET_NAME} DESTINATION pymoab)
endforeach()