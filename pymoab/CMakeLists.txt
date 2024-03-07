# Find necessary packages
find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)

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

    # Create a target for each Cython file
    add_cython_target(${TARGET_NAME} OUTPUT_VAR X CXX)

    # Add the target to the python extension
    add_library(${TARGET_NAME} MODULE ${X})
    python_extension_module(${TARGET_NAME})

    # Link the target with the MOAB library
    target_link_libraries(${TARGET_NAME} ${MOAB_LIBS})

    # Install the target in the pymoab directory
    install(TARGETS ${TARGET_NAME} LIBRARY DESTINATION pymoab)
endforeach()