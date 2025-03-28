#
# FindNETCDF.cmake
#
# This module finds the NETCDF library and headers for C.
#
# It sets the following variables:
#
#   NETCDF_FOUND       - True if NETCDF is found.
#   NETCDF_INCLUDES    - Path to the NETCDF C headers.
#   NETCDF_LIBRARIES   - Path to the NETCDF C library.
#   NETCDF_VERSION     - The found version of NETCDF.
#
set (NETCDF_DIR "/usr" CACHE PATH "Path to search for NETCDF header and library files" )
set (NETCDF_FOUND NO CACHE INTERNAL "Found NETCDF components successfully." )

# Query nc-config script if available
find_program(NC_CONFIG_EXECUTABLE NAMES nc-config)
if (NC_CONFIG_EXECUTABLE)
    execute_process(
        COMMAND ${NC_CONFIG_EXECUTABLE} --prefix
        OUTPUT_VARIABLE NETCDF_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Check if the user specified NETCDF_DIR matches the prefix from nc-config
    # This is a workaround for the fact that nc-config does not support --with-netcdf-dir
    # and we need to ensure that the user is not trying to use a different version of NetCDF
    # than the one that nc-config is pointing to.
    # This is important for cross-compilation scenarios where the user might have
    # a different version of NetCDF installed on the host system than the one
    # that is being used for cross-compilation.
    # This check is not perfect, but it is better than nothing.
    # The user can always override the NETCDF_DIR variable if they know what they are doing.
    if (NOT "${NETCDF_DIR}" STREQUAL "${NETCDF_PREFIX}")
        message(FATAL_ERROR "Error: User specified NETCDF_DIR (${NETCDF_DIR}) does not match the prefix from nc-config (${NETCDF_PREFIX})")
    endif()

    execute_process(
        COMMAND ${NC_CONFIG_EXECUTABLE} --includedir
        OUTPUT_VARIABLE NETCDF_INCLUDES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${NC_CONFIG_EXECUTABLE} --libs
        OUTPUT_VARIABLE NETCDF_LIBRARIES_RAW
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    include(ResolveCompilerPaths)
    # Resolve libraries from the raw string
    RESOLVE_LIBRARIES(NETCDF_LIBRARIES "${NETCDF_LIBRARIES_RAW}")

else(NC_CONFIG_EXECUTABLE)

    find_package(PkgConfig QUIET)
    pkg_check_modules(PC_NETCDF QUIET netcdf)

    find_path(NETCDF_INCLUDES
        NAMES netcdf.h
        HINTS ${PC_NETCDF_INCLUDE_DIRS} ${NETCDF_DIR}/include
    )

    find_library(NETCDF_LIBRARIES
        NAMES netcdf
        HINTS ${PC_NETCDF_LIBRARY_DIRS} ${NETCDF_DIR}/lib
        NO_DEFAULT_PATH
    )

endif(NC_CONFIG_EXECUTABLE)

# Get the version of NETCDF if found
if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
    execute_process(
        COMMAND ${CMAKE_C_COMPILER} -E -dM ${NETCDF_INCLUDES}/netcdf.h
        OUTPUT_VARIABLE _nc_defines
        ERROR_QUIET
    )
    string(REGEX MATCH "#define NC_VERSION_MAJOR ([0-9]+)" _major_match "${_nc_defines}")
    string(REGEX MATCH "#define NC_VERSION_MINOR ([0-9]+)" _minor_match "${_nc_defines}")
    string(REGEX MATCH "#define NC_VERSION_PATCH ([0-9]+)" _patch_match "${_nc_defines}")
    if (_major_match AND _minor_match AND _patch_match)
        set(NETCDF_VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
    endif (_major_match AND _minor_match AND _patch_match)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF
    REQUIRED_VARS NETCDF_INCLUDES NETCDF_LIBRARIES
    VERSION_VAR NETCDF_VERSION
)

mark_as_advanced(NETCDF_INCLUDES NETCDF_LIBRARIES)
