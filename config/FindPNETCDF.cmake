#
# Find PNetCDF include directories and libraries
#
# PNETCDF_INCLUDES            - list of include paths to find netcdf.h
# PNETCDF_LIBRARIES           - list of libraries to link against when using NetCDF
# PNETCDF_FOUND               - Do not attempt to use PNetCDF if "no", "0", or undefined.

set (PNETCDF_DIR "" CACHE PATH "Path to search for PNetCDF header and library files" )
set (PNETCDF_FOUND NO CACHE INTERNAL "Found PNetCDF components successfully." )

IF (MOAB_HAVE_MPI AND ENABLE_PNETCDF)

  find_path( PNETCDF_INCLUDES pnetcdf.h
    ${PNETCDF_DIR}
    ${PNETCDF_DIR}/include
    ENV CPLUS_INCLUDE_PATH
    NO_DEFAULT_PATH
  )

find_library( PNETCDF_LIBRARIES
    NAMES pnetcdf libpnetcdf.a
    HINTS ${PNETCDF_DIR}
    ${PNETCDF_DIR}/lib64
    ${PNETCDF_DIR}/lib
    NO_DEFAULT_PATH
  )

  IF (NOT PNETCDF_FOUND)
    if ( PNETCDF_INCLUDES AND PNETCDF_LIBRARIES )
      set( PNETCDF_FOUND YES )
    else ( PNETCDF_INCLUDES AND PNETCDF_LIBRARIES )
      set( PNETCDF_FOUND NO )
      message("finding PNetCDF failed, please try to set the var PNETCDF_DIR")
    endif ( PNETCDF_INCLUDES AND PNETCDF_LIBRARIES )
  ENDIF (NOT PNETCDF_FOUND)

  mark_as_advanced(
    PNETCDF_INCLUDES
    PNETCDF_LIBRARIES
  )
ELSE (MOAB_HAVE_MPI AND ENABLE_PNETCDF)
  message (STATUS "Not configuring with PNetCDF since MPI installation not specified or explicitly disabled by user")
ENDIF (MOAB_HAVE_MPI AND ENABLE_PNETCDF)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PNETCDF "PNetCDF not found, check the CMake PNETCDF_DIR variable"
  PNETCDF_DIR PNETCDF_INCLUDES PNETCDF_LIBRARIES)
