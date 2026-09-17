#
# Find TempestRemap include directories and libraries
#
# TEMPESTREMAP_INCLUDES            - list of include paths to find TempestRemapAPI.h
# TEMPESTREMAP_LIBRARIES           - list of libraries to link against when using TempestRemap
# TEMPESTREMAP_FOUND               - Do not attempt to use TempestRemap if "no", "0", or undefined.

set (TEMPESTREMAP_DIR "" CACHE PATH "Path to search for TempestRemap header and library files" )
set (TEMPESTREMAP_FOUND NO CACHE INTERNAL "Found TempestRemap components successfully." )

IF (TEMPESTREMAP_DIR)

  find_path( TEMPESTREMAP_INCLUDE_DIR TempestRemapAPI.h
    ${TEMPESTREMAP_DIR}/include
  )

  find_library( TEMPESTREMAP_LIBRARY
    NAMES TempestRemap
    HINTS ${TEMPESTREMAP_DIR}
    ${TEMPESTREMAP_DIR}/lib
    ${TEMPESTREMAP_DIR}/lib64
  )

  if ( TEMPESTREMAP_INCLUDE_DIR AND TEMPESTREMAP_LIBRARY )
    set( TEMPESTREMAP_FOUND YES )
    if (EXISTS ${TEMPESTREMAP_INCLUDE_DIR}/OfflineMap.h)
      SET(TEMPESTREMAP_INCLUDES ${TEMPESTREMAP_INCLUDE_DIR})
      SET(TEMPESTREMAP_LIBRARIES ${TEMPESTREMAP_LIBRARY})
    endif(EXISTS ${TEMPESTREMAP_INCLUDE_DIR}/OfflineMap.h)
  else ( TEMPESTREMAP_INCLUDE_DIR AND TEMPESTREMAP_LIBRARY )
    set( TEMPESTREMAP_FOUND NO )
    message("finding TempestRemap failed, please try to set the var TEMPESTREMAP_DIR")
  endif ( TEMPESTREMAP_INCLUDE_DIR AND TEMPESTREMAP_LIBRARY )
ENDIF (TEMPESTREMAP_DIR)

mark_as_advanced(
  TEMPESTREMAP_DIR
  TEMPESTREMAP_INCLUDES
  TEMPESTREMAP_LIBRARIES
)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (TEMPESTREMAP "TempestRemap not found, check environment variables TEMPESTREMAP_DIR"
  TEMPESTREMAP_DIR TEMPESTREMAP_INCLUDES TEMPESTREMAP_LIBRARIES)

if (TEMPESTREMAP_FOUND)
  include(MOABTPLTargets)
  moab_declare_tpl_target(TempestRemap::TempestRemap TEMPESTREMAP_INCLUDES TEMPESTREMAP_LIBRARIES)
endif (TEMPESTREMAP_FOUND)

