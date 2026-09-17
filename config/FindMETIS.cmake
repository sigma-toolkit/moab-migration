# Find the METIS includes and libraries
#
# METIS is a serial library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/index.html
#
# METIS_INCLUDE_DIR - where to find autopack.h
# METIS_LIBRARIES   - List of fully qualified libraries to link against.
# METIS_FOUND       - Do not attempt to use if "no" or undefined.

set (METIS_DIR "" CACHE PATH "Path to search for Metis header and library files")
set (METIS_FOUND NO CACHE INTERNAL "Found Metis components successfully." )

FIND_LIBRARY(METIS_LIBRARY metis
  HINTS
  ${METIS_DIR}
  ${METIS_DIR}/lib
  ${PARMETIS_DIR}
  ${PARMETIS_DIR}/lib
  )

FIND_PATH(METIS_INCLUDE_DIR metis.h
  HINTS
  ${METIS_DIR}
  ${METIS_DIR}/include
  )

IF (NOT METIS_FOUND)
  if ( METIS_INCLUDE_DIR AND METIS_LIBRARY )
    set( METIS_FOUND YES )
    SET(METIS_INCLUDES ${METIS_INCLUDE_DIR})
    SET(METIS_LIBRARIES ${METIS_LIBRARY})

    # Detect Metis version from metis.h
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_ver_major REGEX "^#define METIS_VER_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_ver_minor REGEX "^#define METIS_VER_MINOR[ \t]+[0-9]+")
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_ver_sub   REGEX "^#define METIS_VER_SUBMINOR[ \t]+[0-9]+")
    if (_metis_ver_major)
      string(REGEX REPLACE "^#define METIS_VER_MAJOR[ \t]+([0-9]+)" "\\1" METIS_VERSION_MAJOR "${_metis_ver_major}")
      string(REGEX REPLACE "^#define METIS_VER_MINOR[ \t]+([0-9]+)" "\\1" METIS_VERSION_MINOR "${_metis_ver_minor}")
      string(REGEX REPLACE "^#define METIS_VER_SUBMINOR[ \t]+([0-9]+)" "\\1" METIS_VERSION_SUBMINOR "${_metis_ver_sub}")
      set(METIS_VERSION "${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_SUBMINOR}")
      message(STATUS "        Version   : ${METIS_VERSION}")
    endif()

    # Metis >= 5.2.0 requires GKlib as an explicit link dependency
    if (METIS_VERSION AND NOT "${METIS_VERSION}" VERSION_LESS "5.2.0")
      find_library(GKLIB_LIBRARY GKlib
        HINTS
        ${METIS_DIR}
        ${METIS_DIR}/lib
      )
      if (GKLIB_LIBRARY)
        list(APPEND METIS_LIBRARIES ${GKLIB_LIBRARY})
        message(STATUS "        GKlib     : ${GKLIB_LIBRARY}")
      else()
        message(WARNING "Metis ${METIS_VERSION} requires GKlib but libGKlib was not found in ${METIS_DIR}/lib")
      endif()
    endif()

  else ( METIS_INCLUDE_DIR AND METIS_LIBRARY )
    set( METIS_FOUND NO )
    message("finding Metis failed, please try to set the var METIS_DIR")
  endif ( METIS_INCLUDE_DIR AND METIS_LIBRARY )
ENDIF (NOT METIS_FOUND)

mark_as_advanced(
  METIS_INCLUDES
  METIS_LIBRARIES
)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (
  METIS "Metis not found, check environment variables METIS_DIR"
  METIS_INCLUDES
  METIS_LIBRARIES
  )

if (METIS_FOUND)
  include(MOABTPLTargets)
  moab_declare_tpl_target(METIS::METIS METIS_INCLUDES METIS_LIBRARIES)
endif (METIS_FOUND)
