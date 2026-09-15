# Component-coupler communication graph.  ParCommGraph records which source
# PEs send to which target PEs between two E3SM components, and is what the
# iMOAB send/receive mesh and tag migration calls are built on.  Only built
# when MPI is available; the parent gates this subdirectory on it.
#
# The public header is declared in src/CMakeLists.txt alongside the rest of
# MOAB_INSTALL_HEADERS, which is where every other component's headers are
# installed from.
include_directories( ${MOAB_BINARY_DIR}/src
                     ${MOAB_SOURCE_DIR}/src
                     ${MOAB_SOURCE_DIR}/src/parallel
                     ${MOAB_BINARY_DIR}/src/parallel )

if( MOAB_HAVE_HDF5 )
  include_directories(${HDF5_INCLUDES})
endif()

if( MOAB_HAVE_ZOLTAN )
  include_directories(${ZOLTAN_INCLUDES})
endif()

set( MOAB_EARTHSYSTEM_PARALLEL_SRCS ParCommGraph.cpp )

add_library( moab_earthsystem_parallel OBJECT ${MOAB_EARTHSYSTEM_PARALLEL_SRCS} )
set_target_properties( moab_earthsystem_parallel PROPERTIES
                       COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
