# Conservative remapping between climate grids, built on TempestRemap.
# Only built when TempestRemap is available.
include_directories( ${MOAB_BINARY_DIR}/src
                     ${MOAB_SOURCE_DIR}/src
                     ${MOAB_SOURCE_DIR}/src/climate/remapping
                     ${MOAB_BINARY_DIR}/src/parallel
                     ${TEMPESTREMAP_INCLUDES} )

if( MOAB_HAVE_HDF5 )
  include_directories(${HDF5_INCLUDES})
endif()
if( MOAB_HAVE_NETCDF )
  include_directories(${NETCDF_INCLUDES})
endif()
if( MOAB_HAVE_PNETCDF )
  include_directories(${PNETCDF_INCLUDES})
endif()
if( MOAB_HAVE_ZOLTAN )
  include_directories(${ZOLTAN_INCLUDES})
endif()

set( MOAB_REMAPPING_SRCS
        ApplyWeights.cpp
        TempestRemapper.cpp
        TempestOnlineMap.cpp
        TempestLinearRemap.cpp
        TempestOnlineMapIO.cpp
        IntegerReprosum.cpp
        netcdfcpp_par.hpp
   )

add_library( moab_climate_remapping OBJECT ${MOAB_REMAPPING_SRCS} )
set_target_properties( moab_climate_remapping PROPERTIES
                       COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
