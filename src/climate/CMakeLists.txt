# Climate/E3SM-specific MOAB sources.
#
# This subtree isolates functionality that exists to support climate
# modeling workflows (E3SM component-coupler mesh migration, conservative
# remapping on the sphere, and the climate NetCDF grid readers/writers) from
# the general-purpose MOAB core and its domain-agnostic features.  Keeping it
# separate makes the climate stack easier to refactor and maintain without
# perturbing the core mesh database.
#
# Each component is gated on the dependency it needs and exports its object
# files to the parent scope, where they are folded into the single libMOAB.

set( MOAB_CLIMATE_OBJECTS )

# Climate NetCDF grid readers/writers.  Gated on NetCDF, matching the
# condition these sources were built under when they lived in src/io.
if( MOAB_HAVE_NETCDF )
  add_subdirectory(io)
  list( APPEND MOAB_CLIMATE_OBJECTS $<TARGET_OBJECTS:moab_climate_io> )
endif()

# Spherical mesh intersection, the geometric kernel underneath conservative
# remapping.  Always built, as it has no external dependencies.
add_subdirectory(intx_mesh)
list( APPEND MOAB_CLIMATE_OBJECTS $<TARGET_OBJECTS:moab_intx_mesh> )

# TempestRemap-based conservative remapping.
if( MOAB_HAVE_TEMPESTREMAP )
  add_subdirectory(remapping)
  list( APPEND MOAB_CLIMATE_OBJECTS $<TARGET_OBJECTS:moab_climate_remapping> )
endif()

set( MOAB_CLIMATE_OBJECTS ${MOAB_CLIMATE_OBJECTS} PARENT_SCOPE )
