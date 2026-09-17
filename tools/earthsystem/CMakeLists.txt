# E3SM remapping tools.  The parent gates this subdirectory on
# MOAB_BUILD_MBTEMPEST, which is only true when MOAB was configured with
# TempestRemap, so nothing here needs to re-test for it.
#
# add_moab_tool is a macro defined in tools/CMakeLists.txt before this
# subdirectory is added, so it is in scope here.

include_directories( ${CMAKE_SOURCE_DIR}/src/earthsystem/intx_mesh
                     ${CMAKE_SOURCE_DIR}/src/earthsystem/remapping )

# MBTEMPEST
add_moab_tool( mbtempest mbtempest.cpp MOAB_HAVE_MBTEMPEST "${TEMPESTREMAP_LIBRARIES}" )

# H5MTOSCRIP - upstream added this to the CMake build (tools/CMakeLists.txt at
# 0b62137e4); registered here because h5mtoscrip.cpp lives in this directory.
if ( MOAB_BUILD_H5MTOSCRIP )
  add_moab_tool( h5mtoscrip h5mtoscrip.cpp MOAB_HAVE_H5MTOSCRIP "${TEMPESTREMAP_LIBRARIES}" )
endif( )
