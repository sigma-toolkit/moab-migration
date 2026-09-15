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

# h5mtoscrip is built here by the autotools build under the same conditional,
# but has never been part of the CMake build.  Left as-is rather than added
# silently; see tools/Makefile.am.
