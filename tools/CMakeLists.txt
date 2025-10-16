# Documentation
# Optional tools
set( MOAB_BUILD_MBCONVERT    TRUE  CACHE INTERNAL  "Build the MOAB mesh converter tool?"           )
set( MOAB_BUILD_HEXMODOPS    TRUE  CACHE INTERNAL  "Build the hexahedral modops tool?"             )
set( MOAB_BUILD_MBSIZE       TRUE  CACHE INTERNAL  "Build the mbsize tool?"                        )
set( MOAB_BUILD_MBMERGE      TRUE  CACHE INTERNAL  "Build the mbmerge tool?"                       )
set( MOAB_BUILD_MBMEM        TRUE  CACHE INTERNAL  "Build the mbmem tool?"                         )
set( MOAB_BUILD_MBSKIN       TRUE  CACHE INTERNAL  "Build the mbskin tool?"                        )
set( MOAB_BUILD_MBDEPTH      TRUE  CACHE INTERNAL  "Build the mbdepth tool?"                       )
set( MOAB_BUILD_MBTAGPROP    TRUE  CACHE INTERNAL  "Build the mbtagprop tool?"                     )
set( MOAB_BUILD_MBGSETS      TRUE  CACHE INTERNAL  "Build the mbgsets tool?"                       )
set( MOAB_BUILD_SPHEREDECOMP TRUE  CACHE INTERNAL  "Build the sphere decomposition tool?"          )
set( MOAB_BUILD_MBSURFPLOT   TRUE  CACHE INTERNAL  "Build the mbsurfplot application?"             )
set( MOAB_BUILD_SPHEREDECOMP TRUE  CACHE INTERNAL  "Build the spheredecomp application?"           )
set( MOAB_BUILD_MBPART       FALSE CACHE INTERNAL  "Build the mbpart (Zoltan-based) partitioner?"  )
set( MOAB_BUILD_MBSLAVEPART  TRUE  CACHE INTERNAL  "Build the slave partitioner tool?"             )
set( MOAB_BUILD_MBCOUPLER    TRUE  CACHE INTERNAL  "Build the mesh coupler tool?"                  )
set( MOAB_BUILD_MBHONODES    TRUE  CACHE INTERNAL  "Build the hex8 to hex27 converter tool?"       )
set( MOAB_BUILD_MBUMR        TRUE  CACHE INTERNAL  "Build the uniform mesh refinement tool?"       )
set( MOAB_BUILD_MBQUALITY    TRUE  CACHE INTERNAL  "Build the quality metrics tool?"               )
set( MOAB_BUILD_MBTEMPEST    FALSE CACHE INTERNAL  "Build the TempestRemap coupler tool?"          )
set( MOAB_BUILD_MBCMPMAPS    FALSE CACHE INTERNAL  "Build the remapping weights comparison tool?"  )
set( MOAB_BUILD_MBCMPFILES   FALSE CACHE INTERNAL  "Build the h5m mesh/data comparison tool?"      )
set( MOAB_BUILD_MBVISUMAP    FALSE CACHE INTERNAL  "Build the remapping weight visualization tool?")
set( MOAB_BUILD_MBCLEANTAGS  FALSE CACHE INTERNAL  "Build the tag modification and cleaning tool?" )

include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/parallel
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/parallel
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (MOAB_HAVE_TEMPESTREMAP)
  include_directories( ${TEMPESTREMAP_INCLUDES} ${EIGEN3_DIR} )
  set(MOAB_BUILD_MBTEMPEST ${MOAB_HAVE_TEMPESTREMAP})
endif (MOAB_HAVE_TEMPESTREMAP)

# Build mbpart if Zoltan or Metis is configured
if (MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS)
  set(MOAB_BUILD_MBPART TRUE)
endif(MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS)

set(MOAB_BUILD_MBCMPFILES ${MOAB_HAVE_HDF5})

set( TOOLS_LINK_LIBS MOAB ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

macro(add_moab_tool target source_file variable)
  add_executable(${target} ${source_file})
  set_source_files_properties(${source_file} COMPILE_FLAGS "${MOAB_DEFINES}")
  target_link_libraries(${target} ${TOOLS_LINK_LIBS} ${ARGN})
  set( ${variable} 1 )
  if(PYMOAB_BINARY_RPATH)
    # Set the install RPATH
    set_target_properties(${target} PROPERTIES INSTALL_RPATH "${PYMOAB_BINARY_RPATH}")
  endif()
  install(
    TARGETS ${target}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  generate_and_install_python_script(${target})
endmacro()

# MBSIZE
if ( MOAB_BUILD_MBSIZE )
  add_moab_tool( mbsize "size.cpp;measure.cpp" MOAB_HAVE_MBSIZE )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbsize)
endif ( )

# MBMERGE
if ( MOAB_BUILD_MBMERGE )
  add_moab_tool( mbmerge merge.cpp MOAB_HAVE_MBMERGE )
endif ( )


# MBMEM
if ( MOAB_BUILD_MBMEM )
  add_moab_tool( mbmem mbmem.cpp MOAB_HAVE_MBMEM )
endif ( )

# MBCONVERT
if ( MOAB_BUILD_MBCONVERT )
  add_moab_tool( mbconvert convert.cpp MOAB_HAVE_MBCONVERT "${TEMPESTREMAP_LIBRARIES}" )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbconvert)
  INSTALL(FILES mbconvert.man DESTINATION ${CMAKE_INSTALL_MANDIR}/man1 RENAME mbconvert.1)
endif ( )

# MBHONODES
if ( MOAB_BUILD_MBHONODES )
  add_moab_tool( mbhonodes honodes.cpp MOAB_HAVE_MBHONODES )
endif ( )

# MBDEPTH
if ( MOAB_BUILD_MBDEPTH )
  add_moab_tool( mbdepth depth.cpp MOAB_HAVE_MBDEPTH )
endif ( )

# MBSKIN
if ( MOAB_BUILD_MBSKIN )
  add_moab_tool( mbskin skin.cpp MOAB_HAVE_MBSKIN )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbskin)
endif ( )

# MBSURFPLOT
if ( MOAB_BUILD_MBSURFPLOT )
  add_moab_tool( mbsurfplot surfplot.cpp MOAB_HAVE_MBSURFPLOT )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbsurfplot)
endif ( )

# MBTAGPROP
if ( MOAB_BUILD_MBTAGPROP )
  add_moab_tool( mbtagprop "parse.cpp;propagate_tags.cpp" MOAB_HAVE_MBTAGPROP )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbtagprop)
endif ( )

# MBGSETS
if ( MOAB_BUILD_MBGSETS )
  add_moab_tool( mbgsets gsets.cpp MOAB_HAVE_MBGSETS )
endif( )

# MBUMR
if ( MOAB_BUILD_MBUMR )
  add_moab_tool( mbumr umr.cpp MOAB_HAVE_MBUMR )
endif( )

# MBSLAVEPART
if ( MOAB_BUILD_MBSLAVEPART )
  add_moab_tool( mbslavepart mbslavepart.cpp MOAB_HAVE_MBSLAVEPART )
endif( )

# MBTEMPEST
if ( MOAB_BUILD_MBTEMPEST )
  add_moab_tool( mbtempest mbtempest.cpp MOAB_HAVE_MBTEMPEST "${TEMPESTREMAP_LIBRARIES}")
endif( )

if (MOAB_HAVE_NETCDF OR MOAB_HAVE_PNETCDF)

  # MOAB_BUILD_MBCMPMAPS
  if ( MOAB_BUILD_MBCMPMAPS )
    add_moab_tool( mbcmpmaps compareMaps.cpp MOAB_HAVE_MBCMPMAPS )
  endif( )

  # MBVISUMAP
  if ( MOAB_BUILD_MBVISUMAP )
    add_moab_tool( mbvisumap visuMapVtk.cpp MOAB_HAVE_MBVISUMAP )
  endif( )

endif (MOAB_HAVE_NETCDF OR MOAB_HAVE_PNETCDF)

# MBCMPFILES
if ( MOAB_BUILD_MBCMPFILES )
  add_moab_tool( mbcmpfiles compareFiles.cpp MOAB_HAVE_MBCMPFILES )
endif( )

# MBQUALITY
if ( MOAB_BUILD_MBQUALITY )
  add_moab_tool( mbquality quality.cpp MOAB_HAVE_MBQUALITY )
endif( )

# MOAB_BUILD_MBCLEANTAGS
if ( MOAB_BUILD_MBCLEANTAGS )
  add_moab_tool( mbcleantags cleanTags.cpp MOAB_BUILD_MBCLEANTAGS )
endif( )


# MBPART
if ( MOAB_BUILD_MBPART )
  set(MOAB_BUILD_MBPART ON)
  # Build the mbpart executable
  add_moab_tool( mbpart mbpart.cpp MOAB_HAVE_MBPART "${ZOLTAN_LIBRARIES}" "${PARMETIS_LIBRARIES}" "${METIS_LIBRARIES}" )
  include_directories( ${ZOLTAN_INCLUDES} ${PARMETIS_INCLUDES} ${METIS_INCLUDES} )
endif ( MOAB_BUILD_MBPART )

# MBCoupler
if ( MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
  add_subdirectory( mbcoupler )
endif ( )

# hexmodops
if( MOAB_BUILD_HEXMODOPS )
  add_moab_tool( hexmodops makeops.cpp MOAB_HAVE_HEXMODOPS )
  SET(MOAB_DOCS ${MOAB_DOCS} README.hexmodops)
endif()

# MBSPHEREDECOMP
if ( MOAB_BUILD_SPHEREDECOMP )
  add_moab_tool( spheredecomp "SphereDecomp.cpp;SphereDecomp.hpp;sphere_decomp.cpp" MOAB_HAVE_SPHEREDECOMP )
  SET(MOAB_DOCS ${MOAB_DOCS} README.spheredecomp)
endif ( )

install(FILES ${MOAB_DOCS} DESTINATION ${MOABDocLocation})
