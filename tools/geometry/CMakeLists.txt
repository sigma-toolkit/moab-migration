# Check that CGM has been enabled

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src )

if ( CGM_FOUND )
  include_directories( ${CGM_INCLUDE_DIRS} )
endif( )


if ( MOAB_HAVE_CGM )
  if ( MOAB_ENABLE_TESTING )
    target_link_libraries( test_geom ${CGM_LIBRARIES} )
  endif()
  option ( MOAB_BUILD_CGM2MOAB    "Build cgm2moab?"  OFF )
  if(MOAB_BUILD_CGM2MOAB)
    add_executable( cgm2moab cgm2moab.cpp
                             obb_analysis.cpp
                             cgm2moab.hpp )
    target_link_libraries( cgm2moab
                           MOAB
                           ${CGM_LIBRARIES} )

    if(SKBUILD)
      # Set the install RPATH
      set_target_properties(cgm2moab PROPERTIES INSTALL_RPATH "${PYMOAB_BINARY_RPATH}")
    endif()
    install(
      TARGETS cgm2moab
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
  endif()
endif( )

