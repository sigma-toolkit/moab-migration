include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_SOURCE_DIR}/src/moab/point_locater/lotte
    ${CMAKE_SOURCE_DIR}/test
)

set( MBCOUPLER_SRCS Coupler.cpp
                    ElemUtil.cpp )

set( MBCOUPLER_INCLUDE_SRCS Coupler.hpp
                            ElemUtil.hpp )

set_source_files_properties( ${MBCOUPLER_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES} -DIS_BUILDING_MB"
)

add_library( mbcoupler
  ${MBCOUPLER_SRCS}
  ${MBCOUPLER_INCLUDE_SRCS} )
add_library( MOAB::mbcoupler ALIAS mbcoupler )

target_include_directories( mbcoupler
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mbcoupler>" )

# Link the same imported targets libMOAB does, so that mbcoupler's exported
# interface is target names too rather than the absolute paths BLAS/LAPACK and
# MPI resolve to on the build machine.  MOAB is PUBLIC because Coupler.hpp
# includes MOAB headers; the numerics are an implementation detail.
target_link_libraries( mbcoupler PUBLIC MOAB )
if ( MOAB_HAVE_LAPACK AND TARGET LAPACK::LAPACK )
  target_link_libraries( mbcoupler PRIVATE LAPACK::LAPACK )
endif ( )
if ( MOAB_HAVE_BLAS AND TARGET BLAS::BLAS )
  target_link_libraries( mbcoupler PRIVATE BLAS::BLAS )
endif ( )
if ( MOAB_HAVE_MPI AND TARGET MPI::MPI_CXX )
  target_link_libraries( mbcoupler PUBLIC MPI::MPI_CXX )
endif ( )

if ( MOAB_HAVE_NETCDF )
  set(netcdf_progs update_coords)
endif()

set( MOAB_HAVE_MBCOUPLER 1 )
install(
  TARGETS mbcoupler
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install(
  FILES ${MBCOUPLER_INCLUDE_SRCS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mbcoupler
  )


if ( MOAB_HAVE_MPI )
  # Exercise the coupler and not just findpt:
  add_executable( mbcoupler_test mbcoupler_test.cpp )
  set_source_files_properties( mbcoupler_test.cpp
    COMPILE_FLAGS "${MOAB_DEFINES} -DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest" )
  target_link_libraries( mbcoupler_test mbcoupler MOAB ${MPI_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( TestMBCoupler ${EXECUTABLE_OUTPUT_PATH}/mbcoupler_test )

  add_executable( ssn_test ssn_test.cpp )
  set_source_files_properties( ssn_test.cpp
    COMPILE_FLAGS "${MOAB_DEFINES} -DTEST -DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest" )
  target_link_libraries( ssn_test mbcoupler MOAB ${MPI_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( TestSSN ${EXECUTABLE_OUTPUT_PATH}/ssn_test )
endif ( MOAB_HAVE_MPI )
