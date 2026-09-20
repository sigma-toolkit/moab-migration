# Earth system tests, mirroring src/earthsystem: intersection, remapping, the
# NetCDF/SCRIP readers and the E3SM coupler workflows.  The parent adds this
# subdirectory after setting TEST_COMP_FLAGS, and MOAB_DEFINES and
# EXECUTABLE_OUTPUT_PATH come from the top level, so all three are in scope.
#
# Only the tests the CMake build already registered are registered here; the
# ones that were autotools-only (commgraph_test, par_intx_sph,
# test_mapweights_bcast, gcrm_par, read_esmf, mpastrvpart, ucdtrvpart,
# read_gcrm_nc, and the parallel NetCDF drivers read_nc_par and write_nc_par)
# stay that way rather than being enabled silently as part of a move.
# See test/Makefile.am and test/parallel/Makefile.am.

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/test
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/parallel )

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

# Serial tests over the intersection module.
set( ES_TESTS arc_intx_tests.cpp
              intx_in_plane_test.cpp
              intx_on_sphere_test.cpp )

if( MOAB_HAVE_NETCDF )
  list( APPEND ES_TESTS intx_rll_cs_sphere_test.cpp
                        intx_orientation_test.cpp )
endif()

foreach( fname ${ES_TESTS} )
  string( REPLACE ".cpp" "" base ${fname} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

# Earth system I/O formats.  These carry the same flags and link line as the
# rest of test/io, which is where they used to be registered.
if( MOAB_HAVE_NETCDF )
  set( ES_IO_TESTS io/read_nc.cpp
                   io/read_ucd_nc.cpp
                   io/read_mpas_nc.cpp
                   io/write_nc.cpp )

  foreach( fname ${ES_IO_TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    string( REPLACE "io/" "" base ${base} )
    add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
    target_link_libraries( ${base} MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif()

# Serial remapping drivers.
if( MOAB_HAVE_TEMPESTREMAP )
  set( ES_TR_TESTS imoab_remap.cpp test_remapping.cpp )

  foreach( fname ${ES_TR_TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( ${base} MOAB ${TEMPESTREMAP_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif()

if( MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_HAVE_HDF5_PARALLEL )
  set( ES_IMOAB_CORE_TESTS parallel/iMOABPtest2.cpp
                             parallel/iMOABSetWithGid.cpp )

  foreach( fname ${ES_IMOAB_CORE_TESTS} )
    get_filename_component( target ${fname} NAME_WE )
    add_executable( ${target} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
    set_target_properties( ${target} PROPERTIES
      COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST"
      CXX_STANDARD 17
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF )
    target_link_libraries( ${target} MOAB )
    add_test( ${target} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${target} ${MPIEXEC_POSTFLAGS} )
  endforeach()

  if ( ENABLE_FORTRAN )
    add_executable( iMOABPtest ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp parallel/iMOABPtest.F )
    set_property( TARGET iMOABPtest PROPERTY LINKER_LANGUAGE Fortran )
    set_target_properties( iMOABPtest PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( iMOABPtest MOAB ${MPI_Fortran_LIBRARIES} )
    add_test( iMOABPtest ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/iMOABPtest ${MPIEXEC_POSTFLAGS} )
  endif()
endif()

# Coupler workflows. C++ sources use the project spelling "iMOAB" throughout;
# executable and CTest names follow the renamed sources.
if( MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_HAVE_TEMPESTREMAP )
  set( ES_PAR_TESTS parallel/iMOABCoupler.cpp parallel/iMOABCouplerBilin.cpp
                    parallel/iMOABCouplerTwoHop.cpp parallel/iMOABDomainCulling.cpp
                    parallel/iMOABMapR2L.cpp parallel/iMOABMapTarget.cpp
                    parallel/iMOABReadComputeMap.cpp parallel/iMOABScrip.cpp
                    parallel/iMOABPhAtmOcnCoupler.cpp parallel/iMOABSpMVBFB.cpp
                    parallel/iMOABDualMapCAAS.cpp )

  # iMOABApg2OlCoupler is the only test here that reads the downloaded E3SM
  # meshes rather than MeshFiles/unittest, and it has no fallback: without
  # them it fails on an unreadable input. Build it only when the meshes exist.
  if ( MOAB_HAVE_EXTMESHES )
    list( APPEND ES_PAR_TESTS parallel/iMOABApg2OlCoupler.cpp )
  else ( MOAB_HAVE_EXTMESHES )
    message( STATUS "Skipping iMOABApg2OlCoupler: the E3SM meshes it reads are not available" )
  endif ( MOAB_HAVE_EXTMESHES )

  foreach( fname ${ES_PAR_TESTS} )
    get_filename_component( target ${fname} NAME_WE )
    add_executable( ${target} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
    set_target_properties( ${target} PROPERTIES
      COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST"
      CXX_STANDARD 17
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF )
    target_link_libraries( ${target} MOAB ${TEMPESTREMAP_LIBRARIES} )
    add_test( ${target} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${target} ${MPIEXEC_POSTFLAGS} )
  endforeach()

  if ( ENABLE_FORTRAN )
    set( ES_PARF_TESTS parallel/iMOABCouplerFortran.F90
                       parallel/iMOABCouplerCAAS.F90 )

    foreach( fname ${ES_PARF_TESTS} )
      get_filename_component( target ${fname} NAME_WE )
      add_executable( ${target} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
      set_property( TARGET ${target} PROPERTY LINKER_LANGUAGE Fortran )
      set_target_properties( ${target} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
      target_link_libraries( ${target} MOAB ${TEMPESTREMAP_LIBRARIES} ${MPI_Fortran_LIBRARIES} )
      add_test( ${target} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${target} ${MPIEXEC_POSTFLAGS} )
    endforeach()
  endif()
endif()
