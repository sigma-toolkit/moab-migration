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

# Coupler workflows.  Executables keep their unprefixed names so the regression
# configs and any existing scripts still find them in ${EXECUTABLE_OUTPUT_PATH}.
if( MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_HAVE_TEMPESTREMAP )
  set( ES_PAR_TESTS parallel/imoab_coupler.cpp parallel/imoab_coupler_bilin.cpp
                    parallel/imoab_coupler_twohop.cpp parallel/imoab_domain_culling.cpp
                    parallel/imoab_map_r2l.cpp parallel/imoab_map_target.cpp
                    parallel/imoab_read_compute_map.cpp parallel/imoab_scrip.cpp
                    parallel/imoab_setwithGid.cpp parallel/imoab_phatm_ocn_coupler.cpp
                    parallel/imoab_spmv_bfb.cpp
                    parallel/imoab_dualmap_caas.cpp )

  # imoab_apg2_ol_coupler is the only test here that reads the downloaded E3SM
  # meshes rather than MeshFiles/unittest, and it has no fallback: without
  # them it fails on an unreadable input.  Build it only when the meshes are
  # actually on disk, which is how test/earthsystem/parallel/Makefile.am has
  # always gated it (if HAVE_EXTMESHES).
  if ( MOAB_HAVE_EXTMESHES )
    set( ES_PAR_TESTS ${ES_PAR_TESTS} parallel/imoab_apg2_ol_coupler.cpp )
  else ( MOAB_HAVE_EXTMESHES )
    message( STATUS "Skipping imoab_apg2_ol_coupler: the E3SM meshes it reads are not available" )
  endif ( MOAB_HAVE_EXTMESHES )

  foreach( fname ${ES_PAR_TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    string( REPLACE "parallel/" "" noparbase ${base} )
    add_executable( ${noparbase} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
    set_target_properties( ${noparbase} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( ${noparbase} MOAB ${TEMPESTREMAP_LIBRARIES} )
    add_test( ${noparbase} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${noparbase} ${MPIEXEC_POSTFLAGS} )
  endforeach()

  # Fortran coupler drivers.  Guarded on ENABLE_FORTRAN: without it the
  # top-level CMakeLists never calls enable_language(Fortran), and handing a
  # .F90 source to add_executable then fails for want of a Fortran compiler.
  if ( ENABLE_FORTRAN )
    set( ES_PARF_TESTS parallel/imoab_coupler_fortran.F90 parallel/imoab_coupler_caas.F90 )

    foreach( fname ${ES_PARF_TESTS} )
      string( REPLACE ".F90" "" base ${fname} )
      string( REPLACE "parallel/" "" noparbase ${base} )
      add_executable( ${noparbase} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname} )
      set_property( TARGET ${noparbase} PROPERTY LINKER_LANGUAGE Fortran )
      set_target_properties( ${noparbase} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
      # MOAB propagates MPI_CXX_LIBRARIES only, which does not define the mpi_*_
      # symbols an F90 unit references, so the MPI Fortran bindings have to be
      # named here as well.
      target_link_libraries( ${noparbase} MOAB ${TEMPESTREMAP_LIBRARIES} ${MPI_Fortran_LIBRARIES} )
      add_test( ${noparbase} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${noparbase} ${MPIEXEC_POSTFLAGS} )
    endforeach()
  endif ( ENABLE_FORTRAN )
endif()
