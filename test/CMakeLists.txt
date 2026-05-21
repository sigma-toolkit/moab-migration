include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_BINARY_DIR}/src/parallel )

add_subdirectory(io)
if( MOAB_HAVE_HDF5 )
  add_subdirectory(h5file)
endif()
add_subdirectory(obb)
add_subdirectory(perf)

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

set( TESTS range_test.cpp
           verdict_test.cpp
           test_Matrix3.cpp
           scdseq_test.cpp
           scd_test_partn.cpp
           test_adj.cpp
           GeomUtilTests.cpp
           oriented_box_test.cpp
           adaptive_kd_tree_tests.cpp
           kd_tree_test.cpp
           bsp_tree_test.cpp
           reorder_test.cpp
           elem_eval_test.cpp
           VarLenTagTest.cpp
           TagTest.cpp
           spatial_locator_test.cpp
           bsp_tree_poly_test.cpp
           test_prog_opt.cpp
           coords_connect_iterate.cpp
           urefine_mesh_test.cpp
           adj_moab_test.cpp
           test_boundbox.cpp
           spherical_area_test.cpp
           arc_intx_tests.cpp
           intx_in_plane_test.cpp
           intx_on_sphere_test.cpp
           dual_test.cpp )

if(MOAB_HAVE_HDF5)
  set( TESTS ${TESTS}
             lloyd_smoother_test.cpp
             mhdf_public_test.cpp
             mergemesh_test.cpp
             mbfacet_test.cpp
             mbground_test.cpp
             gttool_test.cpp
             restore_topo_geom_incl_test.cpp
             crop_vol_test.cpp
	     gqt_simple_test.cpp
	     gqt_rayfire_test.cpp
	     gqt_pointinvol_test.cpp
            gqt_findvol_test.cpp
	     test_geom_gqt.cpp)
  set( HDF5_C_TESTS  imoab_test.c )
endif()

if(MOAB_HAVE_NETCDF)
  set ( TESTS ${TESTS}
              intx_rll_cs_sphere_test.cpp )
endif()

foreach( fname ${HDF5_C_TESTS} )
  string( REPLACE ".c" "" base ${fname} )
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" base ${fname} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

if(MOAB_HAVE_TEMPESTREMAP)
  set( MBTR_TESTS imoab_remap.cpp test_remapping.cpp )
  
  foreach( fname ${MBTR_TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( ${base} MOAB ${TEMPESTREMAP_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()

endif()

if(ENABLE_PYMOAB)
  ADD_SUBDIRECTORY(pymoab)
  find_package(
    Python
    COMPONENTS Interpreter
    REQUIRED
  )
  SET(PYMOAB_TESTS test_core test_hom test_scd test_range test_tag test_topo_util test_skinner)
  FOREACH(pymoab_test ${PYMOAB_TESTS})
    ADD_TEST(NAME pymoab_${pymoab_test}
      COMMAND "env" "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib/:$ENV{LD_LIBRARY_PATH}" "PYTHONPATH=${PROJECT_BINARY_DIR}/Testing:$ENV{PYTHONPATH}" "${Python_EXECUTABLE}" "${pymoab_test}.py"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/pymoab/")
  ENDFOREACH(pymoab_test)
endif()

add_executable( TestTypeSequenceManager ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp TestTypeSequenceManager.cpp)
set_target_properties( TestTypeSequenceManager PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST -DIS_BUILDING_MB" )

set(TestTypeSequenceManager_LIBS MOAB )
if(${Threads_FOUND})
  list(APPEND TestTypeSequenceManager_LIBS Threads::Threads)
endif()
target_link_libraries(TestTypeSequenceManager ${TestTypeSequenceManager_LIBS})

add_test( TestTypeSequenceManager ${EXECUTABLE_OUTPUT_PATH}/TestTypeSequenceManager )

#Test_MBMeshSet.cpp
add_executable( Test_MBMeshSet ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp Test_MBMeshSet.cpp)
set_target_properties( Test_MBMeshSet PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST -DIS_BUILDING_MB" )
target_link_libraries( Test_MBMeshSet MOAB  )
add_test( Test_MBMeshSet ${EXECUTABLE_OUTPUT_PATH}/Test_MBMeshSet )

if(ZCAT_EXE AND NOT WIN32)
  find_program( ZCAT_EXE NAMES zcap ZCAT PATHS /bin /usr/bin )
  add_custom_command( OUTPUT mb_big_test.g
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz > mb_big_test.g
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz )
  add_custom_command( OUTPUT cell1.gen
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz > cell1.gen
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz )
  add_custom_command( OUTPUT cell2.gen
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz > cell2.gen
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz)

  add_custom_target(moab_test_files DEPENDS mb_big_test.g cell1.gen cell2.gen)

  #MBTest.cpp:
  set_source_files_properties( MBTest.cpp
                               COMPILE_FLAGS "-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR} ${TEST_COMP_FLAGS}" )
  add_executable( moab_test MBTest.cpp )
  target_link_libraries( moab_test MOAB  )
  add_dependencies(moab_test moab_test_files)
  add_test( TestMOAB ${EXECUTABLE_OUTPUT_PATH}/moab_test )
endif()

add_executable( var_len_test_no_template VarLenTagTest.cpp )
set_source_files_properties( VarLenTagTest.cpp
  COMPILE_FLAGS "-DMOAB_TEMPLATE_SPECIALIZATION ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( var_len_test_no_template MOAB  )
add_test( var_len_test_no_template ${EXECUTABLE_OUTPUT_PATH}/var_len_test_no_template )

add_executable( mbcn_test mbcn_test.cpp )
set_target_properties( mbcn_test PROPERTIES
  COMPILE_FLAGS "-DMOAB_TEMPLATE_SPECIALIZATION ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( mbcn_test MOAB  )
add_test( mbcn_test ${EXECUTABLE_OUTPUT_PATH}/mbcn_test )

add_executable( file_options_test ${CMAKE_SOURCE_DIR}/src/FileOptions.cpp )
set_source_files_properties( ${CMAKE_SOURCE_DIR}/src/FileOptions.cpp
                             COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( file_options_test MOAB  )
add_test( file_options_test ${EXECUTABLE_OUTPUT_PATH}/file_options_test )

add_executable( homxform_test homxform_test.cpp )
set_source_files_properties( homxform_test.cpp
                             COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( homxform_test MOAB  )
add_test( TestHomXform ${EXECUTABLE_OUTPUT_PATH}/homxform_test )

add_executable( xform_test affinexform_test.cpp )
set_source_files_properties( affinexform_test.cpp
                             COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( xform_test MOAB  )
add_test( xform_test ${EXECUTABLE_OUTPUT_PATH}/xform_test )

if ( MOAB_HAVE_MPI AND MPI_FOUND )

  add_executable ( mbparallelcomm_test parallel/mbparallelcomm_test.cpp )
  target_link_libraries( mbparallelcomm_test MOAB  )
  set_source_files_properties( parallel/mbparallelcomm_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_test( TestParallelComm-BcastDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} 0 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -1 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -2 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-Broadcast
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -3 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )

  if ( MOAB_HAVE_HDF5 )
    include_directories(
      ${HDF5_INCLUDES}
      ${CMAKE_SOURCE_DIR}/src/io/mhdf/include
    )
    set_source_files_properties( parallel/mhdf_parallel.c
      COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
    add_executable( mhdf_parallel parallel/mhdf_parallel.c )
    target_link_libraries( mhdf_parallel MOAB )
    add_test( TestMHDFParallel ${EXECUTABLE_OUTPUT_PATH}/mhdf_parallel )

    set_source_files_properties( parallel/parmerge_test.cpp
      COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
    add_executable( parmerge_test parallel/parmerge_test.cpp )
    target_link_libraries( parmerge_test MOAB )
    add_test( Testparmerge_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/parmerge_test ${MPIEXEC_POSTFLAGS}  )

  endif ( MOAB_HAVE_HDF5 )

  set_source_files_properties( parallel/parallel_unit_tests.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_executable ( parallel_unit_tests parallel/parallel_unit_tests.cpp )
  target_link_libraries( parallel_unit_tests MOAB )
  add_test( TestParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/parallel_unit_tests ${MPIEXEC_POSTFLAGS} )

  set_source_files_properties( parallel/pcomm_serial.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_executable ( pcomm_serial parallel/pcomm_serial.cpp )
  target_link_libraries( pcomm_serial MOAB  )
  add_test( TestPCommSerial
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/pcomm_serial ${MPIEXEC_POSTFLAGS} 1 ${MOAB_SOURCE_DIR}/MeshFiles/unittest/ptest.cub )


  if(MOAB_HAVE_TEMPESTREMAP)
    set( MBTR_TESTS parallel/imoab_coupler.cpp parallel/imoab_coupler_bilin.cpp parallel/imoab_coupler_twohop.cpp
                    parallel/imoab_domain_culling.cpp parallel/imoab_map_r2l.cpp parallel/imoab_map_target.cpp
                    parallel/imoab_ptest2.cpp parallel/imoab_read_compute_map.cpp
                    parallel/imoab_scrip.cpp parallel/imoab_setwithGid.cpp parallel/imoab_phatm_ocn_coupler.cpp
                    parallel/imoab_apg2_ol_coupler.cpp parallel/imoab_spmv_bfb.cpp
                    parallel/imoab_dualmap_caas.cpp )
    
    foreach( fname ${MBTR_TESTS} )
      string( REPLACE ".cpp" "" base ${fname} )
      string( REPLACE "parallel/" "" noparbase ${base} )
      add_executable( ${noparbase} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
      set_target_properties( ${noparbase} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
      target_link_libraries( ${noparbase} MOAB ${TEMPESTREMAP_LIBRARIES} )
      add_test( ${noparbase} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${noparbase} ${MPIEXEC_POSTFLAGS} )
    endforeach()

    # iMOAB BfB / tolerance / exit-code regression harness.
    #
    # Each config under test/parallel/regression/configs/ becomes one CTest
    # test, gated on MOAB_HAVE_TEMPESTREMAP (this whole block) and run when
    # the user invokes `ctest`. Pass/fail comes from the dispatcher's exit
    # code: 0 = harness succeeded, non-zero = a BfB/tolerance/exit-code
    # check failed for that config.
    #
    # Selecting just these from ctest:
    #     ctest -L imoab_bfb
    # Excluding them from a wider ctest run:
    #     ctest -LE imoab_bfb
    # One-shot convenience target (mirrors `make imoabcheck` in autotools):
    #     cmake --build <builddir> --target imoab_bfb_check
    #
    # The dispatcher honors IMOAB_BINDIR to resolve binary locations
    # (CMake puts them in ${EXECUTABLE_OUTPUT_PATH}, autotools in
    # ${MOAB_BUILD}/test/parallel); MBCMPFILES is the h5m comparator the
    # tol-mode configs require. Both are wired through ENVIRONMENT below.
    set(_imoab_bfb_configs
        dualmap_caas spmv_bfb
        coupler coupler_bfb coupler_bilin coupler_twohop
        map_target read_compute_map map_r2l
        phatm_ocn_coupler domain_culling scrip setwithgid ptest2)

    set(_imoab_bfb_dispatcher
        ${CMAKE_SOURCE_DIR}/test/parallel/regression/imoab_test_runner.sh)
    set(_imoab_bfb_env
        "MOAB_BUILD=${CMAKE_BINARY_DIR}"
        "MOAB_MESH_DIR=${CMAKE_SOURCE_DIR}/MeshFiles"
        "MPIRUN=${MPIEXEC_EXECUTABLE}"
        "MBCMPFILES=${EXECUTABLE_OUTPUT_PATH}/mbcmpfiles"
        "IMOAB_BINDIR=${EXECUTABLE_OUTPUT_PATH}")

    foreach(_cfg ${_imoab_bfb_configs})
      add_test(NAME imoab_bfb_${_cfg}
               COMMAND ${_imoab_bfb_dispatcher} ${_cfg})
      set_tests_properties(imoab_bfb_${_cfg} PROPERTIES
        LABELS "imoab_bfb"
        ENVIRONMENT "${_imoab_bfb_env}")
    endforeach()

    add_custom_target(imoab_bfb_check
      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L imoab_bfb
      COMMENT "Running all iMOAB BfB regression configs via ctest")

    # Fortran tests
    set( MBTRF_TESTS parallel/imoab_ptest.F parallel/imoab_coupler_fortran.F90 parallel/imoab_coupler_caas.F90  )
    
    foreach( fname ${MBTRF_TESTS} )
      string( REPLACE ".F90" "" base ${fname} ) # For F90
      string( REPLACE ".F" "" base ${base} ) # For F77
      string( REPLACE "parallel/" "" noparbase ${base} )
      add_executable( ${noparbase} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
      set_property(TARGET ${noparbase} PROPERTY LINKER_LANGUAGE Fortran)
      set_target_properties( ${noparbase} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
      target_link_libraries( ${noparbase} MOAB ${TEMPESTREMAP_LIBRARIES} )
      add_test( ${noparbase} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/${noparbase} ${MPIEXEC_POSTFLAGS} )
    endforeach()


  endif()


endif ( MOAB_HAVE_MPI AND MPI_FOUND )
