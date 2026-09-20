function(_moab_cache_assignment output name)
  if(DEFINED CACHE{${name}})
    get_property(type CACHE "${name}" PROPERTY TYPE)
    get_property(value CACHE "${name}" PROPERTY VALUE)
    if(value MATCHES "]========]")
      message(FATAL_ERROR "Cannot serialize distcheck cache value ${name}")
    endif()
    set(${output} "set(${name} [========[${value}]========] CACHE ${type} \"distcheck parent value\" FORCE)\n" PARENT_SCOPE)
  else()
    set(${output} "" PARENT_SCOPE)
  endif()
endfunction()

function(DISTCHECK_SETUP)
  if(NOT UNIX)
    message(WARNING "The distcheck target currently supports only UNIX-like platforms.")
    return()
  endif()

  set(DISTCHECK_PARALLEL_LEVEL "" CACHE STRING "Parallel build level used by distcheck")
  if(ENABLE_PYMOAB)
    message(FATAL_ERROR "CMake distcheck does not support deprecated ENABLE_PYMOAB installs; install PyMOAB separately")
  endif()
  set(_distcheck_cache_names
    BUILD_DOCUMENTATION BUILD_SHARED_LIBS CMAKE_BUILD_TYPE CMAKE_C_COMPILER
    CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_RELWITHDEBINFO CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_Fortran_COMPILER CMAKE_Fortran_FLAGS
    CMAKE_Fortran_FLAGS_DEBUG CMAKE_Fortran_FLAGS_MINSIZEREL CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS
    CMAKE_MODULE_LINKER_FLAGS CMAKE_STATIC_LINKER_FLAGS CMAKE_INSTALL_BINDIR
    CMAKE_INSTALL_DATADIR CMAKE_INSTALL_DATAROOTDIR CMAKE_INSTALL_DOCDIR
    CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_LIBEXECDIR
    CMAKE_INSTALL_SBINDIR CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
    CMAKE_OSX_SYSROOT CMAKE_PREFIX_PATH CMAKE_TOOLCHAIN_FILE CMAKE_SYSROOT
    CMAKE_FIND_ROOT_PATH ENABLE_BLASLAPACK ENABLE_CGNS ENABLE_CPM
    ENABLE_EIGEN ENABLE_FORTRAN ENABLE_HDF5 ENABLE_METIS ENABLE_MPI ENABLE_NETCDF
    ENABLE_PARMETIS ENABLE_PNETCDF ENABLE_PYMOAB ENABLE_TEMPESTREMAP ENABLE_TESTING
    ENABLE_ZOLTAN EIGEN3_INCLUDE_DIR HDF5_C_COMPILER_EXECUTABLE HDF5_DIR HDF5_ROOT
    METIS_DIR MPIEXEC_EXECUTABLE MPIEXEC_MAX_NUMPROCS MPIEXEC_NUMPROC_FLAG
    MPIEXEC_POSTFLAGS MPIEXEC_PREFLAGS MPI_C_COMPILER MPI_CXX_COMPILER MPI_Fortran_COMPILER
    NC_CONFIG_EXECUTABLE NETCDF_DIR PARMETIS_DIR PNETCDF_DIR TEMPESTREMAP_DIR ZOLTAN_DIR
    MOAB_BUILD_MBTEMPEST_TESTS MOAB_FORCE_32_BIT_HANDLES MOAB_FORCE_64_BIT_HANDLES)
  set(DISTCHECK_INITIAL_CACHE "")
  foreach(name IN LISTS _distcheck_cache_names)
    _moab_cache_assignment(assignment "${name}")
    string(APPEND DISTCHECK_INITIAL_CACHE "${assignment}")
  endforeach()
  string(APPEND DISTCHECK_INITIAL_CACHE
    "set(ENABLE_TESTING ON CACHE BOOL \"Tests are required by distcheck\" FORCE)\n")

  set(DISTCHECK_SIGN_ARG "")
  if(MOAB_DIST_SIGN)
    set(DISTCHECK_SIGN_ARG "-DGPG_EXECUTABLE=${GPG_EXECUTABLE}")
  endif()

  set(DISTCHECK_SOURCE_DIR "${CMAKE_BINARY_DIR}/_distcheck/source")
  set(DISTCHECK_BUILD_DIR "${CMAKE_BINARY_DIR}/_distcheck/build")
  set(DISTCHECK_INSTALL_DIR "${CMAKE_BINARY_DIR}/_distcheck/install")
  set(DISTCHECK_EXAMPLES_BUILD_DIR "${CMAKE_BINARY_DIR}/_distcheck/examples-build")
  configure_file(
    "${CMAKE_SOURCE_DIR}/config/distcheck.cmake"
    "${CMAKE_BINARY_DIR}/distcheck-driver.cmake"
    @ONLY)
  file(WRITE "${CMAKE_BINARY_DIR}/distcheck-cache.cmake" "${DISTCHECK_INITIAL_CACHE}")

  add_custom_target(distcheck
    COMMAND "${CMAKE_COMMAND}" -DMODE=distcheck -P "${CMAKE_BINARY_DIR}/distcheck-driver.cmake"
    COMMENT "Checking the generated source tarball..."
    USES_TERMINAL
    VERBATIM)
endfunction()

if(MODE STREQUAL "distcheck")
  cmake_minimum_required(VERSION 3.20)

  set(source_root "@DISTCHECK_SOURCE_DIR@")
  set(build_dir "@DISTCHECK_BUILD_DIR@")
  set(install_dir "@DISTCHECK_INSTALL_DIR@")
  set(examples_build_dir "@DISTCHECK_EXAMPLES_BUILD_DIR@")
  set(dist_name "@MOAB_DIST_NAME@")
  set(dist_dir "@MOAB_DIST_DIR@")
  set(artifact_file "@MOAB_DIST_ARTIFACT_FILE@")

  function(run_step description)
    execute_process(COMMAND ${ARGN} RESULT_VARIABLE result COMMAND_ECHO STDOUT)
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "${description} failed with exit status ${result}")
    endif()
  endfunction()

  function(assert_no_files directory description)
    file(GLOB_RECURSE remaining LIST_DIRECTORIES FALSE "${directory}/*")
    if(remaining)
      list(JOIN remaining "\n  " formatted)
      message(FATAL_ERROR "${description} left files behind:\n  ${formatted}")
    endif()
  endfunction()

  file(REMOVE_RECURSE "@CMAKE_BINARY_DIR@/_distcheck" "${dist_dir}")
  file(REMOVE "${artifact_file}")
  file(MAKE_DIRECTORY "${source_root}")

  run_step("distribution directory generation"
    "@CMAKE_COMMAND@" --build "@CMAKE_BINARY_DIR@" --target distdir)
  run_step("source tarball generation"
    "@CMAKE_COMMAND@"
    -DMODE=archive
    -DARCHIVE_FORMAT=gz
    -DSOURCE_DIR=@CMAKE_SOURCE_DIR@
    -DBINARY_DIR=@CMAKE_BINARY_DIR@
    -DPROJECT_VERSION=@CMAKE_PROJECT_VERSION@
    -DGIT_EXECUTABLE=@GIT_EXECUTABLE@
    -DPERL_EXECUTABLE=@PERL_EXECUTABLE@
    @DISTCHECK_SIGN_ARG@
    -P "@CMAKE_SOURCE_DIR@/config/dist.cmake")
  if(NOT EXISTS "${artifact_file}")
    message(FATAL_ERROR "Distribution target did not record the generated tarball path")
  endif()
  file(READ "${artifact_file}" tarball)
  string(STRIP "${tarball}" tarball)
  if(NOT EXISTS "${tarball}")
    message(FATAL_ERROR "Distribution tarball does not exist: ${tarball}")
  endif()
  file(REMOVE_RECURSE "${dist_dir}")
  run_step("source tarball extraction"
    "@CMAKE_COMMAND@" -E chdir "${source_root}" "@CMAKE_COMMAND@" -E tar xzf "${tarball}")

  file(RENAME "${source_root}/${dist_name}" "${source_root}/tree")
  set(source_dir "${source_root}/tree")
  file(GLOB_RECURSE source_files RELATIVE "${source_dir}" LIST_DIRECTORIES FALSE "${source_dir}/*")
  foreach(path IN LISTS source_files)
    file(SHA256 "${source_dir}/${path}" digest)
    string(APPEND source_manifest "${digest}  ${path}\n")
  endforeach()

  set(configure_command
    "@CMAKE_COMMAND@" -S "${source_dir}" -B "${build_dir}"
    -G "@CMAKE_GENERATOR@"
    -C "@CMAKE_BINARY_DIR@/distcheck-cache.cmake"
    -DCMAKE_INSTALL_PREFIX:PATH=${install_dir})
  if(NOT "@CMAKE_GENERATOR_PLATFORM@" STREQUAL "")
    list(APPEND configure_command -A "@CMAKE_GENERATOR_PLATFORM@")
  endif()
  if(NOT "@CMAKE_GENERATOR_TOOLSET@" STREQUAL "")
    list(APPEND configure_command -T "@CMAKE_GENERATOR_TOOLSET@")
  endif()
  run_step("tarball configuration" ${configure_command})

  set(build_args "@CMAKE_COMMAND@" --build "${build_dir}")
  if(NOT "@CMAKE_BUILD_TYPE@" STREQUAL "")
    list(APPEND build_args --config "@CMAKE_BUILD_TYPE@")
  endif()
  if(NOT "@DISTCHECK_PARALLEL_LEVEL@" STREQUAL "")
    list(APPEND build_args --parallel "@DISTCHECK_PARALLEL_LEVEL@")
  endif()
  run_step("tarball build" ${build_args})

  set(test_args "@CMAKE_CTEST_COMMAND@" --test-dir "${build_dir}" --output-on-failure)
  if(NOT "@CMAKE_BUILD_TYPE@" STREQUAL "")
    list(APPEND test_args -C "@CMAKE_BUILD_TYPE@")
  endif()
  run_step("tarball test suite" ${test_args})

  set(install_args "@CMAKE_COMMAND@" --install "${build_dir}" --prefix "${install_dir}")
  if(NOT "@CMAKE_BUILD_TYPE@" STREQUAL "")
    list(APPEND install_args --config "@CMAKE_BUILD_TYPE@")
  endif()
  run_step("tarball installation" ${install_args})

  set(examples_configure_command
    "@CMAKE_COMMAND@"
    -S "${install_dir}/@MOABDocLocation@/examples"
    -B "${examples_build_dir}"
    -G "@CMAKE_GENERATOR@"
    -DMOAB_ROOT:PATH=${install_dir})
  if(NOT "@CMAKE_BUILD_TYPE@" STREQUAL "")
    list(APPEND examples_configure_command -DCMAKE_BUILD_TYPE:STRING=@CMAKE_BUILD_TYPE@)
  endif()
  if(NOT "@CMAKE_GENERATOR_PLATFORM@" STREQUAL "")
    list(APPEND examples_configure_command -A "@CMAKE_GENERATOR_PLATFORM@")
  endif()
  if(NOT "@CMAKE_GENERATOR_TOOLSET@" STREQUAL "")
    list(APPEND examples_configure_command -T "@CMAKE_GENERATOR_TOOLSET@")
  endif()
  run_step("installed examples configuration" ${examples_configure_command})

  set(examples_build_args "@CMAKE_COMMAND@" --build "${examples_build_dir}")
  if(NOT "@CMAKE_BUILD_TYPE@" STREQUAL "")
    list(APPEND examples_build_args --config "@CMAKE_BUILD_TYPE@")
  endif()
  if(NOT "@DISTCHECK_PARALLEL_LEVEL@" STREQUAL "")
    list(APPEND examples_build_args --parallel "@DISTCHECK_PARALLEL_LEVEL@")
  endif()
  run_step("installed examples build" ${examples_build_args})

  run_step("tarball uninstall"
    "@CMAKE_COMMAND@" --build "${build_dir}" --target uninstall)
  assert_no_files("${install_dir}" "Uninstall")
  run_step("tarball clean"
    "@CMAKE_COMMAND@" --build "${build_dir}" --target clean)

  set(source_manifest_after "")
  foreach(path IN LISTS source_files)
    if(NOT EXISTS "${source_dir}/${path}")
      message(FATAL_ERROR "Build removed source file: ${path}")
    endif()
    file(SHA256 "${source_dir}/${path}" digest)
    string(APPEND source_manifest_after "${digest}  ${path}\n")
  endforeach()
  file(GLOB_RECURSE source_files_after RELATIVE "${source_dir}" LIST_DIRECTORIES FALSE "${source_dir}/*")
  if(NOT "${source_files}" STREQUAL "${source_files_after}" OR
     NOT "${source_manifest}" STREQUAL "${source_manifest_after}")
    message(FATAL_ERROR "The build modified the extracted source tree")
  endif()

  file(REMOVE_RECURSE "@CMAKE_BINARY_DIR@/_distcheck")
  message(STATUS "${dist_name} is ready for distribution")
endif()
