# ==============================================================================
# mbtempest Regression Test Suite
# ==============================================================================
#
# Comprehensive integration and regression tests for the mbtempest CLI tool.
# Tests all discretization methods, option combinations, and parallel workflows.
#
# Usage:
#   cmake -DMOAB_BUILD_MBTEMPEST_TESTS=ON ..
#   make -j$(nproc)
#   ctest -L mbtempest                    # Run all mbtempest regression tests
#   ctest -L mbtempest_gen                # Mesh generation tests only
#   ctest -L mbtempest_remap              # Remapping tests only
#   ctest -L mbtempest_parallel           # Parallel tests only
#   ctest -L mbtempest_convergence        # Convergence analysis tests only
#   ctest -L mbtempest_options            # Advanced option tests only
#
# Baseline generation:
#   cmake -DMOAB_MBTEMPEST_GENERATE_BASELINES=ON ..
#   ctest -L mbtempest_gen -L mbtempest_remap
#   # Upload generated baselines from ${MBTEMPEST_BASELINE_DIR}
#
# ==============================================================================

# Guard: only build if TempestRemap is available
if(NOT MOAB_HAVE_TEMPESTREMAP)
  return()
endif()

option(MOAB_MBTEMPEST_GENERATE_BASELINES
  "Generate golden reference files instead of comparing against them" OFF)

# ==============================================================================
# Path setup
# ==============================================================================

set(MBTEMPEST_EXE ${EXECUTABLE_OUTPUT_PATH}/mbtempest)
set(MBCMPMAPS_EXE ${EXECUTABLE_OUTPUT_PATH}/mbcmpmaps)
set(MBCMPFILES_EXE ${EXECUTABLE_OUTPUT_PATH}/mbcmpfiles)
set(MBPART_EXE ${EXECUTABLE_OUTPUT_PATH}/mbpart)

# Working directory for test outputs (per-build, not committed)
set(MBTEMPEST_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/work)
file(MAKE_DIRECTORY ${MBTEMPEST_TEST_DIR})

# Baseline directory (golden reference files)
set(MBTEMPEST_BASELINE_DIR ${CMAKE_CURRENT_BINARY_DIR}/baselines)
file(MAKE_DIRECTORY ${MBTEMPEST_BASELINE_DIR})

# Mesh file directory (from unittest data)
set(MESH_DIR ${MOAB_ABSSRC_DIR}/MeshFiles/unittest)

# Pre-existing mesh files in MOAB h5m format (examples/advanced/)
set(EXAMPLE_MESH_DIR ${MOAB_ABSSRC_DIR}/examples/advanced)

# External data directory (already set up by top-level CMakeLists.txt)
if(DEFINED EXT_DATA_DIR)
  set(MBTEMPEST_EXT_MESH_DIR ${EXT_DATA_DIR})
endif()

# Symlink pre-existing mesh files into working directory for test convenience
foreach(meshfile ne30.h5m ico30.h5m icod30.h5m rll90.h5m)
  if(EXISTS ${EXAMPLE_MESH_DIR}/${meshfile})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${EXAMPLE_MESH_DIR}/${meshfile}
      ${MBTEMPEST_TEST_DIR}/${meshfile})
  endif()
endforeach()

# ==============================================================================
# Download golden reference files when not generating baselines
# ==============================================================================

# Placeholder URL — user will upload baselines and set this
set(MBTEMPEST_BASELINE_URL "" CACHE STRING
  "URL prefix for downloading mbtempest regression baselines")

# When baselines are available remotely and not generating, download them
# For now, baselines are generated locally via MOAB_MBTEMPEST_GENERATE_BASELINES=ON
# and then uploaded. The download logic will be enabled once a URL is configured.
#
# if(NOT MOAB_MBTEMPEST_GENERATE_BASELINES AND MBTEMPEST_BASELINE_URL)
#   foreach(baseline_file ${MBTEMPEST_BASELINE_FILES})
#     if(NOT EXISTS "${MBTEMPEST_BASELINE_DIR}/${baseline_file}")
#       file(DOWNLOAD
#         "${MBTEMPEST_BASELINE_URL}/${baseline_file}"
#         "${MBTEMPEST_BASELINE_DIR}/${baseline_file}"
#         STATUS download_status)
#     endif()
#   endforeach()
# endif()

# ==============================================================================
# Helper macros
# ==============================================================================

# Add a serial mbtempest test (runs through mpirun -np 1 when MPI is available,
# since mbtempest is MPI-linked and h5m files require parallel HDF5 reader)
# Usage: add_mbtempest_test(NAME test_name ARGS ... LABELS ... [DEPENDS ...])
macro(add_mbtempest_test)
  cmake_parse_arguments(MT "" "NAME" "ARGS;LABELS;DEPENDS" ${ARGN})
  if(MOAB_HAVE_MPI AND MPI_FOUND)
    add_test(NAME ${MT_NAME}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
              ${MPIEXEC_PREFLAGS} ${MBTEMPEST_EXE} ${MT_ARGS} ${MPIEXEC_POSTFLAGS}
      WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
    )
  else()
    add_test(NAME ${MT_NAME}
      COMMAND ${MBTEMPEST_EXE} ${MT_ARGS}
      WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
    )
  endif()
  set_tests_properties(${MT_NAME} PROPERTIES
    LABELS "mbtempest;${MT_LABELS}"
    TIMEOUT 300
  )
  if(MT_DEPENDS)
    set_tests_properties(${MT_NAME} PROPERTIES
      DEPENDS "${MT_DEPENDS}"
    )
  endif()
endmacro()

# Add a parallel mbtempest test (MPI)
# Usage: add_mbtempest_mpi_test(NAME test_name NPROCS n ARGS ... LABELS ... [DEPENDS ...])
macro(add_mbtempest_mpi_test)
  cmake_parse_arguments(MT "" "NAME;NPROCS" "ARGS;LABELS;DEPENDS" ${ARGN})
  if(MOAB_HAVE_MPI AND MPI_FOUND)
    add_test(NAME ${MT_NAME}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MT_NPROCS}
              ${MPIEXEC_PREFLAGS} ${MBTEMPEST_EXE} ${MT_ARGS} ${MPIEXEC_POSTFLAGS}
      WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
    )
    set_tests_properties(${MT_NAME} PROPERTIES
      LABELS "mbtempest;mbtempest_parallel;${MT_LABELS}"
      TIMEOUT 600
      PROCESSORS ${MT_NPROCS}
    )
    if(MT_DEPENDS)
      set_tests_properties(${MT_NAME} PROPERTIES
        DEPENDS "${MT_DEPENDS}"
      )
    endif()
  endif()
endmacro()

# Add a map comparison test using mbcmpmaps
# Usage: add_map_compare_test(NAME test_name MAP1 file1 MAP2 file2 TOL tolerance [DEPENDS ...])
macro(add_map_compare_test)
  cmake_parse_arguments(MC "" "NAME;MAP1;MAP2;TOL" "LABELS;DEPENDS" ${ARGN})
  add_test(NAME ${MC_NAME}
    COMMAND ${MBCMPMAPS_EXE} -i ${MC_MAP1} -j ${MC_MAP2} -t ${MC_TOL}
    WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
  )
  set_tests_properties(${MC_NAME} PROPERTIES
    LABELS "mbtempest;mbtempest_verify;${MC_LABELS}"
    TIMEOUT 120
  )
  if(MC_DEPENDS)
    set_tests_properties(${MC_NAME} PROPERTIES
      DEPENDS "${MC_DEPENDS}"
    )
  endif()
endmacro()

# Add a mesh/tag comparison test using mbcmpfiles
# Usage: add_file_compare_test(NAME test_name FILE1 f1 FILE2 f2 [TAG tag] [DIM d] [DEPENDS ...])
macro(add_file_compare_test)
  cmake_parse_arguments(FC "" "NAME;FILE1;FILE2;TAG;DIM" "LABELS;DEPENDS" ${ARGN})
  set(_fc_args -i ${FC_FILE1} -j ${FC_FILE2})
  if(FC_TAG)
    list(APPEND _fc_args -n ${FC_TAG})
  endif()
  if(FC_DIM)
    list(APPEND _fc_args -d ${FC_DIM})
  endif()
  add_test(NAME ${FC_NAME}
    COMMAND ${MBCMPFILES_EXE} ${_fc_args}
    WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
  )
  set_tests_properties(${FC_NAME} PROPERTIES
    LABELS "mbtempest;mbtempest_verify;${FC_LABELS}"
    TIMEOUT 120
  )
  if(FC_DEPENDS)
    set_tests_properties(${FC_NAME} PROPERTIES
      DEPENDS "${FC_DEPENDS}"
    )
  endif()
endmacro()

# Conditionally copy output as baseline or compare against baseline
# Usage: add_baseline_test(NAME test_name OUTPUT output_file BASELINE baseline_name TOL tol DEPENDS dep)
macro(add_baseline_test)
  cmake_parse_arguments(BL "" "NAME;OUTPUT;BASELINE;TOL" "LABELS;DEPENDS" ${ARGN})
  if(MOAB_MBTEMPEST_GENERATE_BASELINES)
    # In generate mode: copy the output as the new baseline
    add_test(NAME ${BL_NAME}
      COMMAND ${CMAKE_COMMAND} -E copy
        ${MBTEMPEST_TEST_DIR}/${BL_OUTPUT}
        ${MBTEMPEST_BASELINE_DIR}/${BL_BASELINE}
      WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
    )
    set_tests_properties(${BL_NAME} PROPERTIES
      LABELS "mbtempest;mbtempest_baseline;${BL_LABELS}"
      TIMEOUT 30
    )
  else()
    # In verify mode: compare output against baseline
    add_map_compare_test(NAME ${BL_NAME}
      MAP1 ${MBTEMPEST_BASELINE_DIR}/${BL_BASELINE}
      MAP2 ${MBTEMPEST_TEST_DIR}/${BL_OUTPUT}
      TOL ${BL_TOL}
      LABELS "${BL_LABELS}"
    )
  endif()
  if(BL_DEPENDS)
    set_tests_properties(${BL_NAME} PROPERTIES
      DEPENDS "${BL_DEPENDS}"
    )
  endif()
endmacro()

# Add a mesh partitioning test using mbpart with --inferred for co-partitioning
# Usage: add_partition_test(NAME test_name NPROCS n TARGET tgt.h5m SOURCE src.h5m
#                           OUTPUT_TARGET tgt_pN.h5m [DEPENDS ...])
macro(add_partition_test)
  cmake_parse_arguments(PT "" "NAME;NPROCS;TARGET;SOURCE;OUTPUT_TARGET" "LABELS;DEPENDS" ${ARGN})
  if(MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_HAVE_ZOLTAN)
    add_test(NAME ${PT_NAME}
      COMMAND ${MBPART_EXE} ${PT_NPROCS} -z RCB
              ${PT_TARGET} ${PT_OUTPUT_TARGET}
              --inferred ${PT_SOURCE}
              --project_on_sphere 2 --recompute_rcb_box --scale_sphere
      WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
    )
    set_tests_properties(${PT_NAME} PROPERTIES
      LABELS "mbtempest;mbtempest_partition;${PT_LABELS}"
      TIMEOUT 120
    )
    if(PT_DEPENDS)
      set_tests_properties(${PT_NAME} PROPERTIES
        DEPENDS "${PT_DEPENDS}"
      )
    endif()
  endif()
endmacro()

# ==============================================================================
# Category 1: Mesh Generation Tests
# ==============================================================================

# Cubed-Sphere meshes
add_mbtempest_test(NAME mbtempest_gen_cs5
  ARGS -t 0 -r 5 -f cs5.h5m
  LABELS "mbtempest_gen")

add_mbtempest_test(NAME mbtempest_gen_cs10
  ARGS -t 0 -r 10 -f cs10.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_cs5)

add_mbtempest_test(NAME mbtempest_gen_cs25
  ARGS -t 0 -r 25 -f cs25.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_cs10)

# RLL (regular lat-lon) meshes
add_mbtempest_test(NAME mbtempest_gen_rll5
  ARGS -t 1 -r 5 -f rll5.h5m
  LABELS "mbtempest_gen")

add_mbtempest_test(NAME mbtempest_gen_rll10
  ARGS -t 1 -r 10 -f rll10.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_rll5)

add_mbtempest_test(NAME mbtempest_gen_rll25
  ARGS -t 1 -r 25 -f rll25.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_rll10)

# Icosahedral meshes
add_mbtempest_test(NAME mbtempest_gen_ico5
  ARGS -t 2 -r 5 -f ico5.h5m
  LABELS "mbtempest_gen")

add_mbtempest_test(NAME mbtempest_gen_ico5_dual
  ARGS -t 2 -r 5 -d -f icod5.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_ico5)

add_mbtempest_test(NAME mbtempest_gen_ico10
  ARGS -t 2 -r 10 -f ico10.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_ico5_dual)

add_mbtempest_test(NAME mbtempest_gen_ico10_dual
  ARGS -t 2 -r 10 -d -f icod10.h5m
  LABELS "mbtempest_gen"
  DEPENDS mbtempest_gen_ico10)

# ==============================================================================
# SE Mesh Preparation: generate CS meshes with GLL DOF tags for SE tests
# ==============================================================================
# Step 1: Generate TempestRemap .g meshes (these are Exodus format)
# Step 2: Convert to MOAB h5m with GLOBAL_DOFS tag via mbconvert -B
#
# Note: mbconvert -B with -r N only works for CS meshes (GLL DOF numbering).
# DGLL tests are currently unsupported because DGLL needs per-element DOF
# numbering which mbconvert does not yet generate.
# ==============================================================================

set(MBCONVERT_EXE ${EXECUTABLE_OUTPUT_PATH}/mbconvert)

# Generate CS30 .g mesh for SE source/target
add_mbtempest_test(NAME mbtempest_se_gen_cs30
  ARGS -t 0 -r 30 -f ne30_se.g
  LABELS "mbtempest_gen;mbtempest_se_prep")

# Generate CS25 .g mesh for SE-SE target
add_mbtempest_test(NAME mbtempest_se_gen_cs25
  ARGS -t 0 -r 25 -f ne25_se.g
  LABELS "mbtempest_gen;mbtempest_se_prep"
  DEPENDS mbtempest_se_gen_cs30)

# Convert CS30 with order 4 DOFs
add_test(NAME mbtempest_se_convert_cs30_o4
  COMMAND ${MBCONVERT_EXE} -B ne30_se.g ne30_se_o4.h5m -i GLOBAL_DOFS -r 4
  WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR})
set_tests_properties(mbtempest_se_convert_cs30_o4 PROPERTIES
  LABELS "mbtempest;mbtempest_se_prep"
  TIMEOUT 60
  DEPENDS mbtempest_se_gen_cs30)

# Convert CS30 with order 2 DOFs
add_test(NAME mbtempest_se_convert_cs30_o2
  COMMAND ${MBCONVERT_EXE} -B ne30_se.g ne30_se_o2.h5m -i GLOBAL_DOFS -r 2
  WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR})
set_tests_properties(mbtempest_se_convert_cs30_o2 PROPERTIES
  LABELS "mbtempest;mbtempest_se_prep"
  TIMEOUT 60
  DEPENDS mbtempest_se_gen_cs30)

# Convert CS25 with order 4 DOFs (for SE-SE tests)
add_test(NAME mbtempest_se_convert_cs25_o4
  COMMAND ${MBCONVERT_EXE} -B ne25_se.g ne25_se_o4.h5m -i GLOBAL_DOFS -r 4
  WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR})
set_tests_properties(mbtempest_se_convert_cs25_o4 PROPERTIES
  LABELS "mbtempest;mbtempest_se_prep"
  TIMEOUT 60
  DEPENDS mbtempest_se_gen_cs25)

# ==============================================================================
# Category 2: FV-FV Remapping (serial + parallel)
#   Analytical function: SH (slow harmonic) for order 1
# ==============================================================================

# --- CS → RLL ---
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_o1
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_rll90_fv_fv_o1.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_o1
  OUTPUT map_ne30_rll90_fv_fv_o1.nc
  BASELINE gold_ne30_rll90_fv_fv_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_o1)

# Partition meshes for 2-rank parallel tests
add_partition_test(NAME mbtempest_partition_ne30_rll90_np2
  NPROCS 2
  TARGET rll90.h5m SOURCE ne30.h5m
  OUTPUT_TARGET rll90_p2.h5m
  LABELS "mbtempest_remap")

# Parallel: 2 ranks (uses partitioned meshes: ne30_inferred.h5m and rll90_p2.h5m)
add_mbtempest_mpi_test(NAME mbtempest_remap_ne30_rll90_fv_fv_o1_np2
  NPROCS 2
  ARGS -t 5 -l ne30_inferred.h5m -l rll90_p2.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_rll90_fv_fv_o1_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ne30_rll90_np2")

# Serial vs parallel consistency
add_map_compare_test(NAME mbtempest_compare_ne30_rll90_fv_fv_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ne30_rll90_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ne30_rll90_fv_fv_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ne30_rll90_fv_fv_o1;mbtempest_remap_ne30_rll90_fv_fv_o1_np2")

# Partition meshes for 4-rank parallel tests
add_partition_test(NAME mbtempest_partition_ne30_rll90_np4
  NPROCS 4
  TARGET rll90.h5m SOURCE ne30.h5m
  OUTPUT_TARGET rll90_p4.h5m
  LABELS "mbtempest_remap")

# Parallel: 4 ranks
add_mbtempest_mpi_test(NAME mbtempest_remap_ne30_rll90_fv_fv_o1_np4
  NPROCS 4
  ARGS -t 5 -l ne30_inferred.h5m -l rll90_p4.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_rll90_fv_fv_o1_np4.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ne30_rll90_np4")

add_map_compare_test(NAME mbtempest_compare_ne30_rll90_fv_fv_o1_serial_vs_np4
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ne30_rll90_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ne30_rll90_fv_fv_o1_np4.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ne30_rll90_fv_fv_o1;mbtempest_remap_ne30_rll90_fv_fv_o1_np4")

# --- CS → ICO ---
add_mbtempest_test(NAME mbtempest_remap_ne30_ico30_fv_fv_o1
  ARGS -t 5 -l ne30.h5m -l ico30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_ico30_fv_fv_o1.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_ico30_fv_fv_o1
  OUTPUT map_ne30_ico30_fv_fv_o1.nc
  BASELINE gold_ne30_ico30_fv_fv_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ne30_ico30_fv_fv_o1)

add_partition_test(NAME mbtempest_partition_ne30_ico30_np2
  NPROCS 2
  TARGET ico30.h5m SOURCE ne30.h5m
  OUTPUT_TARGET ico30_p2.h5m
  LABELS "mbtempest_remap")

add_mbtempest_mpi_test(NAME mbtempest_remap_ne30_ico30_fv_fv_o1_np2
  NPROCS 2
  ARGS -t 5 -l ne30_inferred.h5m -l ico30_p2.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_ico30_fv_fv_o1_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ne30_ico30_np2")

add_map_compare_test(NAME mbtempest_compare_ne30_ico30_fv_fv_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ne30_ico30_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ne30_ico30_fv_fv_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ne30_ico30_fv_fv_o1;mbtempest_remap_ne30_ico30_fv_fv_o1_np2")

# --- CS → ICOD (dual) ---
add_mbtempest_test(NAME mbtempest_remap_ne30_icod30_fv_fv_o1
  ARGS -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_icod30_fv_fv_o1.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_icod30_fv_fv_o1
  OUTPUT map_ne30_icod30_fv_fv_o1.nc
  BASELINE gold_ne30_icod30_fv_fv_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ne30_icod30_fv_fv_o1)

add_partition_test(NAME mbtempest_partition_ne30_icod30_np2
  NPROCS 2
  TARGET icod30.h5m SOURCE ne30.h5m
  OUTPUT_TARGET icod30_p2.h5m
  LABELS "mbtempest_remap")

add_mbtempest_mpi_test(NAME mbtempest_remap_ne30_icod30_fv_fv_o1_np2
  NPROCS 2
  ARGS -t 5 -l ne30_inferred.h5m -l icod30_p2.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ne30_icod30_fv_fv_o1_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ne30_icod30_np2")

add_map_compare_test(NAME mbtempest_compare_ne30_icod30_fv_fv_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ne30_icod30_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ne30_icod30_fv_fv_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ne30_icod30_fv_fv_o1;mbtempest_remap_ne30_icod30_fv_fv_o1_np2")

# --- RLL → CS ---
add_mbtempest_test(NAME mbtempest_remap_rll90_ne30_fv_fv_o1
  ARGS -t 5 -l rll90.h5m -l ne30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_rll90_ne30_fv_fv_o1.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_rll90_ne30_fv_fv_o1
  OUTPUT map_rll90_ne30_fv_fv_o1.nc
  BASELINE gold_rll90_ne30_fv_fv_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_rll90_ne30_fv_fv_o1)

add_partition_test(NAME mbtempest_partition_rll90_ne30_np2
  NPROCS 2
  TARGET ne30.h5m SOURCE rll90.h5m
  OUTPUT_TARGET ne30_p2.h5m
  LABELS "mbtempest_remap")

add_mbtempest_mpi_test(NAME mbtempest_remap_rll90_ne30_fv_fv_o1_np2
  NPROCS 2
  ARGS -t 5 -l rll90_inferred.h5m -l ne30_p2.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_rll90_ne30_fv_fv_o1_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_rll90_ne30_np2")

add_map_compare_test(NAME mbtempest_compare_rll90_ne30_fv_fv_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_rll90_ne30_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_rll90_ne30_fv_fv_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_rll90_ne30_fv_fv_o1;mbtempest_remap_rll90_ne30_fv_fv_o1_np2")

# --- ICO → CS ---
add_mbtempest_test(NAME mbtempest_remap_ico30_ne30_fv_fv_o1
  ARGS -t 5 -l ico30.h5m -l ne30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ico30_ne30_fv_fv_o1.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ico30_ne30_fv_fv_o1
  OUTPUT map_ico30_ne30_fv_fv_o1.nc
  BASELINE gold_ico30_ne30_fv_fv_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ico30_ne30_fv_fv_o1)

add_partition_test(NAME mbtempest_partition_ico30_ne30_np2
  NPROCS 2
  TARGET ne30.h5m SOURCE ico30.h5m
  OUTPUT_TARGET ne30_p2_ico.h5m
  LABELS "mbtempest_remap")

add_mbtempest_mpi_test(NAME mbtempest_remap_ico30_ne30_fv_fv_o1_np2
  NPROCS 2
  ARGS -t 5 -l ico30_inferred.h5m -l ne30_p2_ico.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH
       -f map_ico30_ne30_fv_fv_o1_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ico30_ne30_np2")

add_map_compare_test(NAME mbtempest_compare_ico30_ne30_fv_fv_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ico30_ne30_fv_fv_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ico30_ne30_fv_fv_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ico30_ne30_fv_fv_o1;mbtempest_remap_ico30_ne30_fv_fv_o1_np2")

# --- Higher order FV-FV (CS10 → RLL10) ---
add_mbtempest_test(NAME mbtempest_remap_ne30_icod30_fv_fv_o2
  ARGS -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 2 -o 2 --verify --var SH
       -f map_ne30_icod30_fv_fv_o2.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_icod30_fv_fv_o2
  OUTPUT map_ne30_icod30_fv_fv_o2.nc
  BASELINE gold_ne30_icod30_fv_fv_o2.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ne30_icod30_fv_fv_o2)

add_mbtempest_mpi_test(NAME mbtempest_remap_ne30_icod30_fv_fv_o2_np2
  NPROCS 2
  ARGS -t 5 -l ne30_inferred.h5m -l icod30_p2.h5m -w -m fv -m fv -o 2 -o 2 --verify --var SH
       -f map_ne30_icod30_fv_fv_o2_np2.nc
  LABELS "mbtempest_remap"
  DEPENDS "mbtempest_partition_ne30_icod30_np2")

add_map_compare_test(NAME mbtempest_compare_ne30_icod30_fv_fv_o2_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_ne30_icod30_fv_fv_o2.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_ne30_icod30_fv_fv_o2_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel"
  DEPENDS "mbtempest_remap_ne30_icod30_fv_fv_o2;mbtempest_remap_ne30_icod30_fv_fv_o2_np2")

add_mbtempest_test(NAME mbtempest_remap_ne30_icod30_fv_fv_o3
  ARGS -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 3 -o 3 --verify --var SH
       -f map_ne30_icod30_fv_fv_o3.nc
  LABELS "mbtempest_remap"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_icod30_fv_fv_o3
  OUTPUT map_ne30_icod30_fv_fv_o3.nc
  BASELINE gold_ne30_icod30_fv_fv_o3.nc
  TOL 1e-12
  LABELS "mbtempest_remap"
  DEPENDS mbtempest_remap_ne30_icod30_fv_fv_o3)

# ==============================================================================
# Category 3: CGLL-FV Remapping (serial only for now)
#   Uses SE-tagged CS meshes (ne30_se_o*.h5m) with GLOBAL_DOFS tag
#   High-order uses SV (stationary vortex) for analytical test
#
#   Note: DGLL tests are currently disabled. DGLL requires per-element DOF
#   numbering which mbconvert -B does not yet generate. Only CGLL is supported.
# ==============================================================================

# CGLL → FV, order 2-1
add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o2_o1
  ARGS -t 5 -l ne30_se_o2.h5m -l rll90.h5m -w -m cgll -m fv -o 2 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV
       -f map_cgll_fv_o2_o1.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o2")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o2_o1
  OUTPUT map_cgll_fv_o2_o1.nc
  BASELINE gold_cgll_fv_o2_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS mbtempest_remap_cgll_fv_o2_o1)

# CGLL → FV, order 4-1 (high-order, vortex field)
add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o4_o1
  ARGS -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV
       -f map_cgll_fv_o4_o1.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o4_o1
  OUTPUT map_cgll_fv_o4_o1.nc
  BASELINE gold_cgll_fv_o4_o1.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS mbtempest_remap_cgll_fv_o4_o1)

# CGLL → FV, order 4-1 with monotonicity
add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o4_o1_mono1
  ARGS -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 1 --verify --var SV
       -f map_cgll_fv_o4_o1_mono1.nc
  LABELS "mbtempest_remap;mbtempest_se;mbtempest_options"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o4_o1_mono1
  OUTPUT map_cgll_fv_o4_o1_mono1.nc
  BASELINE gold_cgll_fv_o4_o1_mono1.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se;mbtempest_options"
  DEPENDS mbtempest_remap_cgll_fv_o4_o1_mono1)

# CGLL → FV, order 4-1 with nobubble
add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o4_o1_nobubble
  ARGS -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --nobubble --verify --var SV
       -f map_cgll_fv_o4_o1_nobubble.nc
  LABELS "mbtempest_remap;mbtempest_se;mbtempest_options"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o4_o1_nobubble
  OUTPUT map_cgll_fv_o4_o1_nobubble.nc
  BASELINE gold_cgll_fv_o4_o1_nobubble.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se;mbtempest_options"
  DEPENDS mbtempest_remap_cgll_fv_o4_o1_nobubble)

# --- Parallel CGLL-FV tests ---
# Partition SE source mesh and FV target mesh for 2-rank runs
add_partition_test(NAME mbtempest_partition_se_cgll_fv_np2
  NPROCS 2
  TARGET rll90.h5m SOURCE ne30_se_o4.h5m
  OUTPUT_TARGET rll90_p2_se.h5m
  LABELS "mbtempest_se")

# CGLL→FV order 4, parallel (2 ranks)
add_mbtempest_mpi_test(NAME mbtempest_remap_cgll_fv_o4_o1_np2
  NPROCS 2
  ARGS -t 5 -l ne30_se_o4_inferred.h5m -l rll90_p2_se.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV
       -f map_cgll_fv_o4_o1_np2.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4;mbtempest_partition_se_cgll_fv_np2")

# Serial vs parallel consistency for CGLL→FV
add_map_compare_test(NAME mbtempest_compare_cgll_fv_o4_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_cgll_fv_o4_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_cgll_fv_o4_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel;mbtempest_se"
  DEPENDS "mbtempest_remap_cgll_fv_o4_o1;mbtempest_remap_cgll_fv_o4_o1_np2")

# Partition SE order 2 source mesh
add_partition_test(NAME mbtempest_partition_se_cgll_o2_fv_np2
  NPROCS 2
  TARGET rll90.h5m SOURCE ne30_se_o2.h5m
  OUTPUT_TARGET rll90_p2_se_o2.h5m
  LABELS "mbtempest_se")

# CGLL→FV order 2, parallel (2 ranks)
add_mbtempest_mpi_test(NAME mbtempest_remap_cgll_fv_o2_o1_np2
  NPROCS 2
  ARGS -t 5 -l ne30_se_o2_inferred.h5m -l rll90_p2_se_o2.h5m -w -m cgll -m fv -o 2 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV
       -f map_cgll_fv_o2_o1_np2.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o2;mbtempest_partition_se_cgll_o2_fv_np2")

add_map_compare_test(NAME mbtempest_compare_cgll_fv_o2_o1_serial_vs_np2
  MAP1 ${MBTEMPEST_TEST_DIR}/map_cgll_fv_o2_o1.nc
  MAP2 ${MBTEMPEST_TEST_DIR}/map_cgll_fv_o2_o1_np2.nc
  TOL 1e-12
  LABELS "mbtempest_parallel;mbtempest_se"
  DEPENDS "mbtempest_remap_cgll_fv_o2_o1;mbtempest_remap_cgll_fv_o2_o1_np2")

# ==============================================================================
# Category 4: FV-CGLL Remapping (serial only)
#   Target is a CS mesh with GLOBAL_DOFS tag
# ==============================================================================

# FV → CGLL, order 1-4
add_mbtempest_test(NAME mbtempest_remap_fv_cgll_o1_o4
  ARGS -t 5 -l rll90.h5m -l ne30_se_o4.h5m -w -m fv -m cgll -o 1 -o 4
       -g GLOBAL_ID -g GLOBAL_DOFS --verify --var SV
       -f map_fv_cgll_o1_o4.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_fv_cgll_o1_o4
  OUTPUT map_fv_cgll_o1_o4.nc
  BASELINE gold_fv_cgll_o1_o4.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS mbtempest_remap_fv_cgll_o1_o4)

# ==============================================================================
# Category 5: CGLL-CGLL Remapping (serial only)
#   Both source and target are CS meshes with GLOBAL_DOFS tag
# ==============================================================================

# CGLL → CGLL, order 4-4 (CS30 → CS25)
add_mbtempest_test(NAME mbtempest_remap_cgll_cgll_o4_o4
  ARGS -t 5 -l ne30_se_o4.h5m -l ne25_se_o4.h5m -w -m cgll -m cgll -o 4 -o 4
       -g GLOBAL_DOFS -g GLOBAL_DOFS --verify --var SV
       -f map_cgll_cgll_o4_o4.nc
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4;mbtempest_se_convert_cs25_o4")

add_baseline_test(NAME mbtempest_gold_cgll_cgll_o4_o4
  OUTPUT map_cgll_cgll_o4_o4.nc
  BASELINE gold_cgll_cgll_o4_o4.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS mbtempest_remap_cgll_cgll_o4_o4)

# ==============================================================================
# Category 6: FV Sub-methods (--fvmethod)
# ==============================================================================

# Bilinear (skips intersection)
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_bilin
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod bilin
       --verify --var SH -f map_ne30_rll90_fv_fv_bilin.nc
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_bilin
  OUTPUT map_ne30_rll90_fv_fv_bilin.nc
  BASELINE gold_ne30_rll90_fv_fv_bilin.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_bilin)

# Delaunay (skips intersection)
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_delaunay
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod delaunay
       --verify --var SH -f map_ne30_rll90_fv_fv_delaunay.nc
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_delaunay
  OUTPUT map_ne30_rll90_fv_fv_delaunay.nc
  BASELINE gold_ne30_rll90_fv_fv_delaunay.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_delaunay)

# Integrated bilinear
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_intbilin
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod intbilin
       --verify --var SH -f map_ne30_rll90_fv_fv_intbilin.nc
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_intbilin
  OUTPUT map_ne30_rll90_fv_fv_intbilin.nc
  BASELINE gold_ne30_rll90_fv_fv_intbilin.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_intbilin)

# Inverse distance
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_invdist
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod invdist
       --verify --var SH -f map_ne30_rll90_fv_fv_invdist.nc
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_invdist
  OUTPUT map_ne30_rll90_fv_fv_invdist.nc
  BASELINE gold_ne30_rll90_fv_fv_invdist.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_fvmethod"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_invdist)

# ==============================================================================
# Category 7: Monotonicity Tests
# ==============================================================================

add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_mono1
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --monotonicity 1
       --verify --var SH -f map_ne30_rll90_fv_fv_mono1.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_mono1
  OUTPUT map_ne30_rll90_fv_fv_mono1.nc
  BASELINE gold_ne30_rll90_fv_fv_mono1.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_mono1)

# Note: --monotonicity 2 and 3 are only valid for SE (CGLL/DGLL) discretizations,
# not for FV-FV. These tests use CGLL source with SE-tagged mesh.

add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o4_o1_mono2
  ARGS -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 2 --verify --var SV
       -f map_cgll_fv_o4_o1_mono2.nc
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o4_o1_mono2
  OUTPUT map_cgll_fv_o4_o1_mono2.nc
  BASELINE gold_cgll_fv_o4_o1_mono2.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS mbtempest_remap_cgll_fv_o4_o1_mono2)

add_mbtempest_test(NAME mbtempest_remap_cgll_fv_o4_o1_mono3
  ARGS -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1
       -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 3 --verify --var SV
       -f map_cgll_fv_o4_o1_mono3.nc
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_cgll_fv_o4_o1_mono3
  OUTPUT map_cgll_fv_o4_o1_mono3.nc
  BASELINE gold_cgll_fv_o4_o1_mono3.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS mbtempest_remap_cgll_fv_o4_o1_mono3)

# ==============================================================================
# Category 8: CAAS Limiter Tests
# ==============================================================================

add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_caas1
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --limiter 1
       --verify --var SH -f map_ne30_rll90_fv_fv_caas1.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_caas1
  OUTPUT map_ne30_rll90_fv_fv_caas1.nc
  BASELINE gold_ne30_rll90_fv_fv_caas1.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_caas1)

add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_caas2
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --limiter 2
       --verify --var SH -f map_ne30_rll90_fv_fv_caas2.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_caas2
  OUTPUT map_ne30_rll90_fv_fv_caas2.nc
  BASELINE gold_ne30_rll90_fv_fv_caas2.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_caas2)

add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_caas3
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --limiter 3
       --verify --var SH -f map_ne30_rll90_fv_fv_caas3.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_caas3
  OUTPUT map_ne30_rll90_fv_fv_caas3.nc
  BASELINE gold_ne30_rll90_fv_fv_caas3.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_caas3)

add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_fv_fv_caas4
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --limiter 4
       --verify --var SH -f map_ne30_rll90_fv_fv_caas4.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_fv_fv_caas4
  OUTPUT map_ne30_rll90_fv_fv_caas4.nc
  BASELINE gold_ne30_rll90_fv_fv_caas4.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_fv_fv_caas4)

# ==============================================================================
# Category 9: Advanced Options
# ==============================================================================

# Gnomonic projection
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_gnomonic
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --gnomonic
       --verify --var SH -f map_ne30_rll90_gnomonic.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_gnomonic
  OUTPUT map_ne30_rll90_gnomonic.nc
  BASELINE gold_ne30_rll90_gnomonic.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_gnomonic)

# Advancing front algorithm
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_advfront
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --advfront
       --verify --var SH -f map_ne30_rll90_advfront.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_advfront
  OUTPUT map_ne30_rll90_advfront.nc
  BASELINE gold_ne30_rll90_advfront.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_advfront)

# No conservation
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_noconserve
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --noconserve
       --verify --var SH -f map_ne30_rll90_noconserve.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_noconserve
  OUTPUT map_ne30_rll90_noconserve.nc
  BASELINE gold_ne30_rll90_noconserve.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_noconserve)

# Enforce convexity
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_convexity
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --enforce_convexity
       --verify --var SH -f map_ne30_rll90_convexity.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_convexity
  OUTPUT map_ne30_rll90_convexity.nc
  BASELINE gold_ne30_rll90_convexity.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_convexity)

# No bubble (SE method)
# nobubble test moved to Category 3 (SE tests) with proper DOF tags

# Sparse constraints
add_mbtempest_test(NAME mbtempest_remap_ne30_rll90_sparse
  ARGS -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --sparseconstraints
       --verify --var SH -f map_ne30_rll90_sparse.nc
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS "")

add_baseline_test(NAME mbtempest_gold_ne30_rll90_sparse
  OUTPUT map_ne30_rll90_sparse.nc
  BASELINE gold_ne30_rll90_sparse.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options"
  DEPENDS mbtempest_remap_ne30_rll90_sparse)

# Volumetric projection (only valid for FV→CGLL or FV→DGLL, not FV→FV)
add_mbtempest_test(NAME mbtempest_remap_fv_cgll_volumetric
  ARGS -t 5 -l rll90.h5m -l ne30_se_o4.h5m -w -m fv -m cgll -o 1 -o 4
       -g GLOBAL_ID -g GLOBAL_DOFS --volumetric --verify --var SV
       -f map_fv_cgll_volumetric.nc
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4")

add_baseline_test(NAME mbtempest_gold_fv_cgll_volumetric
  OUTPUT map_fv_cgll_volumetric.nc
  BASELINE gold_fv_cgll_volumetric.nc
  TOL 1e-12
  LABELS "mbtempest_remap;mbtempest_options;mbtempest_se"
  DEPENDS mbtempest_remap_fv_cgll_volumetric)

# ==============================================================================
# Category 10: Convergence Analysis
#   Each convergence test runs mbtempest at two different mesh resolutions (or
#   discretization orders) and verifies that the L_inf error decreases.
#   The test is a single shell command that runs both cases and compares.
# ==============================================================================

set(CONV_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/run_convergence_test.sh)

# Helper macro: run mbtempest at two configurations and verify error decreases.
# The bash script runs both, extracts L_inf error, and compares.
# Usage: add_convergence_test(NAME name
#          COARSE_ARGS "arg1;arg2;..." FINE_ARGS "arg1;arg2;..." [DEPENDS ...])
macro(add_convergence_test)
  cmake_parse_arguments(CV "" "NAME" "COARSE_ARGS;FINE_ARGS;LABELS;DEPENDS" ${ARGN})
  if(MOAB_HAVE_MPI AND MPI_FOUND)
    set(_mpi_prefix "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
  else()
    set(_mpi_prefix "")
  endif()
  add_test(NAME ${CV_NAME}
    COMMAND bash ${CONV_SCRIPT} ${CV_NAME} ${MBTEMPEST_EXE} "${_mpi_prefix}"
            -- ${CV_COARSE_ARGS} -- ${CV_FINE_ARGS}
    WORKING_DIRECTORY ${MBTEMPEST_TEST_DIR}
  )
  set_tests_properties(${CV_NAME} PROPERTIES
    LABELS "mbtempest;mbtempest_convergence;${CV_LABELS}"
    TIMEOUT 600
  )
  if(CV_DEPENDS)
    set_tests_properties(${CV_NAME} PROPERTIES
      DEPENDS "${CV_DEPENDS}"
    )
  endif()
endmacro()

# --- FV-FV order 1 vs order 2: higher order should have lower error ---
add_convergence_test(NAME mbtempest_conv_fv_o1_vs_o2
  COARSE_ARGS "-t;5;-l;ne30.h5m;-l;rll90.h5m;-w;-m;fv;-m;fv;-o;1;-o;1;--verify;--var;SH;-f;conv_fv_o1.nc"
  FINE_ARGS   "-t;5;-l;ne30.h5m;-l;rll90.h5m;-w;-m;fv;-m;fv;-o;2;-o;2;--verify;--var;SH;-f;conv_fv_o2.nc"
  LABELS "mbtempest_remap")

# --- FV-FV: order 2 vs order 3 (monotonic improvement with order) ---
add_convergence_test(NAME mbtempest_conv_fv_o2_vs_o3
  COARSE_ARGS "-t;5;-l;ne30.h5m;-l;icod30.h5m;-w;-m;fv;-m;fv;-o;2;-o;2;--verify;--var;SH;-f;conv_fv_icod30_o2.nc"
  FINE_ARGS   "-t;5;-l;ne30.h5m;-l;icod30.h5m;-w;-m;fv;-m;fv;-o;3;-o;3;--verify;--var;SH;-f;conv_fv_icod30_o3.nc"
  LABELS "mbtempest_remap")

# --- CGLL-FV order 2 vs order 4: higher order should have lower error ---
add_convergence_test(NAME mbtempest_conv_cgll_o2_vs_o4
  COARSE_ARGS "-t;5;-l;ne30_se_o2.h5m;-l;rll90.h5m;-w;-m;cgll;-m;fv;-o;2;-o;1;-g;GLOBAL_DOFS;-g;GLOBAL_ID;--verify;--var;SV;-f;conv_cgll_o2.nc"
  FINE_ARGS   "-t;5;-l;ne30_se_o4.h5m;-l;rll90.h5m;-w;-m;cgll;-m;fv;-o;4;-o;1;-g;GLOBAL_DOFS;-g;GLOBAL_ID;--verify;--var;SV;-f;conv_cgll_o4.nc"
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o2;mbtempest_se_convert_cs30_o4")

# --- CGLL-FV vs CGLL-CGLL: SE-SE should be more accurate than SE-FV ---
add_convergence_test(NAME mbtempest_conv_cgll_fv_vs_cgll_cgll
  COARSE_ARGS "-t;5;-l;ne30_se_o4.h5m;-l;rll90.h5m;-w;-m;cgll;-m;fv;-o;4;-o;1;-g;GLOBAL_DOFS;-g;GLOBAL_ID;--verify;--var;SV;-f;conv_cgll_fv.nc"
  FINE_ARGS   "-t;5;-l;ne30_se_o4.h5m;-l;ne25_se_o4.h5m;-w;-m;cgll;-m;cgll;-o;4;-o;4;-g;GLOBAL_DOFS;-g;GLOBAL_DOFS;--verify;--var;SV;-f;conv_cgll_cgll.nc"
  LABELS "mbtempest_remap;mbtempest_se"
  DEPENDS "mbtempest_se_convert_cs30_o4;mbtempest_se_convert_cs25_o4")

# ==============================================================================
# Summary
# ==============================================================================
#
# Test count breakdown:
#   Mesh generation:          10 tests
#   FV-FV serial remap:        7 tests + 7 gold comparisons
#   FV-FV parallel remap:     ~8 tests + ~5 serial-vs-parallel comparisons
#   SE-FV serial + parallel:   8 tests + 4 gold + 4 serial-vs-parallel
#   FV-SE serial only:         2 tests + 2 gold
#   SE-SE serial only:         4 tests + 4 gold
#   FV sub-methods:            4 tests + 4 gold
#   Monotonicity:              4 tests + 4 gold
#   CAAS limiter:              4 tests + 4 gold
#   Advanced options:          7 tests + 7 gold
#   Convergence:               3 additional tests + 1 gold
#   -----------------------------------------------
#   Total:                    ~95 test targets
#
