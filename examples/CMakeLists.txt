cmake_minimum_required( VERSION 3.5 )

set ( MOAB_ROOT "" CACHE PATH "Path to MOAB library installation" )

# find_package(MOAB) has to run before project() below, because the compilers
# MOAB was built with are only known afterwards and project() bakes the compiler
# in permanently.  But MOABTargets.cmake declares the MOAB library with
# add_library(MOAB SHARED IMPORTED), and add_library() rejects SHARED while the
# global TARGET_SUPPORTS_SHARED_LIBS property is unset - which it is until a
# language is enabled - with "the target platform does not support dynamic
# linking".  The property is meaningless for an IMPORTED target since nothing is
# built from it, so seed it here to get through the declaration; project() then
# sets it authoritatively from the platform module and this value is discarded.
set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

find_package(MOAB REQUIRED
  HINTS ${MOAB_ROOT} ${MOAB_ROOT}/lib/cmake/MOAB)

# Reuse the compilers and flags MOAB itself was built with.  The names exported
# by MOABConfig.cmake are MOAB_CFLAGS/MOAB_CXXFLAGS/MOAB_FCFLAGS - with no
# underscore before FLAGS - so they have to be spelled exactly that way here.
set ( CMAKE_C_COMPILER ${MOAB_CC} )
set ( CMAKE_CXX_COMPILER ${MOAB_CXX} )
set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )
set ( CMAKE_C_FLAGS "${MOAB_CFLAGS}" )
set ( CMAKE_CXX_FLAGS "${MOAB_CXXFLAGS}" )
set ( CMAKE_Fortran_FLAGS "${MOAB_FCFLAGS}" )

message ( "++ Found MOAB configuration details at: ${MOAB_DIR}" )

project( MOAB_Examples )

# The add_test() calls below were dead without this: CTest never generated a
# test file for the examples tree, so `ctest` here reported "No tests found".
enable_testing()

# MigrateMesh.F90 only calls the iMOAB Fortran entry points, which live in
# libMOAB itself, so it needs a Fortran compiler.  Enable the
# language only when a compiler is actually available, so that a MOAB built
# without Fortran still configures.  This has to happen before find_package(MPI)
# below, otherwise FindMPI never probes the Fortran bindings.
include ( CheckLanguage )
check_language ( Fortran )
if ( CMAKE_Fortran_COMPILER )
  enable_language ( Fortran )
  set ( MOAB_EXAMPLES_HAVE_FORTRAN TRUE )
else ( CMAKE_Fortran_COMPILER )
  set ( MOAB_EXAMPLES_HAVE_FORTRAN FALSE )
  message ( STATUS "No Fortran compiler available; skipping the Fortran examples" )
endif ( CMAKE_Fortran_COMPILER )

# MOAB may have been configured against MPI through FindMPI rather than through
# the compiler wrappers, in which case MOAB_CC/MOAB_CXX are the plain compilers
# and nothing above supplies the MPI headers.  Locate MPI ourselves so that
# moab_mpi.h can find mpi.h either way.
if ( MOAB_USE_MPI )
  find_package( MPI REQUIRED )
endif ( MOAB_USE_MPI )

set ( BEGINNER_EXAMPLES
  beginner/mbex1.cpp
  beginner/mbex2.cpp
  beginner/mbex3.cpp
  beginner/mbex4.cpp
)

set ( BASIC_EXAMPLES
  basic/ErrorHandlingSimulation.cpp
  basic/GetEntities.cpp
  basic/HelloMOAB.cpp
  basic/ProjectToSphere.cpp
  basic/ReadWriteTest.cpp
  basic/SetsNTags.cpp
  basic/StructuredMeshSimple.cpp
  basic/TestErrorHandling.cpp
  basic/TestExodusII.cpp
  basic/UniformRefinement.cpp
)

if ( MOAB_USE_MPI )
  set ( BASIC_EXAMPLES ${BASIC_EXAMPLES}
                basic/HelloParMOAB.cpp
                basic/TestErrorHandlingPar.cpp
    )
endif ( MOAB_USE_MPI )

set ( INTERMEDIATE_EXAMPLES
  intermediate/DirectAccessNoHoles.cpp
  intermediate/DirectAccessWithHoles.cpp
  intermediate/LoadPartial.cpp
  intermediate/PointInElementSearch.cpp
)
if ( MOAB_USE_MPI )
  set ( INTERMEDIATE_EXAMPLES ${INTERMEDIATE_EXAMPLES}
                              intermediate/ReduceExchangeTags.cpp
    )
endif ( MOAB_USE_MPI )

set ( ADVANCED_EXAMPLES
  advanced/CrystalRouterExample.cpp
  advanced/ExtrudePoly.cpp
  # LloydRelaxation now lives in advanced/smoothing with the other smoothers,
  # but it stays in this unguarded list rather than in the MPI+Eigen3 block
  # below: it has #ifdef MOAB_HAVE_MPI guards and builds serially, which the
  # two AutoDiff optimization examples do not.
  advanced/smoothing/LloydRelaxation.cpp
  advanced/VisTags.cpp
)
if ( MOAB_USE_MPI )
  set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES}
                          advanced/GenLargeMesh.cpp
    )
  # advanced/smoothing was not referenced here at all, so an installed CMake
  # examples tree built none of it.  The two optimization examples take their
  # gradients from Eigen's unsupported AutoDiff module, so Eigen3 is a hard
  # requirement rather than an optimization.
  if ( MOAB_USE_EIGEN )
    set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES}
                            advanced/smoothing/MeshOptimization.cpp
                            advanced/smoothing/SphericalCVTOptimization.cpp
      )
  else ( MOAB_USE_EIGEN )
    message ( STATUS "Skipping the mesh optimization examples: MOAB was built without Eigen3" )
  endif ( MOAB_USE_EIGEN )
  # LaplacianSmoother is left to advanced/smoothing/makefile: it needs an input
  # mesh as well as HDF5, so it has nothing to run as an argument-free CTest.
  # CoupleMGen is guarded by MOAB_HAVE_MBCOUPLER and links against the separate
  # mbcoupler library, so it can only be built when MOAB was configured with the
  # mesh coupler.  "mbcoupler" on its own is not a target here, it has to be
  # resolved to the installed library.
  if ( MOAB_USE_MBCOUPLER )
    find_library ( MBCOUPLER_LIBRARY NAMES mbcoupler
                   HINTS ${MOAB_LIBRARY_DIRS} ${MOAB_ROOT}/lib )
    if ( MBCOUPLER_LIBRARY )
      set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES} advanced/CoupleMGen.cpp )
      set ( CoupleMGen_DEPS ${MBCOUPLER_LIBRARY} )
    else ( MBCOUPLER_LIBRARY )
      message ( STATUS "Skipping CoupleMGen: libmbcoupler was not found" )
    endif ( MBCOUPLER_LIBRARY )
  else ( MOAB_USE_MBCOUPLER )
    message ( STATUS "Skipping CoupleMGen: MOAB is not configured with the mesh coupler" )
  endif ( MOAB_USE_MBCOUPLER )
endif ( MOAB_USE_MPI )

# Fortran examples.  These all go through the iMOAB Fortran module, which lives
# in libMOAB itself.
set ( FORTRAN_EXAMPLES )
if ( MOAB_EXAMPLES_HAVE_FORTRAN )
  # DirectAccessF90 builds its mesh in memory and needs neither MPI nor an
  # input file, so it is gated on the Fortran compiler alone.
  set ( FORTRAN_EXAMPLES ${FORTRAN_EXAMPLES}
          fortran/DirectAccessF90.F90
    )
  if ( MOAB_USE_MPI )
    set ( FORTRAN_EXAMPLES ${FORTRAN_EXAMPLES}
            fortran/MigrateMesh.F90
      )
  endif ( MOAB_USE_MPI )
endif ( MOAB_EXAMPLES_HAVE_FORTRAN )

# Every source in earthsystem/earth is listed here, so that what CMake compiles
# matches what earthsystem/earth/makefile builds and what both build systems
# install.  The four E3SM post-processing tools take their meshes on the command
# line; run with no arguments - as CTest does below - they print their usage and
# exit 0, which is the same smoke test the other examples get.
set ( APPS_EXAMPLES
  earthsystem/earth/ContinentsOnGlobe.cpp
  earthsystem/earth/BoundaryContinents.cpp
  earthsystem/earth/addPCdata.cpp
  earthsystem/earth/AddFieldtoPG2.cpp
  earthsystem/earth/copyPartition.cpp
  earthsystem/earth/ExtractLand.cpp
  earthsystem/earth/addConnec.cpp
)

# MOAB's own headers, plus the third-party headers its public headers include
# (Eigen from Matrix3.hpp, NetCDF/TempestRemap from the Earth system headers, ...).
include_directories(${MOAB_INCLUDE_DIRS} ${MOAB_PACKAGE_INCLUDES})
if ( MOAB_USE_MPI )
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif ( MOAB_USE_MPI )

set ( ALLCXXEXAMPLES ${BEGINNER_EXAMPLES} ${BASIC_EXAMPLES} ${INTERMEDIATE_EXAMPLES} ${ADVANCED_EXAMPLES} ${APPS_EXAMPLES} )

# A few examples require a mesh file on the command line rather than falling
# back to a default, so give them one when CTest runs them. Same "<base>_..."
# convention as the per-example dependency lists above.
set ( UniformRefinement_TEST_ARGS "${MOAB_MESH_DIR}/32hex_ef.h5m" )

# ProjectToSphere requires -i; with no arguments it tries to load "" and exits
# non-zero.  sphere_tris_5.vtk is used rather than an .h5m so the test does not
# silently depend on HDF5.
set ( ProjectToSphere_TEST_ARGS "-i" "${MOAB_MESH_DIR}/sphere_tris_5.vtk"
                                "-o" "ProjectToSphere_out.vtk" "-R" "1.0" )

# MeshOptimization's default run would optimize a 10x10 mesh; --selftest is both
# faster and an actual check, since it verifies the metric against closed-form
# values and against Verdict, and the AutoDiff gradient against finite
# differences.  Keep the spherical case small so it stays a smoke test.
set ( MeshOptimization_TEST_ARGS "--selftest" )
set ( SphericalCVTOptimization_TEST_ARGS "-N" "120" "-l" "10" "-n" "20" )

set ( EXAMPLE_LIBRARIES ${MOAB_LIBRARIES} )
if ( MOAB_USE_MPI )
  list ( APPEND EXAMPLE_LIBRARIES ${MPI_CXX_LIBRARIES} )
endif ( MOAB_USE_MPI )

foreach( fname ${ALLCXXEXAMPLES} )
  string( REPLACE ".cpp" "" base ${fname} )
  get_filename_component(base ${base} NAME_WE)
  message("-- Driver [${base}] (C++ sources) :: ${fname}")
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
  target_link_libraries( ${base} LINK_PUBLIC ${${base}_DEPS} ${EXAMPLE_LIBRARIES} )
  add_test( NAME ${base} COMMAND $<TARGET_FILE:${base}> ${${base}_TEST_ARGS} )
endforeach()

# The Fortran examples link against libMOAB like the C++ examples, but are
# compiled and linked as Fortran, so they additionally need the MPI Fortran
# bindings - MPI_CXX_LIBRARIES does not provide mpif.h or the mpi_*_ symbols.
foreach( fname ${FORTRAN_EXAMPLES} )
  get_filename_component( base ${fname} NAME_WE )
  message("-- Driver [${base}] (F90 sources) :: ${fname}")
  add_executable( ${base} ${fname} )
  set_target_properties( ${base} PROPERTIES
    COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}"
    LINKER_LANGUAGE Fortran )
  target_include_directories( ${base} PRIVATE ${MPI_Fortran_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_PATH} )
  # The link is driven by the Fortran compiler, which does not pull in the C++
  # runtime on its own, and libMOAB is C++.  MOABTargets.cmake asks for
  # "stdc++" unconditionally, which is wrong for a libc++ toolchain, so use the
  # runtime CMake actually detected for this compiler instead.
  target_link_libraries( ${base} LINK_PUBLIC ${EXAMPLE_LIBRARIES} ${MPI_Fortran_LIBRARIES}
                         ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} )
endforeach()

if ( TARGET MigrateMesh )
  # MigrateMesh reads spherecube.h5m from the working directory and writes its
  # output beside it, so stage the input in the build tree rather than running
  # the test inside the source tree.
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/fortran/spherecube.h5m
                  ${CMAKE_CURRENT_BINARY_DIR}/spherecube.h5m COPYONLY )
  add_test( NAME MigrateMesh
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:MigrateMesh> ${MPIEXEC_POSTFLAGS}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
endif ( TARGET MigrateMesh )

if ( TARGET DirectAccessF90 )
  add_test( NAME DirectAccessF90 COMMAND $<TARGET_FILE:DirectAccessF90> )
endif ( TARGET DirectAccessF90 )

