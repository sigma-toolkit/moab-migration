cmake_minimum_required( VERSION 3.5 )

set ( MOAB_ROOT "" CACHE PATH "Path to MOAB library installation" )

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
# libMOAB itself, so it needs a Fortran compiler but not iMesh.  Enable the
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
  advanced/LloydRelaxation.cpp
  advanced/VisTags.cpp
)
if ( MOAB_USE_MPI )
  set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES}
                          advanced/GenLargeMesh.cpp
    )
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

# Fortran examples that go through iMesh.  iMesh is only built by the autotools
# build of MOAB, so these are gated on ENABLE_IMESH further down.
set ( FORTRAN_EXAMPLES
  fortran/DirectAccessNoHolesF90.F90
)
if ( MOAB_USE_MPI )
  set ( FORTRAN_EXAMPLES ${FORTRAN_EXAMPLES}
          fortran/PushParMeshIntoMoabF90.F90
    )
endif ( MOAB_USE_MPI )

# Fortran examples that only use iMOAB, and so link against libMOAB directly.
set ( IMOAB_FORTRAN_EXAMPLES )
if ( MOAB_USE_MPI AND MOAB_EXAMPLES_HAVE_FORTRAN )
  set ( IMOAB_FORTRAN_EXAMPLES ${IMOAB_FORTRAN_EXAMPLES}
          fortran/MigrateMesh.F90
    )
endif ( MOAB_USE_MPI AND MOAB_EXAMPLES_HAVE_FORTRAN )

set ( APPS_EXAMPLES
  earthsystem/earth/ContinentsOnGlobe.cpp
  earthsystem/earth/BoundaryContinents.cpp
  earthsystem/earth/addPCdata.cpp
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

#set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )
if ( ENABLE_IMESH )
  foreach( fname ${FORTRAN_EXAMPLES} )
    string( REPLACE ".F90" "" base ${fname} )
    get_filename_component(base ${base} NAME_WE)
    message("-- Driver [${base}] (F90 sources) :: ${fname}")
    add_executable( ${base} ${fname})
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" LINKER_LANGUAGE Fortran )
    target_link_libraries( ${base} LINK_PUBLIC ${IMESH_LIBRARIES} )
    add_test( NAME ${base} COMMAND $<TARGET_FILE:${base}> ${${base}_TEST_ARGS} )
  endforeach()
endif ( ENABLE_IMESH )

# iMOAB Fortran examples.  These link against libMOAB like the C++ examples,
# but are compiled and linked as Fortran, so they additionally need the MPI
# Fortran bindings - MPI_CXX_LIBRARIES does not provide mpif.h or the
# mpi_*_ symbols.
foreach( fname ${IMOAB_FORTRAN_EXAMPLES} )
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

