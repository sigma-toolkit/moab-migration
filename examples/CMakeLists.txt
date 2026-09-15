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

if ( ENABLE_IMESH )
  enable_language(Fortran)
endif ( ENABLE_IMESH )

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

set ( FORTRAN_EXAMPLES
  fortran/DirectAccessNoHolesF90.F90
)
if ( MOAB_USE_MPI )
  set ( FORTRAN_EXAMPLES ${FORTRAN_EXAMPLES}
          fortran/PushParMeshIntoMoabF90.F90
          fortran/MigrateMesh.F90
    )
endif ( MOAB_USE_MPI )

set ( APPS_EXAMPLES
  apps/earth/ContinentsOnGlobe.cpp
  apps/earth/BoundaryContinents.cpp
  apps/earth/addPCdata.cpp
)

# MOAB's own headers, plus the third-party headers its public headers include
# (Eigen from Matrix3.hpp, NetCDF/TempestRemap from the climate headers, ...).
include_directories(${MOAB_INCLUDE_DIRS} ${MOAB_PACKAGE_INCLUDES})
if ( MOAB_USE_MPI )
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif ( MOAB_USE_MPI )

set ( ALLCXXEXAMPLES ${BEGINNER_EXAMPLES} ${BASIC_EXAMPLES} ${INTERMEDIATE_EXAMPLES} ${ADVANCED_EXAMPLES} ${APPS_EXAMPLES} )

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
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
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
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif ( ENABLE_IMESH )

