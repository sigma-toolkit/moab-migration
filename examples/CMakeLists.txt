cmake_minimum_required( VERSION 3.20 )

set ( MOAB_ROOT "" CACHE PATH "Path to MOAB library installation" )
project( MOAB_Examples NONE )

# Imported shared libraries require platform properties that are normally set
# when project() enables a language. This property is harmless for imported
# targets and allows the first package pass to discover MOAB's compilers.
set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
find_package(MOAB REQUIRED
  HINTS ${MOAB_ROOT} ${MOAB_ROOT}/lib/cmake/MOAB)

set ( CMAKE_C_COMPILER ${MOAB_CC} )
set ( CMAKE_CXX_COMPILER ${MOAB_CXX} )
set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )
set ( CMAKE_C_FLAGS ${MOAB_C_FLAGS} )
set ( CMAKE_CXX_FLAGS ${MOAB_CXX_FLAGS} )
set ( CMAKE_Fortran_FLAGS ${MOAB_FC_FLAGS} )
set ( CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} )

enable_language(C)
enable_language(CXX)
if ( ENABLE_IMESH )
  enable_language(Fortran)
endif ( ENABLE_IMESH )

# Re-run dependency discovery after the relevant languages are enabled.
find_package(MOAB REQUIRED
  HINTS ${MOAB_ROOT} ${MOAB_ROOT}/lib/cmake/MOAB)

message ( "++ Found MOAB configuration details at: ${MOAB_DIR}" )

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
  find_path( MBCOUPLER_INCLUDE_DIR mbcoupler/Coupler.hpp HINTS ${MOAB_ROOT}/include )
  find_library( MBCOUPLER_LIBRARY mbcoupler HINTS ${MOAB_ROOT}/lib )
  if ( MBCOUPLER_INCLUDE_DIR AND MBCOUPLER_LIBRARY )
    set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES} advanced/CoupleMGen.cpp )
    set ( CoupleMGen_DEPS ${MBCOUPLER_LIBRARY} )
    set ( CoupleMGen_DEFINES -DMOAB_HAVE_MBCOUPLER )
  endif ()
endif ( MOAB_USE_MPI )

if ( MOAB_USE_MPI AND MOAB_USE_TEMPESTREMAP )
  set ( ADVANCED_EXAMPLES ${ADVANCED_EXAMPLES}
                          advanced/ParallelRemapTemplate.cpp
    )
endif ( MOAB_USE_MPI AND MOAB_USE_TEMPESTREMAP )

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

include_directories(${MOAB_INCLUDE_DIRS})

set ( ALLCXXEXAMPLES ${BEGINNER_EXAMPLES} ${BASIC_EXAMPLES} ${INTERMEDIATE_EXAMPLES} ${ADVANCED_EXAMPLES} ${APPS_EXAMPLES} )

set ( CMAKE_LINKER ${MOAB_CXX} )
foreach( fname ${ALLCXXEXAMPLES} )
  string( REPLACE ".cpp" "" base ${fname} )
  get_filename_component(base ${base} NAME_WE)
  message("-- Driver [${base}] (C++ sources) :: ${fname}")
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES} ${${base}_DEFINES}" )
  if ( base STREQUAL "ParallelRemapTemplate" )
    set_target_properties( ${base} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF )
  endif()
  target_link_libraries( ${base} LINK_PUBLIC ${${base}_DEPS} ${MOAB_LIBRARIES} )
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
