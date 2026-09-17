SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR}
                    parallel
                    io io/mhdf
                    earthsystem/io
                    ${MOAB_BINARY_DIR}/src/parallel )

if(MOAB_HAVE_MPI AND MPI_FOUND)
  add_subdirectory(parallel)
  set( MOAB_MPI moabpar )
  set( MOABPAR_OBJECTS $<TARGET_OBJECTS:moabpar> )
endif()

#Create the moab/EntityHandle.hpp file
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/moab/EntityHandle.hpp.in infile2)
string( REGEX REPLACE "#undef" "#cmakedefine" outstring2 "${infile2}" )
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp.in" "${outstring2}")
configure_file(${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp)

set(moab_link_lib)
add_subdirectory(io)
add_subdirectory(earthsystem)

if( MOAB_HAVE_HDF5 )
  include_directories(${HDF5_INCLUDES})
  set( MOAB_HDF_LIBS ${HDF5_LIBRARIES} )
endif()
if( MOAB_HAVE_NETCDF )
  include_directories(${NETCDF_INCLUDES})
  if ( MOAB_HAVE_MPI )
    set( MOAB_NETCDF_LIBS ${PNETCDF_LIBRARIES} ${NETCDF_LIBRARIES} )
  else ( MOAB_HAVE_MPI )
    set( MOAB_NETCDF_LIBS ${NETCDF_LIBRARIES} )
  endif ( MOAB_HAVE_MPI )
endif()


set(MOAB_LIB_SOURCES
        AdaptiveKDTree.cpp
        AEntityFactory.hpp AEntityFactory.cpp
        AffineXform.hpp    AffineXform.cpp
        AxisBox.hpp        AxisBox.cpp
        BitPage.hpp        BitPage.cpp
        BitTag.hpp         BitTag.cpp
        BoundBox.cpp
        BSPTree.cpp
        BSPTreePoly.cpp
        BVHTree.cpp
        CN.cpp
        CartVect.cpp
        Core.cpp
        DebugOutput.hpp DebugOutput.cpp
        DenseTag.hpp    DenseTag.cpp
        DualTool.cpp
        ElementSequence.hpp
        EntitySequence.hpp  EntitySequence.cpp
        ErrorHandler.cpp
        ErrorOutput.hpp  ErrorOutput.cpp
        Factory.cpp
        FBEngine.cpp
        FileOptions.cpp
        GeomUtil.cpp
        GeomTopoTool.cpp
        GeomQueryTool.cpp
        HalfFacetRep.cpp
        HigherOrderFactory.cpp
        HomXform.cpp
        iMOAB.cpp
        Internals.hpp
        LloydSmoother.cpp
        MBCNArrays.hpp
        MergeMesh.cpp
        MeshSet.hpp         MeshSet.cpp
        MeshSetSequence.hpp MeshSetSequence.cpp
        MeshTag.hpp         MeshTag.cpp
        MeshTopoUtil.cpp
        OrientedBox.cpp
        OrientedBoxTreeTool.cpp
        lotte/poly.c
        lotte/findpt.c
        lotte/errmem.c
        lotte/tensor.c
        local_discretization/ElemEvaluator.cpp
        local_discretization/LinearHex.cpp
        local_discretization/LinearQuad.cpp
        local_discretization/LinearTet.cpp
        local_discretization/LinearTri.cpp
        local_discretization/QuadraticHex.cpp
        MeshGeneration.cpp
        NestedRefine.cpp
        PolyElementSeq.hpp  PolyElementSeq.cpp
        ProgOptions.cpp
        Range.cpp
        RangeSeqIntersectIter.hpp RangeSeqIntersectIter.cpp
        ReadUtil.hpp              ReadUtil.cpp
        ReaderWriterSet.cpp
        ReorderTool.cpp
        ScdElementData.hpp        ScdElementData.cpp
        ScdInterface.cpp
        ScdVertexData.hpp         ScdVertexData.cpp
        SequenceData.hpp          SequenceData.cpp
        SequenceManager.hpp       SequenceManager.cpp
        SetIterator.cpp
        Skinner.cpp
        SmoothCurve.hpp           SmoothCurve.cpp
        SmoothFace.hpp            SmoothFace.cpp
        SparseTag.hpp             SparseTag.cpp
        SpatialLocator.cpp
        SpectralMeshTool.cpp
        StructuredElementSeq.hpp  StructuredElementSeq.cpp
        SweptElementData.hpp      SweptElementData.cpp
        SweptElementSeq.hpp       SweptElementSeq.cpp
        SweptVertexData.hpp       SweptVertexData.cpp
        SysUtil.hpp               SysUtil.cpp
        TagCompare.hpp
        TagInfo.hpp               TagInfo.cpp
        TupleList.cpp
        Tree.cpp
        Types.cpp
        TypeSequenceManager.hpp   TypeSequenceManager.cpp
        UnstructuredElemSeq.hpp   UnstructuredElemSeq.cpp
        Util.cpp
        VarLenDenseTag.hpp        VarLenDenseTag.cpp
        VarLenSparseTag.hpp       VarLenSparseTag.cpp
        VarLenTag.hpp
        VertexSequence.hpp        VertexSequence.cpp
        WriteUtil.hpp             WriteUtil.cpp
        verdict/v_vector.h
        verdict/V_EdgeMetric.cpp
        verdict/V_GaussIntegration.cpp
        verdict/V_HexMetric.cpp
        verdict/V_KnifeMetric.cpp
        verdict/V_PyramidMetric.cpp
        verdict/V_QuadMetric.cpp
        verdict/V_TetMetric.cpp
        verdict/V_TriMetric.cpp
        verdict/V_WedgeMetric.cpp
        verdict/VerdictVector.cpp
        verdict/VerdictWrapper.cpp
        verdict/verdict_defines.hpp
        verdict/VerdictVector.hpp
        verdict/V_GaussIntegration.hpp
        moab_mpe.h )

set( MOAB_INSTALL_HEADERS
        moab/AdaptiveKDTree.hpp
        moab/BoundBox.hpp
        moab/BSPTree.hpp
        moab/BSPTreePoly.hpp
        moab/BVHTree.hpp
        moab/CN.hpp
        moab/CartVect.hpp
        moab/Compiler.hpp
        moab/Core.hpp
        moab/CpuTimer.hpp
        moab/DualTool.hpp
        moab/Error.hpp
        moab/GeomTopoTool.hpp
        moab/GeomQueryTool.hpp
        moab/HalfFacetRep.hpp
        moab/HigherOrderFactory.hpp
        moab/HomXform.hpp
        moab/EntityType.hpp
        moab/ErrorHandler.hpp
        moab/FBEngine.hpp
        moab/FileOptions.hpp
        moab/FindPtFuncs.h
        moab/Forward.hpp
        moab/GeomUtil.hpp
        moab/iMOAB.h
        moab/Interface.hpp
        moab/earthsystem/intx_mesh/Intx2Mesh.hpp
        moab/earthsystem/intx_mesh/IntxUtils.hpp
        moab/earthsystem/intx_mesh/Intx2MeshInPlane.hpp
        moab/earthsystem/intx_mesh/Intx2MeshOnSphere.hpp
        moab/earthsystem/intx_mesh/Intx2MeshEdges.hpp
        moab/earthsystem/intx_mesh/IntxRllCssphere.hpp
        moab/LloydSmoother.hpp
        moab/point_locater/tree/common_tree.hpp
        moab/point_locater/tree/element_tree.hpp
        moab/point_locater/tree/bvh_tree.hpp
        moab/point_locater/io.hpp
        moab/point_locater/element_maps/linear_hex_map.hpp
        moab/point_locater/element_maps/linear_tet_map.hpp
        moab/point_locater/element_maps/spectral_hex_map.hpp
        moab/point_locater/element_maps/quadratic_hex_map.hpp
        moab/point_locater/point_locater.hpp
        moab/point_locater/parametrizer.hpp
        moab/discrete_geometry/DGMSolver.hpp
        moab/discrete_geometry/HiReconstruction.hpp
        moab/ElemEvaluator.hpp
        moab/local_discretization/LinearHex.hpp
        moab/local_discretization/LinearQuad.hpp
        moab/local_discretization/LinearTet.hpp
        moab/local_discretization/LinearTri.hpp
        moab/local_discretization/QuadraticHex.hpp
        moab/Matrix3.hpp
        moab/MergeMesh.hpp
        moab/MeshTopoUtil.hpp
        moab/MeshGeneration.hpp
        moab/OrientedBox.hpp
        moab/OrientedBoxTreeTool.hpp
        moab/PartitionerBase.hpp
        moab/ProgOptions.hpp
        moab/Range.hpp
        moab/RangeMap.hpp
        moab/ReadUtilIface.hpp
        moab/ReaderIface.hpp
        moab/ReaderWriterSet.hpp
        moab/ReorderTool.hpp
        moab/NestedRefine.hpp
        moab/NestedRefineTemplates.hpp
        moab/ScdInterface.hpp
        moab/SetIterator.hpp
        moab/Skinner.hpp
        moab/SpatialLocator.hpp
        moab/SpatialLocatorTimes.hpp
        moab/SpectralMeshTool.hpp
        moab/Tree.hpp
        moab/TreeStats.hpp
        moab/TupleList.hpp
        moab/Types.hpp
        moab/mhdf_public.h
        moab/UnknownInterface.hpp
        moab/Util.hpp
        moab/WriteUtilIface.hpp
        moab/WriterIface.hpp
        moab/verdict.h
        moab/verdict/VerdictWrapper.hpp
        moab/win32_config.h
        MBEntityType.h
        MBCN.h
        MBCN_protos.h
        MBTagConventions.hpp )


if ( ENABLE_FORTRAN )
  set(MOAB_LIB_SOURCES ${MOAB_LIB_SOURCES} iMOABF.F90 )
  set_source_files_properties(iMOABF.F90  PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/imoab.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/)
endif ( ENABLE_FORTRAN )

set( MOAB_GENERATED_INSTALL_HEADERS
     ${CMAKE_CURRENT_BINARY_DIR}/moab/EntityHandle.hpp
     ${CMAKE_CURRENT_BINARY_DIR}/moab/MOABConfig.h
     ${CMAKE_CURRENT_BINARY_DIR}/moab_export.h
  )

if (MOAB_HAVE_ZOLTAN)
  set( MOAB_LIB_SOURCES ${MOAB_LIB_SOURCES}
          ZoltanPartitioner.cpp )
  set( MOAB_INSTALL_HEADERS ${MOAB_INSTALL_HEADERS}
          moab/ZoltanPartitioner.hpp )
  include_directories( ${ZOLTAN_INCLUDES} )
endif()

if (MOAB_HAVE_METIS)
  set( MOAB_LIB_SOURCES ${MOAB_LIB_SOURCES}
        MetisPartitioner.cpp )
  set( MOAB_INSTALL_HEADERS ${MOAB_INSTALL_HEADERS}
        moab/MetisPartitioner.hpp )

endif()

if (MOAB_HAVE_MPI AND MPI_FOUND)
  # ParCommGraph builds in src/earthsystem/parallel; only the public header is
  # declared here so it is installed with the rest.
  set( MOAB_INSTALL_HEADERS ${MOAB_INSTALL_HEADERS}
        moab/earthsystem/parallel/ParCommGraph.hpp )
endif()

if (MOAB_HAVE_TEMPESTREMAP)
  # Remapping sources build in src/earthsystem/remapping; only the public
  # headers are declared here so they are installed with the rest.
  set( MOAB_INSTALL_HEADERS ${MOAB_INSTALL_HEADERS}
        moab/earthsystem/remapping/Remapper.hpp
        moab/earthsystem/remapping/TempestRemapper.hpp
        moab/earthsystem/remapping/TempestOnlineMap.hpp
        moab/earthsystem/remapping/IntegerReprosum.hpp )

endif (MOAB_HAVE_TEMPESTREMAP)

add_library(
  MOAB
  ${MOAB_LIB_SOURCES} ${MOAB_INSTALL_HEADERS} ${MOAB_GENERATED_INSTALL_HEADERS}
  $<TARGET_OBJECTS:moabio>
  ${MOAB_EARTHSYSTEM_OBJECTS}
  ${MOABPAR_OBJECTS}
  ${MDHF_OBJ}
  )

include(CMakeLibraryFuncs)
SciFixupLibrary(MOAB MOAB MOAB)

# MOAB::MOAB is the name find_package() gives consumers, so offer it here too.
# A project that pulls MOAB in with add_subdirectory() then links the same name
# as one that finds an installed MOAB.
add_library(MOAB::MOAB ALIAS MOAB)

# Public headers are reached as <moab/Core.hpp> and <MBEntityType.h>, so the
# include roots are this directory and its binary counterpart, which is where
# MOABConfig.h and the generate_export_header() output land.  Only the
# INSTALL_INTERFACE was set before, which left the exported build-tree target
# with no usable include path at all.
target_include_directories(MOAB
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

generate_export_header(MOAB)

if(PYMOAB_LIBRARY_RPATH)
  set_target_properties(
    MOAB
    PROPERTIES
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
    CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}"
    INSTALL_RPATH "${PYMOAB_LIBRARY_RPATH}"
    )
else()
  set_target_properties(
    MOAB
    PROPERTIES
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
    CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}"
    VERSION ${MOAB_VERSION}
    SOVERSION ${MOAB_VERSION_MAJOR}
    )
endif()

# MOAB's C++11 use is visible in its installed headers, so record it as a usage
# requirement rather than relying on the consumer to have set CMAKE_CXX_STANDARD
# the way MOAB's own build did.
target_compile_features(MOAB PUBLIC cxx_std_11)

# Third-party dependencies are linked as imported targets, not as the absolute
# paths they resolve to on this machine.  An absolute path is baked verbatim
# into the exported link interface - including SDK-internal .tbd stubs on macOS
# - so an installed MOAB could only ever be consumed from the machine it was
# built on.  A target name is re-resolved by the consumer: MOABConfig.cmake runs
# the matching find_dependency() calls to recreate exactly these names.
#
# PUBLIC vs PRIVATE follows what MOAB's *installed* headers include:
#
#   moab/mhdf_public.h                       <H5Tpublic.h>        HDF5
#   moab/ZoltanPartitioner.hpp               "zoltan_cpp.h"       Zoltan
#   moab/MetisPartitioner.hpp                "metis.h"            Metis
#   moab/Matrix3.hpp                         "Eigen/Dense"        Eigen3
#   .../remapping/TempestRemapper.hpp        "TempestRemapAPI.h"  TempestRemap
#   .../remapping/TempestRemapper.hpp        "netcdfcpp.h"        NetCDF
#   moab_mpi.h                               <mpi.h>              MPI
#
# Everything else is reached only from MOAB's own .cpp files.  PRIVATE still
# puts those in the exported interface - libMOAB is normally static, so the
# consumer has to link them - but wrapped in $<LINK_ONLY:>, so their headers do
# not leak into the consumer's include path.
set(MOAB_PUBLIC_LIBS "")
set(MOAB_PRIVATE_LIBS "")

if (MOAB_HAVE_TEMPESTREMAP AND TARGET TempestRemap::TempestRemap)
  list(APPEND MOAB_PUBLIC_LIBS TempestRemap::TempestRemap)
endif()
if (MOAB_HAVE_NETCDF AND TARGET NetCDF::NetCDF)
  list(APPEND MOAB_PUBLIC_LIBS NetCDF::NetCDF)
endif()
if (MOAB_HAVE_PNETCDF AND TARGET PNetCDF::PNetCDF)
  list(APPEND MOAB_PRIVATE_LIBS PNetCDF::PNetCDF)
endif()
if (MOAB_HAVE_HDF5 AND TARGET HDF5::HDF5)
  list(APPEND MOAB_PUBLIC_LIBS HDF5::HDF5)
endif()
if (CURL_FOUND AND TARGET CURL::libcurl)
  list(APPEND MOAB_PRIVATE_LIBS CURL::libcurl)
endif()
if (MOAB_HAVE_ZOLTAN AND TARGET Zoltan::Zoltan)
  list(APPEND MOAB_PUBLIC_LIBS Zoltan::Zoltan)
endif()
if (MOAB_HAVE_PARMETIS AND TARGET ParMETIS::ParMETIS)
  list(APPEND MOAB_PRIVATE_LIBS ParMETIS::ParMETIS)
endif()
if (MOAB_HAVE_METIS AND TARGET METIS::METIS)
  list(APPEND MOAB_PUBLIC_LIBS METIS::METIS)
endif()
if (MOAB_HAVE_EIGEN3 AND TARGET Eigen3::Eigen)
  list(APPEND MOAB_PUBLIC_LIBS Eigen3::Eigen)
endif()
if (MOAB_HAVE_MPI AND TARGET MPI::MPI_CXX)
  list(APPEND MOAB_PUBLIC_LIBS MPI::MPI_CXX)
endif()
if (MOAB_HAVE_LAPACK AND TARGET LAPACK::LAPACK)
  list(APPEND MOAB_PRIVATE_LIBS LAPACK::LAPACK)
endif()
if (MOAB_HAVE_BLAS AND TARGET BLAS::BLAS)
  list(APPEND MOAB_PRIVATE_LIBS BLAS::BLAS)
endif()
if (TARGET Threads::Threads)
  list(APPEND MOAB_PRIVATE_LIBS Threads::Threads)
endif()

# When Fortran-driven executables link against libMOAB (a C++ static library),
# the Fortran linker does not pull in the C++ standard library on its own.
# CMAKE_CXX_IMPLICIT_LINK_LIBRARIES is what the C++ compiler actually links, so
# it names the right runtime - "stdc++" was hard-coded here before and was wrong
# on every libc++ toolchain, yet the -lstdc++ compile probe still succeeded on
# macOS because the driver silently maps it to libc++.  Keep an explicit stdc++
# only for toolchains that report nothing (Cray wrappers, some NVHPC setups).
#
# Only the runtime is wanted, not the whole implicit set: with an MPI compiler
# wrapper as CXX that set also names mpicxx, mpi and pmpi, which already arrive
# through MPI::MPI_CXX and which, as bare names in the exported interface, the
# consumer's linker would look for in its default library path.
if (ENABLE_FORTRAN)
  set(MOAB_CXX_RUNTIME_LIBS "")
  foreach (lib ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
    if (lib MATCHES "^(c\\+\\+|stdc\\+\\+|c\\+\\+abi|supc\\+\\+)$")
      list(APPEND MOAB_CXX_RUNTIME_LIBS ${lib})
    endif()
  endforeach()
  if (MOAB_CXX_RUNTIME_LIBS)
    list(APPEND MOAB_PRIVATE_LIBS ${MOAB_CXX_RUNTIME_LIBS})
  else()
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES "stdc++")
    check_cxx_source_compiles("int main() { return 0; }" MOAB_HAVE_LIBSTDCXX)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if (MOAB_HAVE_LIBSTDCXX)
      list(APPEND MOAB_PRIVATE_LIBS stdc++)
    endif()
  endif()
endif()

# Filter the two lists: remove bare directory paths and unrecognized linker
# flags (e.g. -fortranlibs from NVHPC LAPACK) that CMake cannot handle in
# target_link_libraries.  Platform-specific linker syntax such as Apple
# "-framework Accelerate", GNU "-Wl,...", "-Xlinker ...", "-pthread" is
# preserved via the shared keep-pattern table in ResolveCompilerPaths.cmake.
# Add new exceptions there, not here.
include(ResolveCompilerPaths)
filter_link_libraries(MOAB_PUBLIC_LIBS)
filter_link_libraries(MOAB_PRIVATE_LIBS)

# clean_library_list() used to run here.  It exists to preserve the meaning of
# the "debug"/"optimized"/"general" keywords while de-duplicating, which it does
# by rewriting every unkeyworded entry into a debug/optimized pair - so each
# library came out of the export twice, wrapped in $<$<CONFIG:DEBUG>:...> and
# $<$<NOT:$<CONFIG:DEBUG>>:...>.  Neither list below uses those keywords, so a
# plain de-duplication is both correct and far easier to read in MOABTargets.
list(REMOVE_DUPLICATES MOAB_PUBLIC_LIBS)
list(REMOVE_DUPLICATES MOAB_PRIVATE_LIBS)

target_link_libraries(MOAB
  PUBLIC ${MOAB_PUBLIC_LIBS}
  PRIVATE ${MOAB_PRIVATE_LIBS})

# If using clang-tidy, set the target property so that
# reformatting code does not use full fledged C++11 semantics
# set_property(TARGET MOAB PROPERTY CXX_STANDARD 98)

install(
  TARGETS MOAB
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${SKBUILD_SUBDIR}${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${SKBUILD_SUBDIR}${CMAKE_INSTALL_LIBDIR}
  )

list(APPEND MOAB_INSTALL_TARGETS  MOAB)
include(AutoconfHeader)
moab_install_headers(${MOAB_INSTALL_HEADERS})
foreach (header ${MOAB_GENERATED_INSTALL_HEADERS})
    STRING(REGEX REPLACE "^${CMAKE_CURRENT_BINARY_DIR}/" "" relative ${header})
    string(REGEX MATCH "(.*)/" DIR ${relative})
    install(
      FILES ${header}
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${DIR}
      )
endforeach()
