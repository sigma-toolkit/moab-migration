# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

# Things to build
lib_LTLIBRARIES = libMOAB.la
AM_CPPFLAGS =
libMOAB_la_LDFLAGS = -version-info @LIBTOOL_VERSION@

# Subdirectories to build
SUBDIRS = io earthsystem .
libMOAB_la_LIBADD = io/libmoabio.la earthsystem/io/libmoabearthsystemio.la
EXTRA_DIST =

# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB \
               -DLOCDIR=src/ \
               -I. \
               -I$(srcdir)/io \
               -I$(srcdir)/io/mhdf/include \
               -I$(srcdir)/earthsystem/io \
               -I$(srcdir)/parallel \
               -I$(srcdir)/local_discretization \
               -I$(srcdir)/verdict \
               -I$(srcdir)/discrete_geometry \
               -I$(srcdir)/earthsystem/intx_mesh

# The directory in which to install headers
libMOAB_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMOAB_la_SOURCES = \
  AdaptiveKDTree.cpp \
  AEntityFactory.cpp \
  AEntityFactory.hpp \
  AffineXform.cpp \
  AffineXform.hpp \
  AxisBox.cpp \
  AxisBox.hpp \
  BitPage.cpp \
  BitPage.hpp \
  BitTag.cpp \
  BitTag.hpp \
  BoundBox.cpp \
  BSPTree.cpp \
  BSPTreePoly.cpp \
  BVHTree.cpp \
  CN.cpp \
  CartVect.cpp \
  Core.cpp \
  DebugOutput.cpp \
  DebugOutput.hpp \
  DenseTag.hpp \
  DenseTag.cpp \
  DualTool.cpp \
  discrete_geometry/DGMSolver.cpp \
  discrete_geometry/HiReconstruction.cpp \
  ElementSequence.hpp \
  EntitySequence.cpp \
  EntitySequence.hpp \
  ErrorHandler.cpp \
  ErrorOutput.cpp \
  ErrorOutput.hpp \
  Factory.cpp \
  FBEngine.cpp \
  FileOptions.cpp \
  GeomUtil.cpp \
  GeomQueryTool.cpp \
  GeomTopoTool.cpp \
  HalfFacetRep.cpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  iMOAB.cpp \
  Internals.hpp \
  earthsystem/intx_mesh/Intx2Mesh.cpp  \
  earthsystem/intx_mesh/Intx2MeshInPlane.cpp  \
  earthsystem/intx_mesh/Intx2MeshOnSphere.cpp \
  earthsystem/intx_mesh/Intx2MeshEdges.cpp \
  earthsystem/intx_mesh/IntxRllCssphere.cpp \
  earthsystem/intx_mesh/IntxUtils.cpp \
  LloydSmoother.cpp \
  MBCNArrays.hpp \
  MergeMesh.cpp \
  MeshSet.cpp \
  MeshSet.hpp \
  MeshSetSequence.cpp \
  MeshSetSequence.hpp \
  MeshTag.hpp \
  MeshTag.cpp \
  MeshTopoUtil.cpp \
  NestedRefine.cpp \
  OrientedBox.cpp \
  OrientedBoxTreeTool.cpp \
  lotte/poly.c \
  lotte/findpt.c \
  lotte/errmem.c \
  lotte/tensor.c \
  local_discretization/ElemEvaluator.cpp \
  local_discretization/LinearHex.cpp \
  local_discretization/LinearQuad.cpp \
  local_discretization/LinearTet.cpp \
  local_discretization/LinearTri.cpp \
  local_discretization/QuadraticHex.cpp \
  MeshGeneration.cpp \
  PolyElementSeq.cpp \
  PolyElementSeq.hpp \
  ProgOptions.cpp \
  Range.cpp \
  RangeSeqIntersectIter.cpp \
  RangeSeqIntersectIter.hpp \
  ReadUtil.cpp \
  ReadUtil.hpp \
  ReaderWriterSet.cpp \
  ReorderTool.cpp \
  ScdElementData.cpp \
  ScdElementData.hpp \
  ScdInterface.cpp \
  ScdVertexData.cpp \
  ScdVertexData.hpp \
  SequenceData.hpp \
  SequenceData.cpp \
  SequenceManager.cpp \
  SequenceManager.hpp \
  SetIterator.cpp \
  Skinner.cpp \
  SmoothCurve.cpp \
  SmoothCurve.hpp \
  SmoothFace.cpp \
  SmoothFace.hpp \
  SparseTag.cpp \
  SparseTag.hpp \
  SpatialLocator.cpp \
  SpectralMeshTool.cpp \
  StructuredElementSeq.cpp \
  StructuredElementSeq.hpp \
  SweptElementData.cpp \
  SweptElementData.hpp \
  SweptElementSeq.cpp \
  SweptElementSeq.hpp \
  SweptVertexData.cpp \
  SweptVertexData.hpp \
  SysUtil.cpp \
  SysUtil.hpp \
  TagCompare.hpp \
  TagInfo.cpp \
  TagInfo.hpp \
  Tree.cpp \
  TupleList.cpp \
  Types.cpp \
  TypeSequenceManager.cpp \
  TypeSequenceManager.hpp \
  UnstructuredElemSeq.cpp \
  UnstructuredElemSeq.hpp \
  Util.cpp \
  VarLenDenseTag.hpp \
  VarLenDenseTag.cpp \
  VarLenSparseTag.hpp \
  VarLenSparseTag.cpp \
  VarLenTag.hpp \
  VertexSequence.hpp \
  VertexSequence.cpp \
  verdict/V_EdgeMetric.cpp \
  verdict/V_GaussIntegration.cpp \
  verdict/V_HexMetric.cpp \
  verdict/V_KnifeMetric.cpp \
  verdict/V_PyramidMetric.cpp \
  verdict/V_QuadMetric.cpp \
  verdict/V_TetMetric.cpp \
  verdict/V_TriMetric.cpp \
  verdict/V_WedgeMetric.cpp \
  verdict/V_GaussIntegration.hpp \
  verdict/verdict_defines.hpp \
  verdict/v_vector.h \
  verdict/VerdictVector.cpp \
  verdict/VerdictVector.hpp \
  verdict/VerdictWrapper.cpp \
  WriteUtil.cpp \
  WriteUtil.hpp \
  moab_mpe.h

nobase_includedir = $(includedir)/moab
nodist_nobase_include_HEADERS = moab/MOABConfig.h

# The list of header files which are to be installed
nobase_libMOAB_la_include_HEADERS = \
  moab/AdaptiveKDTree.hpp \
  moab/BoundBox.hpp \
  moab/BSPTree.hpp \
  moab/BSPTreePoly.hpp \
  moab/BVHTree.hpp \
  moab/CN.hpp \
  moab/CartVect.hpp \
  moab/Compiler.hpp \
  moab/Core.hpp \
  moab/CpuTimer.hpp \
  moab/DualTool.hpp \
  moab/Error.hpp \
  moab/GeomQueryTool.hpp \
  moab/GeomTopoTool.hpp \
  moab/HalfFacetRep.hpp \
  moab/HigherOrderFactory.hpp \
  moab/HomXform.hpp \
  moab/EntityType.hpp \
  moab/EntityHandle.hpp \
  moab/ErrorHandler.hpp \
  moab/FBEngine.hpp \
  moab/FileOptions.hpp \
  moab/FindPtFuncs.h \
  moab/Forward.hpp \
  moab/GeomUtil.hpp \
  moab/Interface.hpp \
  moab/iMOAB.h \
  moab/LloydSmoother.hpp \
  moab/NestedRefine.hpp \
  moab/NestedRefineTemplates.hpp \
  moab/Matrix3.hpp \
  moab/MergeMesh.hpp \
  moab/MeshTopoUtil.hpp \
  moab/MeshGeneration.hpp \
  moab/OrientedBox.hpp \
  moab/OrientedBoxTreeTool.hpp \
  moab/PartitionerBase.hpp \
  moab/ProgOptions.hpp \
  moab/Range.hpp \
  moab/RangeMap.hpp \
  moab/ReadUtilIface.hpp \
  moab/ReaderIface.hpp \
  moab/ReaderWriterSet.hpp \
  moab/ReorderTool.hpp \
  moab/ScdInterface.hpp \
  moab/SetIterator.hpp \
  moab/Skinner.hpp \
  moab/SpatialLocator.hpp \
  moab/SpatialLocatorTimes.hpp \
  moab/SpectralMeshTool.hpp \
  moab/Tree.hpp \
  moab/TreeStats.hpp \
  moab/TupleList.hpp \
  moab/Types.hpp \
  moab/UnknownInterface.hpp \
  moab/Util.hpp \
  moab/verdict.h \
  moab/WriteUtilIface.hpp \
  moab/WriterIface.hpp \
  moab/ElemEvaluator.hpp \
  moab/win32_config.h \
  MBEntityType.h \
  MBCN.h \
  MBCN_protos.h \
  MBTagConventions.hpp 

verdictincludedir = $(includedir)/moab/verdict
verdictinclude_HEADERS = \
  moab/verdict/VerdictWrapper.hpp

discretegeomincludedir = $(includedir)/moab/discrete_geometry
discretegeominclude_HEADERS = \
	moab/discrete_geometry/DGMSolver.hpp \
	moab/discrete_geometry/HiReconstruction.hpp

localdiscincludedir = $(includedir)/moab/local_discretization
localdiscinclude_HEADERS = \
  moab/local_discretization/LinearHex.hpp \
  moab/local_discretization/LinearQuad.hpp \
  moab/local_discretization/LinearTet.hpp \
  moab/local_discretization/LinearTri.hpp \
  moab/local_discretization/QuadraticHex.hpp

pointlocincludedir = $(includedir)/moab/point_locater
pointlocinclude_HEADERS = \
  moab/point_locater/io.hpp  \
  moab/point_locater/point_locater.hpp \
  moab/point_locater/parametrizer.hpp
pointloctreeincludedir = $(pointlocincludedir)/tree
pointloctreeinclude_HEADERS = \
  moab/point_locater/tree/common_tree.hpp \
  moab/point_locater/tree/element_tree.hpp \
  moab/point_locater/tree/bvh_tree.hpp \
  moab/point_locater/tree/initialize_tree.hpp

pointlocelementmapsincludedir = $(pointlocincludedir)/element_maps
pointlocelementmapsinclude_HEADERS = \
  moab/point_locater/element_maps/linear_hex_map.hpp  \
  moab/point_locater/element_maps/linear_tet_map.hpp  \
  moab/point_locater/element_maps/spectral_hex_map.hpp  \
  moab/point_locater/element_maps/quadratic_hex_map.hpp  

intxmeshincludedir = $(includedir)/moab/earthsystem/intx_mesh
intxmeshinclude_HEADERS = \
	moab/earthsystem/intx_mesh/Intx2Mesh.hpp \
	moab/earthsystem/intx_mesh/IntxUtils.hpp \
  moab/earthsystem/intx_mesh/Intx2MeshInPlane.hpp \
  moab/earthsystem/intx_mesh/Intx2MeshOnSphere.hpp \
  moab/earthsystem/intx_mesh/Intx2MeshEdges.hpp \
  moab/earthsystem/intx_mesh/IntxRllCssphere.hpp

# Conditional sources
if ENABLE_FORTRAN
libMOAB_la_SOURCES += iMOABF.F90
imoabincludedir = $(includedir)

if HAVE_LCASE_FCMOD
nodist_imoabinclude_HEADERS = imoab.$(FCMODEXT)
imoab.$(FCMODEXT): moab/iMOAB.h iMOABF.F90
iMOAB.o: imoab.$(FCMODEXT)
else
nodist_imoabinclude_HEADERS = IMOAB.$(FCMODEXT)
IMOAB.$(FCMODEXT): moab/iMOAB.h iMOABF.F90
iMOAB.o: IMOAB.$(FCMODEXT)
endif

endif # ENABLE_FORTRAN

if ENABLE_mbtempest
MOAB_REMAPPING_SRCS =

MOAB_REMAPPING_HDRS = moab/earthsystem/remapping/Remapper.hpp

AM_CPPFLAGS += -I$(srcdir)/earthsystem/remapping

EXTRA_DIST += earthsystem/remapping/netcdfcpp_par.hpp

if HAVE_TEMPESTREMAP
  AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
  libMOAB_la_LIBADD += $(TEMPESTREMAP_LIBS)

  MOAB_REMAPPING_SRCS += \
     earthsystem/remapping/ApplyWeights.cpp \
     earthsystem/remapping/TempestRemapper.cpp  \
     earthsystem/remapping/TempestLinearRemap.cpp \
     earthsystem/remapping/TempestOnlineMap.cpp \
     earthsystem/remapping/TempestOnlineMapIO.cpp \
     earthsystem/remapping/IntegerReprosum.cpp \
     earthsystem/remapping/netcdfcpp_par.hpp

  MOAB_REMAPPING_HDRS += \
     moab/earthsystem/remapping/TempestRemapper.hpp \
     moab/earthsystem/remapping/TempestOnlineMap.hpp \
     moab/earthsystem/remapping/IntegerReprosum.hpp

endif # if HAVE_TEMPESTREMAP

# Add the remapping sources to the library
libMOAB_la_SOURCES += $(MOAB_REMAPPING_SRCS)
# Add the include headers to installation paths
remappingincludedir = $(includedir)/moab/earthsystem/remapping
remappinginclude_HEADERS = $(MOAB_REMAPPING_HDRS)

endif # if ENABLE_mbtempest

if HAVE_HDF5
  nobase_libMOAB_la_include_HEADERS += moab/mhdf_public.h
endif

if PARALLEL
  libMOAB_la_SOURCES += \
     parallel/ParallelComm.cpp \
     parallel/ProcConfig.cpp \
     parallel/ParallelData.cpp \
     parallel/ParallelMergeMesh.cpp \
     parallel/ReadParallel.cpp \
     parallel/ReadParallel.hpp \
     parallel/SharedSetData.cpp \
     parallel/SharedSetData.hpp \
     parallel/ParCommGraph.cpp \
     parallel/gs.cpp

if HAVE_HDF5_PARALLEL
  libMOAB_la_SOURCES += parallel/WriteHDF5Parallel.cpp \
                        parallel/WriteHDF5Parallel.hpp
endif

mpiparallelincludedir = $(includedir)
mpiparallelinclude_HEADERS = \
     parallel/moab_mpi.h \
     parallel/moab_mpi_config.h \
     parallel/MBParallelConventions.h

moabparallelincludedir = $(includedir)/moab
moabparallelinclude_HEADERS = \
     parallel/moab/gs.hpp \
     parallel/moab/ParallelComm.hpp \
     parallel/moab/ParallelMergeMesh.hpp \
     parallel/moab/ProcConfig.hpp \
     parallel/moab/ParallelData.hpp \
     parallel/moab/ParCommGraph.hpp

AM_CPPFLAGS += -I$(srcdir)/parallel

if HAVE_HYPRE

AM_CPPFLAGS += $(HYPRE_INCLUDES) -I$(srcdir)/hypre/moab
AM_LDFLAGS += $(HYPRE_LIBS)

MOAB_HYPRE_SRCS = \
		hypre/HypreParVector.cpp  \
		hypre/hypre_parcsr.cpp \
		hypre/HypreParMatrix.cpp \
		hypre/HypreSolver.cpp
libMOAB_la_SOURCES += $(MOAB_HYPRE_SRCS)

moabhypreincludedir = $(includedir)/moab/hypre
moabhypreinclude_HEADERS = \
		hypre/moab/HypreParVector.hpp \
		hypre/moab/hypre_parcsr.hpp \
		hypre/moab/HypreParMatrix.hpp \
		hypre/moab/HypreSolver.hpp

endif # if HAVE_HYPRE

endif

if HAVE_ZOLTAN
  libMOAB_la_SOURCES += ZoltanPartitioner.cpp
  nobase_libMOAB_la_include_HEADERS += moab/ZoltanPartitioner.hpp
  AM_CPPFLAGS += $(ZOLTAN_INC_FLAGS)
  # not sure yet if this is needed; will see when linking io reader (MPAS so far)
  # the tool mbpart needed it, and it just includes ZoltanPartitioner.hpp;
  # so it needs ways to find the libraries zoltan depends on
  libMOAB_la_LIBADD += $(ZOLTAN_LIBS) $(ZOLTAN_LIB_FLAGS)
endif

if ENABLE_parmetis
  AM_CPPFLAGS += $(PARMETIS_INCLUDES)
  # not sure yet if this is needed; will see when linking io reader (MPAS so far)
  # the tool mbpart needed it, and it just includes MetisPartitioner.hpp;
  # so it needs ways to find the libraries zoltan depends on
  libMOAB_la_LIBADD += $(PARMETIS_LIBS)
endif

if ENABLE_metis
  libMOAB_la_SOURCES += MetisPartitioner.cpp
  nobase_libMOAB_la_include_HEADERS += moab/MetisPartitioner.hpp
  AM_CPPFLAGS += $(METIS_INCLUDES)
  # not sure yet if this is needed; will see when linking io reader (MPAS so far)
  # the tool mbpart needed it, and it just includes MetisPartitioner.hpp;
  # so it needs ways to find the libraries zoltan depends on
  libMOAB_la_LIBADD += $(METIS_LIBS)
endif


CLEANFILES = a.out imoab.mod
EXTRA_DIST += CMakeLists.txt parallel/CMakeLists.txt earthsystem/intx_mesh/CMakeLists.txt

## Rules for Fortran sources compilation
FORTDEFS = $(FCDEFS)
PPFCCOMPILE = $(FC) $(FORTDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
    $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
    $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
    $(AM_FCFLAGS) $(FCFLAGS)
PPF77COMPILE = $(F77) $(FORTDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
   $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FFLAGS) $(FFLAGS)
LTPPF77COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=F77 $(AM_LIBTOOLFLAGS) \
   $(LIBTOOLFLAGS) --mode=compile $(F77) $(FCDEFS) \
   $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
   $(AM_FFLAGS) $(FFLAGS)
