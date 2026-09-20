# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

LDADD = $(top_builddir)/src/libMOAB.la

# The Intel Fortran runtime links in LLVM MSan when compiled with -C (runtime
# checks). Since third-party libraries (MPI, HDF5, hwloc) are not MSan-
# instrumented, all reads from their output produce false positives. Allow
# the test to report warnings but complete successfully.
AM_TESTS_ENVIRONMENT = MSAN_OPTIONS="halt_on_error=0:exit_code=0"

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest
RMESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

DEFS += -DIS_BUILDING_MB \
        -DSRCDIR=$(srcdir) \
        -DMESHDIR=$(MESHDIR)

FCDEFS += $(FC_DEFINE)IS_BUILDING_MB

# -I$(srcdir)/../.. is for TestUtil.hpp in test/.
AM_CPPFLAGS  = -I. \
               -I../.. -I$(srcdir)/../.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/earthsystem/intx_mesh \
               -I$(top_srcdir)/src/io/mhdf/include

# initialize
HDF5_TESTS =
NETCDF_TESTS =
IMOAB_TESTS =
MBCSLAM_TESTS =

# Run parallel tests in parallel
if USE_MPIEXEC
  LOG_COMPILER = ${MPIEXEC}
  AM_LOG_FLAGS = ${MPIEXEC_NP} ${NP}
endif

TESTS = $(NETCDF_TESTS) $(HDF5_TESTS) $(MBCSLAM_TESTS) $(IMOAB_TESTS)

if HAVE_HDF5_PARALLEL
  HDF5_TESTS += commgraph_test
  MBCSLAM_TESTS += par_intx_sph
  IMOAB_TESTS += iMOABSetWithGid iMOABPtest2

if HAVE_TEMPESTREMAP
AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
commgraph_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
commgraph_test_LDADD = $(LDADD) $(FCLIBS)
endif

if HAVE_TEMPESTREMAP
  HDF5_TESTS += iMOABCoupler iMOABCouplerTwoHop iMOABPhAtmOcnCoupler test_mapweights_bcast iMOABScrip iMOABDualMapCAAS iMOABSpMVBFB
  iMOABCoupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABCouplerTwoHop_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABPhAtmOcnCoupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABMapTarget_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABDomainCulling_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABReadComputeMap_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABDualMapCAAS_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  iMOABSpMVBFB_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
if HAVE_NETCDF
  HDF5_TESTS += iMOABMapTarget iMOABCouplerBilin iMOABDomainCulling iMOABReadComputeMap
endif

if HAVE_PNETCDF
  NETCDF_TESTS += iMOABMapR2L
endif

if HAVE_EXTMESHES
  HDF5_TESTS += iMOABApg2OlCoupler
  DEFS += -DEXT_DATA_DIR="$(builddir)/../../../sandbox/MeshFiles/"
  iMOABApg2OlCoupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
endif

  test_mapweights_bcast_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
endif # tempestremap

endif # HAVE_HDF5_PARALLEL

# Parallel readers for the Earth system grid formats.  read_nc_par and
# write_nc_par build the serial sources next door in earthsystem/io with an
# MPI-enabled MOAB; that is how they were configured in test/parallel before
# this subtree existed.
if HAVE_NETCDF
  NETCDF_TESTS += read_nc_par ucdtrvpart mpastrvpart gcrm_par read_esmf
endif

if HAVE_PNETCDF
if !HAVE_NETCDF
  NETCDF_TESTS += read_nc_par ucdtrvpart mpastrvpart write_nc_par read_esmf
else
  NETCDF_TESTS += write_nc_par
endif
endif

if HAVE_PNETCDF
  write_nc_par_LDADD = $(LDADD) $(PNETCDF_LIBS)
endif

if ENABLE_FORTRAN
  IMOAB_TESTS += iMOABPtest
  iMOABPtest_LDADD = $(LDADD) $(FCLIBS)
if HAVE_TEMPESTREMAP
  IMOAB_TESTS += iMOABCouplerFortran iMOABCouplerCAAS
  iMOABCouplerFortran_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
  iMOABCouplerCAAS_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
endif
endif

check_PROGRAMS = $(TESTS)

read_nc_par_SOURCES = ../io/read_nc.cpp
write_nc_par_SOURCES = ../io/write_nc.cpp
ucdtrvpart_SOURCES = ucdtrvpart.cpp
mpastrvpart_SOURCES = mpastrvpart.cpp
read_esmf_SOURCES = read_esmf.cpp
gcrm_par_SOURCES = gcrm_par.cpp
commgraph_test_SOURCES = commgraph_test.cpp
test_mapweights_bcast_SOURCES = test_mapweights_bcast.cpp
par_intx_sph_SOURCES = par_intx_sph.cpp
par_intx_sph_LDADD = $(LDADD)

iMOABCoupler_SOURCES = iMOABCoupler.cpp
iMOABScrip_SOURCES = iMOABScrip.cpp
iMOABCouplerTwoHop_SOURCES = iMOABCouplerTwoHop.cpp
iMOABPhAtmOcnCoupler_SOURCES = iMOABPhAtmOcnCoupler.cpp
iMOABApg2OlCoupler_SOURCES = iMOABApg2OlCoupler.cpp
iMOABMapTarget_SOURCES = iMOABMapTarget.cpp
iMOABReadComputeMap_SOURCES = iMOABReadComputeMap.cpp
iMOABSetWithGid_SOURCES = iMOABSetWithGid.cpp
iMOABPtest2_SOURCES = iMOABPtest2.cpp
iMOABCouplerBilin_SOURCES = iMOABCouplerBilin.cpp
iMOABDomainCulling_SOURCES = iMOABDomainCulling.cpp
iMOABMapR2L_SOURCES = iMOABMapR2L.cpp
iMOABDualMapCAAS_SOURCES = iMOABDualMapCAAS.cpp
iMOABSpMVBFB_SOURCES = iMOABSpMVBFB.cpp
iMOABPtest_SOURCES = iMOABPtest.F
iMOABCouplerFortran_SOURCES = iMOABCouplerFortran.F90
iMOABCouplerCAAS_SOURCES = iMOABCouplerCAAS.F90

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = *.h5m *.nc *.vtk *.txt

# The regression harness stays in test/parallel and drives these renamed tests.
EXTRA_DIST = iMOABCouplerUtils.hpp iMOABDualMapCAASBFB.sh

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

DISTCLEANFILES = a.out
