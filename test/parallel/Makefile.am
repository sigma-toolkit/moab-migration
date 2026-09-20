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

AM_CPPFLAGS  = -I. \
               -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/earthsystem/intx_mesh \
               -I$(top_srcdir)/src/io/mhdf/include

# initialize
HDF5_TESTS =
NETCDF_TESTS =

if ENABLE_mbcoupler
  AM_CPPFLAGS += -I$(top_srcdir)/tools/mbcoupler
endif

# Run parallel tests in parallel
if USE_MPIEXEC
  LOG_COMPILER = ${MPIEXEC}
  AM_LOG_FLAGS = ${MPIEXEC_NP} ${NP}
endif

TESTS = pcomm_unit \
        scdtest \
        pcomm_serial \
        par_spatial_locator_test \
        parallel_unit_tests \
        uber_parallel_test \
        hireconst_test_parallel \
        $(NETCDF_TESTS) \
        $(HDF5_TESTS) \
        $(IMOAB_TESTS)

if HAVE_HDF5_PARALLEL
  HDF5_TESTS += parallel_hdf5_test mhdf_parallel parallel_write_test \
         parallel_adj parmerge_test augment_with_ghosts uref_hirec_test_parallel migrate_test \
         ghost_thin_layers

if HAVE_ZOLTAN
  HDF5_TESTS += migrate_nontrivial
endif

if HAVE_TEMPESTREMAP
AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
migrate_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
migrate_nontrivial_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
migrate_test_LDADD = $(LDADD) $(FCLIBS)
migrate_nontrivial_LDADD = $(LDADD) $(FCLIBS)
endif

# The coupler drivers, commgraph_test, par_intx_sph and test_mapweights_bcast
# moved to test/earthsystem/parallel.

endif # HAVE_HDF5_PARALLEL

# iMOAB BfB/tolerance regression harness — run all or one config.
# Uses the configure-generated wrapper which has build dir, mesh dir, and
# MPI launcher baked in.  Requires TempestRemap (the binaries won't exist
# otherwise), but .PHONY and the pattern rule must live outside any
# Automake conditional block.
#
# Examples:
#   make imoabcheck                                    — run all configs
#   make imoabcheck-coupler_bilin                      — run one config
#   make imoabcheck IMOAB_TEST_FLAGS='--ranks "1 2 4"' — custom rank sweep
#   make imoabcheck-dualmap_caas IMOAB_TEST_FLAGS='--workdir /tmp/run --keep'
#
# IMOAB_TEST_FLAGS is forwarded verbatim to run_imoab_tests.sh, which in turn
# forwards every recognised flag (--ranks, --workdir, --keep, --mbcmp, --tol,
# --extra) to the per-config dispatcher. See run_imoab_tests.sh -h for the
# full flag list.
.PHONY: imoabcheck

imoabcheck:
	$(SHELL) $(top_builddir)/test/parallel/regression/run_imoab_tests.sh $(IMOAB_TEST_FLAGS)

imoabcheck-%:
	$(SHELL) $(top_builddir)/test/parallel/regression/run_imoab_tests.sh $(IMOAB_TEST_FLAGS) $*

# The parallel Earth system format readers - read_nc_par, write_nc_par,
# ucdtrvpart, mpastrvpart, gcrm_par and read_esmf - moved to
# test/earthsystem/parallel.  scdpart is a structured-mesh partition test and
# stays here.
if HAVE_NETCDF
  NETCDF_TESTS += scdpart
endif

if HAVE_PNETCDF
if !HAVE_NETCDF
  NETCDF_TESTS += scdpart
endif
endif

if HAVE_HYPRE
  MBHYPRE_TESTS = hypre_test
else
  MBHYPRE_TESTS =
endif

IMOAB_TESTS =

if HAVE_HYPRE
AM_CPPFLAGS += $(HYPRE_INCLUDES) -I$(top_srcdir)/src/hypre/moab
AM_LDFLAGS += $(HYPRE_LIBS)

hypre_test_SOURCES = hypre_test.cpp
TESTS += hypre_test
endif

check_PROGRAMS = $(TESTS) mbparallelcomm_test partcheck structured3
if HAVE_HDF5_PARALLEL
  check_PROGRAMS += parmerge
endif


pcomm_unit_SOURCES = pcomm_unit.cpp
parallel_hdf5_test_SOURCES = parallel_hdf5_test.cpp
mhdf_parallel_SOURCES = mhdf_parallel.c
mhdf_parallel_LDADD = $(LDADD) $(HDF5_LIBS)
parallel_unit_tests_SOURCES = parallel_unit_tests.cpp
parallel_write_test_SOURCES = parallel_write_test.cpp
uber_parallel_test_SOURCES = uber_parallel_test.cpp
pcomm_serial_SOURCES = pcomm_serial.cpp
mbparallelcomm_test_SOURCES = mbparallelcomm_test.cpp
scdtest_SOURCES = scdtest.cpp
partcheck_SOURCES = partcheck.cpp
structured3_SOURCES = structured3.cpp
parmerge_SOURCES = parmerge.cpp
scdpart_SOURCES = scdpart.cpp
par_spatial_locator_test_SOURCES = par_spatial_locator_test.cpp
parallel_adj_SOURCES = ../adj_moab_test.cpp
parmerge_test_SOURCES = parmerge_test.cpp
augment_with_ghosts_SOURCES = augment_with_ghosts.cpp
hireconst_test_parallel_SOURCES = hireconst_test_parallel.cpp
uref_hirec_test_parallel_SOURCES = uref_hirec_test_parallel.cpp ../geomObject.cpp
migrate_test_SOURCES = migrate_test.cpp
migrate_nontrivial_SOURCES = migrate_nontrivial.cpp
ghost_thin_layers_SOURCES = ghost_thin_layers.cpp


# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = *.h5m *.nc *.vtk *.txt

# iMOABCouplerUtils.hpp and iMOABDualMapCAASBfb.sh moved with the coupler
# drivers. The regression harness stays here and drives the renamed C++ tests.
EXTRA_DIST = regression/imoab_test_runner.sh \
             regression/lib_common.sh \
             regression/lib_bfb_harness.sh \
             regression/lib_tol_harness.sh \
             regression/lib_exitcode_harness.sh \
             regression/README.md \
             regression/configs/coupler.cfg \
             regression/configs/coupler_bfb.cfg \
             regression/configs/coupler_bilin.cfg \
             regression/configs/coupler_twohop.cfg \
             regression/configs/domain_culling.cfg \
             regression/configs/dualmap_caas.cfg \
             regression/configs/map_r2l.cfg \
             regression/configs/map_target.cfg \
             regression/configs/phatm_ocn_coupler.cfg \
             regression/configs/ptest2.cfg \
             regression/configs/read_compute_map.cfg \
             regression/configs/scrip.cfg \
             regression/configs/setwithgid.cfg \
             regression/configs/spmv_bfb.cfg \
             regression/run_imoab_tests.sh.in
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

