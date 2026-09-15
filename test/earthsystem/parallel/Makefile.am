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
  IMOAB_TESTS += imoab_setwithGid

if HAVE_TEMPESTREMAP
AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
commgraph_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
commgraph_test_LDADD = $(LDADD) $(FCLIBS)
endif

if HAVE_TEMPESTREMAP
  HDF5_TESTS += imoab_coupler imoab_coupler_twohop imoab_phatm_ocn_coupler test_mapweights_bcast imoab_scrip imoab_dualmap_caas imoab_spmv_bfb
  imoab_coupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_coupler_twohop_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_phatm_ocn_coupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_map_target_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_domain_culling_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_read_compute_map_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_dualmap_caas_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_spmv_bfb_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
if HAVE_NETCDF
  HDF5_TESTS += imoab_map_target imoab_coupler_bilin imoab_domain_culling imoab_read_compute_map
endif

if HAVE_PNETCDF
  NETCDF_TESTS += imoab_map_r2l
  imoab_map_r2l_SOURCES = imoab_map_r2l.cpp
endif

if HAVE_EXTMESHES
  HDF5_TESTS +=  imoab_apg2_ol_coupler
  DEFS += -DEXT_DATA_DIR="$(builddir)/../../../sandbox/MeshFiles/"
  imoab_apg2_ol_coupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
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
if HAVE_TEMPESTREMAP
IMOAB_TESTS += imoab_coupler_fortran imoab_coupler_caas
imoab_coupler_fortran_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
imoab_coupler_caas_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
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

imoab_coupler_SOURCES = imoab_coupler.cpp
imoab_scrip_SOURCES = imoab_scrip.cpp
imoab_coupler_twohop_SOURCES = imoab_coupler_twohop.cpp
imoab_phatm_ocn_coupler_SOURCES = imoab_phatm_ocn_coupler.cpp
imoab_apg2_ol_coupler_SOURCES = imoab_apg2_ol_coupler.cpp
imoab_map_target_SOURCES = imoab_map_target.cpp
imoab_read_compute_map_SOURCES = imoab_read_compute_map.cpp
imoab_setwithGid_SOURCES = imoab_setwithGid.cpp
imoab_coupler_bilin_SOURCES = imoab_coupler_bilin.cpp
imoab_domain_culling_SOURCES = imoab_domain_culling.cpp
imoab_dualmap_caas_SOURCES = imoab_dualmap_caas.cpp
imoab_spmv_bfb_SOURCES = imoab_spmv_bfb.cpp
imoab_coupler_fortran_SOURCES = imoab_coupler_fortran.F90
imoab_coupler_caas_SOURCES = imoab_coupler_caas.F90

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = *.h5m *.nc *.vtk *.txt

# The regression harness stays in test/parallel: it also drives the generic
# imoab_ptest2, and test/parallel/Makefile.am distributes it.
EXTRA_DIST = imoab_coupler_utils.hpp imoab_dualmap_caas_bfb.sh

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
