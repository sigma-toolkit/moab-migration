# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

LDADD = $(top_builddir)/src/libMOAB.la

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
               -I$(top_srcdir)/src/IntxMesh \
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
        $(MBCSLAM_TESTS) $(IMOAB_TESTS)

if HAVE_HDF5_PARALLEL
  HDF5_TESTS += parallel_hdf5_test mhdf_parallel parallel_write_test \
         parallel_adj parmerge_test augment_with_ghosts uref_hirec_test_parallel migrate_test \
         ghost_thin_layers commgraph_test

if HAVE_ZOLTAN
  HDF5_TESTS += migrate_nontrivial
endif

if HAVE_TEMPESTREMAP
AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE)
migrate_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
migrate_nontrivial_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
commgraph_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
migrate_test_LDADD = $(LDADD) $(FCLIBS)
migrate_nontrivial_LDADD = $(LDADD) $(FCLIBS)
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
  DEFS += -DEXT_DATA_DIR="$(builddir)/../../sandbox/MeshFiles/"
  imoab_apg2_ol_coupler_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
endif

  test_mapweights_bcast_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
endif # tempestremap

endif # HAVE_HDF5_PARALLEL

if HAVE_NETCDF
  NETCDF_TESTS += scdpart read_nc_par ucdtrvpart mpastrvpart gcrm_par read_esmf
endif

if HAVE_PNETCDF
if !HAVE_NETCDF
  NETCDF_TESTS += scdpart read_nc_par ucdtrvpart mpastrvpart write_nc_par read_esmf
else
  NETCDF_TESTS += write_nc_par
endif
endif

if HAVE_PNETCDF
  write_nc_par_LDADD = $(LDADD) $(PNETCDF_LIBS)
endif

if HAVE_HYPRE
  MBHYPRE_TESTS = hypre_test
else
  MBHYPRE_TESTS =
endif

if HAVE_HDF5_PARALLEL
  MBCSLAM_TESTS = par_intx_sph
else
  MBCSLAM_TESTS =
endif

IMOAB_TESTS =
if ENABLE_FORTRAN
   IMOAB_TESTS += imoab_ptest imoab_ptest2
if HAVE_TEMPESTREMAP
IMOAB_TESTS += imoab_coupler_fortran imoab_coupler_caas
imoab_ptest_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
imoab_ptest2_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
imoab_coupler_fortran_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
imoab_coupler_caas_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
imoab_ptest_LDADD = $(LDADD) $(FCLIBS)
imoab_ptest2_LDADD = $(LDADD) $(FCLIBS)
endif
endif

if HAVE_HDF5_PARALLEL
 IMOAB_TESTS += imoab_setwithGid
endif

if HAVE_HYPRE
AM_CPPFLAGS += $(HYPRE_INCLUDES) -I$(top_srcdir)/src/HYPRE/moab
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
read_nc_par_SOURCES = ../io/read_nc.cpp
ucdtrvpart_SOURCES = ucdtrvpart.cpp
mpastrvpart_SOURCES = mpastrvpart.cpp
read_esmf_SOURCES = read_esmf.cpp
gcrm_par_SOURCES = gcrm_par.cpp
write_nc_par_SOURCES = ../io/write_nc.cpp
par_spatial_locator_test_SOURCES = par_spatial_locator_test.cpp
parallel_adj_SOURCES = ../adj_moab_test.cpp
parmerge_test_SOURCES = parmerge_test.cpp
augment_with_ghosts_SOURCES = augment_with_ghosts.cpp
imoab_ptest_SOURCES = imoab_ptest.F
imoab_coupler_fortran_SOURCES = imoab_coupler_fortran.F90
imoab_coupler_caas_SOURCES = imoab_coupler_caas.F90
imoab_ptest2_SOURCES = imoab_ptest2.cpp
hireconst_test_parallel_SOURCES = hireconst_test_parallel.cpp
uref_hirec_test_parallel_SOURCES = uref_hirec_test_parallel.cpp ../geomObject.cpp
migrate_test_SOURCES = migrate_test.cpp
imoab_coupler_SOURCES = imoab_coupler.cpp
imoab_scrip_SOURCES = imoab_scrip.cpp
imoab_coupler_twohop_SOURCES = imoab_coupler_twohop.cpp
imoab_phatm_ocn_coupler_SOURCES = imoab_phatm_ocn_coupler.cpp
imoab_apg2_ol_coupler_SOURCES = imoab_apg2_ol_coupler.cpp
imoab_map_target_SOURCES = imoab_map_target.cpp
imoab_read_compute_map_SOURCES = imoab_read_compute_map.cpp
imoab_setwithGid_SOURCES = imoab_setwithGid.cpp
migrate_nontrivial_SOURCES = migrate_nontrivial.cpp
ghost_thin_layers_SOURCES = ghost_thin_layers.cpp
commgraph_test_SOURCES = commgraph_test.cpp
test_mapweights_bcast_SOURCES = test_mapweights_bcast.cpp
imoab_coupler_bilin_SOURCES = imoab_coupler_bilin.cpp
imoab_domain_culling_SOURCES = imoab_domain_culling.cpp
imoab_dualmap_caas_SOURCES = imoab_dualmap_caas.cpp
imoab_spmv_bfb_SOURCES = imoab_spmv_bfb.cpp

par_intx_sph_SOURCES = par_intx_sph.cpp
par_intx_sph_LDADD = $(LDADD)

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = *.h5m *.nc *.vtk *.txt

EXTRA_DIST = imoab_coupler_utils.hpp
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

