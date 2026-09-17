# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects no-dependencies

SUBDIRS = . io
if PARALLEL
  SUBDIRS += parallel
endif

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

# The Intel Fortran runtime links in LLVM MSan when compiled with -C (runtime
# checks). Since third-party libraries (MPI, HDF5, hwloc) are not MSan-
# instrumented, all reads from their output produce false positives. Allow
# the test to report warnings but complete successfully.
AM_TESTS_ENVIRONMENT = MSAN_OPTIONS="halt_on_error=0:exit_code=0"

# -I$(srcdir)/.. is for TestUtil.hpp, which stays in test/ because the generic
# tests use it too.
AM_CPPFLAGS  = -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/earthsystem/intx_mesh

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

TESTS = arc_intx_tests \
        intx_in_plane_test \
        intx_on_sphere_test

if HAVE_NETCDF
  TESTS += intx_rll_cs_sphere_test \
           intx_orientation_test
endif

arc_intx_tests_SOURCES = arc_intx_tests.cpp
intx_in_plane_test_SOURCES = intx_in_plane_test.cpp
intx_on_sphere_test_SOURCES = intx_on_sphere_test.cpp
intx_rll_cs_sphere_test_SOURCES = intx_rll_cs_sphere_test.cpp
intx_orientation_test_SOURCES = intx_orientation_test.cpp

# Tests for remapping
if ENABLE_mbtempest
if HAVE_TEMPESTREMAP
  test_remapping_SOURCES = test_remapping.cpp
  test_remapping_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/earthsystem/remapping -I$(top_srcdir)/src/earthsystem/intx_mesh $(TEMPESTREMAP_INCLUDE)
  test_remapping_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_remapping_SOURCES = imoab_remap.cpp
  imoab_remapping_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/earthsystem/remapping -I$(top_srcdir)/src/earthsystem/intx_mesh $(TEMPESTREMAP_INCLUDE)
  imoab_remapping_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  TESTS += test_remapping imoab_remapping
endif
endif

check_PROGRAMS = $(TESTS)

EXTRA_DIST = CMakeLists.txt

CLEANFILES = *.g *.gen *.nc *.h5m *.vtk

DISTCLEANFILES = a.out
