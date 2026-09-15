# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign no-dependencies

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

# -I$(srcdir)/../.. is for TestUtil.hpp in test/.
AM_CPPFLAGS  = $(BGDEF)-DSRCDIR=$(srcdir) \
               $(BGDEF)-DMESHDIR=$(MESHDIR) \
               -I$(builddir) -I../.. -I$(srcdir)/../.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/earthsystem/io \
               -I$(top_srcdir)/src/parallel

# These readers and the writer are only built when NetCDF is available, which
# is the same gate they had in test/io/Makefile.am.
if HAVE_NETCDF
  NC_TEST = read_nc read_ucd_nc read_mpas_nc write_nc read_gcrm_nc
  LDADD += $(NETCDF_LIBS) $(PNETCDF_LIBS) $(HDF5_LIBS)
else
  NC_TEST =
endif

TESTS = $(NC_TEST)

check_PROGRAMS = $(TESTS)

read_nc_SOURCES = $(srcdir)/../../TestUtil.hpp read_nc.cpp
read_ucd_nc_SOURCES = $(srcdir)/../../TestUtil.hpp read_ucd_nc.cpp
read_mpas_nc_SOURCES = $(srcdir)/../../TestUtil.hpp read_mpas_nc.cpp
read_gcrm_nc_SOURCES = $(srcdir)/../../TestUtil.hpp read_gcrm_nc.cpp
write_nc_SOURCES = $(srcdir)/../../TestUtil.hpp write_nc.cpp

# No CMakeLists.txt here: the CMake build registers these from
# test/earthsystem/CMakeLists.txt, which distributes itself.

MOSTLYCLEANFILES = *.nc *.h5m *.vtk

DISTCLEANFILES = a.out
