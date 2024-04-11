AM_CPPFLAGS = -I$(top_srcdir)/src \
              -I$(top_builddir)/src \
	      -I$(top_srcdir)/tools

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

bin_PROGRAMS = cgm2moab
noinst_PROGRAMS = quads_to_tris

cgm2moab_SOURCES = cgm2moab.cpp obb_analysis.cpp cgm2moab.hpp
cgm2moab_LDADD = $(top_builddir)/src/libMOAB.la

quads_to_tris_SOURCES = quads_to_tris_driver.cpp quads_to_tris.cpp
quads_to_tris_LDADD = $(top_builddir)/src/libMOAB.la

EXTRA_DIST = CMakeLists.txt
DISTCLEANFILES = a.out

