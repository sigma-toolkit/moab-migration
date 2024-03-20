SUBDIRS = point_location
MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = -I$(top_srcdir)/src -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -DMESHDIR=$(MESHDIR) \
               -DIS_BUILDING_MB

LDADD = $(top_builddir)/src/libMOAB.la

check_PROGRAMS = perf seqperf adj_time perftool adj_mem_time umr_perf
noinst_PROGRAMS =

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel -I$(top_srcdir)/src/parallel/moab
endif

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
adj_time_SOURCES = adj_time.cpp
perftool_SOURCES = perftool.cpp
adj_mem_time_SOURCES = adj_mem_time_test.cpp
umr_perf_SOURCES = umr_perf.cpp

DISTCLEANFILES = a.out
