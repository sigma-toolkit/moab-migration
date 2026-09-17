# subdir-objects is required now that mbtempest and h5mtoscrip build from
# earthsystem/; without it automake warns and puts the objects in this
# directory, where two sources of the same basename would collide.
AUTOMAKE_OPTIONS = foreign subdir-objects
AM_CPPFLAGS =  -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               $(EIGEN3_INCLUDES)

LDADD = $(top_builddir)/src/libMOAB.la

if HAVE_TEMPESTREMAP
  LDADD += $(TEMPESTREMAP_LIBS)
endif

# For old (pre 1.10.x) versions of Automake
docdir = @docdir@

bin_PROGRAMS =
dist_man1_MANS =
# Always install the tools doc index + the build-integration overview,
# regardless of which individual tools are enabled below. Per-tool docs
# are appended inside the corresponding ENABLE_<tool> blocks; this
# parallels the unconditional MOAB_DOCS += in tools/CMakeLists.txt.
doc_DATA = docs/README.md docs/tools.md

EXTRA_DIST = docs/README.md \
             docs/tools.md \
             docs/mbconvert.md \
             docs/mbtempest.md \
             docs/mbIntxCheck.md \
             docs/compareMaps.md \
             docs/compareFiles.md \
             docs/measure.md \
             docs/hexmodops.md \
             docs/mbsize.md \
             docs/mbskin.md \
             docs/mbsurfplot.md \
             docs/mbtagprop.md \
             docs/spheredecomp.md \
             CMakeLists.txt \
             earthsystem/CMakeLists.txt

if HAVE_HDF5
  bin_PROGRAMS += mbcmpfiles mbintxchk mbcleantags
  doc_DATA += docs/compareFiles.md docs/mbIntxCheck.md
if HAVE_NETCDF
  bin_PROGRAMS += mbaddnc mbaddchunk mbedgesplit
if HAVE_EIGEN3
  bin_PROGRAMS += mbcmpmaps mbvisumap
  doc_DATA += docs/compareMaps.md
endif
endif
endif

if ENABLE_mbconvert
  mbconvert_SOURCES = convert.cpp
  bin_PROGRAMS += mbconvert
  doc_DATA += docs/mbconvert.md
endif

if ENABLE_hexmodops
  bin_PROGRAMS += hexmodops
  doc_DATA += docs/hexmodops.md
endif

if ENABLE_mbsize
  bin_PROGRAMS += mbsize
  doc_DATA += docs/mbsize.md docs/measure.md
endif

if ENABLE_mbskin
  bin_PROGRAMS += mbskin
  doc_DATA += docs/mbskin.md
endif

if ENABLE_mbdepth
  bin_PROGRAMS += mbdepth
endif

if ENABLE_mbtagprop
  bin_PROGRAMS += mbtagprop
  doc_DATA += docs/mbtagprop.md
endif

if ENABLE_mbmem
  bin_PROGRAMS += mbmem
endif

if ENABLE_spheredecomp
  bin_PROGRAMS += spheredecomp
  doc_DATA += docs/spheredecomp.md
endif

if ENABLE_mbsurfplot
  bin_PROGRAMS += mbsurfplot
  doc_DATA += docs/mbsurfplot.md
endif

if ENABLE_gsets
  bin_PROGRAMS += mbgsets
endif

if ENABLE_mbcoupler
  mbcoupler_DIR = mbcoupler
else
  mbcoupler_DIR=
endif

if ENABLE_mbpart
  bin_PROGRAMS += mbpart
  AM_CPPFLAGS += $(ZOLTAN_INC_FLAGS) 
  AM_CPPFLAGS += $(METIS_INCLUDES) 
  AM_LDFLAGS += $(ZOLTAN_LIB_FLAGS)
  mbpart_SOURCES = mbpart.cpp
  mbpart_LDADD = $(top_builddir)/src/libMOAB.la $(ZOLTAN_LIBS) $(PARMETIS_LIBS) $(METIS_LIBS)
else
if ENABLE_metis
  bin_PROGRAMS += mbpart
  AM_CPPFLAGS += $(METIS_INCLUDES)
  AM_LDFLAGS += $(METIS_LIB_FLAGS)
  mbpart_SOURCES = mbpart.cpp
  mbpart_LDADD = $(top_builddir)/src/libMOAB.la
endif
endif

if ENABLE_mcnpmit
  mcnpmit_DIR = mcnpmit
else
  mcnpmit_DIR =
endif

if PARALLEL
  parread = parread
else
  parread =
endif

if ENABLE_mbquality
  AM_CPPFLAGS += -I$(top_srcdir)/src/verdict
  bin_PROGRAMS += mbquality
mbquality_SOURCES = quality.cpp
endif

if ENABLE_mbumr
  AM_CPPFLAGS += -I$(top_srcdir)/src/RefineMesh
  bin_PROGRAMS += mbumr
mbumr_SOURCES = umr.cpp
endif

if HAVE_HYPRE

AM_CPPFLAGS += $(HYPRE_INCLUDES) -I$(top_srcdir)/src/hypre/moab
AM_LDFLAGS += $(HYPRE_LIBS)

endif

geometry_DIR = geometry

if ENABLE_mbtempest
  AM_CPPFLAGS += $(TEMPESTREMAP_INCLUDE) -I$(top_srcdir)/src/earthsystem/remapping -I$(top_srcdir)/src/earthsystem/intx_mesh
  bin_PROGRAMS += mbtempest h5mtoscrip
mbtempest_SOURCES = earthsystem/mbtempest.cpp
h5mtoscrip_SOURCES = earthsystem/h5mtoscrip.cpp
  doc_DATA += docs/mbtempest.md
endif

if ENABLE_mbmerge
  bin_PROGRAMS += mbmerge
mbmerge_SOURCES = merge.cpp
endif

if ENABLE_mbhonodes
  bin_PROGRAMS += mbhonodes
endif

noinst_PROGRAMS = $(parread)

SUBDIRS = . \
          $(mbcoupler_DIR) \
          $(mcnpmit_DIR) \
          $(geometry_DIR)

mbgsets_SOURCES = gsets.cpp
hexmodops_SOURCES = makeops.cpp
mbsize_SOURCES = measure.cpp measure.hpp size.cpp
mbskin_SOURCES = skin.cpp
mbdepth_SOURCES = depth.cpp
spheredecomp_SOURCES = SphereDecomp.hpp SphereDecomp.cpp sphere_decomp.cpp
mbsurfplot_SOURCES = surfplot.cpp
mbtagprop_SOURCES = parse.cpp parse.hpp propagate_tags.cpp
mbmem_SOURCES = mbmem.cpp
parread_SOURCES = parread.cpp
mbhonodes_SOURCES = honodes.cpp
mbaddnc_SOURCES = addncdata.cpp
mbaddchunk_SOURCES = mbaddchunk.cpp
mbcmpfiles_SOURCES = compareFiles.cpp
mbintxchk_SOURCES = mbIntxCheck.cpp
mbcmpmaps_SOURCES = compareMaps.cpp
mbcleantags_SOURCES = cleanTags.cpp
mbvisumap_SOURCES = visuMap.cpp
mbedgesplit_SOURCES = mbedgesplit.cpp

bin_PROGRAMS += mbslavepart
mbslavepart_SOURCES = mbslavepart.cpp

DISTCLEANFILES = a.out
