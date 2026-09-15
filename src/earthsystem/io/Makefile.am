# Climate NetCDF readers/writers (CAM Euler/FV/HOMME, MPAS, GCRM, SCRIP,
# ESMF and domain grids).  These are the climate-specific half of the old
# src/io directory; the ExodusII NetCDF reader/writer (ReadNCDF/WriteNCDF)
# is general-purpose and stays in src/io.

AUTOMAKE_OPTIONS = foreign subdir-objects

noinst_LTLIBRARIES = libmoabclimateio.la
libmoabclimateio_la_LIBADD = $(NETCDF_LIBS) $(PNETCDF_LIBS) $(HDF5_LIBS)

AM_CPPFLAGS  = -DIS_BUILDING_MB \
               -DLOCDIR=src/climate/io/ \
               -I../.. -I$(srcdir)/../.. \
               -I$(srcdir)/../../moab \
               -I$(srcdir)/../../io \
               -I$(srcdir)/../../parallel \
               $(EIGEN3_INCLUDES)

if HAVE_HDF5
  AM_CPPFLAGS += $(HDF5_CPPFLAGS)
endif

if ENABLE_mbpart
  AM_CPPFLAGS += $(ZOLTAN_INC_FLAGS)
endif

# The climate NetCDF sources build when NetCDF is available, and also when
# only PNetCDF is available -- mirroring the gating these files had while
# they lived in src/io/Makefile.am.
MOAB_CLIMATE_NETCDF_SRCS = ReadNC.cpp ReadNC.hpp \
                           WriteNC.cpp WriteNC.hpp \
                           NCHelper.cpp NCHelper.hpp \
                           NCHelperEuler.cpp NCHelperEuler.hpp \
                           NCHelperFV.cpp NCHelperFV.hpp \
                           NCHelperDomain.cpp NCHelperDomain.hpp \
                           NCHelperESMF.cpp NCHelperESMF.hpp \
                           NCHelperScrip.hpp NCHelperScrip.cpp \
                           NCHelperHOMME.cpp NCHelperHOMME.hpp \
                           NCHelperMPAS.cpp NCHelperMPAS.hpp \
                           NCHelperGCRM.cpp NCHelperGCRM.hpp \
                           NCWriteHelper.cpp NCWriteHelper.hpp \
                           NCWriteEuler.cpp NCWriteEuler.hpp \
                           NCWriteFV.cpp NCWriteFV.hpp \
                           NCWriteHOMME.cpp NCWriteHOMME.hpp \
                           NCWriteMPAS.cpp NCWriteMPAS.hpp \
                           NCWriteGCRM.cpp NCWriteGCRM.hpp

if HAVE_CLIMATE_NETCDF
  libmoabclimateio_la_SOURCES = $(MOAB_CLIMATE_NETCDF_SRCS)
else
# No NetCDF/PNetCDF: build an empty convenience library so that the link in
# src/Makefile.am stays valid.  A dummy translation unit keeps automake from
# having to infer a linker language from an empty source list.
  libmoabclimateio_la_SOURCES = climate_io_dummy.cpp
endif

EXTRA_DIST = CMakeLists.txt
DISTCLEANFILES = a.out
