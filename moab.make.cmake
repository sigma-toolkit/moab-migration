# The values below are for a copy of MOAB used directly
# from its build directory. These values will be overridden below
# for installed copies of MOAB.

# Feature list
MOAB_MPI_ENABLED = @MOAB_HAVE_MPI@
MOAB_FORTRAN_ENABLED = @ENABLE_FORTRAN@
MOAB_HDF5_ENABLED = @MOAB_HAVE_HDF5@
MOAB_NETCDF_ENABLED = @MOAB_HAVE_NETCDF@
MOAB_PNETCDF_ENABLED = @MOAB_HAVE_PNETCDF@
MOAB_TEMPESTREMAP_ENABLED = @MOAB_HAVE_TEMPESTREMAP@
MOAB_METIS_ENABLED = @MOAB_HAVE_METIS@
MOAB_PARMETIS_ENABLED = @MOAB_HAVE_PARMETIS@
MOAB_ZOLTAN_ENABLED = @MOAB_HAVE_ZOLTAN@

# Library and Include paths
MOAB_LIBDIR = @abs_builddir@/lib
MOAB_INCLUDES = -I@abs_srcdir@/src \
                -I@abs_builddir@/src \
                -I@abs_srcdir@/src/oldinc \
                -I@abs_srcdir@/src/verdict \
                -I@abs_srcdir@/src/parallel \
                -I@abs_builddir@/src/parallel \
                -I@abs_srcdir@/src/LocalDiscretization \
                -I@abs_srcdir@/src/RefineMesh

MOAB_INCLUDES += -I@ZOLTAN_INCLUDES@ -I@NETCDF_INCLUDES@ -I@PNETCDF_INCLUDES@ -I@HDF5_INCLUDES@

MOAB_CPPFLAGS = @CPPFLAGS@ 
MOAB_CXXFLAGS = @CXXFLAGS@ 
MOAB_CFLAGS = @CFLAGS@ 
MOAB_FFLAGS = @FFLAGS@
MOAB_FCFLAGS = @FCFLAGS@
MOAB_LDFLAGS = @EXPORT_LDFLAGS@ @CXX_LDFLAGS@ @LDFLAGS@

# missing support for DAMSEL, CCMIO
# PROBABLY FORMATTED INCORRECTLY
MOAB_EXT_LIBS = @NETCDF_LIBRARIES@ @PNETCDF_LIBRARIES@ @CGNS_LIBRARIESS@ @HDF5_LIBRARIES@ \
                @TEMPESTREMAP_LIBRARIES@ @ZOLTAN_LIBRARIES@ @PARMETIS_LIBRARIES@ \
                @METIS_LIBRARIES@ @HYPRE_LIBRARIES@ @LAPACK_LIBRARIES@ @BLAS_LIBRARIES@ @LIBS@ 
MOAB_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB $(MOAB_EXT_LIBS) 
DAGMC_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} @DAGMC_LIBS@ -lMOAB $(MOAB_EXT_LIBS)

MOAB_CXX = @CXX@
MOAB_CC  = @CC@
MOAB_FC  = @FC@
MOAB_F77  = @F77@

# Override MOAB_LIBDIR and MOAB_INCLUDES from above with the correct
# values for the installed MOAB.

# NEED TO ADD SOMETHING TO MODIFY THIS FILE AT INSTALL TIME (OR FIX ISSUE #30)
