# The values below are for a copy of MOAB used directly
# from its build directory. These values will be overridden below
# for installed copies of MOAB.

# Feature list.  These are consumed by the example makefiles, which compare them
# against the literal string "yes", so they are substituted from the yes/no
# variables CMakeLists.txt derives just before configuring this file - not from
# the ON/OFF values of MOAB_HAVE_<PKG> directly.
MOAB_MPI_ENABLED = @MOAB_MAKE_MPI_ENABLED@
MOAB_FORTRAN_ENABLED = @MOAB_MAKE_FORTRAN_ENABLED@
# The iTAPS/iMesh bindings are not part of this source tree - there is no
# itaps/ directory and no --enable-imesh - so this is fixed at "no" rather than
# substituted.  It still has to be *defined*: examples/fortran/makefile tests
# it, and against an undefined variable that test silently compared against the
# empty string.  Substitute it properly if iMesh support ever returns.
MOAB_IMESH_ENABLED = no
MOAB_HDF5_ENABLED = @MOAB_MAKE_HDF5_ENABLED@
MOAB_NETCDF_ENABLED = @MOAB_MAKE_NETCDF_ENABLED@
MOAB_PNETCDF_ENABLED = @MOAB_MAKE_PNETCDF_ENABLED@
MOAB_TEMPESTREMAP_ENABLED = @MOAB_MAKE_TEMPESTREMAP_ENABLED@
MOAB_METIS_ENABLED = @MOAB_MAKE_METIS_ENABLED@
MOAB_PARMETIS_ENABLED = @MOAB_MAKE_PARMETIS_ENABLED@
MOAB_ZOLTAN_ENABLED = @MOAB_MAKE_ZOLTAN_ENABLED@
MOAB_EIGEN3_ENABLED = @MOAB_MAKE_EIGEN3_ENABLED@

# Library and Include paths
MOAB_LIBDIR = @abs_builddir@/lib
MOAB_INCLUDES = -I@abs_srcdir@/src \
                -I@abs_builddir@/src \
                -I@abs_srcdir@/src/oldinc \
                -I@abs_srcdir@/src/verdict \
                -I@abs_srcdir@/src/parallel \
                -I@abs_builddir@/src/parallel \
                -I@abs_srcdir@/src/local_discretization \
                -I@abs_srcdir@/src/RefineMesh

MOAB_INCLUDES += @MOAB_MAKE_TPL_INCLUDES@

MOAB_CPPFLAGS = @CPPFLAGS@ 
MOAB_CXXFLAGS = @CXXFLAGS@ 
MOAB_CFLAGS = @CFLAGS@ 
MOAB_FFLAGS = @FFLAGS@
MOAB_FCFLAGS = @FCFLAGS@
MOAB_LDFLAGS = @EXPORT_LDFLAGS@ @CXX_LDFLAGS@ @LDFLAGS@

# missing support for DAMSEL, CCMIO
# Assembled in CMakeLists.txt: the individual @..._LIBRARIES@ are CMake lists and
# cannot be substituted directly, or their separating semicolons end up in the
# makefile.
MOAB_EXT_LIBS = @MOAB_MAKE_EXT_LIBS@
# The C++ runtime, for the Fortran examples: libMOAB is C++ but they are linked
# by the Fortran driver, which does not pull it in on its own.
MOAB_CXX_RUNTIME_LIBS = @MOAB_MAKE_CXX_RUNTIME@
# MPI's Fortran bindings, likewise: MOAB_EXT_LIBS carries the C/C++ MPI libraries
# only, and those do not define the mpi_*_ symbols an F90 unit references.
MOAB_FC_MPI_LIBS = @MOAB_MAKE_FC_MPI_LIBS@
MOAB_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB $(MOAB_EXT_LIBS)

MOAB_CXX = @CXX@
MOAB_CC  = @CC@
MOAB_FC  = @FC@
MOAB_F77  = @F77@

# Override MOAB_LIBDIR and MOAB_INCLUDES from above with the correct
# values for the installed MOAB.

# NEED TO ADD SOMETHING TO MODIFY THIS FILE AT INSTALL TIME (OR FIX ISSUE #30)
