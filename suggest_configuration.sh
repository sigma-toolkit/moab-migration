#!/bin/bash
# Configuration command builder

# Build specifications
CMAKE=yes
ENABLE_DEBUG=no
ENABLE_OPTIMIZE=yes
ENABLE_SHARED=no
ENABLE_STATIC=yes
# Compiler specifications
ENABLE_MPI=yes
ENABLE_FORTRAN=yes

# Installation specifications
PREFIX_INSTALL_PATH=$HOME/install/MOAB

ENABLE_EIGEN3=no
ENABLE_TEMPESTREMAP=NO
#####################
### DO NOT MODIFY ###
#####################
HOSTNAME=`hostname`
INTERNAL_OPTIONS=""
if (test "x$CMAKE" != "xyes"); then
  if (test "x$ENABLE_DEBUG" != "xno"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --enable-debug"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --disable-debug"
  fi
  if (test "x$ENABLE_OPTIMIZE" != "xno"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --enable-optimize"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --disable-optimize"
  fi
  if (test "x$ENABLE_SHARED" != "xno"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --enable-shared"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --disable-shared"
  fi
  if (test "x$ENABLE_STATIC" != "xno"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --enable-static"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --disable-static"
  fi

  if (test "x$PREFIX_INSTALL_PATH" != "x"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --prefix=$PREFIX_INSTALL_PATH"
  fi
else
  if (test "x$ENABLE_OPTIMIZE" != "xyes"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_BUILD_TYPE=Debug"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_BUILD_TYPE=RelWithDebInfo"
  fi
  if (test "x$ENABLE_SHARED" != "xno" || test "x$ENABLE_STATIC" != "xyes"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DBUILD_SHARED_LIBS=ON"
  else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DBUILD_SHARED_LIBS=OFF"
  fi

  if (test "x$PREFIX_INSTALL_PATH" != "x"); then
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_INSTALL_PREFIX=$PREFIX_INSTALL_PATH"
  fi
fi

# PRESET CONFIGURATION COMMANDS
CONFIGURE_CMD="$PWD/configure"

# Variables to set
MBHOSTSYS=""
MBCC=""
MBCXX=""
MBFC=""
MBF77=""
test "x`which mpicc 2>&-`" != "x" && MBCC=`which mpicc 2>&-`
test "x$MBCC" != "x" || MBCC=`which cc 2>&-`
test "x$MBCC" != "x" || MBCC=`which gcc 2>&-`
test "x$MBCC" != "x" || MBCC=`which icc 2>&-`
test "x$MBCC" != "x" || MBCC=`which clang 2>&-`
test "x`which mpicxx 2>&-`" != "x" && MBCXX=`which mpicxx 2>&-`
test "x$MBCXX" != "x" || MBCXX=`which CC 2>&-`
test "x$MBCXX" != "x" || MBCXX=`which g++ 2>&-`
test "x$MBCXX" != "x" || MBCXX=`which icpc 2>&-`
test "x$MBCXX" != "x" || MBCXX=`which clang++ 2>&-`
test "`which mpif90 2>&-`" != "x" && MBFC=`which mpif90 2>&-`
test "x$MBFC" != "x" || MBFC=`which ftn 2>&-`
test "x$MBFC" != "x" || MBFC=`which gfortran 2>&-`
test "x$MBFC" != "x" || MBFC=`which ifort 2>&-`
test "x`which mpif77 2>&-`" != "x" && MBF77=`which mpif77 2>&-`
test "x$MBF77" != "x" || MBF77=`which ftn 2>&-`
test "x$MBF77" != "x" || MBF77=`which gfortran 2>&-`
test "x$MBF77" != "x" || MBF77=`which ifort 2>&-`
MBMPI_DIR=""
MBHDF5_DIR=""
MBSZIP_DIR=""
MBZLIB_DIR=""
MBNETCDF_DIR=""
MBPNETCDF_DIR=""
MBMETIS_DIR=""
MBPARMETIS_DIR=""
MBZOLTAN_DIR=""
MBSCOTCH_DIR=""
MBPTSCOTCH_DIR=""
MBVTK_DIR=""
MBCGM_DIR=""
CROSSCOMPILE="no"

MBDWLD_OPTIONS=""
case "$HOSTNAME" in
  *vesta* | *mira*)
    MBCC="mpixlc_r"
    MBCXX="mpixlcxx_r"
    MBFC="mpixlf90_r"
    MBF77="mpixlf77_r"
    MBNMPICC="xlc_r"
    MBNMPICXX="xlc++_r"
    MBNMPIFC="xlf_r"
    MBNMPIF77="xlf_r"
    MBHOSTSYS="powerpc64-bgq-linux"
    MBHDF5_DIR="/soft/libraries/hdf5/1.8.17/cnk-xl/current"
    MBZLIB_DIR="/soft/libraries/alcf/current/xl/ZLIB"
    MBPNETCDF_DIR="/soft/libraries/pnetcdf/current/cnk-xl/current"
    MBNETCDF_DIR="/soft/libraries/netcdf/4.3.3-f4.4.1/cnk-xl/current"
    MBPARMETIS_DIR="/soft/libraries/alcf/current/xl/PARMETIS"
    MBMETIS_DIR="/soft/libraries/alcf/current/xl/METIS"
    MBDWLD_OPTIONS="--with-zlib=/soft/libraries/alcf/current/xl/ZLIB"
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS --enable-static --enable-all-static"
    PREREQ="soft add +mpiwrapper-xl.legacy.ndebug"
    LIBS="-L/soft/compilers/ibmcmp-may2016/xlmass/bg/7.3/bglib64 -lmassv -lmass -L/soft/libraries/alcf/current/xl/LAPACK/lib -llapack -L/soft/libraries/essl/current/essl/5.1/lib64 -lesslsmpbg -L/soft/compilers/ibmcmp-may2016/xlf/bg/14.1/bglib64 -lxlf90_r -L/soft/compilers/ibmcmp-may2016/xlsmp/bg/3.1/bglib64 -lxlsmp -lxlopt -lxlfmath -lxl -Wl,--allow-multiple-definition"
    ;;
  *blogin*)
    MBNMPICC="icc"
    MBNMPICXX="icpc"
    MBNMPIFC="ifort"
    MBNMPIF77="ifort"
    MBMPI_DIR="/software/mvapich2-intel-psm-1.9.5"
    MBHDF5_DIR="/soft/hdf5/1.8.12-parallel/intel-13.1/mvapich2-1.9"
    MBSZIP_DIR="/soft/szip/2.1/intel-13.1"
    MBNETCDF_DIR="/soft/netcdf/4.3.1-parallel/intel-13.1/mvapich2-1.9"
    MBPNETCDF_DIR="/soft/pnetcdf/1.6.1-gnu4.4-mvapich2"
    MBMETIS_DIR="/soft/metis/5.0.3"
    PREREQ="soft add +intel-15.0 +mvapich2-2.2b-intel-15.0"
    ;;
  *theta*)
    CROSSCOMPILE="yes"
    MBCONFARCH="haswell"
    MBCC="cc"
    MBCXX="CC"
    MBFC="ftn"
    MBF77="ftn"
    MBNMPICC="cc"
    MBNMPICXX="CC"
    MBNMPIFC="ftn"
    MBNMPIF77="ftn"
    MBMPI_DIR=""
    MBHDF5_DIR="/opt/cray/pe/hdf5-parallel/default/intel/51"
    MBNETCDF_DIR="/opt/cray/pe/netcdf/default/intel/51"
    MBMETIS_DIR="/opt/cray/pe/tpsl/default/intel/51/$MBCONFARCH"
    LIBS="-L/opt/cray/pe/libsci/default/intel/51/$MBCONFARCH/lib -lsci_gnu_mp"
    PREREQ="module load gcc/6.3.0 craype-haswell"
    ;;
  *cori*)
    CROSSCOMPILE="yes"
    MBCONFARCH="haswell"
    MBCC="cc"
    MBCXX="CC"
    MBFC="ftn"
    MBF77="ftn"
    MBNMPICC="cc"
    MBNMPICXX="CC"
    MBNMPIFC="ftn"
    MBNMPIF77="ftn"
    MBMPI_DIR=""
    MBHDF5_DIR="/opt/cray/pe/hdf5-parallel/default/intel/150"
    MBNETCDF_DIR="/opt/cray/pe/netcdf-hdf5parallel/default/intel/150"
    MBPNETCDF_DIR="/opt/cray/pe/parallel-netcdf/default/intel/150"
    MBBLASLAPACK_LIBS="/opt/cray/pe/libsci/default/intel/150/$MBCONFARCH/lib/libsci_intel.a"
    MBMETIS_DIR="/opt/cray/pe/tpsl/default/INTEL/150/$MBCONFARCH"
    PREREQ="module load craype-haswell"
    ;;
  *edison* )
    CROSSCOMPILE="yes"
    MBCONFARCH="haswell"
    MBCC="cc"
    MBCXX="CC"
    MBFC="ftn"
    MBF77="ftn"
    MBNMPICC="cc"
    MBNMPICXX="CC"
    MBNMPIFC="ftn"
    MBNMPIF77="ftn"
    MBMPI_DIR=""
    MBHDF5_DIR="/opt/cray/hdf5-parallel/default/intel/15.0"
    MBNETCDF_DIR="/opt/cray/netcdf-hdf5parallel/4.4.0/intel/15.0"
    MBPNETCDF_DIR="/opt/cray/parallel-netcdf/1.7.0/intel/15.0"
    MBBLASLAPACK_LIBS="/opt/cray/libsci/default/intel/15.0/haswell/lib/libsci_intel.a"
    MBMETIS_DIR="/opt/cray/tpsl/16.07.1/INTEL/15.0/haswell"
    PREREQ="module load craype-haswell"
    ;;
  *)  # Nothing to do
    ;;
esac

# Finalize the configuration command
if (test "x$MBHOSTSYS" != "x"); then
  INTERNAL_OPTIONS="$INTERNAL_OPTIONS --host=$MBHOSTSYS"
fi

if (test "x$ENABLE_MPI" != "xno"); then
  if (test "x$MBMPI_DIR" != "x"); then
    if (test "x$CMAKE" != "xyes"); then
      INTERNAL_OPTIONS="$INTERNAL_OPTIONS --with-mpi=$MBMPI_DIR"
    else
      INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DENABLE_MPI=ON -DMPI_HOME=$MBMPI_DIR"
    fi
  else
    if (test "x$CMAKE" != "xyes"); then
      INTERNAL_OPTIONS="$INTERNAL_OPTIONS --with-mpi"
    fi
  fi
else
  MBCC=$MBNMPICC
  MBCXX=$MBNMPICXX
  MBF77=$MBNMPIF77
  MBFC=$MBNMPIFC
fi

if (test "x$CMAKE" != "xyes"); then
  INTERNAL_OPTIONS="$INTERNAL_OPTIONS CC=$MBCC CXX=$MBCXX"
else
  INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_C_COMPILER=$MBCC"
  INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_CXX_COMPILER=$MBCXX"
fi

if (test "x$ENABLE_FORTRAN" != "xno" && test "x$MBFC" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
	  INTERNAL_OPTIONS="$INTERNAL_OPTIONS FC=$MBFC F77=$MBF77"
	else
    INTERNAL_OPTIONS="$INTERNAL_OPTIONS -DCMAKE_Fortran_COMPILER=$MBFC"
  fi
else
	INTERNAL_OPTIONS="$INTERNAL_OPTIONS --disable-fortran"
fi

if (test "x$CMAKE" != "xyes" && test "x$CROSSCOMPILE" != "xno"); then
  INTERNAL_OPTIONS="$INTERNAL_OPTIONS cross_compiling=yes"
fi

LDFLAGSLIBS=""
DEPENDENCY_OPTIONS=""
if (test "x$MBBLASLAPACK_LIBS" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-blas=\"$MBBLASLAPACK_LIBS\" --with-lapack=\"$MBBLASLAPACK_LIBS\""
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_BLASLAPACK=ON LDFLAGSLIBS=\"$LDFLAGSLIBS $MBBLASLAPACK_LIBS\""
  fi
fi

if (test "x$MBHDF5_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-hdf5=$MBHDF5_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_HDF5=ON -DHDF5_ROOT\"$MBHDF5_DIR\""
  fi
fi

if (test "x$MBZLIB_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-zlib=$MBZLIB_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_ZLIB=ON -DZLIB_ROOT=$MBZLIB_DIR"
  fi
fi

if (test "x$MBSZIP_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-szip=$MBSZIP_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_SZIP=ON -DSZIP_ROOT=$MBSZIP_DIR"
  fi
fi

if (test "x$MBNETCDF_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-netcdf=$MBNETCDF_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_NETCDF=ON -DNETCDF_ROOT=$MBNETCDF_DIR"
  fi
fi

if (test "x$MBPNETCDF_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-pnetcdf=$MBPNETCDF_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_PNETCDF=ON -DPNETCDF_ROOT=$MBPNETCDF_DIR"
  fi
fi

if (test "x$MBMETIS_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-metis=$MBMETIS_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_METIS=ON -DMETIS_ROOT=$MBMETIS_DIR"
  fi
fi

if (test "x$MBPARMETIS_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-parmetis=$MBPARMETIS_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_PARMETIS=ON -DPARMETIS_ROOT=$MBPARMETIS_DIR"
  fi
fi

if (test "x$MBZOLTAN_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-zoltan=$MBZOLTAN_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DENABLE_ZOLTAN=ON -DZOLTAN_ROOT=$MBZOLTAN_DIR"
  fi
fi

if (test "x$MBSCOTCH_DIR" != "x"); then
  DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-scotch=$MBSCOTCH_DIR"
fi

if (test "x$MBPTSCOTCH_DIR" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-ptscotch=$MBPTSCOTCH_DIR"
  else
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS --with-ptscotch=$MBPTSCOTCH_DIR"
  fi
fi

if (test "x$LIBS" != "x"); then
  if (test "x$CMAKE" != "xyes"); then
    DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS LIBS=\"$LIBS\""
    MBDWLD_OPTIONS="$MBDWLD_OPTIONS LIBS=\"$LIBS\""
  else
    LDFLAGSLIBS=\"$LDFLAGSLIBS $LIBS\"
  fi
fi

if (test "x$CMAKE" != "xno"); then
  DEPENDENCY_OPTIONS="$DEPENDENCY_OPTIONS -DCMAKE_EXE_LINKER_FLAGS=\"$LDFLAGSLIBS\""
fi

# Put them all together
if (test "x$CMAKE" != "xyes"); then
  DWLD_CONFIGURE_CMD="$CONFIGURE_CMD $INTERNAL_OPTIONS $MBDWLD_OPTIONS --with-pic=1 --enable-tools --download-hdf5 --download-netcdf --download-metis"
fi
CONFIGURE_CMD="$CONFIGURE_CMD $INTERNAL_OPTIONS $DEPENDENCY_OPTIONS"

# PRINT OUT INFORMATION
echo "###########################################"
echo "   Hostname: $HOSTNAME"
echo "###########################################"

if (test "x$PREREQ" != "x"); then
  echo "MOAB Prerequisites     = $PREREQ"
fi
echo "MOAB Install path      = $PREFIX_INSTALL_PATH"
echo "Enable debug info      = $ENABLE_DEBUG"
echo "Enable optimization    = $ENABLE_OPTIMIZE"
echo "Enable shared build    = $ENABLE_SHARED"
echo "Enable static build    = $ENABLE_STATIC"
echo "Enable MPI parallelism = $ENABLE_MPI"
echo "Enable Fortran support = $ENABLE_FORTRAN"
echo -n "Compilers: CC=$MBCC, CXX=$MBCXX"
if (test "x$ENABLE_FORTRAN" != "xno" && test "x$MBFC" != "x"); then
  echo ", FC=$MBFC, F77=$MBF77"
else
  echo ""
fi
echo ""
echo "Configure command to use:"
echo "-------------------------"
echo "$CONFIGURE_CMD"
echo ""
if (test "x$CMAKE" != "xyes"); then
  echo "         OR"
  echo ""
  echo "$DWLD_CONFIGURE_CMD"
fi

# Done.
