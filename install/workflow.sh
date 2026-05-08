# Usage commands on machines
#
# GCE
MPICH_DIR=/nfs/gce/projects/climate/software/linux-ubuntu22.04-x86_64/mpich/4.1.2/gcc-12.1.0
./install-moab-e3sm.sh --mpi-root=$MPICH_DIR \
  --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_PATH --pnetcdf-root=$PNETCDF_PATH
# Perlmutter: works the same for GNU and Intel envs
./install-moab-e3sm.sh --mpi-root=$MPICH_DIR --cc=cc --cxx=CC --fc=ftn --f77=ftn --prefix=$PWD/installs \
  --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_PATH --pnetcdf-root=$PNETCDF_PATH  
# Bebop
./install-moab-e3sm.sh --mpi-root=/software/software/custom-built/openmpi/4.1.8/gcc/13.2.0 \
  --prefix=$PWD/installs --hdf5-root=/lcrc/group/e3sm/soft/bebop/hdf5/1.12.3/gcc-13.2.0/openmpi-4.1.8 \
  --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH  \
  --extra-tempestremap='--with-blas="/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/libblas.a -lgfortran" --with-lapack="/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/liblapack.a -lm"' \
  --extra='--with-blas="/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/libblas.a -lgfortran" --with-lapack="/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/liblapack.a -lm"'

