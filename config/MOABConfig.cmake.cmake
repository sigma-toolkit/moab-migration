# Config file for MOAB; use the CMake find_package() function to pull this into
# your own CMakeLists.txt file.
#
# This file defines the following variables:
# MOAB_FOUND        - boolean indicating that MOAB is found
# MOAB_VERSION      - version of MOAB
# MOAB_INCLUDE_DIRS - include directories from which to pick up MOAB includes
# MOAB_LIBRARIES    - libraries need to link to MOAB; use this in target_link_libraries for MOAB-dependent targets
# MOAB_CXX, MOAB_CC, MOAB_F77, MOAB_FC - compilers used to compile MOAB
# MOAB_CXXFLAGS, MOAB_CCFLAGS, MOAB_FFLAGS, MOAB_FCFLAGS - compiler flags used to compile MOAB; possibly need to use these in add_definitions or CMAKE_<LANG>_FLAGS_<MODE> 

set(MOAB_FOUND 1)

set(MOAB_VERSION @PACKAGE_VERSION@)

# Check whether the requested PACKAGE_FIND_VERSION is compatible
if("${MOAB_VERSION}" VERSION_LESS "${PACKAGE_FIND_VERSION}")
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if("${MOAB_VERSION}" VERSION_EQUAL "${PACKAGE_FIND_VERSION}")
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()

set(MOAB_CC "@CMAKE_C_COMPILER@")
set(MOAB_CXX "@CMAKE_CXX_COMPILER@")
set(MOAB_FC "@CMAKE_Fortran_COMPILER@")
set(MOAB_F77 "@CMAKE_Fortran_COMPILER@")
# Compiler flags used by MOAB
set(MOAB_CFLAGS "@CFLAGS@ @CPPFLAGS@")
set(MOAB_CXXFLAGS "@CXXFLAGS@ @CPPFLAGS@")
set(MOAB_FCFLAGS "@FFLAGS@")
set(MOAB_FFLAGS "@FFLAGS@")

set(MOAB_BUILT_SHARED @BUILD_SHARED_LIBS@)
set(MOAB_USE_MPI @MOAB_HAVE_MPI@)
set(MPI_DIR "@MPI_HOME@")
set(MOAB_USE_HDF5 @MOAB_HAVE_HDF5@)
set(MOAB_USE_HDF5_PARALLEL @MOAB_HAVE_HDF5_PARALLEL@)
set(HDF5_DIR "@HDF5_ROOT@")
set(MOAB_USE_ZLIB @MOAB_HAVE_ZLIB@)
set(ZLIB_DIR "@ZLIB_DIR@")
set(MOAB_USE_SZIP @MOAB_HAVE_SZIP@)
set(SZIP_DIR "@SZIP_DIR@")
set(MOAB_USE_NETCDF @MOAB_HAVE_NETCDF@)
set(NETCDF_DIR "@NETCDF_DIR@")
set(MOAB_USE_PNETCDF @MOAB_HAVE_PNETCDF@)
set(PNETCDF_DIR "@PNETCDF_DIR@")
set(MOAB_USE_METIS @MOAB_HAVE_METIS@)
set(METIS_DIR "@METIS_DIR@")
set(MOAB_USE_PARMETIS @MOAB_HAVE_PARMETIS@)
set(PARMETIS_DIR "@PARMETIS_DIR@")
set(MOAB_USE_ZOLTAN @MOAB_HAVE_ZOLTAN@)
set(ZOLTAN_DIR "@ZOLTAN_DIR@")
set(MOAB_USE_BLAS @MOAB_HAVE_BLAS@)
set(BLAS_LIBRARIES "@BLAS_LIBRARIES@")
set(MOAB_USE_LAPACK @MOAB_HAVE_LAPACK@)
set(LAPACK_LIBRARIES "@LAPACK_LIBRARIES@")
set(MOAB_USE_EIGEN @MOAB_HAVE_EIGEN3@)
set(EIGEN3_DIR "@EIGEN3_DIR@")
set(TEMPESTREMAP_DIR "@TEMPESTREMAP_DIR@")
set(MOAB_USE_TEMPESTREMAP @MOAB_HAVE_TEMPESTREMAP@)

set(MOAB_MESH_DIR "@CMAKE_SOURCE_DIR@/MeshFiles/unittest")

# Library and include defs
get_filename_component(MOAB_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

include (${MOAB_CMAKE_DIR}/ResolveCompilerPaths.cmake)

if(NOT @SKBUILD@)
  # missing support for DAMSEL, CCMIO
  set (MOAB_PACKAGE_LIBS "@ZOLTAN_LIBRARIES@ @TEMPESTREMAP_LIBRARIES@ @PNETCDF_LIBRARIES@ @NETCDF_LIBRARIES@ @CGNS_LIBRARIES@ @HDF5_LIBRARIES@ @PARMETIS_LIBRARIES@ @METIS_LIBRARIES@ @LAPACK_LIBRARIES@ @BLAS_LIBRARIES@ @MPI_CXX_LIBRARIES@" )
  string(STRIP "${MOAB_PACKAGE_LIBS}" MOAB_PACKAGE_LIBS)
  set(MOAB_PACKAGE_LIBS_LIST ${MOAB_PACKAGE_LIBS})
  separate_arguments(MOAB_PACKAGE_LIBS_LIST)
  list(REMOVE_DUPLICATES MOAB_PACKAGE_LIBS_LIST)
  set(MOAB_PACKAGE_LIBS "${MOAB_PACKAGE_LIBS_LIST}")

  set (MOAB_PACKAGE_INCLUDES_LIST "@ZOLTAN_INCLUDES@ @PNETCDF_INCLUDES@ @NETCDF_INCLUDES@ @HDF5_INCLUDES@ @PARMETIS_INCLUDES@ @METIS_INCLUDES@ @TEMPESTREMAP_INCLUDES@ @EIGEN3_INCLUDES@" )
  string(STRIP "${MOAB_PACKAGE_INCLUDES_LIST}" MOAB_PACKAGE_INCLUDES_LIST)
  RESOLVE_INCLUDES(MOAB_PACKAGE_INCLUDES "${MOAB_PACKAGE_INCLUDES_LIST}")
  separate_arguments(MOAB_PACKAGE_INCLUDES)
  list(REMOVE_DUPLICATES MOAB_PACKAGE_INCLUDES)

  if(NOT TARGET MOAB AND NOT MOAB_BINARY_DIR)
    include("${MOAB_CMAKE_DIR}/MOABTargets.cmake")
  endif()

  # Target information
  if(MOAB_USE_HDF5)
    if(EXISTS "@HDF5_ROOT@/share/cmake/hdf5/hdf5-config.cmake")
      include(@HDF5_ROOT@/share/cmake/hdf5/hdf5-config.cmake)
    endif()
  endif()
  set(MOAB_LIBRARY_DIRS "@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_LIBDIR@")
  set(MOAB_INCLUDE_DIRS "@CMAKE_INSTALL_PREFIX@/include" ${MOAB_PACKAGE_INCLUDES})
  set(MOAB_LIBS "-lMOAB")
  set(MOAB_LIBRARIES "-L@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_LIBDIR@ ${MOAB_LIBS} ${MOAB_PACKAGE_LIBS}")
else()
  find_package(Python COMPONENTS Interpreter REQUIRED)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import pymoab; print(pymoab.include_path)"
    OUTPUT_VARIABLE MOAB_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import pymoab; print(pymoab.lib_path)"
    OUTPUT_VARIABLE MOAB_LIBRARY_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import pymoab; print(' '.join(pymoab.extra_lib))"
    OUTPUT_VARIABLE MOAB_EXTRA_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(MOAB_EXTRA_LIBRARIES)
    message(FATAL_ERROR
        "This build of MOAB is not supported. "
        "It appears that the wheel was repaired using tools like auditwheel or delocate, "
        "that modifies the shared libraries, which may cause problems.\n"
        "MOAB_EXTRA_LIBRARIES is not empty: ${MOAB_EXTRA_LIBRARIES}.\n"
        "To resolve this, please build MOAB from scratch. "
        "For more information, visit: https://bitbucket.org/fathomteam/moab\n"
    )
  endif()
endif()
