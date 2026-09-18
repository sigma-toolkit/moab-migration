# Config file for MOAB; use the CMake find_package() function to pull this into
# your own CMakeLists.txt file.
#
# The supported way to consume MOAB is the imported target:
#
#   find_package(MOAB REQUIRED)
#   target_link_libraries(myapp PRIVATE MOAB::MOAB)
#
# which carries MOAB's include directories, compile features and link
# dependencies with it.  MOAB::mbcoupler is available when MOAB was built with
# the mesh coupler.
#
# The variables below are still defined, in terms of those targets, so that
# existing projects keep building.  Prefer the target in new code.
#
# MOAB_FOUND        - boolean indicating that MOAB is found
# MOAB_VERSION      - version of MOAB
# MOAB_INCLUDE_DIRS - include directories from which to pick up MOAB includes
# MOAB_LIBRARIES    - what to link against; this is now the MOAB::MOAB target
# MOAB_CXX, MOAB_CC, MOAB_F77, MOAB_FC - compilers used to compile MOAB
# MOAB_CXXFLAGS, MOAB_CCFLAGS, MOAB_FFLAGS, MOAB_FCFLAGS - compiler flags used to compile MOAB; possibly need to use these in add_definitions or CMAKE_<LANG>_FLAGS_<MODE>

@PACKAGE_INIT@

# @PACKAGE_VERSION@ was used here, but that variable is only ever set much later
# in the top-level CMakeLists, so MOAB_VERSION came out empty in every generated
# config and find_package(MOAB 5.6) could not check anything.
set(MOAB_VERSION @MOAB_VERSION_STRING@)

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
set(MPI_DIR "@MPI_ROOT@")
set(MOAB_USE_HDF5 @MOAB_HAVE_HDF5@)
set(MOAB_USE_HDF5_PARALLEL @MOAB_HAVE_HDF5_PARALLEL@)
set(HDF5_DIR "@HDF5_DIR@")
# HDF5_DIR is only set when HDF5 was found through its own config package; a
# plain FindHDF5 search records the prefix in HDF5_ROOT instead, and one of the
# two is what find_dependency(HDF5) below needs as a hint.
set(HDF5_ROOT "@HDF5_ROOT@")
# The six <PKG>_DIR hints below are cache entries rather than plain variables,
# and that is load-bearing.  Each MOAB Find module opens with
#
#   set(<PKG>_DIR "<default>" CACHE PATH "...")
#
# and under CMP0126 OLD - which is what any consumer requiring less than CMake
# 3.21 gets - that call *deletes* a normal variable of the same name.  A hint
# passed as a normal variable therefore disappeared at the very moment the
# module started to use it, and the module's own default won instead:
# find_dependency(NETCDF) aborted with "User specified NETCDF_DIR (/usr) does
# not match the prefix from nc-config".  Seeding the cache gives the intended
# precedence for free, since set(CACHE) without FORCE leaves an existing entry
# alone: a consumer's own -D<PKG>_DIR still wins over what MOAB recorded here.
set(MOAB_USE_NETCDF @MOAB_HAVE_NETCDF@)
set(NETCDF_DIR "@NETCDF_DIR@" CACHE PATH "Path to search for NETCDF header and library files")
set(MOAB_USE_PNETCDF @MOAB_HAVE_PNETCDF@)
set(PNETCDF_DIR "@PNETCDF_DIR@" CACHE PATH "Path to search for PNetCDF header and library files")
set(MOAB_USE_METIS @MOAB_HAVE_METIS@)
set(METIS_DIR "@METIS_DIR@" CACHE PATH "Path to search for Metis header and library files")
set(MOAB_USE_PARMETIS @MOAB_HAVE_PARMETIS@)
set(PARMETIS_DIR "@PARMETIS_DIR@" CACHE PATH "Path to search for ParMetis header and library files")
set(MOAB_USE_ZOLTAN @MOAB_HAVE_ZOLTAN@)
set(ZOLTAN_DIR "@ZOLTAN_DIR@" CACHE PATH "Path to search for Zoltan header and library files")
set(MOAB_USE_BLAS @MOAB_HAVE_BLAS@)
set(BLAS_LIBRARIES "@BLAS_LIBRARIES@")
set(MOAB_USE_LAPACK @MOAB_HAVE_LAPACK@)
set(LAPACK_LIBRARIES "@LAPACK_LIBRARIES@")
set(MOAB_USE_EIGEN @MOAB_HAVE_EIGEN3@)
set(EIGEN3_DIR "@EIGEN3_DIR@")
set(TEMPESTREMAP_DIR "@TEMPESTREMAP_DIR@" CACHE PATH "Path to search for TempestRemap header and library files")
set(MOAB_USE_TEMPESTREMAP @MOAB_HAVE_TEMPESTREMAP@)
set(MOAB_USE_MBCOUPLER @MOAB_HAVE_MBCOUPLER@)
set(MOAB_USE_SKBUILD @SKBUILD@)
set(MOAB_MESH_DIR "@CMAKE_SOURCE_DIR@/MeshFiles/unittest")

# Library and include defs
get_filename_component(MOAB_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

include (${MOAB_CMAKE_DIR}/ResolveCompilerPaths.cmake)

# missing support for DAMSEL, CCMIO
#set (MOAB_PACKAGE_LIBS @ZOLTAN_LIBRARIES@ @TEMPESTREMAP_LIBRARIES@ @PNETCDF_LIBRARIES@ @NETCDF_LIBRARIES@ @CGNS_LIBRARIES@ @HDF5_LIBRARIES@ @PARMETIS_LIBRARIES@ @METIS_LIBRARIES@ @LAPACK_LIBRARIES@ @BLAS_LIBRARIES@ @MPI_CXX_LIBRARIES@ )
set (MOAB_PACKAGE_LIBS @MOAB_LIBRARIES@)
string(STRIP "${MOAB_PACKAGE_LIBS}" MOAB_PACKAGE_LIBS)
set(MOAB_PACKAGE_LIBS_LIST ${MOAB_PACKAGE_LIBS})
separate_arguments(MOAB_PACKAGE_LIBS_LIST)
# separate_arguments() splits on whitespace, which tears a macOS framework
# reference such as "-framework Accelerate" into two list items.  Consumers
# then pass "-framework" and "Accelerate" to target_link_libraries()
# separately and CMake turns the latter into "-lAccelerate", so the link
# fails with "framework '-lAccelerate' not found".  Re-pair each flag with
# the argument that follows it before de-duplicating, so the pair stays a
# single item and dedups as a unit.
set(_moab_pkg_libs "")
set(_moab_pending "")
foreach(_moab_item IN LISTS MOAB_PACKAGE_LIBS_LIST)
  if(_moab_pending)
    list(APPEND _moab_pkg_libs "${_moab_pending} ${_moab_item}")
    set(_moab_pending "")
  elseif(_moab_item MATCHES "^(-framework|-Xlinker)$")
    set(_moab_pending "${_moab_item}")
  else()
    list(APPEND _moab_pkg_libs "${_moab_item}")
  endif()
endforeach()
if(_moab_pending)
  list(APPEND _moab_pkg_libs "${_moab_pending}")
endif()
set(MOAB_PACKAGE_LIBS_LIST ${_moab_pkg_libs})
unset(_moab_pkg_libs)
unset(_moab_pending)
unset(_moab_item)
list(REMOVE_DUPLICATES MOAB_PACKAGE_LIBS_LIST)
set(MOAB_PACKAGE_LIBS "${MOAB_PACKAGE_LIBS_LIST}")

set (MOAB_PACKAGE_INCLUDES_LIST "@ZOLTAN_INCLUDES@ @PNETCDF_INCLUDES@ @NETCDF_INCLUDES@ @HDF5_INCLUDES@ @PARMETIS_INCLUDES@ @METIS_INCLUDES@ @TEMPESTREMAP_INCLUDES@ @EIGEN3_INCLUDES@" )
string(STRIP "${MOAB_PACKAGE_INCLUDES_LIST}" MOAB_PACKAGE_INCLUDES_LIST)
RESOLVE_INCLUDES(MOAB_PACKAGE_INCLUDES "${MOAB_PACKAGE_INCLUDES_LIST}")
separate_arguments(MOAB_PACKAGE_INCLUDES)
list(REMOVE_DUPLICATES MOAB_PACKAGE_INCLUDES)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
#
# libMOAB's exported link interface names imported targets - MPI::MPI_CXX,
# HDF5::HDF5, NetCDF::NetCDF and so on - rather than the absolute paths those
# resolved to on the machine MOAB was built on.  Recreating those names here is
# what makes an installed MOAB usable elsewhere, so this block has to run before
# MOABTargets.cmake is included below.
#
# The MOAB-specific Find modules are installed next to this file, which is also
# where MOABTPLTargets.cmake (the helper that declares the targets) lives.
include(CMakeFindDependencyMacro)
if(@MOAB_CONFIG_IS_BUILD_TREE@)
  # Nothing is installed yet, so read them straight out of the source tree.
  list(APPEND CMAKE_MODULE_PATH "@PROJECT_SOURCE_DIR@/config")
else()
  list(APPEND CMAKE_MODULE_PATH "${MOAB_CMAKE_DIR}")
endif()

# Every module below probes with the C or C++ compiler, and FindMPI fails
# outright when asked for a component whose language is not enabled.
# find_package(MOAB) before project() is a legitimate thing to do - it is how a
# project reuses MOAB_CXX as its own compiler - so in that case define the
# variables and skip the targets.  Such a caller must call find_package(MOAB)
# again once it has enabled its languages, before it links MOAB::MOAB.
get_property(_moab_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if("CXX" IN_LIST _moab_languages)

  # moab_mpi.h is a public header and includes <mpi.h>, so a consumer of a
  # parallel MOAB needs the MPI include path whether or not it uses MPI itself.
  if(MOAB_USE_MPI AND NOT TARGET MPI::MPI_CXX)
    find_dependency(MPI COMPONENTS CXX)
  endif()

  if(MOAB_USE_HDF5 AND NOT TARGET HDF5::HDF5)
    if(MOAB_USE_HDF5_PARALLEL)
      set(HDF5_PREFER_PARALLEL TRUE)
    endif()
    # A "<var>-NOTFOUND" left in HDF5_DIR would send config mode looking in a
    # directory that does not exist; clear it and let FindHDF5 use HDF5_ROOT.
    if(HDF5_DIR MATCHES "NOTFOUND$")
      unset(HDF5_DIR)
    endif()
    find_dependency(HDF5 COMPONENTS C HL)
    # HDF5 can be built against libcurl (the ROS3 virtual file driver); MOAB
    # adds it to the link line when it is present, so it has to come back.
    find_dependency(CURL)
  endif()

  # The <PKG>_DIR hints these modules search were seeded near the top of this
  # file, where the reason they have to be cache entries is spelled out.
  if(MOAB_USE_NETCDF AND NOT TARGET NetCDF::NetCDF)
    find_dependency(NETCDF)
  endif()

  if(MOAB_USE_PNETCDF AND NOT TARGET PNetCDF::PNetCDF)
    # FindPNETCDF does its whole search inside "if(MOAB_HAVE_MPI AND
    # ENABLE_PNETCDF)".  Those are knobs from MOAB's own configure and mean
    # nothing in a consumer's project, so the module skipped the search and
    # then failed on the empty result.  Reaching this line at all means MOAB
    # was built with parallel NetCDF, so both are true by construction.
    set(MOAB_HAVE_MPI TRUE)
    set(ENABLE_PNETCDF TRUE)
    find_dependency(PNETCDF)
  endif()

  if(MOAB_USE_METIS AND NOT TARGET METIS::METIS)
    find_dependency(METIS)
  endif()

  if(MOAB_USE_PARMETIS AND NOT TARGET ParMETIS::ParMETIS)
    find_dependency(PARMETIS)
  endif()

  if(MOAB_USE_ZOLTAN AND NOT TARGET Zoltan::Zoltan)
    find_dependency(ZOLTAN)
  endif()

  if(MOAB_USE_TEMPESTREMAP AND NOT TARGET TempestRemap::TempestRemap)
    find_dependency(TEMPESTREMAP)
  endif()

  if(MOAB_USE_EIGEN AND NOT TARGET Eigen3::Eigen)
    if(NOT EIGEN3_INCLUDE_DIR)
      set(EIGEN3_INCLUDE_DIR "@EIGEN3_INCLUDE_DIR@")
    endif()
    # MOAB may have been built against an unpacked Eigen source tree, which ships
    # no Eigen3Config.cmake.  find_dependency() cannot see such a tree, so it
    # would either fail outright or silently bind this consumer to a different
    # system-wide Eigen than MOAB was compiled against.  Recreate the
    # header-only target from the recorded path instead; this mirrors the
    # EIGEN3_INCLUDE_DIR branch in MOAB's own top-level CMakeLists.txt.
    if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Eigen")
      add_library(Eigen3::Eigen INTERFACE IMPORTED)
      set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    else()
      find_dependency(Eigen3)
    endif()
  endif()

  if(MOAB_USE_LAPACK AND NOT TARGET LAPACK::LAPACK)
    find_dependency(LAPACK)
  endif()
  if(MOAB_USE_BLAS AND NOT TARGET BLAS::BLAS)
    find_dependency(BLAS)
  endif()

  if(NOT TARGET Threads::Threads)
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    find_dependency(Threads)
  endif()

endif()
unset(_moab_languages)

if(MOAB_USE_SKBUILD)
  # Find the Python interpreter and ensure it's available.
  find_package(Python COMPONENTS Interpreter REQUIRED)

  # Function to run Python commands and validate their execution.
  function(run_python_command output_var command)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -c "${command}"
      OUTPUT_VARIABLE ${output_var}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE result
      )
    # Check if the command was successful
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "Failed to run Python command: ${command}")
    else()
      # Add the output variable to the parent scope
      set(${output_var} "${${output_var}}" PARENT_SCOPE)
    endif()
  endfunction()

  # Extract MOAB include paths, library paths, and extra libraries
  run_python_command(MOAB_INCLUDE_DIRS "import pymoab; print(pymoab.include_path[0])")
  run_python_command(MOAB_LIBRARY_DIRS "import pymoab; print(pymoab.lib_path[0])")
  run_python_command(MOAB_EXTRA_LIBRARIES "import pymoab; print(' '.join(pymoab.extra_lib))")

  # Check if the wheel was repaired using auditwheel or delocate
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

  # Add MOAB targets
  file(TO_CMAKE_PATH "${MOAB_LIBRARY_DIRS}/cmake/MOAB/MOABTargets.cmake" MOAB_TARGETS_FILE)
  include(${MOAB_TARGETS_FILE})

  set(MOAB_INCLUDE_DIRS ${MOAB_INCLUDE_DIRS} ${MOAB_PACKAGE_INCLUDES})
elseif(@MOAB_CONFIG_IS_BUILD_TREE@)
  if(NOT TARGET MOAB::MOAB AND NOT MOAB_BINARY_DIR)
    include("${MOAB_CMAKE_DIR}/MOABTargets.cmake")
  endif()
  set(MOAB_LIBRARY_DIRS "@PROJECT_BINARY_DIR@/lib")
  set(MOAB_INCLUDE_DIRS "@CMAKE_SOURCE_DIR@/src" "@PROJECT_BINARY_DIR@/src" ${MOAB_PACKAGE_INCLUDES})
else()
  if(NOT TARGET MOAB::MOAB AND NOT MOAB_BINARY_DIR)
    include("${MOAB_CMAKE_DIR}/MOABTargets.cmake")
  endif()
  # PACKAGE_PREFIX_DIR is set by the package-init preamble above and is derived
  # from where this file actually sits, so an install tree that has been moved
  # or relocated still resolves.  These used to be the configure-time prefix.
  # (Do not name the init placeholder in a comment here: it is a substitution
  # token, so configure_package_config_file() would expand it a second time.)
  set_and_check(MOAB_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")
  set_and_check(MOAB_LIBRARY_DIRS "${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_LIBDIR@")
  set(MOAB_INCLUDE_DIRS "${MOAB_INCLUDE_DIR}" ${MOAB_PACKAGE_INCLUDES})
endif()

# MOAB_LIBRARIES was a raw "-L<prefix>/lib -lMOAB <third party libs>" string.
# CMake rejects a link item with leading or trailing whitespace outright, and
# even when it did not, the string bypassed the imported target: none of MOAB's
# include directories, compile features or transitive dependencies came with it.
# Name the target instead - target_link_libraries(${MOAB_LIBRARIES}) keeps
# working for every existing consumer, and now carries the usage requirements.
set(MOAB_LIBRARIES MOAB::MOAB)

# Projects written against older MOAB releases link the un-namespaced "MOAB"
# name directly.  The export is namespaced now, so give them an INTERFACE
# target that forwards.  It is IMPORTED, so nothing is built for it.
if(NOT TARGET MOAB)
  add_library(MOAB INTERFACE IMPORTED)
  set_target_properties(MOAB PROPERTIES INTERFACE_LINK_LIBRARIES MOAB::MOAB)
endif()

# Include standard argument handling for finding packages
include(FindPackageHandleStandardArgs)

# Validates that the necessary variables are set
find_package_handle_standard_args(MOAB
  REQUIRED_VARS MOAB_LIBRARIES MOAB_INCLUDE_DIRS
  VERSION_VAR MOAB_VERSION
  )

# Defined by the package-init preamble.  MOAB exports no optional components,
# so this only reports an error if a caller asks for one that does not exist.
check_required_components(MOAB)
