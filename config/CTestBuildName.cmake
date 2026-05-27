# Derive CTEST_SITE and CTEST_BUILD_NAME from the resolved CMake state.
#
# Include this module from the top-level CMakeLists.txt AFTER all ENABLE_*
# options and find_package() calls have run, so that the feature suffix
# reflects what is actually enabled (not just what the user requested).
#
# Naming scheme:
#   <os>-<arch>-<compiler><major>-<buildtype>[+tag1+tag2...]
# Examples:
#   macOS-arm64-clang17-Rel+mpi+h5+nc+pnc+zlt+tr+eig
#   linux-x86_64-gcc13-Dbg+mpi+h5
#
# Overrides:
#   CTEST_SITE        env or -DCTEST_SITE=...        wins over auto-detection
#   CTEST_BUILD_NAME  env or -DCTEST_BUILD_NAME=...  wins over auto-detection
#
# Add new feature tags by appending "<OPTION>;<short-tag>" pairs to
# _moab_feature_map below.  Keep the order fixed so a given configuration
# always produces the same string (so CDash collapses re-runs into one row).

include_guard(GLOBAL)

# --- Site --------------------------------------------------------------------
if(DEFINED ENV{CTEST_SITE} AND NOT "$ENV{CTEST_SITE}" STREQUAL "")
    set(_moab_site "$ENV{CTEST_SITE}")
else()
    cmake_host_system_information(RESULT _moab_site QUERY HOSTNAME)
    # Trim trailing ".local" that macOS adds to Bonjour names.
    string(REGEX REPLACE "\\.local$" "" _moab_site "${_moab_site}")
endif()
set(CTEST_SITE "${_moab_site}" CACHE STRING "CDash site name" FORCE)

# --- OS / arch / compiler / build type --------------------------------------
if(APPLE)
    set(_moab_os "macOS")
elseif(WIN32)
    set(_moab_os "windows")
else()
    string(TOLOWER "${CMAKE_SYSTEM_NAME}" _moab_os)
endif()

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _moab_arch)
if(_moab_arch STREQUAL "")
    set(_moab_arch "unknown")
endif()

if(CMAKE_CXX_COMPILER_ID)
    string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" _moab_cc)
elseif(CMAKE_C_COMPILER_ID)
    string(TOLOWER "${CMAKE_C_COMPILER_ID}" _moab_cc)
else()
    set(_moab_cc "cc")
endif()

set(_moab_cc_ver "")
if(CMAKE_CXX_COMPILER_VERSION)
    string(REGEX REPLACE "^([0-9]+).*" "\\1"
           _moab_cc_ver "${CMAKE_CXX_COMPILER_VERSION}")
elseif(CMAKE_C_COMPILER_VERSION)
    string(REGEX REPLACE "^([0-9]+).*" "\\1"
           _moab_cc_ver "${CMAKE_C_COMPILER_VERSION}")
endif()
set(_moab_compiler "${_moab_cc}${_moab_cc_ver}")

if(CMAKE_BUILD_TYPE)
    string(SUBSTRING "${CMAKE_BUILD_TYPE}" 0 3 _moab_buildtype)
else()
    set(_moab_buildtype "Rel")
endif()

# --- Feature tags -----------------------------------------------------------
# Ordered (OPTION, tag) pairs. Append, do not rearrange — order is the
# stable contract that makes re-runs collapse to one CDash row.
set(_moab_feature_map
    "ENABLE_MPI;mpi"
    "ENABLE_HDF5;h5"
    "ENABLE_NETCDF;nc"
    "ENABLE_PNETCDF;pnc"
    "ENABLE_EIGEN3;eig"
    "ENABLE_ZOLTAN;zlt"
    "ENABLE_METIS;mts"
    "ENABLE_PARMETIS;pmts"
    "ENABLE_TEMPESTREMAP;tr"
    "ENABLE_BLASLAPACK;bla"
    "ENABLE_CGNS;cgns"
    "ENABLE_FORTRAN;f"
    "ENABLE_PYMOAB;py"
    )

set(_moab_tags "")
list(LENGTH _moab_feature_map _moab_n)
math(EXPR _moab_last "${_moab_n} - 1")
foreach(_i RANGE 0 ${_moab_last} 2)
    math(EXPR _j "${_i} + 1")
    list(GET _moab_feature_map ${_i} _moab_opt)
    list(GET _moab_feature_map ${_j} _moab_tag)
    if(${_moab_opt})
        set(_moab_tags "${_moab_tags}+${_moab_tag}")
    endif()
endforeach()

# --- Assemble (respecting any user override) --------------------------------
if(DEFINED ENV{CTEST_BUILD_NAME} AND NOT "$ENV{CTEST_BUILD_NAME}" STREQUAL "")
    set(_moab_buildname "$ENV{CTEST_BUILD_NAME}")
else()
    set(_moab_buildname
        "${_moab_os}-${_moab_arch}-${_moab_compiler}-${_moab_buildtype}${_moab_tags}")
endif()
set(CTEST_BUILD_NAME "${_moab_buildname}" CACHE STRING "CDash build name" FORCE)

# Mirror to BUILDNAME/SITE because some CTest paths still read the legacy
# names when generating DartConfiguration.tcl.
set(BUILDNAME "${CTEST_BUILD_NAME}" CACHE STRING "CDash build name" FORCE)
set(SITE      "${CTEST_SITE}"       CACHE STRING "CDash site name"  FORCE)

message(STATUS "CDash submission identity:")
message(STATUS "  Site       : ${CTEST_SITE}")
message(STATUS "  Build name : ${CTEST_BUILD_NAME}")
