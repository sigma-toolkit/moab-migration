cmake_minimum_required( VERSION 3.15 )

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT WIN32)
  #This all breaks on windows.
  set (CMAKE_Fortran_COMPILER_INIT ${CMAKE_GENERATOR_FC})
  set (CMAKE_CXX_FLAGS_INIT "-fPIC -DPIC")
  set (CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  set (CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  set (CMAKE_C_FLAGS_INIT "-fPIC -DPIC")
  set (CMAKE_C_FLAGS_DEBUG_INIT "-O0 -g")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  set (CMAKE_C_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_INIT "-fPIC")
  set (CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2")
  set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS_INIT "")
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS_INIT "")
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel.")
endif()

project( MOAB 
         VERSION 5.5.1
         DESCRIPTION "Mesh-Oriented datABase (MOAB) library"
         LANGUAGES C CXX Fortran
       )
set (PACKAGE_NAME "MOAB")
set (PACKAGE_VERSION "5.5.1")

# Include standard CMake macros
include(GNUInstallDirs)
include(FeatureSummary)
include(ExternalData)
# Include local CMake macros
include(config/logging.cmake)
include(config/dist.cmake)
include(config/distcheck.cmake)
#include(config/release.cmake)

LOGGING_INITIALIZE()
_SETUP_PROJECT_DIST()
DISTCHECK_SETUP()
#RELEASE_SETUP()
message(STATUS "MOAB: Configured 'distcheck' target")

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})
################################################################################
# Set up version info
################################################################################
include (config/GetAcInitVersion.cmake)
get_ac_init_version()
set ( MOAB_VERSION_STRING "${PACKAGE_NAME} ${VERSION_STRING}" )
set ( MOAB_VERSION        "${VERSION_STRING}" )
set ( MOAB_VERSION_MAJOR  ${MAJOR_VERSION}  )
set ( MOAB_VERSION_MINOR  ${MINOR_VERSION}  )
if ( DEFINED PATCH_VERSION )
  set ( MOAB_VERSION_PATCH "${PATCH_VERSION}" )
else ( DEFINED PATCH_VERSION )
  if ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (alpha)" )
  else ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (beta)" )
  endif ( MOAB_VERSION_MINOR EQUAL 99 )
endif ( DEFINED PATCH_VERSION )

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0020 NEW)
  #cmake_policy(SET CMP0042 OLD)
ENDIF ()
if (NOT (CMAKE_VERSION VERSION_LESS "3.12.0"))
# Allow use of <PackageName>_ROOT as a cmake variable
  cmake_policy(SET CMP0074 NEW)
endif ()


set ( abs_srcdir   ${CMAKE_SOURCE_DIR} )
set ( abs_builddir ${PROJECT_BINARY_DIR} )

################################################################################
# Install Related Settings
################################################################################

get_filename_component( MOAB_ABSSRC_DIR2 moab.make.in REALPATH )
execute_process( COMMAND dirname ${MOAB_ABSSRC_DIR2}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE MOAB_ABSSRC_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

## Set the directory where the binaries will be stored
set( EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Directory where all executable will be stored"
)

## Set the directory where the libraries will be stored
set( LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Directory where all the libraries will be stored"
)
mark_as_advanced(
  MOAB_ABSSRC_DIR
  EXECUTABLE_OUTPUT_PATH
  LIBRARY_OUTPUT_PATH)

include ( CheckIncludeFile )
include ( CheckFunctionExists )
include ( CheckTypeSize )

# Compiler defines... this should really be in a config file.
set( MOAB_DEFINES "" )
set( MOAB_INCLUDES "" )
set( MOAB_LIBRARIES "" )
set( MOAB_INSTALL_TARGETS "" )
set( MOAB_MESH_DIR        "${CMAKE_SOURCE_DIR}/MeshFiles/" )

# TODO: Need to update this to GNUInstallDirs
set(MOABDocLocation share/doc/moab)

include ( config/GenerateScript.cmake )
include ( config/ResolveCompilerPaths.cmake )

# Macro to automatically enable a package if its directory variable is set
macro(AUTO_ENABLE_PACKAGE PACKAGE_NAME)
    if (DEFINED ${PACKAGE_NAME}_DIR)
        set(ENABLE_${PACKAGE_NAME} ON CACHE BOOL "Include ${PACKAGE_NAME} support in the build?" FORCE)
    endif()
endmacro()

################################################################################
# Options that the user controls
################################################################################
option ( BUILD_SHARED_LIBS   "Should shared or static libraries be created?"   ON  )
option ( ENABLE_CPM       "Should build with C++/CMake Package Manager (CPM) support?"  OFF )
option ( ENABLE_CGNS       "Should build with CGNS support?"                 OFF )
option ( ENABLE_MPI        "Should MOAB be compiled with MPI support?"       OFF )
option ( ENABLE_HDF5       "Include HDF I/O interfaces in the build?"                   OFF )
option ( ENABLE_NETCDF     "Include NetCDF (ExodusII) interfaces in the build?" OFF )
option ( ENABLE_PNETCDF    "Include PNetCDF interfaces in the build (Requires NetCDF) ?" OFF )
option ( ENABLE_METIS      "Include Metis support for partitioning algorithms?" OFF )
option ( ENABLE_PARMETIS   "Include ParMetis support for partitioning algorithms?" OFF )
option ( ENABLE_ZOLTAN     "Include Zoltan support for partitioning algorithms?" OFF )
option ( ENABLE_TESTING    "Enable Testing"                                  ON  )
option ( MOAB_FORCE_64_BIT_HANDLES "Force MBEntityHandle to be 64 bits (uint64_t)" OFF )
option ( MOAB_FORCE_32_BIT_HANDLES "Force MBEntityHandle to be 32 bits (uint32_t)" OFF )
option ( ENABLE_PYMOAB       "Build python interface to MOAB?" OFF )
option ( ENABLE_BLASLAPACK   "Should require BLAS/LAPACK?" ON )
option ( ENABLE_EIGEN3       "Should require Eigen3 header-only library?" OFF )
option ( ENABLE_FORTRAN      "Should enable FORTRAN?" ON )
option ( ENABLE_TEMPESTREMAP "Should enable TempestRemap dependency?" OFF )

add_feature_info(SharedLibs BUILD_SHARED_LIBS "Should shared or static libraries be created?")
add_feature_info(EnableFortran ENABLE_FORTRAN "Should FORTRAN support be enabled?")

# Automatically enable NetCDF if NETCDF_DIR is set
auto_enable_package(MPI)
auto_enable_package(HDF5)
auto_enable_package(NETCDF)
auto_enable_package(PNETCDF)
auto_enable_package(METIS)
auto_enable_package(PARMETIS)
auto_enable_package(ZOLTAN)
auto_enable_package(EIGEN3)
auto_enable_package(TEMPESTREMAP)

find_package(PkgConfig QUIET)

mark_as_advanced(
  MOAB_FORCE_64_BIT_HANDLES
  MOAB_FORCE_32_BIT_HANDLES
)

if ( ENABLE_FORTRAN )
  include (${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
endif ()

# patch added to ensure libs have a full-path ID set on Mac
if(APPLE)
  option(CMAKE_MACOSX_RPATH "Use MacOSX RPATH" OFF)
  if(SKBUILD)
    set(CMAKE_MACOSX_RPATH ON)
  endif()
  message(STATUS "CMAKE_MACOSX_RPATH = ${CMAKE_MACOSX_RPATH}")
  if(NOT CMAKE_INSTALL_NAME_DIR AND NOT CMAKE_MACOSX_RPATH)
    set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
  endif()
  message(STATUS "CMAKE_INSTALL_NAME_DIR = ${CMAKE_INSTALL_NAME_DIR}")
endif()

################################################################################
# Check for system include files
################################################################################
check_include_file( inttypes.h   MOAB_HAVE_INTTYPES_H )
check_include_file( stdint.h     MOAB_HAVE_STDINT_H )
check_include_file( stddef.h     MOAB_HAVE_STDDEF_H )
check_include_file( stdlib.h     MOAB_HAVE_STDLIB_H )
check_include_file( sys/types.h  MOAB_HAVE_SYS_TYPES_H )
set( HAVE_SYS_TYPES_H ${MOAB_HAVE_SYS_TYPES_H} )
set( HAVE_STDDEF_H    ${MOAB_HAVE_STDDEF_H} )
set( HAVE_STDINT_H    ${MOAB_HAVE_STDINT_H} )
set( HAVE_INTTYPES_H    ${MOAB_HAVE_INTTYPES_H} )
set( HAVE_STDLIB_H    ${MOAB_HAVE_STDLIB_H} )
check_include_file( memory.h     MOAB_HAVE_MEMORY_H )

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)
################################################################################
# Integer size Related Settings
################################################################################
if ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )
  message( FATAL_ERROR
      "You may not turn both MOAB_FORCE_64_BIT_HANDLES and MOAB_FORCE_32_BIT_HANDLES on. Turn one off to continue."
    )
endif ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )

if ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )
  if ( MOAB_HAVE_INTTYPES_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};inttypes.h" )
  endif ( MOAB_HAVE_INTTYPES_H )
  if ( MOAB_HAVE_STDLIB_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};stdlib.h" )
    #set ( CHECK_TYPE_SIZE_PREMAIN "${CHECK_TYPE_SIZE_PREMAIN}\n#include <stdlib.h>\n" )
  endif ( MOAB_HAVE_STDLIB_H )
  check_type_size(  size_t       MOAB_HAVE_SIZE_T )
  check_type_size(  ptrdiff_t    MOAB_HAVE_PTRDIFF_T )
  set ( HAVE_SIZE_T ${MOAB_HAVE_SIZE_T} )
  set ( HAVE_PTRDIFF_T ${MOAB_HAVE_PTRDIFF_T} )
endif ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )

# Check of presence of finite or isfinite macros
include (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(std::isfinite MOAB_HAVE_STDISFINITE)
CHECK_FUNCTION_EXISTS(isfinite MOAB_HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(finite MOAB_HAVE_FINITE)

################################################################################
# Find packages
################################################################################

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

include(config/CMakeLibraryFuncs.cmake)
include(GenerateExportHeader)
if (BUILD_SHARED_LIBS)
  if (WIN32)
    add_definitions(/DMOAB_DLL)
    if (NOT DEFINED CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS)
# https://cmake.org/cmake/help/v3.4/prop_tgt/WINDOWS_EXPORT_ALL_SYMBOLS.html
# Export all symbols
      set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_USE_MATH_DEFINES")
    endif ()
  else ()
    # add_definitions(-DMOAB_DLL) # Not needed?
  endif ()
endif ()


# Macro to automatically enable a package if its directory variable is set
macro(SET_VARS_IF_FOUND PACKAGE_NAME PACKAGE_TYPE PRINT_NAME DESCRIPTION PACKAGE_INCLUDES PACKAGE_LIBRARIES)
    if ( ${PACKAGE_NAME}_FOUND )
      set (MOAB_HAVE_${PACKAGE_NAME} ON CACHE INTERNAL "Found necessary ${PACKAGE_NAME} components. Configure MOAB with it.")

      # print some Metis configuration details
      message (STATUS "   ${PRINT_NAME} Configuration ::")
      if (DEFINED ${PACKAGE_NAME}_DIR AND NOT "${${PACKAGE_NAME}_DIR}" STREQUAL "")
        message (STATUS "        Directory : ${${PACKAGE_NAME}_DIR}")
      endif (DEFINED ${PACKAGE_NAME}_DIR AND NOT "${${PACKAGE_NAME}_DIR}" STREQUAL "")
      if (DEFINED ${PACKAGE_INCLUDES} AND NOT "${${PACKAGE_INCLUDES}}" STREQUAL "")
        message (STATUS "        Includes  : ${${PACKAGE_INCLUDES}}")
        include_directories( ${${PACKAGE_INCLUDES}} )
        #RESOLVE_INCLUDES(_tmp_resolved_pkgincs "${${PACKAGE_INCLUDES}}")
        set(MOAB_INCLUDES ${MOAB_INCLUDES} "${${PACKAGE_INCLUDES}}")
        #list( APPEND MOAB_INCLUDES "${_tmp_resolved_pkgincs}" )
      endif (DEFINED ${PACKAGE_INCLUDES} AND NOT "${${PACKAGE_INCLUDES}}" STREQUAL "")
      if (DEFINED ${PACKAGE_LIBRARIES} AND NOT "${${PACKAGE_LIBRARIES}}" STREQUAL "")
        message (STATUS "        Libraries : ${${PACKAGE_LIBRARIES}}")
        RESOLVE_LIBRARIES(_tmp_resolved_pkglibs "${${PACKAGE_LIBRARIES}}")
        # set( MOAB_LIBRARIES ${${PACKAGE_LIBRARIES}} ${MOAB_LIBRARIES} )
        # list( APPEND MOAB_LIBRARIES ${_tmp_resolved_pkglibs} )
        list( INSERT MOAB_LIBRARIES 0 "${_tmp_resolved_pkglibs}" )
      endif (DEFINED ${PACKAGE_LIBRARIES} AND NOT "${${PACKAGE_LIBRARIES}}" STREQUAL "")

      set_package_properties(${PACKAGE_NAME} PROPERTIES
                        TYPE ${PACKAGE_TYPE}
                        PURPOSE ${DESCRIPTION})

    else ( ${PACKAGE_NAME}_FOUND )
      set (MOAB_HAVE_${PACKAGE_NAME} OFF CACHE INTERNAL "Did not find necessary ${PACKAGE_NAME} components. Configure MOAB without it.")
    endif ( ${PACKAGE_NAME}_FOUND )
endmacro()

# ---- TPL Dependencies ----
if (ENABLE_CPM)
  set(FETCHCONTENT_BASE_DIR "${PROJECT_BINARY_DIR}/sandbox")
  file(MAKE_DIRECTORY ${FETCHCONTENT_BASE_DIR})
  set(CPM_SOURCE_CACHE "${PROJECT_BINARY_DIR}/sandbox/archives")
  file(MAKE_DIRECTORY ${CPM_SOURCE_CACHE})
  include(config/CPM.cmake)

  CPMAddPackage("gh:fmtlib/fmt#7.1.3")
endif (ENABLE_CPM)
# ---------------------------

# Check for BLAS and LAPACK availability
if ( ENABLE_BLASLAPACK  )

  # Find BLAS
  FIND_PACKAGE(BLAS REQUIRED)
  SET_VARS_IF_FOUND(BLAS RECOMMENDED "BLAS library" "Native BLAS implementations" BLAS_INCLUDES BLAS_LIBRARIES )

  # Find LAPACK
  FIND_PACKAGE(LAPACK REQUIRED)
  if (LAPACK_FOUND)
    set(ENABLE_FORTRAN ON)
  endif(LAPACK_FOUND)

  SET_VARS_IF_FOUND(LAPACK RECOMMENDED "LAPACK library" "Native LAPACK implementations" LAPACK_INCLUDES LAPACK_LIBRARIES )
  list(APPEND MOAB_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

  LIST(APPEND LOGGING_WATCHED_VARIABLES
    BLAS_LIBRARIES
    LAPACK_LIBRARIES )

endif (ENABLE_BLASLAPACK)

# Check for Eigen3 availability
if ( ENABLE_EIGEN3 )

  # Eigen3 is optional if LAPACK has been found
  FIND_PACKAGE(Eigen3 OPTIONAL_COMPONENTS)
  
  ## Use CPM to auto-download Eigen3 directly
  if ( DOWNLOAD_EIGEN3 AND ENABLE_CPM )
    CPMAddPackage(
      NAME Eigen
      VERSION 3.4.0
      URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
      # Eigen's CMakelists are not intended for library use
      DOWNLOAD_ONLY YES
    )
    
    if(Eigen_ADDED)
      add_library(Eigen INTERFACE IMPORTED)
      target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
    endif()
  endif ( DOWNLOAD_EIGEN3 AND ENABLE_CPM )

  if ( EIGEN3_FOUND )
    INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )
    SET(MOAB_DEFINES "${MOAB_DEFINES} -I${EIGEN3_INCLUDE_DIR}")
  endif (EIGEN3_FOUND)

  SET_VARS_IF_FOUND(EIGEN3 RECOMMENDED "Eigen3" "Eigen3 headers as supplement to BLAS/LAPACK implementations" EIGEN3_INCLUDE_DIR EIGEN3_LIBRARIES )

  LIST(APPEND LOGGING_WATCHED_VARIABLES
    EIGEN3_INCLUDES
    EIGEN3_LIBRARIES )

endif ( ENABLE_EIGEN3 )

# Check for Fortran language support
if ( ENABLE_FORTRAN )
  enable_language(Fortran)

  # Include CMake's Fortran mangling evaluation macros
  include(FortranCInterface)

  if ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "NAME")
  else ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "name")
  endif ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
  if (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL_SUFFIX "## ${FortranCInterface_GLOBAL_SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
  if (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL__SUFFIX "## ${FortranCInterface_GLOBAL__SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")

  set(MOAB_F77_FUNC  "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_F77_FUNC_ "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
  set(MOAB_FC_FUNC   "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_FC_FUNC_  "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
endif ( ENABLE_FORTRAN )

# Check for MPI compiler support and set the necessary flags
if ( ENABLE_MPI )
  set (MPI_C_COMPILER ${CC} )
  set (MPI_CXX_COMPILER ${CXX} )
  set (MPI_Fortran_COMPILER ${FC} )
  find_package( MPI REQUIRED )
  if ( MPI_FOUND )

    set(MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_OPTIONS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    # find_package(MPI) sets MPI_LIBARY, 
    set(MPI_LIBRARIES ${MPI_CXX_LIBRARIES} )
    set(MPI_INCLUDES ${MPI_CXX_HEADER_DIR} ${MPI_CXX_ADDITIONAL_INCLUDE_DIRS} )
    if ( ENABLE_FORTRAN )
      set(MPI_INCLUDES ${MPI_INCLUDES} ${MPI_Fortran_MODULE_DIR} )
    endif ( ENABLE_FORTRAN )

    # Print some MPI configuration details
    SET_VARS_IF_FOUND(MPI RECOMMENDED "MPI" "MPI library to enable node-level parallelism" MPI_INCLUDES MPI_LIBRARIES )
    # Print more MPI flags and variables
    foreach (var MPI_C_COMPILER MPI_CXX_COMPILER MPI_Fortran_COMPILER MPIEXEC_EXECUTABLE MPIEXEC_NUMPROC_FLAG MPI_INCLUDES MPI_LIBRARIES )
      message (STATUS "        ${var}   : ${${var}}")
    endforeach ()

    #set( MOAB_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES} ${MOAB_LIBRARIES} )
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      MPI_DIR
      MPI_C_COMPILER
      MPI_CXX_COMPILER
      MPI_Fortran_COMPILER
      MPIEXEC_EXECUTABLE
      MPIEXEC_NUMPROC_FLAG
      MPIEXEC_MAX_NUMPROCS )

  endif ( MPI_FOUND )
endif ( ENABLE_MPI )

function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
        endif()
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endfunction()

include (config/CheckCompilerFlags.cmake)

# Check for HDF5 availability
# set (HDF5_DIR "$ENV{HDF5_DIR}" CACHE PATH "Path to search for HDF5 header and library files")
set (MOAB_HAVE_HDF5_PARALLEL OFF CACHE INTERNAL "Found necessary parallel HDF5 components. Configure MOAB with it." )
if ( ENABLE_HDF5 )
  if ( MOAB_HAVE_MPI )
    set(HDF5_PREFER_PARALLEL true)
  endif( MOAB_HAVE_MPI )
  set(HDF5_FIND_DEBUG false)

  # Call CMake-native FindHDF5 module and specify the components to check
  find_package( HDF5 REQUIRED COMPONENTS C HL )

  if ( HDF5_FOUND )
    if ( HDF5_IS_PARALLEL )
      set (MOAB_HAVE_HDF5_PARALLEL ON)
      if ( NOT MPI_FOUND )
        message( FATAL_ERROR
            "HDF5 was found with parallel support, but MPI was not found. Please provide MPI_DIR path and reconfigure."
          )
      endif ( NOT MPI_FOUND )
    endif ( HDF5_IS_PARALLEL )

    set(HDF5_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})

    # check for CURL as well
    find_package( CURL )
    if ( CURL_FOUND )
      set( HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CURL_LIBRARIES} )
    endif ( CURL_FOUND )

    include_directories( src/io/mhdf/include )

    # Print some HDF5 configuration details
    SET_VARS_IF_FOUND(HDF5 RECOMMENDED "HDF5" "HDF5 library interfaces for MOAB native file format support" HDF5_INCLUDE_DIRS HDF5_LIBRARIES )
    message (STATUS "        Parallel     : ${HDF5_IS_PARALLEL}")

    LIST(APPEND LOGGING_WATCHED_VARIABLES
      HDF5_DIR
      HDF5_INCLUDE_DIRS
      HDF5_LIBRARIES )

  endif ( HDF5_FOUND )

endif ( ENABLE_HDF5 )

# Check for PNetCDF availability
if ( ENABLE_PNETCDF )
	find_package( PNETCDF REQUIRED )

  # Print some PNetCDF configuration details
  SET_VARS_IF_FOUND(PNETCDF RECOMMENDED "PNetCDF" "PNetCDF library interfaces for MOAB nc-file format support" PNETCDF_INCLUDES PNETCDF_LIBRARIES )

  if (MOAB_HAVE_PNETCDF)
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      PNETCDF_DIR
      PNETCDF_INCLUDES
      PNETCDF_LIBRARIES )
  endif (MOAB_HAVE_PNETCDF)

endif ( ENABLE_PNETCDF )

# Check for NetCDF availability
if ( ENABLE_NETCDF )
  find_package( NETCDF REQUIRED )

  if (MOAB_HAVE_NETCDF)
    set (MOAB_HAVE_NETCDF_H 1)
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      NETCDF_DIR
      NETCDF_INCLUDES
      NETCDF_LIBRARIES )
  endif (MOAB_HAVE_NETCDF)

  # Print some NetCDF configuration details
  SET_VARS_IF_FOUND(NETCDF RECOMMENDED "NetCDF" "NetCDF library interfaces for MOAB ExodusII file format support" NETCDF_INCLUDES NETCDF_LIBRARIES )

endif ( ENABLE_NETCDF )

# Check for Metis partitioner availability
if ( ENABLE_METIS )
  find_package( METIS REQUIRED )

  # Print some Metis configuration details
  SET_VARS_IF_FOUND(METIS RECOMMENDED "Metis" "Metis library interfaces for MOAB offline (serial) mesh partitioning support" METIS_INCLUDES METIS_LIBRARIES )

  if (MOAB_HAVE_METIS)
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      METIS_DIR
      METIS_INCLUDES
      METIS_LIBRARIES )
  endif (MOAB_HAVE_METIS)
endif (ENABLE_METIS )

# Check for ParMetis partitioner availability
if ( ENABLE_PARMETIS AND METIS_FOUND )
  find_package( PARMETIS REQUIRED )

  # Print some ParMetis configuration details
  SET_VARS_IF_FOUND(PARMETIS OPTIONAL "ParMetis" "ParMetis library interfaces for MOAB offline and online parallel mesh partitioning support" PARMETIS_INCLUDES PARMETIS_LIBRARIES )

  if (MOAB_HAVE_PARMETIS)
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      PARMETIS_DIR
      PARMETIS_INCLUDES
      PARMETIS_LIBRARIES )
  endif (MOAB_HAVE_PARMETIS)
endif ( ENABLE_PARMETIS AND METIS_FOUND )

# Check for Zoltan partitioner availability
if ( ENABLE_ZOLTAN )
  find_package( ZOLTAN REQUIRED )

  # Print some Zoltan configuration details
  SET_VARS_IF_FOUND(ZOLTAN RECOMMENDED "Zoltan" "Zoltan library interfaces for MOAB offline/online parallel mesh partitioning support" ZOLTAN_INCLUDES ZOLTAN_LIBRARIES )

  if (MOAB_HAVE_ZOLTAN)
    LIST(APPEND LOGGING_WATCHED_VARIABLES
      ZOLTAN_DIR
      ZOLTAN_INCLUDES
      ZOLTAN_LIBRARIES )
  endif (MOAB_HAVE_ZOLTAN)
endif (ENABLE_ZOLTAN )

# Check for TempestRemap availability if Eigen3 headers have been found already
if ( ENABLE_TEMPESTREMAP )
  if (NOT EIGEN3_FOUND)
    message("WARNING: TempestRemap cannot be configured without Eigen3 support for MOAB. Skipping...")
  else (NOT EIGEN3_FOUND)
    find_package( TEMPESTREMAP REQUIRED )

    # Print some TempestRemap configuration details
    SET_VARS_IF_FOUND(TEMPESTREMAP OPTIONAL "TempestRemap" "TempestRemap library interfaces for MOAB climate remapping support" TEMPESTREMAP_INCLUDES TEMPESTREMAP_LIBRARIES )

    if (MOAB_HAVE_TEMPESTREMAP)
      LIST(APPEND LOGGING_WATCHED_VARIABLES
        TEMPESTREMAP_DIR
        TEMPESTREMAP_INCLUDES
        TEMPESTREMAP_LIBRARIES )

      # Make sure the directory exists
      set(EXT_DATA_DIR "${CMAKE_BINARY_DIR}/sandbox/MeshFiles")
      set(EXT_E3SM_DATA_DIR "${CMAKE_BINARY_DIR}/sandbox/MeshFiles/e3sm")
      file(MAKE_DIRECTORY "${EXT_E3SM_DATA_DIR}")
      file(MAKE_DIRECTORY "${EXT_E3SM_DATA_DIR}/ne4pg2_o240")
      file(MAKE_DIRECTORY "${EXT_E3SM_DATA_DIR}/ne30pg2_r05_oECv3_ICG")
      set(MOAB_DEFINES "${MOAB_DEFINES} -DEXT_DATA_DIR=${EXT_DATA_DIR}/")

      # case: ne4pg2_o240
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne4pg2_o240/AtmPhys_pg2.h5m "${EXT_E3SM_DATA_DIR}/ne4pg2_o240/AtmPhys_pg2.h5m" )
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne4pg2_o240/land_p8.h5m "${EXT_E3SM_DATA_DIR}/ne4pg2_o240/land_p8.h5m" )
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne4pg2_o240/ne4pg2_p8.h5m "${EXT_E3SM_DATA_DIR}/ne4pg2_o240/ne4pg2_p8.h5m" )

      # case: ne30pg2_r05_oECv3_ICG
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne30pg2_r05_oECv3_ICG/AtmPhys_01.h5m "${EXT_E3SM_DATA_DIR}/ne30pg2_r05_oECv3_ICG/AtmPhys_01.h5m" )
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne30pg2_r05_oECv3_ICG/domain_clean_p16.h5m "${EXT_E3SM_DATA_DIR}/ne30pg2_r05_oECv3_ICG/domain_clean_p16.h5m" )
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne30pg2_r05_oECv3_ICG/outCS30_pg2_p8.h5m "${EXT_E3SM_DATA_DIR}/ne30pg2_r05_oECv3_ICG/outCS30_pg2_p8.h5m" )
      file(DOWNLOAD https://web.cels.anl.gov/projects/sigma/downloads/MeshFiles/e3sm/ne30pg2_r05_oECv3_ICG/wholeMPAS.h5m "${EXT_E3SM_DATA_DIR}/ne30pg2_r05_oECv3_ICG/wholeMPAS.h5m" )
    endif(MOAB_HAVE_TEMPESTREMAP)
  endif (NOT EIGEN3_FOUND)
endif( ENABLE_TEMPESTREMAP )

if (ENABLE_CGNS)
  set( MOABIO_LIBS ${MOABIO_LIBS} ${CGNS_LIBRARIES} )
  set (MOAB_HAVE_CGNS ON)
endif(ENABLE_CGNS)

set(tmp_libs ${MOAB_LIBRARIES})
string(REPLACE ";" " " tmp_libs "${MOAB_LIBRARIES}")
RESOLVE_LIBRARIES(MOAB_LIBRARIES "${tmp_libs}")
list (REVERSE MOAB_LIBRARIES)
list (REMOVE_DUPLICATES MOAB_LIBRARIES)
list (REVERSE MOAB_LIBRARIES)

set(tmp_incs ${MOAB_INCLUDES})
string(REPLACE ";" " " tmp_incs "${MOAB_INCLUDES}")
RESOLVE_INCLUDES(MOAB_INCLUDES "${tmp_incs}")
list (REMOVE_DUPLICATES MOAB_INCLUDES)

# Strip our custom flags
string(STRIP "${CMAKE_EXE_LINKER_FLAGS} " CMAKE_EXE_LINKER_FLAGS)
string(STRIP "${MOAB_CXX_FLAGS} " MOAB_CXX_FLAGS)
string(STRIP "${MOAB_DEFINES} " MOAB_DEFINES)

###############################################################################
# Python Module Build
###############################################################################
if(SKBUILD OR ENABLE_PYMOAB)

  # Find necessary packages
  find_package(
      Python
      COMPONENTS Interpreter Development.Module NumPy
      REQUIRED
  )
  if(SKBUILD)
    # By default, scikit-build will install everything to ${SKBUILD_PLATLIB_DIR}/pymoab
    # Install include and bin directories to ${SKBUILD_PLATLIB_DIR}/pymoab/$SKBUILD_SUBDIR
    set(SKBUILD_SUBDIR core/)
    set(CMAKE_INSTALL_BINDIR ${SKBUILD_SUBDIR}${CMAKE_INSTALL_BINDIR})
    set(CMAKE_INSTALL_INCLUDEDIR ${SKBUILD_SUBDIR}${CMAKE_INSTALL_INCLUDEDIR})
    set(CMAKE_INSTALL_DATADIR ${SKBUILD_DATA_DIR}/${CMAKE_INSTALL_DATADIR})
    set(CMAKE_INSTALL_DOCDIR ${SKBUILD_DATA_DIR}/${CMAKE_INSTALL_DOCDIR})
    set(CMAKE_INSTALL_MANDIR ${SKBUILD_DATA_DIR}/${CMAKE_INSTALL_MANDIR})
    set(MOABDocLocation ${SKBUILD_DATA_DIR}/${MOABDocLocation})

    # Set the RPATH
    if(APPLE)
      set(PYMOAB_LIBRARY_RPATH "@loader_path")
      set(PYMOAB_PYTHON_MODULE_RPATH "@loader_path/${SKBUILD_SUBDIR}${CMAKE_INSTALL_LIBDIR}")
      set(PYMOAB_BINARY_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    elseif(UNIX)
      set(PYMOAB_LIBRARY_RPATH "$ORIGIN")
      set(PYMOAB_PYTHON_MODULE_RPATH "$ORIGIN/${SKBUILD_SUBDIR}${CMAKE_INSTALL_LIBDIR}")
      set(PYMOAB_BINARY_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()

    # Add the pymoab subdirectory
    add_subdirectory(pymoab)
  else()
    # Get a list of all Python files in the pymoab directory
    file(GLOB PY_FILES "${CMAKE_SOURCE_DIR}/pymoab/*.py")
    file(GLOB PXD_FILES "${CMAKE_SOURCE_DIR}/pymoab/*.pxd")
    file(GLOB PYX_FILES "${CMAKE_SOURCE_DIR}/pymoab/*.pyx")

    # Set the PyMOAB prefix
    if (NOT PYMOAB_PREFIX)
        set(PYMOAB_PREFIX ${CMAKE_INSTALL_PREFIX})
    endif()

    # Show PyMOAB prefix
    message(STATUS "PyMOAB prefix: ${PYMOAB_PREFIX}")

    # Copy the Python files to the current binary directory
    file( COPY ${PY_FILES} ${PXD_FILES} ${PYX_FILES} DESTINATION ${CMAKE_BINARY_DIR}/pymoab)

    # Copy the pyproject.toml file to the binary directory
    file( COPY ${CMAKE_SOURCE_DIR}/deprecated/pymoab/pyproject.toml DESTINATION ${CMAKE_BINARY_DIR})

    # Configure the setup.py file
    configure_file(${CMAKE_SOURCE_DIR}/deprecated/pymoab/setup.py.cmake ${CMAKE_BINARY_DIR}/setup.py)

    # Install the PyMOAB package for testing
    add_custom_target(
      pymoab-test-install
      DEPENDS MOAB
    )
    add_custom_command(
      TARGET pymoab-test-install
      COMMAND ${Python_EXECUTABLE} -m pip install . --target ${CMAKE_BINARY_DIR}/Testing
      COMMENT "Installing PyMOAB package for testing..."
    )
    add_test(NAME pymoab_test COMMAND ${CMAKE_COMMAND} --build . --target pymoab-test-install)

    # Install the PyMOAB package
    install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} -m pip install . --prefix ${PYMOAB_PREFIX})")

    # Print a deprecation warning
    message(WARNING "ENABLE_PYMOAB and PYMOAB_PREFIX is deprecated and will be removed in the next release. "
      "Consider installing PyMOAB via pip instead. "
      "See README.md for more information."
      )
  endif()
endif()

################################################################################
# Add Directories
################################################################################
add_subdirectory( src )
if (NOT WIN32 OR MSYS OR MINGW) # Need further work to prepare for windows
  add_subdirectory( tools )
endif ()

################################################################################
# Testing Related Settings
################################################################################
#turn on ctest if we want testing
if ( ENABLE_TESTING )
  enable_testing()
  add_subdirectory( test )
endif()

###############################################################################
# Generate the MOABTargets.cmake file
###############################################################################
export(
  TARGETS ${MOAB_INSTALL_TARGETS}
  FILE "${PROJECT_BINARY_DIR}/MOABTargets.cmake"
  )
set(MOABTargetLocation ${SKBUILD_SUBDIR}${CMAKE_INSTALL_LIBDIR}/cmake/MOAB)
install(
  EXPORT MOABTargets
  DESTINATION ${MOABTargetLocation}
)

# For consistency with autoconf, create a "make check" equivalent to "ctest"
ADD_CUSTOM_TARGET( check COMMAND ctest DEPENDS )

################################
# PkgConfig
################################
find_package(PkgConfig QUIET)

# For consistency with autoconf, create "make dist" and "make distcheck" targets
# INCLUDE( CMakeDistCheck )
# DEFINE_DISTCHECK ( )
################################
# MakeDist
################################
# Enable 'dist and distcheck'.
# File adapted from http://ensc.de/cmake/FindMakeDist.cmake
#find_package(MakeDist)
# End 'enable dist and distcheck'


################################################################################
# Generate the MOABConfig.cmake file
################################################################################
set(CXX ${CMAKE_CXX_COMPILER})
set(CC ${CMAKE_C_COMPILER})
if (ENABLE_FORTRAN)
  set(F77 ${CMAKE_Fortran_COMPILER} )
  set(FC ${CMAKE_Fortran_COMPILER} )
  LIST(APPEND LOGGING_WATCHED_VARIABLES
    CMAKE_Fortran_COMPILER )
endif (ENABLE_FORTRAN)

set(CPPFLAGS "${CMAKE_CXX_FLAGS} ${MOAB_DEFINES}")
IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
  set(CFLAGS ${CMAKE_C_FLAGS_DEBUG})
  if (ENABLE_FORTRAN)
    set(FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
  endif (ENABLE_FORTRAN)
ELSE (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
  set(CFLAGS ${CMAKE_C_FLAGS_RELEASE})
  if (ENABLE_FORTRAN)
    set(FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  endif (ENABLE_FORTRAN)
ENDIF (CMAKE_BUILD_TYPE MATCHES "Debug")

LIST(APPEND LOGGING_WATCHED_VARIABLES
    MOAB_DEFINES
    MOAB_CXX_FLAGS
    MOAB_INCLUDES
    MOAB_LIBRARIES
    MOAB_INSTALL_TARGETS
    )

# CMake does not define uninstall target. Look at install manifest and remove
# all the files recursively.
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
    )

## Find and use clang-format if available
find_program( CLANG_FORMAT_EXE NAMES "clang-format" DOC "Path to clang-format executable" )
if(CLANG_FORMAT_EXE)
	message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
	set(DO_CLANG_FORMAT "${CLANG_FORMAT_EXE}" "-i -style=file")
endif()

# get all project files
file(GLOB_RECURSE MOAB_PROJECT_SOURCE_FILES
			RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
			LIST_DIRECTORIES false
			src/*.cpp src/*.c src/*.hpp
			tools/*.cpp tools/*.hpp
			test/*.cpp
			examples/*.cpp
		)

function(prepend var prefix)
  set(listVar "")
  foreach(f ${ARGN})
    list(APPEND listVar "${prefix}/${f}")
  endforeach()
  set(${var} "${listVar}" PARENT_SCOPE)
endfunction()

if(CLANG_FORMAT_EXE)
	prepend(FILES_TO_FORMAT ${CMAKE_CURRENT_SOURCE_DIR} ${MOAB_PROJECT_SOURCE_FILES})
	add_custom_target( clang-format-project-files COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${FILES_TO_FORMAT} )
endif()

## Check for clang-tidy as well for performing pre-commit checks
find_program(CLANG_TIDY_EXE NAMES "clang-tidy" DOC "Path to clang-tidy executable")
if(CLANG_TIDY_EXE)
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  SET(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE}
    -format-style=file
    -checks=clan*,mpi-*,port*,bug*,llvm-namespace-comment,perf*,modernize-deprecated-headers,modernize-avoid-bind,modernize-use-bool-literals,misc-noexcept-move-constructor,-readability-braces-around-statements,-misc-unused-parameters
    -header-filter="\"moab/*\""
    -fix -fix-errors
    -p=${CMAKE_CURRENT_BINARY_DIR} )

  find_program(RUN_CLANG_TIDY_BIN /usr/local/Cellar/llvm/10.0.0_3/share/clang/run-clang-tidy.py)
  find_program(RUN_CLANG_APPLY_REPLACEMENTS_BIN clang-apply-replacements PATH /usr/local/bin /usr/local/opt/llvm/bin)
  list(APPEND RUN_CLANG_TIDY_BIN_ARGS
        -clang-tidy-binary ${CLANG_TIDY_EXE}
        -clang-apply-replacements-binary ${RUN_CLANG_APPLY_REPLACEMENTS_BIN}
        -header-filter="\"moab/*\""
        -checks=clan*,mpi-*,port*,bug*,llvm-namespace-comment,perf*,modernize-deprecated-headers,modernize-avoid-bind,modernize-use-bool-literals,misc-noexcept-move-constructor,-readability-braces-around-statements,-misc-unused-parameters
        -fix
        )

	add_custom_target(
		tidy
		COMMAND ${RUN_CLANG_TIDY_BIN} ${RUN_CLANG_TIDY_BIN_ARGS}
		COMMENT "running clang tidy"
		)
	add_custom_target( clang-tidy-project-files COMMAND ${CMAKE_CXX_CLANG_TIDY} ${FILES_TO_FORMAT} )

endif(CLANG_TIDY_EXE)

# Generate the CMake input header file and then pass it through configuration
## execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -e "s/#undef /#cmakedefine MOAB_/g"
#execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -r "s/#undef ([A-Z].*)/#cmakedefine MOAB_\\1 @MOAB_\\1@/"
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/config
#  OUTPUT_FILE MOABConfig.cmake2.h.in
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#  )
CONFIGURE_FILE(config/MOABConfig.cmake.h.in "${PROJECT_BINARY_DIR}/src/moab/MOABConfig.h" @ONLY)
CONFIGURE_FILE(config/MOABConfig.cmake.cmake "${PROJECT_BINARY_DIR}/MOABConfig.cmake" @ONLY)
CONFIGURE_FILE(config/MOABConfigVersion.cmake.cmake "${PROJECT_BINARY_DIR}/MOABConfigVersion.cmake" @ONLY)
CONFIGURE_FILE(moab.make.cmake "${PROJECT_BINARY_DIR}/lib/moab.make" @ONLY)
CONFIGURE_FILE(config/MOAB.pc.cmake.in "${PROJECT_BINARY_DIR}/MOAB.pc" @ONLY)

# Install the CMake configuration files
set(MOABMakeLocation ${CMAKE_INSTALL_LIBDIR})
set(MOABConfigLocation ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB)
set(MOABPKGCONFIG_LOCATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Skip installing the make and pkgconfig files if SKBUILD is set
if(NOT SKBUILD)
  install(FILES "${PROJECT_BINARY_DIR}/lib/moab.make"
    DESTINATION ${MOABMakeLocation}
    )
  install(FILES ${PROJECT_BINARY_DIR}/MOAB.pc
    DESTINATION ${MOABPKGCONFIG_LOCATION}
    )
  install(CODE
    "file(APPEND  ${CMAKE_INSTALL_LIBDIR}/moab.make MOAB_LIBDIR=${CMAKE_INSTALL_FULL_LIBDIR}\\n )"
    )
  install(CODE
    "file(APPEND  ${CMAKE_INSTALL_LIBDIR}/moab.make MOAB_INCLUDES=-I${CMAKE_INSTALL_FULL_INCLUDEDIR}\\n )"
    )
else()
  set(MOABConfigLocation ${SKBUILD_DATA_DIR}/${MOABConfigLocation})
endif()

install(
  FILES "${PROJECT_BINARY_DIR}/MOABConfig.cmake" "${PROJECT_BINARY_DIR}/MOABConfigVersion.cmake"
  DESTINATION ${MOABConfigLocation}
  )
install(
  FILES "${PROJECT_SOURCE_DIR}/config/ResolveCompilerPaths.cmake"
  DESTINATION ${MOABConfigLocation}
  )

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )
###################
# INSTALL documentation and examples
####################

SET(MOAB_DOCS ${MOAB_DOCS}
  LICENSE README.md RELEASE_NOTES.md TODO.md
)
install(FILES ${MOAB_DOCS} DESTINATION ${MOABDocLocation})

CONFIGURE_FILE(examples/makefile.config.cmake "${PROJECT_BINARY_DIR}/share/examples/makefile.config" @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/share/examples/makefile.config
  DESTINATION ${CMAKE_INSTALL_DATADIR}/examples)

SET(MOAB_EXAMPLES
  examples/makefile
  examples/makefile.config.in
  examples/CMakeLists.txt )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/beginner/README
  examples/beginner/makefile
  examples/beginner/mbex1.cpp
  examples/beginner/mbex2.cpp
  examples/beginner/mbex3.cpp
  examples/beginner/mbex4.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/basic/makefile
  examples/basic/ErrorHandlingSimulation.cpp
  examples/basic/GetEntities.cpp
  examples/basic/HelloMOAB.cpp
  examples/basic/HelloParMOAB.cpp
  examples/basic/ReadWriteTest.cpp
  examples/basic/SetsNTags.cpp
  examples/basic/StructuredMeshSimple.cpp
  examples/basic/TestErrorHandling.cpp
  examples/basic/TestErrorHandlingPar.cpp
  examples/basic/TestExodusII.cpp
  examples/basic/UniformRefinement.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/intermediate/makefile
  examples/intermediate/DirectAccessNoHoles.cpp
  examples/intermediate/DirectAccessWithHoles.cpp
  examples/intermediate/LoadPartial.cpp
  examples/intermediate/PointInElementSearch.cpp
  examples/intermediate/ReduceExchangeTags.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/advanced/makefile
  examples/advanced/CrystalRouterExample.cpp
  examples/advanced/ExtrudePoly.cpp
  examples/advanced/GenLargeMesh.cpp
  examples/advanced/LloydRelaxation.cpp
  examples/advanced/VisTags.cpp
  examples/advanced/CoupleMGen.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/fortran/makefile
  examples/fortran/DirectAccessNoHolesF90.F90
  examples/fortran/PushParMeshIntoMoabF90.F90
  examples/fortran/MigrateMesh.F90 )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/python/makefile
  examples/python/laplaciansmoother.py )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/apps/makefile )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/apps/earth/makefile
  examples/apps/earth/ContinentsOnGlobe.cpp
  examples/apps/earth/BoundaryContinents.cpp
  examples/apps/earth/addPCdata.cpp )

include(AutoconfHeader)
moab_install_docs(${MOAB_EXAMPLES})

feature_summary(WHAT ALL #ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND
                DESCRIPTION "MOAB Configuration Summary"
                INCLUDE_QUIET_PACKAGES
                FILENAME ${CMAKE_BINARY_DIR}/summary.log
                )
MESSAGE("-- MOAB Configuration summary has been written to ${CMAKE_BINARY_DIR}/summary.log.")

LOGGING_FINALIZE()
# All done.

# dump_cmake_variables("MPI_")
# dump_cmake_variables("HDF5")
