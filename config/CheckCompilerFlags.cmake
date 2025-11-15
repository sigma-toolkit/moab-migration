
# Include ability to add flags unless already present
include (config/ForceAddFlags.cmake)
include(CheckFortranCompilerFlag)

SET(MOAB_CXX_FLAGS "")

#
# Set -pedantic if the compiler supports it.
#
IF(NOT (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
        CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.4"))
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-pedantic")
ENDIF()

ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fp-model=precise")
if (ENABLE_FORTRAN)
  ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS "-fp-model source")
endif (ENABLE_FORTRAN)
# Check for compiler types and add flags accordingly
if ( CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )

  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fpic")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wall")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-pipe")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-long-long")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wextra")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wshadow")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-cast-align")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wsign-compare")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wpointer-arith")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wformat")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wformat-security")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wunused-parameter")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fstack-protector-all")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fpermissive")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fsignaling-nans")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-ftrapping-math")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-fnon-call-exceptions")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-ignored-attributes")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-variadic-macros")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-deprecated-declarations")
  ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-shadow")
  # Need to enable or check for this only if user asks for C++11 support
  # ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-c++11-long-long")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8)
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-unused-local-typedefs")
  endif()
  if(APPLE) # Clang / Mac OS only
    # Required on OSX to compile c++11
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-stdlib=libc++")
  endif()

  if (ENABLE_FORTRAN)
    # gfortran
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops;-fno-f2c")
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS_DEBUG "-ffpe-trap=invalid,zero,overflow") 
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS "-ffree-line-length-none")
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS "-ffixed-line-length-none")
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")
    ENABLE_IF_SUPPORTED_FC(CMAKE_Fortran_FLAGS "-fconvert=big-endian")
  endif (ENABLE_FORTRAN)
else ( CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-ansi")  # standard compliant
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-w2")    # verbose warnings
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wall")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-Wno-shadow")
    # disable some warnings -- consistent with autoconf
    # -wd981 -wd279 -wd1418 -wd383 -wd1572 -wd2259
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd981")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd279")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd1418")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd1572")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd383")
    ENABLE_IF_SUPPORTED(MOAB_CXX_FLAGS "-wd2259")

    if (ENABLE_FORTRAN)
      # ifort (untested)
      set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl")
      set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl")
    endif (ENABLE_FORTRAN)
  endif()
endif ( CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )

# Debug targets
IF (CMAKE_BUILD_TYPE MATCHES "Debug")

  if(NOT WIN32)

    ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "-ggdb")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-ggdb")
    ENABLE_IF_SUPPORTED(CMAKE_EXE_LINKER_FLAGS "-ggdb")
    #
    # If -ggdb is not available, fall back to -g:
    #
    IF(NOT CMAKE_HAVE_FLAG_ggdb)
      ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "-g")
      ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-g")
      if (ENABLE_FORTRAN)
        FORCE_ADD_FLAGS(CMAKE_Fortran_FLAGS "-g")
      endif (ENABLE_FORTRAN)
      ENABLE_IF_SUPPORTED(CMAKE_EXE_LINKER_FLAGS "-g")
    ENDIF()

    IF(CMAKE_SETUP_COVERAGE)
      #
      # Enable test coverage
      #
      ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fno-elide-constructors")
      ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-ftest-coverage -fprofile-arcs")
      ENABLE_IF_SUPPORTED(CMAKE_EXE_LINKER_FLAGS "-ftest-coverage -fprofile-arcs")
    ENDIF()
  else(NOT WIN32)
    ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "/Od /w")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "/Od /w")
  endif(NOT WIN32)

ENDIF()

FORCE_ADD_FLAGS(CMAKE_CXX_FLAGS "${MOAB_CXX_FLAGS}")

# Force C99 standard support
FORCE_ADD_FLAGS(CMAKE_C_FLAGS "-std=c99")

# Release targets
IF (CMAKE_BUILD_TYPE MATCHES "Release")
  #
  # General optimization flags:
  #
  if(NOT WIN32)

    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-ip")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-funroll-loops")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-funroll-all-loops")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-fstrict-aliasing")

    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused")
  else(NOT WIN32)
    ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "/O2 /w")
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "/O2 /w")
  endif(NOT WIN32)

ENDIF()

mark_as_advanced(CMAKE_Fortran_FLAGS CMAKE_C_FLAGS CMAKE_CXX_FLAGS MOAB_CXX_FLAGS)
