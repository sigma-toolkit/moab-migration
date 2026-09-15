# EarthSystem NetCDF readers/writers (CAM Euler/FV/HOMME, MPAS, GCRM, SCRIP,
# ESMF and domain grids).  These are the Earth-system-specific half of the old
# src/io directory; the ExodusII NetCDF reader/writer (ReadNCDF/WriteNCDF)
# is general-purpose and stays in src/io.
#
# Built only when NetCDF (or PNetCDF alone) is available, mirroring the
# gating that previously lived in src/io/CMakeLists.txt.
include_directories( ${MOAB_BINARY_DIR}/src
                     ${MOAB_SOURCE_DIR}/src
                     ${MOAB_SOURCE_DIR}/src/io
                     ${MOAB_SOURCE_DIR}/src/earthsystem/io
                     ${MOAB_BINARY_DIR}/src/parallel
                    )

set( MOAB_EARTHSYSTEM_IO_SRCS
        ReadNC.hpp        ReadNC.cpp
        WriteNC.hpp       WriteNC.cpp
        NCHelper.hpp      NCHelper.cpp
        NCWriteHelper.hpp NCWriteHelper.cpp
        NCHelperEuler.hpp NCHelperEuler.cpp
        NCWriteEuler.hpp  NCWriteEuler.cpp
        NCHelperFV.hpp    NCHelperFV.cpp
        NCWriteFV.hpp     NCWriteFV.cpp
        NCHelperDomain.hpp NCHelperDomain.cpp
        NCHelperScrip.hpp NCHelperScrip.cpp
        NCHelperESMF.hpp  NCHelperESMF.cpp
        NCHelperGCRM.hpp  NCHelperGCRM.cpp
        NCWriteGCRM.hpp   NCWriteGCRM.cpp
        NCHelperHOMME.hpp NCHelperHOMME.cpp
        NCWriteHOMME.hpp  NCWriteHOMME.cpp
        NCHelperMPAS.hpp  NCHelperMPAS.cpp
        NCWriteMPAS.hpp   NCWriteMPAS.cpp
   )

if( MOAB_HAVE_NETCDF )
  include_directories(${NETCDF_INCLUDES})
endif()
if( MOAB_HAVE_PNETCDF )
  include_directories(${PNETCDF_INCLUDES})
endif()
if( MOAB_HAVE_HDF5 )
  include_directories(${HDF5_INCLUDES})
endif()
if( MOAB_HAVE_ZOLTAN )
  include_directories(${ZOLTAN_INCLUDES})
endif()

add_library( moab_earthsystem_io OBJECT ${MOAB_EARTHSYSTEM_IO_SRCS} )
set_target_properties( moab_earthsystem_io PROPERTIES
                       COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
