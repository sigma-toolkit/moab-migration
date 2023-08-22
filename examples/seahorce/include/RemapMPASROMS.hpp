#ifndef __remap_mpas_roms_hpp__
#define __remap_mpas_roms_hpp__

#include "moab/MOABConfig.h"

// 3D settings
constexpr int mpas_zreflevels = 60;
constexpr int mpas_zlevels    = 60;
constexpr int roms_zlevels    = 60;
constexpr int nvars           = 2;

// tag name data
const char* mpas_twod_tagnames[nvars]       = { "salinity", "temperature" };
const char* mpas_threed_cum_tagnames[nvars] = { "salinity_3d", "temperature_3d" };
const char* mpas_threed_tagnames[nvars]     = { "Salinity3d", "Temperature3d" };
const char* roms_twod_tagnames[nvars]       = { "Salinity2DROMS", "Temperature2DROMS" };
const char* roms_threed_tagnames[nvars]     = { "Salinity3dROMS", "Temperature3dROMS" };

// Error routines for use with MPI API
#define MPICHKERR( CODE, MSG )                 \
    do                                         \
    {                                          \
        if( 0 != ( CODE ) )                    \
        {                                      \
            std::cout << ( MSG ) << std::endl; \
            MPI_Finalize();                    \
        }                                      \
    } while( false )

#define dbgprint( MSG )                \
    do                                 \
    {                                  \
        std::cout << MSG << std::endl; \
    } while( false )

#define dbgprintall( MSG )                                     \
    do                                                         \
    {                                                          \
        std::cout << "[" << rank << "]: " << MSG << std::endl; \
    } while( false )

#endif  // __remap_mpas_roms_hpp__