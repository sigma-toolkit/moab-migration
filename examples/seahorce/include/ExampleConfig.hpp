#ifndef __example_config_hpp__
#define __example_config_hpp__

// MOAB config include
#include "moab/MOABConfig.h"

// 3D settings
constexpr int mpas_zreflevels = 80;
constexpr int mpas_zlevels    = 80;
constexpr int roms_zlevels    = 100;
constexpr int nvars           = 4;
// constexpr double axial_scaling = 100000.0;
// constexpr double axial_scaling = 6371220.0;
constexpr double axial_scaling = 1.0;

// Fields that are currently projected
enum Fields
{
    Bathymetry,
    Salinity,
    Temperature,
    SeaSurfaceHeight,
    VelocityX,
    VelocityY
};

//
const char* mpas_twod_standardtagnames[2] = { "bottomDepth", "timeDaily_avg_ssh" };
const char* roms_twod_standardtagnames[2] = { "Bathymetry", "SSH" };

// tag name data
// const char* mpas_twod_tagnames[4] = { "salinity", "temperature", "VX", "VY" };
// const char* mpas_threed_cum_tagnames[4] = { "salinity_3d", "temperature_3d", "velocityX",
//                                                 "velocityY" };
const char* mpas_tagnames[4] = { "timeDaily_avg_activeTracers_salinity_3d",
                                 "timeDaily_avg_activeTracers_temperature_3d", "timeDaily_avg_velocityMeridional_3d",
                                 "timeDaily_avg_velocityZonal_3d" };
const char* mpas_ele_tagnames[4] = { "MPAS_Salinity", "MPAS_Temperature", "MPAS_VelMeridional", "MPAS_VelZonal" };
const char* roms_tagnames[4]     = { "ROMS_Salinity", "ROMS_Temperature", "ROMS_VelMeridional", "ROMS_VelZonal" };

// write the map file to disk; comment out to just compute in-memory
#define VERTICAL_INTERPOLATION
// #define VERTICAL_INTERPOLANT_LINEAR
#define WRITE_MAP_FILE
// ROMS stretching function vs constant delta
#define USE_STRETCHING_FUNCTION
// Use the centroid of the polyhedra to generate tetrahedra
#define USE_DUAL_TETS
// Maximum number of vertex to elemetn adjacencies (ifndef USE_DUAL_TETS)
#define N_MAX_V2EADJACENCIES 128
// Output debug files
#define VERBOSE_OUTPUT

#endif  // __example_config_hpp__