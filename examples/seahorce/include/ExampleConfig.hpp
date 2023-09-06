#ifndef __example_config_hpp__
#define __example_config_hpp__

// MOAB config include
#include "moab/MOABConfig.h"

// 3D settings
constexpr int mpas_zreflevels = 60;
constexpr int mpas_zlevels    = 60;
constexpr int roms_zlevels    = 120;
constexpr int nvars           = 2;
// constexpr double axial_scaling = 100000.0;
// constexpr double axial_scaling = 6371220.0;
constexpr double axial_scaling = 1.0;

// Fields that are currently projected
enum Fields
{
    Bathymetry,
    Salinity,
    Temperature
};

// tag name data
const char* mpas_twod_tagnames[nvars]       = { "salinity", "temperature" };
const char* mpas_threed_cum_tagnames[nvars] = { "salinity_3d", "temperature_3d" };
const char* mpas_threed_tagnames[nvars]     = { "Salinity3d", "Temperature3d" };
const char* roms_twod_tagnames[nvars]       = { "Salinity", "Temperature" };
const char* roms_threed_tagnames[nvars]     = { "Salinity", "Temperature" };

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