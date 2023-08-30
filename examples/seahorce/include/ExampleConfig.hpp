#ifndef __example_config_hpp__
#define __example_config_hpp__

// MOAB config include
#include "moab/MOABConfig.h"

// 3D settings
constexpr int mpas_zreflevels  = 60;
constexpr int mpas_zlevels     = 60;
constexpr int roms_zlevels     = 100;
constexpr int nvars            = 2;
constexpr double axial_scaling = 1.0E4;

// tag name data
const char* mpas_twod_tagnames[nvars]       = { "salinity", "temperature" };
const char* mpas_threed_cum_tagnames[nvars] = { "salinity_3d", "temperature_3d" };
const char* mpas_threed_tagnames[nvars]     = { "Salinity3d", "Temperature3d" };
const char* roms_twod_tagnames[nvars]       = { "Salinity2DROMS", "Temperature2DROMS" };
const char* roms_threed_tagnames[nvars]     = { "Salinity3dROMS", "Temperature3dROMS" };

// write the map file to disk; comment out to just compute in-memory
#define VERTICAL_INTERPOLATION
// #define VERTICAL_INTERPOLANT_LINEAR
#define WRITE_MAP_FILE
// #define VERBOSE_OUTPUT

#endif // __example_config_hpp__