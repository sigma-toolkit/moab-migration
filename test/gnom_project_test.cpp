/*
 * gnom_project_test.cpp
 * will test new global gnomonic projection method, to be used by zoltan for partitioning
 *
 */
#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/climate/intx_mesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"
#include "TestUtil.hpp"

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{
    string filein  = STRINGIFY( MESHDIR ) "/mbcslam/eulerHomme.vtk";
    string fileout = "project.vtk";

    ProgOptions opts;
    opts.addOpt< std::string >( "model,m", "input file ", &filein );

    opts.addOpt< std::string >( "output,o", "output filename", &fileout );

    opts.parseCommandLine( argc, argv );

    Core moab;
    Interface* mb = &moab;
    EntityHandle sf;
    MB_CHK_ERR( mb->create_meshset( MESHSET_SET, sf ) );

    MB_CHK_ERR( mb->load_file( filein.c_str(), &sf ) );

    double R = 1.;  // should be input
    EntityHandle outSet;
    MB_CHK_ERR( mb->create_meshset( MESHSET_SET, outSet ) );

    MB_CHK_ERR( IntxUtils::global_gnomonic_projection( mb, sf, R, false, outSet ) );
    MB_CHK_ERR( mb->write_file( fileout.c_str(), 0, 0, &outSet, 1 ) );
    // skip the test for 32 bit builds, CHECK_ARRAYS_EQUAL macro is not working correctly
#ifndef MOAB_FORCE_32_BIT_HANDLES
    // check first cell position
    Range cells;
    MB_CHK_ERR( mb->get_entities_by_dimension( outSet, 2, cells ) );
    EntityHandle firstCell = cells[0];
    double coords[3];
    MB_CHK_ERR( mb->get_coords( &firstCell, 1, coords ) );
    // check position
    double values[3] = { -0.78867513420303437, 0.78867513420303437, 0 };
    CHECK_ARRAYS_EQUAL( coords, 3, values, 3 );
#endif

    return 0;
}
