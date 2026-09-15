/*
 * gnom_project_test.cpp
 * will test new global gnomonic projection method, to be used by zoltan for partitioning
 *
 */
#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/earthsystem/intx_mesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{
    string filein  = "target_1.h5m";
    string fileout = "project.h5m";

    ProgOptions opts;
    opts.addOpt< std::string >( "model,m", "input file ", &filein );

    opts.addOpt< std::string >( "output,o", "output filename", &fileout );

    opts.parseCommandLine( argc, argv );

    Core moab;
    Interface* mb = &moab;
    EntityHandle sf;
    ErrorCode rval = mb->create_meshset( MESHSET_SET, sf );MB_CHK_ERR( rval );

    rval = mb->load_file( filein.c_str(), &sf );MB_CHK_ERR( rval );

    EntityHandle outSet;
    rval = mb->create_meshset( MESHSET_SET, outSet );MB_CHK_ERR( rval );
    // get the coords of the first vertex entity handle
    EntityHandle v1 = 1;  // we know it must exist; we could get a specific one too :)

    CartVect P;
    rval = mb->get_coords( &v1, 1, P.array() );MB_CHK_ERR( rval );

    rval = IntxUtils::global_gnomonic_projection_general( mb, sf, P, outSet );MB_CHK_ERR( rval );
    std::cout << "writing output file " << fileout.c_str() << "\n";
    rval = mb->write_file( fileout.c_str(), 0, 0, &outSet, 1 );MB_CHK_ERR( rval );

    return 0;
}
