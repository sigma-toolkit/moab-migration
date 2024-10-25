/** \brief This example shows how to change coords for mali land ice.
 *
 * MPAS file has z coords 0, so we want to use lat and lon on Vertex tags to change it
 * to unit sphere; compare then with the h5m converted from scrip file (by Vijay)
 * If they are on top of each other, we should be fine
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"


using namespace moab;
using namespace std;


int main( int argc, char** argv )
{
    ProgOptions opts;

    string inputFile;
    opts.addOpt< string >( "inFile,i", "Specify the input file name string ", &inputFile );

    string outFile = "out.h5m";
    opts.addOpt< string >( "outFile,o", "Specify the output file name string ", &outFile );

    opts.parseCommandLine( argc, argv );


    // Instantiate
    Core mb;

    ErrorCode rval = mb.load_file( inputFile.c_str() );MB_CHK_SET_ERR( rval, "Error loading file" );

    // Get all vertices in the file set, and latVertex and lonVertex tags
    Range verts;
    rval = mb.get_entities_by_dimension( 0, 0, verts );MB_CHK_SET_ERR( rval, "Error getting vertices" );

    Tag latv, lonv;
    rval = mb.tag_get_handle("latVertex", latv);MB_CHK_SET_ERR( rval, "Error getting lat tag" );
    rval = mb.tag_get_handle("lonVertex", lonv);MB_CHK_SET_ERR( rval, "Error getting lon tag" );
    vector<double> latVerts(verts.size()), lonVerts(verts.size());
    rval = mb.tag_get_data(latv, verts, &latVerts[0]);MB_CHK_SET_ERR( rval, "Error getting lat tag values" );
    rval = mb.tag_get_data(lonv, verts, &lonVerts[0]);MB_CHK_SET_ERR( rval, "Error getting lon tag values" );

    double *x_ptr, *y_ptr, *z_ptr;
    int count;
    rval = mb.coords_iterate( verts.begin(), verts.end(), x_ptr, y_ptr, z_ptr, count );MB_CHK_SET_ERR( rval, "Error getting coords" );
    // change now coordinates, based on spherical coordinates
    for (int i=0; i<count; i++)
    {
    	struct IntxUtils::SphereCoords spP;
    	spP.R=1.0; spP.lon = lonVerts[i]; spP.lat = latVerts[i];
    	CartVect cart = IntxUtils::spherical_to_cart(spP);
    	*x_ptr=cart[0];
    	*y_ptr=cart[1];
    	*z_ptr=cart[2];
    	x_ptr++; y_ptr++; z_ptr++;
    }

    rval = mb.write_file( outFile.c_str());MB_CHK_SET_ERR( rval, "Error writing output file" );

    return 0;
}
