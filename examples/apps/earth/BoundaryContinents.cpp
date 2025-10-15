/**
 * @file BoundaryContinents.cpp
 * @brief Example demonstrating creation of continent boundary edge meshes
 *
 * This example shows how to:
 * - Load continent boundary point data and loop information
 * - Create edge meshes from boundary point sequences
 * - Convert 3D spherical coordinates to 2D lat-lon coordinates
 * - Handle periodic boundary conditions for longitude wrapping
 * - Remove long edges that cross the date line
 * - Create mesh sets for boundary visualization
 * - Write boundary meshes in VTK format
 *
 * This tool is useful for creating continent boundary meshes
 * for climate and geophysical applications.
 *
 * @author MOAB Development Team
 * @date 2024
 *
 * @example BoundaryContinents.cpp
 * Description: read boundary points and loops that form islands and continents
      and create an edge mesh file \n
 *
 *    BoundaryContinents  <boundary_points.dat> \c SaveLoopCounts
 * (default values can run if users don't specify input files)
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/CartVect.hpp"
#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;

string bd_name = string( MESH_DIR ) + string( "/../examples/earth/boundary_points.dat" );
string loops   = string( MESH_DIR ) + string( "/../examples/earth/SaveLoopCounts" );

double getLat( CartVect p )
{
    p.normalize();
    return asin( p[2] );
}
double getLon( CartVect p )
{
    p.normalize();
    double lon;

    lon = atan2( p[1], p[0] );
    if( lon < -2.95 )
    {  // this is Behring Strait :)
        return 2.0 * M_PI + lon;
    }
    else
    {
        return lon;
    }
}

int main( int argc, char** argv )
{
    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    if( argc > 1 )
    {
        // User has input file for boundary describing continents
        bd_name = argv[1];
    }
    if( argc > 2 )
    {
        // User has input file for loops for continents/islands
        loops = argv[2];
    }

    std::vector< double > coords;
    ifstream bdFile( bd_name.c_str() );
    if( !bdFile )
    {
        cout << endl << "Failed to open file " << bd_name << endl;
        ;
        return 1;
    }
    coords.resize( 0 );
    while( !bdFile.eof() )
    {
        double co;
        bdFile >> co;
        coords.push_back( co );
    }
    cout << " number of boundary points:" << coords.size() / 3 << "\n";
    /// get loops for boundaries
    vector< int > loopsindx;
    ifstream loopFile( loops.c_str() );
    while( !loopFile.eof() )
    {
        int indx;
        loopFile >> indx;
        loopsindx.push_back( indx );
    }
    cout << "number of loops: " << loopsindx.size() / 2 << "\n";
    ReadUtilIface* readMeshIface;
    mb->query_interface( readMeshIface );
    Range verts;

    MB_CHK_SET_ERR( mb->create_vertices( &coords[0], coords.size() / 3, verts ), "do not create boundary vertices" );

    Tag gid = mb->globalId_tag();

    int num_edges = 0;
    for( size_t i = 0; i < loopsindx.size() / 2; i++ )
    {
        num_edges += ( loopsindx[2 * i + 1] - loopsindx[2 * i] + 1 );
    }
    EntityHandle handle;
    EntityHandle* conn_array;
    MB_CHK_SET_ERR( readMeshIface->get_element_connect( num_edges, 2, MBEDGE, 1, handle, conn_array ), "do not elems" );
    int i1 = 0;  // index in edge connectivity
    for( size_t i = 0; i < loopsindx.size() / 2; i++ )
    {
        for( int j = loopsindx[2 * i]; j <= loopsindx[2 * i + 1]; j++ )
        {
            int j2;
            j2 = j;
            if( j == loopsindx[2 * i + 1] ) j2 = loopsindx[2 * i] - 1;  // close the loop
            conn_array[2 * i1]     = verts[j - 1];                      // vertices for sure start at 1
            conn_array[2 * i1 + 1] = verts[j2];                         // vertices for sure start at 1
            i1++;
        }
        // the last vertex is first one in the loop
    }

    Range bedges( handle, handle + num_edges - 1 );
    EntityHandle bound_set;
    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, bound_set ), "Can't create boundary edges set" );

    MB_CHK_SET_ERR( mb->add_entities( bound_set, bedges ), "Can't add edges to boundary set" );

    // set global ids for vertices and edges
    vector< int > gids;
    gids.resize( verts.size() );
    for( int j = 0; j < (int)verts.size(); j++ )
    {
        gids[j] = j + 1;
    }
    MB_CHK_SET_ERR( mb->tag_set_data( gid, verts, &gids[0] ), "Can't set global ids on verts" );
    // we have less edges than vertices, we can reuse the array for edges too
    MB_CHK_SET_ERR( mb->tag_set_data( gid, bedges, &gids[0] ), "Can't set global ids on edges" );

    mb->write_file( "bound.vtk", 0, 0, &bound_set, 1 );

    // now project the mesh in 2d plane
    // get all the vertices coordinates
    std::vector< CartVect > co3;
    co3.resize( verts.size() );
    MB_CHK_SET_ERR( mb->get_coords( verts, &( co3[0][0] ) ), "Can't get vertex coords" );
    for( size_t i = 0; i < verts.size(); i++ )
    {
        CartVect p  = co3[i];
        double lat1 = getLat( p );
        double lon1 = getLon( p );
        co3[i]      = CartVect( lon1, lat1, 0. );
    }
    MB_CHK_SET_ERR( mb->set_coords( verts, &( co3[0][0] ) ), "Can't set new vertex coords" );
    // remove edges in 2d that are too long (longer than 6; they are on the other side..., periodic)

    Range longEdges;
    for( Range::iterator eit = bedges.begin(); eit != bedges.end(); ++eit )
    {
        EntityHandle eh = *eit;
        const EntityHandle* conn;
        int nconn;
        MB_CHK_ERR( mb->get_connectivity( eh, conn, nconn ) );
        // get coordinates of nodes
        CartVect c2[2];
        MB_CHK_ERR( mb->get_coords( conn, 2, &( c2[0][0] ) ) );
        double dist = ( c2[1] - c2[0] ).length_squared();
        if( dist > 36. )  // too long edge in 2d, remove it
            longEdges.insert( eh );
    }
    MB_CHK_ERR( mb->delete_entities( longEdges ) );

    mb->write_file( "bound2d.vtk" );
    delete mb;

    return 0;
}
