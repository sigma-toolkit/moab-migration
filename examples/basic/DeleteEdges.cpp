/** @example DeleteEdges.cpp
 * This example demonstrates edge deletion from a mesh.
 * It shows how to load a mesh from a file,
 * retrieve all edges (1D entities) from the mesh,
 * delete all edges from the mesh,
 * and write the modified mesh to a new file.
 *
 * The resulting mesh will have vertices and higher-dimensional entities
 * (faces, volumes) but no edges, which can be useful for certain
 * mesh processing workflows.
 *
 * To run: ./DeleteEdges [meshfile] [outfile]
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string( MESH_DIR ) + string( "/hex01.vtk" );
string out_file       = string( "outFile.h5m" );

int main( int argc, char** argv )
{
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }
    if( argc > 2 )
    {
        // User has specified an output file
        out_file = argv[2];
    }

    // Instantiate & load a mesh from a file
    Core* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;
    MB_CHK_ERR( mb->load_mesh( test_file_name.c_str() ) );

    Range edges;
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 1, edges ) );
    MB_CHK_ERR( mb->delete_entities( edges ) );

    MB_CHK_ERR( mb->write_file( out_file.c_str() ) );
    delete mb;

    return 0;
}
