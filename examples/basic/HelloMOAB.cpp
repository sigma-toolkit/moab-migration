/** @example HelloMOAB.cpp
 * A simple example demonstrating basic MOAB functionality.
 * This example demonstrates the fundamental operations in MOAB:
 * - Creating a MOAB instance
 * - Loading a mesh from a file
 * - Querying entities by type and dimension
 * - Basic error handling
 *
 * The program reads a mesh file (default: 3k-tri-sphere.vtk) and reports
 * the number of vertices, edges, faces, and elements in the mesh.
 *
 * To run: ./HelloMOAB [meshfile]
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string( MESH_DIR ) + string( "/3k-tri-sphere.vtk" );

int main( int argc, char** argv )
{
    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    // Need option handling here for input filename
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }

    // Load the mesh from vtk file
    MB_CHK_ERR( mb->load_mesh( test_file_name.c_str() ) );

    // Get verts entities, by type
    Range verts;
    MB_CHK_ERR( mb->get_entities_by_type( 0, MBVERTEX, verts ) );

    // Get edge entities, by type
    Range edges;
    MB_CHK_ERR( mb->get_entities_by_type( 0, MBEDGE, edges ) );

    // Get faces, by dimension, so we stay generic to entity type
    Range faces;
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 2, faces ) );

    // Get regions, by dimension, so we stay generic to entity type
    Range elems;
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 3, elems ) );

    // Output the number of entities
    cout << "Number of vertices is " << verts.size() << endl;
    cout << "Number of edges is " << edges.size() << endl;
    cout << "Number of faces is " << faces.size() << endl;
    cout << "Number of elements is " << elems.size() << endl;

    delete mb;

    return 0;
}
