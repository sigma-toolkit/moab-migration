/** @example GetEntities.cpp
 * Description: Get entities and report non-vertex entity connectivity and vertex adjacencies.\n
 * This example shows how to get connectivity and adjacencies.\n
 *
 * To run: ./GetEntities [meshfile]\n
 * (default values can run if users don't specify a mesh file)
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

int main( int argc, char** argv )
{
    int dimension = 3;
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
        dimension = atoi(argv[2]);
    }

    // Instantiate & load a mesh from a file
    Core* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;
    ErrorCode rval = mb->load_mesh( test_file_name.c_str() );MB_CHK_ERR( rval );

    moab::EntityHandle dimset;
    rval = mb->create_meshset( moab::MESHSET_SET, dimset );MB_CHK_SET_ERR( rval, "Can't create new set" );

    Range ents;
    // Get all entities in the database
    rval = mb->get_entities_by_dimension( 0, dimension, ents );MB_CHK_ERR( rval );

    rval = mb->add_entities( dimset, ents );MB_CHK_ERR( rval );

    rval = mb->write_file( "converted_mesh_file.vtk", "VTK", "", &dimset, 1 );MB_CHK_ERR( rval );

    delete mb;

    return 0;
}
