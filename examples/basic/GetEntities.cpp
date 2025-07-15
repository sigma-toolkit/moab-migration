/** @example GetEntities.cpp
 * This example demonstrates entity querying and connectivity access.
 * It shows how to get all entities in the mesh database,
 * access entity connectivity (vertex connectivity for elements),
 * query vertex adjacencies (elements connected to vertices),
 * and iterate through entities and examine their properties.
 *
 * The program reads a mesh file and reports connectivity information
 * for all non-vertex entities and adjacency information for vertices.
 *
 * To run: ./GetEntities [meshfile]
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
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }

    // Instantiate & load a mesh from a file
    Core* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;
    MB_CHK_ERR( mb->load_mesh( test_file_name.c_str() ) );

    Range ents;

    // Get all entities in the database
    MB_CHK_ERR( mb->get_entities_by_handle( 0, ents ) );

    for( Range::iterator it = ents.begin(); it != ents.end(); ++it )
    {
        if( MBVERTEX == mb->type_from_handle( *it ) )
        {
            Range adjs;
            MB_CHK_ERR( mb->get_adjacencies( &( *it ), 1, 3, false, adjs ) );
            cout << "Vertex " << mb->id_from_handle( *it ) << " adjacencies:" << endl;
            adjs.print();
        }
        else if( mb->type_from_handle( *it ) < MBENTITYSET )
        {
            const EntityHandle* connect;
            int num_connect;
            MB_CHK_ERR( mb->get_connectivity( *it, connect, num_connect ) );
            cout << CN::EntityTypeName( mb->type_from_handle( *it ) ) << " " << mb->id_from_handle( *it )
                 << " vertex connectivity is: ";
            for( int i = 0; i < num_connect; i++ )
                cout << mb->id_from_handle( connect[i] ) << " ";
            cout << endl;
        }
    }

    delete mb;

    return 0;
}
