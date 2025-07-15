/**
 * @file SetsNTags.cpp
 * @brief Example demonstrating entity sets and tag operations
 *
 * @details This example shows how to:
 * - Work with entity sets (material sets, boundary condition sets)
 * - Access conventional tags from MBTagConventions.hpp
 * - Query entities by type and tag
 * - Retrieve set contents recursively
 * - Access tag data on sets
 *
 * The program reads a mesh file and identifies material sets, Dirichlet
 * boundary condition sets, and Neumann boundary condition sets, then
 * reports their contents.
 *
 * @author MOAB Team
 * @date 2024
 *
 * @param[in] argc Number of command line arguments
 * @param[in] argv Command line arguments array
 * @param[in] argv[1] Optional mesh file path (default: hex01.vtk)
 *
 * @return 0 on success, 1 on failure
 *
 * @par Usage:
 * @code
 * ./SetsNTags [meshfile]
 * @endcode
 *
 * @par Example:
 * @code
 * ./SetsNTags my_mesh.vtk
 * @endcode
 *
 * @see Core, Interface, Range, Tag, MBTagConventions
 *

 * Description: Get the sets representing materials and Dirichlet/Neumann boundary conditions and
 * list their contents.\n This example shows how to get entity sets, and tags on those sets.
 *
 * To run: ./SetsNTags [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string( MESH_DIR ) + string( "/hex01.vtk" );

// Tag names for these conventional tags come from MBTagConventions.hpp
const char* tag_nms[] = { MATERIAL_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME };

int main( int argc, char** argv )
{
    // Get the material set tag handle
    Tag mtag;
    ErrorCode rval;
    Range sets, set_ents;

    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    // Need option handling here for input filename
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }

    // Load a file
    MB_CHK_ERR( mb->load_file( test_file_name.c_str() ) );

    // Loop over set types
    for( int i = 0; i < 3; i++ )
    {
        // Get the tag handle for this tag name; tag should already exist (it was created during
        // file read)
        MB_CHK_ERR( mb->tag_get_handle( tag_nms[i], 1, MB_TYPE_INTEGER, mtag ) );

        // Get all the sets having that tag (with any value for that tag)
        sets.clear();
        MB_CHK_ERR( mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &mtag, NULL, 1, sets ) );

        // Iterate over each set, getting the entities and printing them
        Range::iterator set_it;
        for( set_it = sets.begin(); set_it != sets.end(); ++set_it )
        {
            // Get the id for this set
            int set_id;
            MB_CHK_ERR( mb->tag_get_data( mtag, &( *set_it ), 1, &set_id ) );

            // Get the entities in the set, recursively
            MB_CHK_ERR( mb->get_entities_by_handle( *set_it, set_ents, true ) );

            cout << tag_nms[i] << " " << set_id << " has " << set_ents.size() << " entities:" << endl;
            set_ents.print( "   " );
            set_ents.clear();
        }
    }

    delete mb;

    return 0;
}
