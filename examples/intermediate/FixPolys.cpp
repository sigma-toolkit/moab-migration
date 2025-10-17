/**
 * @file FixPolys.cpp
 * @brief Example demonstrating how to fix polygons with duplicated vertices
 *
 * This example shows how to:
 * - Load a mesh file containing polygons
 * - Identify polygons with duplicated vertices
 * - Remove consecutive duplicate vertices from polygon connectivity
 * - Create new polygons with cleaned connectivity
 * - Preserve global ID tags during polygon reconstruction
 * - Write the corrected mesh to a new file
 *
 * The program is useful for cleaning up meshes that have been created
 * with padded vertices to reduce data sequences, but now need to be
 * optimized by removing redundant vertex references.
 *
 * @author MOAB Development Team
 * @date 2024
 *

 * \brief This test shows how to fix polygons that have duplicated vertices
 *
 * We sometimes use padded vertices option, to reduce the number of data sequences
 * We identify first the polygons that have padded vertices
 * then we set the new ones with reduced number of vertices, but with the same global id tag
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;
using namespace std;

#ifdef MOAB_HAVE_HDF5
string test_file_name = string( MESH_DIR ) + string( "/64bricks_512hex_256part.h5m" );
#endif
int main( int argc, char** argv )
{
    ProgOptions opts;

    string inputFile = test_file_name;
    opts.addOpt< string >( "inFile,i", "Specify the input file name string ", &inputFile );

    string outFile = "out.h5m";
    opts.addOpt< string >( "outFile,o", "Specify the output file name string ", &outFile );

    opts.parseCommandLine( argc, argv );

    // Instantiate
    Core mb;

    MB_CHK_SET_ERR( mb.load_file( inputFile.c_str() ), "Error loading file" );

    cout << " reading file " << inputFile << "\n";

    // Get all 2d elements in the file set
    Range elems;
    MB_CHK_SET_ERR( mb.get_entities_by_dimension( 0, 2, elems ), "Error getting 2d elements" );

    cout << "number of cells: " << elems.size() << "\n";

    Tag gidTag = mb.globalId_tag();
    Range OldCells;
    for( Range::iterator it = elems.begin(); it != elems.end(); ++it )
    {
        EntityHandle cell = *it;
        const EntityHandle* conn;
        int number_nodes;
        MB_CHK_SET_ERR( mb.get_connectivity( cell, conn, number_nodes ), "Error getting connectivity" );
        // now check if we have consecutive duplicated vertices, and if so, create a new cell
        std::vector< EntityHandle > new_verts;
        // push to it, if we do not have duplicates

        EntityHandle current = conn[0];
        // new_verts.push_back(current);
        for( int i = 1; i <= number_nodes; i++ )
        {
            EntityHandle nextV;
            if( i < number_nodes )
                nextV = conn[i];
            else
                nextV = conn[0];  // first vertex
            if( current != nextV )
            {
                new_verts.push_back( current );
                current = nextV;
            }
        }
        if( number_nodes > (int)new_verts.size() )
        {
            // create a new poly, and put this in a list to be removed
            int gid;
            MB_CHK_SET_ERR( mb.tag_get_data( gidTag, &cell, 1, &gid ), "Error getting global id tag" );
            EntityHandle newCell;
            MB_CHK_SET_ERR( mb.create_element( MBPOLYGON, &new_verts[0], (int)new_verts.size(), newCell ),
                            "Error creating new polygon " );
            MB_CHK_SET_ERR( mb.tag_set_data( gidTag, &newCell, 1, &gid ), "Error setting global id tag" );
            cout << "delete old cell " << cell << " with num_nodes vertices: " << number_nodes
                 << " and with global id: " << gid << "\n";
            for( int i = 0; i < number_nodes; i++ )
            {
                cout << " " << conn[i];
            }
            cout << "\n";
            OldCells.insert( cell );
        }
        mb.delete_entities( OldCells );
    }
    MB_CHK_SET_ERR( mb.write_file( outFile.c_str() ), "Error writing file" );

    return 0;
}
