/** @example ReadPartFile.cpp
 * This example demonstrates how to read partition files created by Zoltan processes
 * and apply them to MOAB meshes. It shows how to:
 * - Load a mesh file and a partition file
 * - Remove existing partition sets from the mesh
 * - Create new partition sets based on the partition file
 * - Assign entities to appropriate partition sets
 * - Handle parallel partitioning with global IDs
 * - Write the partitioned mesh to a new file
 *
 * The partition file contains entity-to-partition assignments where entities
 * are identified by their global IDs. This is useful for load balancing
 * and parallel mesh processing workflows.
 *
 * Usage: ./ReadPartFile <input_mesh> <partition_file> <num_parts> <output_file>
 *
 * Example:
 * ./ReadPartFile mesh.h5m partition.txt 4 partitioned_mesh.h5m
 *
 * The partition file should contain one integer per entity (in order)
 * indicating which partition (0 to num_parts-1) each entity belongs to.
 */

#include "moab/Core.hpp"
#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string( MESH_DIR ) + string( "/3k-tri-sphere.vtk" );
string part_file_name;
int nparts;

int main( int argc, char** argv )
{
    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    // Need option handling here for input filename
    if( argc > 4 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
        part_file_name = argv[2];
        nparts         = atoi( argv[3] );
    }
    else
    {
        cerr << " usage is " << argv[0] << " <input file> <part file> <#parts> <output file> \n";
        exit( 0 );
    }

    ifstream inFile;
    inFile.open( part_file_name.c_str() );
    if( !inFile )
    {
        cerr << "Unable to open file " << part_file_name << "\n";
        exit( 1 );  // call system to stop
    }

    // Load the mesh from file
    MB_CHK_ERR( mb->load_mesh( test_file_name.c_str() ) );

    // Get sets entities, by type
    Range sets;
    MB_CHK_ERR( mb->get_entities_by_type( 0, MBENTITYSET, sets ) );

    // Output the number of sets
    cout << "Number of sets is " << sets.size() << endl;

    // remove the sets that have a PARALLEL_PARTITION tag
    Tag tag;
    MB_CHK_ERR( mb->tag_get_handle( "PARALLEL_PARTITION", tag ) );

    int i                = 0;
    int num_deleted_sets = 0;
    ErrorCode rval;
    for( Range::iterator it = sets.begin(); it != sets.end(); it++, i++ )
    {
        EntityHandle eh = *it;
        // cout << " set :" << mb->id_from_handle(eh) <<"\n";
        int val = -1;
        rval    = mb->tag_get_data( tag, &eh, 1, &val );
        if( val != -1 )
        {
            num_deleted_sets++;
            rval = mb->delete_entities( &eh, 1 );  // delete the set, we will have a new partition soon
        }
    }
    if( num_deleted_sets ) cout << "delete " << num_deleted_sets << " existing  partition sets, and create new ones \n";

    Range cells;  // get them by dimension 2!
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 2, cells ) );
    EntityHandle* psets = new EntityHandle[nparts];
    for( int i = 0; i < nparts; i++ )
    {
        MB_CHK_ERR( mb->create_meshset( MESHSET_SET, psets[i] ) );
        MB_CHK_ERR( mb->tag_set_data( tag, &( psets[i] ), 1, &i ) );
    }

    for( Range::iterator it = cells.begin(); it != cells.end(); it++ )
    {
        int part;
        EntityHandle eh = *it;
        inFile >> part;
        MB_CHK_ERR( mb->add_entities( psets[part], &eh, 1 ) );
    }

    MB_CHK_ERR( mb->write_file( argv[4] ) );

    delete[] psets;
    delete mb;

    return 0;
}
