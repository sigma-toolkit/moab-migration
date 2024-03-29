/** @example WritePartitionFile.cpp \n
 * \brief After partitioning a file with mbpart, write the
 * partition file as needed by mpas framework.
 */

#include <iostream>
#include <fstream>
#include <vector>

// Include header for MOAB instance and tag conventions for
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

int main( int argc, char** argv )
{
    moab::ErrorCode rval;
    ProgOptions opts;

    // Get MOAB instance
    moab::Interface* mb = new( std::nothrow ) moab::Core;
    if( NULL == mb ) return 1;

    std::string inputFile = std::string( MESH_DIR ) + "/globalmpas_deltri.h5m";
    opts.addOpt< std::string >( "inFile,i", "Specify the input file name string ", &inputFile );

    std::string outputFile = "Partfile.part";

    opts.addOpt< std::string >( "outFile,o", "Specify the output file name string ", &outputFile );

    opts.parseCommandLine( argc, argv );

    // This file is in the mesh files directory

    rval = mb->load_file( inputFile.c_str() );MB_CHK_SET_ERR( rval, "Failed to read" );

    // get all cells of dimension 2;
    moab::Tag part_set_tag;

    rval = mb->tag_get_handle( "PARALLEL_PARTITION", part_set_tag );MB_CHK_SET_ERR( rval, "Failed to get partition tag" );
    moab::Range tagged_sets;
    rval = mb->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &part_set_tag, NULL, 1, tagged_sets,
                                             moab::Interface::UNION );MB_CHK_SET_ERR( rval, "Failed to get partition sets" );
    if( tagged_sets.empty() )
    {
        MB_CHK_SET_ERR( moab::MB_FAILURE, "no partition sets" );
    }
    int num_sets = (int)tagged_sets.size();
    moab::Tag idtag;
    rval = mb->tag_get_handle( "GLOBAL_ID", idtag );MB_CHK_SET_ERR( rval, "Failed to read" );

    std::vector< int > part_tags_vals;
    part_tags_vals.resize( num_sets );
    rval = mb->tag_get_data( part_set_tag, tagged_sets, &part_tags_vals[0] );MB_CHK_SET_ERR( rval, "Failed to get tag vals" );

    std::map< int, int > mapid;
    int minNbCells, maxNbCells;
    for( int i = 0; i < num_sets; i++ )
    {
        moab::EntityHandle pset = tagged_sets[i];
        int val                 = part_tags_vals[i];
        moab::Range cells;
        rval = mb->get_entities_by_handle( pset, cells );MB_CHK_SET_ERR( rval, "Failed to get cells in set " );
        int nbCells = (int)cells.size();
        if( 0 == i )
        {
            minNbCells = maxNbCells = nbCells;
        }
        else
        {
            if( minNbCells > nbCells ) minNbCells = nbCells;
            if( maxNbCells < nbCells ) maxNbCells = nbCells;
        }
        std::vector< int > global_ids;
        global_ids.resize( cells.size() );
        rval = mb->tag_get_data( idtag, cells, &global_ids[0] );MB_CHK_SET_ERR( rval, "Failed to get global ids " );
        for( int j = 0; j < (int)cells.size(); j++ )
        {
            mapid[global_ids[j]] = val;
        }
    }
    std::cout << " number of cells: " << mapid.size() << " number of parts:" << num_sets << "\n";
    std::cout << " min/max nb cells per partition " << minNbCells << " / " << maxNbCells << "\n";
    std::fstream file;
    file.open( outputFile.c_str(), std::ios_base::out );
    // assume global ids start at 1
    for( int i = 0; i < (int)mapid.size(); i++ )
    {
        file << mapid[i + 1] << "\n";
    }
    file.close();

    delete mb;

    return 0;
}
