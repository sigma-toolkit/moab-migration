/*
 * proj1.cpp
 *
 *  project on a sphere of radius R, delete sets if needed, and delete edges between parts
 *  (created by resolve shared ents)
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>
#include <cmath>

#include "moab/climate/intx_mesh/IntxUtils.hpp"
#include <cassert>
using namespace moab;

double radius = 1.;  // in m:  6371220.

int main( int argc, char** argv )
{
    std::string file1;
    std::string file2;
    ProgOptions opts;

    opts.addOpt< std::string >( "input,i", "input file", &file1 );
    opts.addOpt< std::string >( "output,o", "output file", &file2 );

    double radius = 1.0;
    opts.addOpt< double >( std::string( "radius,R" ), std::string( "project to radius" ), &radius );

    opts.addOpt< void >( "deletePartitionSets,D", "delete partition sets from output file" );
    opts.addOpt< void >( "deleteEdges,E", "delete edges from output file" );

    opts.parseCommandLine( argc, argv );

    bool delete_partition_sets = opts.numOptSet( "deletePartitionSets" ) > 0;
    bool delete_edges          = opts.numOptSet( "deleteEdges" ) > 0;

    Core moab;
    Interface& mb = moab;

    ErrorCode rval = mb.load_mesh( file1.c_str() );MB_CHK_SET_ERR( rval, "can't read input file" );

    std::cout << "project to radius " << radius << " this input: " << file1 << " to output: " << file2 << "\n";

    Range verts;
    rval = mb.get_entities_by_dimension( 0, 0, verts );MB_CHK_SET_ERR( rval, "can't get vertices" );

    double *x_ptr, *y_ptr, *z_ptr;
    int count;
    rval = mb.coords_iterate( verts.begin(), verts.end(), x_ptr, y_ptr, z_ptr, count );MB_CHK_SET_ERR( rval, "can't coords iterate" );

    assert( count == (int)verts.size() );  // should end up with just one contiguous chunk of vertices

    for( int v = 0; v < count; v++ )
    {
        // EntityHandle v = verts[v];
        CartVect pos( x_ptr[v], y_ptr[v], z_ptr[v] );
        pos      = pos / pos.length();
        pos      = radius * pos;
        x_ptr[v] = pos[0];
        y_ptr[v] = pos[1];
        z_ptr[v] = pos[2];
    }

    if( delete_edges )
    {
        Range edges;
        rval = mb.get_entities_by_dimension( 0, 1, edges );
        if( MB_SUCCESS != rval ) return 1;
        // write edges to a new set, and after that, write the set, delete the edges and the set
        EntityHandle sf1;
        rval = mb.create_meshset( MESHSET_SET, sf1 );MB_CHK_SET_ERR( rval, "can't create edges set" );
        rval = mb.add_entities( sf1, edges );MB_CHK_SET_ERR( rval, "can't add edges to new set" );
        rval = mb.write_mesh( "edgesOnly.h5m", &sf1, 1 );MB_CHK_SET_ERR( rval, "can't write edges only set" );
        rval = mb.delete_entities( &sf1, 1 );MB_CHK_SET_ERR( rval, "can't delete edge set from database" );
        mb.delete_entities( edges );
    }

    if( delete_partition_sets )
    {
        Tag par_tag;
        rval = mb.tag_get_handle( "PARALLEL_PARTITION", par_tag );
        if( MB_SUCCESS == rval )

        {
            Range par_sets;
            rval =
                mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &par_tag, NULL, 1, par_sets, moab::Interface::UNION );
            if( !par_sets.empty() ) mb.delete_entities( par_sets );
            mb.tag_delete( par_tag );
        }
    }

    mb.write_file( file2.c_str() );

    return 0;
}
