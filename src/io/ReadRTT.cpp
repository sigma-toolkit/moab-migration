/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

/**
 * \class ReadRTT
 * \brief ReadRTT based on ReadNASTRAN
 *
 * See:
 *
 * \author Andrew Davis
 */

#include "ReadRTT.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "FileTokenizer.hpp"
#include "Internals.hpp"  // for MB_START_ID
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include "moab/ErrorHandler.hpp"
#include "moab/FileOptions.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ReadUtilIface.hpp"

namespace moab
{

ReaderIface* ReadRTT::factory( Interface* iface )
{
    return new ReadRTT( iface );
}

// constructor
ReadRTT::ReadRTT( Interface* impl )
    : MBI( impl ), geom_tag( 0 ), id_tag( 0 ), name_tag( 0 ), category_tag( 0 ), faceting_tol_tag( 0 )
{
    assert( NULL != impl );
    myGeomTool = new GeomTopoTool( impl );
    MBI->query_interface( readMeshIface );
    assert( NULL != readMeshIface );

    // this section copied from ReadCGM initalisation
    int negone  = -1;
    double zero = 0.;
    ErrorCode rval;
    rval = MBI->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag, MB_TAG_SPARSE | MB_TAG_CREAT,
                                &negone );
    MB_CHK_ERR_CONT( rval );
    id_tag = MBI->globalId_tag();
    rval = MBI->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE, name_tag, MB_TAG_SPARSE | MB_TAG_CREAT );
    MB_CHK_ERR_CONT( rval );
    rval = MBI->tag_get_handle( CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, MB_TYPE_OPAQUE, category_tag,
                                MB_TAG_SPARSE | MB_TAG_CREAT );
    MB_CHK_ERR_CONT( rval );
    rval =
        MBI->tag_get_handle( "FACETING_TOL", 1, MB_TYPE_DOUBLE, faceting_tol_tag, MB_TAG_SPARSE | MB_TAG_CREAT, &zero );
    MB_CHK_ERR_CONT( rval );
}

// destructor
ReadRTT::~ReadRTT()
{
    if( readMeshIface )
    {
        MBI->release_interface( readMeshIface );
        readMeshIface = 0;
    }

    delete myGeomTool;
}

ErrorCode ReadRTT::read_tag_values( const char* /*file_name*/,
                                    const char* /*tag_name*/,
                                    const FileOptions& /*opts*/,
                                    std::vector< int >& /*tag_values_out*/,
                                    const SubsetList* /*subset_list*/ )
{
    return MB_NOT_IMPLEMENTED;
}

// load the file as called by the Interface function
ErrorCode ReadRTT::load_file( const char* filename,
                              const EntityHandle*,
                              const FileOptions&,
                              const ReaderIface::SubsetList* subset_list,
                              const Tag* /*file_id_tag*/ )
{
    ErrorCode rval;

    // at this time there is no support for reading a subset of the file
    if( subset_list )
    {
        std::cout << "Subset reading not supported for RTT meshes" << std::endl;
        return MB_UNSUPPORTED_OPERATION;
    }

    // test to see if file exists
    FILE* file = NULL;
    file       = fopen( filename, "r" );
    if( file == NULL ) return MB_FILE_DOES_NOT_EXIST;
    // otherwise close the file
    fclose( file );

    // read the header
    rval = ReadRTT::read_header( filename );
    if( rval != MB_SUCCESS ) return rval;

    // read the side_flag data
    rtt_flags side_flags;
    rval = ReadRTT::read_side_flags( filename, side_flags );
    if( rval != MB_SUCCESS ) return rval;

    //process sides
    std::vector< side > side_data;
    rval = ReadRTT::side_process_faces( side_flags, side_data );
    if( rval != MB_SUCCESS ) return rval;

    // read the cell data
    rtt_flags cell_flags;
    rval = ReadRTT::read_cell_flags( filename, cell_flags );
    if( rval != MB_SUCCESS ) return rval;

    //process REGIONS
    rval = ReadRTT::cell_process_flag( cell_flags, "REGIONS" , regions_data );
    if( rval != MB_SUCCESS ) return rval;

    // process the ABAQUS_PARTS
    rval = ReadRTT::cell_process_flag( cell_flags, "ABAQUS_PARTS" , abaqus_parts_data );
    if( rval != MB_SUCCESS ) return rval; 

    // process the MCNP_PSEUDO-CELLS
    // Give priority to MCNOP_PSEUDO-CELLS for cell ID.
    rval = ReadRTT::cell_process_flag( cell_flags, "MCNP_PSEUDO-CELLS" , mcnp_pseudo_cells_data );
    if( rval != MB_SUCCESS ) return rval;
  
    // build the index of the cell flags data
    build_idx();

    // store the cell data definitions, MCNP_PSEUDO-CELLS if present
    // otherwise use the REGIONS data
    if( mcnp_pseudo_cells_data.size() > 0 )
    {
        cell_data = mcnp_pseudo_cells_data;
        cell_data_idx = mcnp_pseudo_cells_idx;
    }
    else
    {
        cell_data = regions_data;
        cell_data_idx = regions_idx;
    }

    // read the node data
    std::vector< node > node_data;
    rval = ReadRTT::read_nodes( filename, node_data );
    if( rval != MB_SUCCESS ) return rval;

    // read the facet data
    std::vector< facet > facet_data;
    rval = ReadRTT::read_facets( filename, facet_data );
    if( rval != MB_SUCCESS ) return rval;

    // read the tetrahedra data
    std::vector< tet > tet_data;
    rval = ReadRTT::read_tets( filename, tet_data );
    if( rval != MB_SUCCESS ) return rval;

    // make the map of surface number in the rttmesh to the surface meshset
    std::map< int, EntityHandle > surface_map;  // corrsespondance of surface number to entity handle
    rval = ReadRTT::generate_topology( side_data, cell_data, surface_map );
    if( rval != MB_SUCCESS ) return rval;

    // generate the rest of the database, triangles to surface meshsets etc
    rval = ReadRTT::build_moab( node_data, facet_data, tet_data, surface_map );
    if( rval != MB_SUCCESS ) return rval;

    return MB_SUCCESS;
}


ErrorCode ReadRTT::build_idx()
{
    // build the index of the cell data
    for( size_t i = 0; i < regions_data.size(); ++i )
    {
        cell tmp = regions_data[i];
        regions_idx[tmp.id] = i;
    }

    // build the index of the abaqus parts
    for ( size_t i = 0; i < abaqus_parts_data.size(); ++i )
    {
        cell tmp = abaqus_parts_data[i];
        abaqus_parts_idx[tmp.id] = i;
    }

    // build the index of the mcnp pseudo cells
    for ( size_t i = 0; i < mcnp_pseudo_cells_data.size(); ++i )
    {
        cell tmp = mcnp_pseudo_cells_data[i];
        mcnp_pseudo_cells_idx[tmp.id] = i;
    }
    
    return MB_SUCCESS;
}

/*
 * builds the topology of the problem
 */
ErrorCode ReadRTT::generate_topology( std::vector< side > side_data,
                                      std::vector< cell > cell_data,
                                      std::map< int, EntityHandle >& surface_map )
{

    ErrorCode rval;
    std::vector< EntityHandle > entmap[4];
    int num_ents[4];  // number of entities in each dimension

    const char geom_categories[][CATEGORY_TAG_SIZE] = { "Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0" };

    std::vector< int > surface_numbers;  // the surface numbers in the problem

    // corresponds to number of cad like surfaces and cad like volumes
    num_ents[2] = side_data.size();
    num_ents[3] = cell_data.size();

    // loop over surfaces & volumes
    for( int dim = 2; dim <= 3; dim++ )
    {
        for( int i = 0; i != num_ents[dim]; i++ )
        {
            EntityHandle handle;
            // create a meshset for each entity surface/volume
            rval = MBI->create_meshset( dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle );
            // if failure
            if( rval != MB_SUCCESS ) return rval;

            // collect the entity handles into an
            entmap[dim].push_back( handle );

            // set the dimension tag
            rval = MBI->tag_set_data( geom_tag, &handle, 1, &dim );
            // if fail
            if( MB_SUCCESS != rval ) return rval;
            // if we are a surface
            if( dim == 2 )
            {
                // tag the id onto the surface meshset
                rval = MBI->tag_set_data( id_tag, &handle, 1, &side_data[i].id );
                // inesert entity into the map
                surface_map[side_data[i].id] = handle;
            }
            else
            {
                // otherwise we set the volume tag data, loop is only 2 & 3 dim
                rval = MBI->tag_set_data( id_tag, &handle, 1, &cell_data[i].id );
            }
            // if fail
            if( MB_SUCCESS != rval ) return rval;
            // set the category tag
            rval = MBI->tag_set_data( category_tag, &handle, 1, &geom_categories[dim] );
            if( MB_SUCCESS != rval ) return rval;
        }
    }

    // generate parent child links
    // best to loop over the surfaces and assign them to volumes, we can then
    // assign facets to
    // to each surface
    generate_parent_child_links( num_ents, entmap, side_data, cell_data );

    // set the surface senses
    set_surface_senses( num_ents, entmap, side_data, cell_data );

    // set the group data
    rval = setup_group_data( entmap );

    return MB_SUCCESS;
}

/*
 * builds the moab representation of the mesh
 */
ErrorCode ReadRTT::build_moab( std::vector< node > node_data,
                               std::vector< facet > facet_data,
                               std::vector< tet > tet_data,
                               std::map< int, EntityHandle > surface_map )
{

    ErrorCode rval;         // reusable return value
    EntityHandle file_set;  // the file handle
    // create the file set
    rval = MBI->create_meshset( MESHSET_SET, file_set );
    if( MB_SUCCESS != rval ) return rval;

    // create the vertices
    EntityHandle handle;
    std::vector< node >::iterator it;  // iterate over the nodes
    Range mb_coords;                   // range of coordinates
    for( it = node_data.begin(); it != node_data.end(); ++it )
    {
        node tmp         = *it;
        double coords[3] = { tmp.x, tmp.y, tmp.z };
        rval             = MBI->create_vertex( coords, handle );
        if( MB_SUCCESS != rval ) return rval;
        mb_coords.insert( handle );  // inesert handle into the coordinate range
    }

    // add verts to set
    rval = MBI->add_entities( file_set, mb_coords );

    // create sense tag
    Tag side_id_tag, surface_number_tag;
    //  int zero = 0;
    rval = MBI->tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER, side_id_tag, MB_TAG_SPARSE | MB_TAG_CREAT );
    rval =
        MBI->tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER, surface_number_tag, MB_TAG_SPARSE | MB_TAG_CREAT );

    // create the facets
    EntityHandle triangle;
    std::vector< facet >::iterator it_f;
    // range of triangles
    Range mb_tris;
    // loop over the facet data
    for( it_f = facet_data.begin(); it_f != facet_data.end(); ++it_f )
    {
        facet tmp = *it_f;
        // get the nodes for the triangle
        EntityHandle tri_nodes[3] = { mb_coords[tmp.connectivity[0] - 1], mb_coords[tmp.connectivity[1] - 1],
                                      mb_coords[tmp.connectivity[2] - 1] };
        // create a triangle element
        rval = MBI->create_element( MBTRI, tri_nodes, 3, triangle );
        // tag in side id on the triangle
        rval = MBI->tag_set_data( side_id_tag, &triangle, 1, &tmp.side_id );
        // tag the surface number on the triangle
        rval = MBI->tag_set_data( surface_number_tag, &triangle, 1, &tmp.surface_number );
        // insert vertices and triangles into the appropriate surface meshset
        EntityHandle meshset_handle = surface_map[tmp.surface_number];
        // also set surface tag
        rval = MBI->tag_set_data( side_id_tag, &meshset_handle, 1, &tmp.side_id );
        rval = MBI->tag_set_data( surface_number_tag, &meshset_handle, 1, &tmp.surface_number );
        // add vertices to the mesh
        rval = MBI->add_entities( meshset_handle, &( *tri_nodes ), 3 );
        // add triangles to the meshset
        rval = MBI->add_entities( meshset_handle, &triangle, 1 );
        // ineter triangles into large run
        mb_tris.insert( triangle );
    }
    // add tris to set to fileset
    rval = MBI->add_entities( file_set, mb_tris );

    // create material number tag
    Tag mat_num_tag;
    //  int zero = 0;
    rval = MBI->tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER, mat_num_tag, MB_TAG_SPARSE | MB_TAG_CREAT );

    Tag mat_name_tag;
    rval = MBI->tag_get_handle( "MATERIAL_NAME", 1, MB_TYPE_OPAQUE, mat_name_tag, MB_TAG_SPARSE | MB_TAG_CREAT );


    // create the tets
    EntityHandle tetra;  // handle for a specific tet
    std::vector< tet >::iterator it_t;
    Range mb_tets;
    // loop over all tets
    for( it_t = tet_data.begin(); it_t != tet_data.end(); ++it_t )
    {
        tet tmp = *it_t;
        // get the handles for the tet
        EntityHandle tet_nodes[4] = { mb_coords[tmp.connectivity[0] - 1], mb_coords[tmp.connectivity[1] - 1],
                                      mb_coords[tmp.connectivity[2] - 1], mb_coords[tmp.connectivity[3] - 1] };
        // create the tet
        rval           = MBI->create_element( MBTET, tet_nodes, 4, tetra );
        int mat_number = tmp.material_number;
        // tag the tet with the material number
        rval = MBI->tag_set_data( mat_num_tag, &tetra, 1, &mat_number );

        // // Add material name tag
        auto idx = cell_data_idx.find(mat_number);
        if(idx != cell_data_idx.end()) 
        {
            const char* mat_name = cell_data[idx->second].name.c_str();
            // set the material name tag
            rval = MBI->tag_set_data(mat_name_tag, &tetra, 1, mat_name);
        }



        // set the tag data
        mb_tets.insert( tetra );
    }
    // add tris to set
    rval = MBI->add_entities( file_set, mb_tets );

    rval = add_metadata( file_set );

    return MB_SUCCESS;
}

moab::ErrorCode ReadRTT::add_metadata( EntityHandle file_set )
{
    moab::ErrorCode rval = MB_FAILURE;

    // Create CONTIGUITY tag and set its value
    Tag contiguity_tag;
    const char* contiguity_value = header_data.contiguity.c_str();
    rval = MBI->tag_get_handle( "CONTIGUITY", strlen( contiguity_value ) + 1, MB_TYPE_OPAQUE, contiguity_tag,
                                MB_TAG_SPARSE | MB_TAG_CREAT );
    if( rval != MB_SUCCESS ) return rval;
    rval = MBI->tag_set_data( contiguity_tag, &file_set, 1, contiguity_value );

    return rval;
}

/*
 * read the header data from the filename pointed to
 */
ErrorCode ReadRTT::read_header( const char* filename )
{
    std::ifstream input_file( filename );  // filename for rtt file
    // file ok?
    if( !input_file.good() )
    {
        std::cout << "Problems reading file = " << filename << std::endl;
        return MB_FAILURE;
    }

    // if it works
    std::string line;
    moab::ErrorCode rval = MB_FAILURE;
    if( input_file.is_open() )
    {
        while( std::getline( input_file, line ) )
        {
            if( line.compare( "header" ) == 0 )
            {
                rval = get_header_data( input_file );
            }
            else if( line.compare( "dims" ) == 0 )
            {
                rval = parse_dims( input_file );
            }
            else if( line.compare( "cell_defs" ) == 0 )
            {
                rval = ReadRTT::read_cell_defs( input_file );
            }
        }
        input_file.close();
    }
    return rval;
}

/*
 * reads the side data from the filename pointed to
 */
ErrorCode ReadRTT::read_side_flags( const char* filename, rtt_flags& side_flags )
{
    ErrorCode rval = MB_FAILURE;
    // read all the side data
    rval = read_all_flags( filename, dim_data.nside_flags, "side", side_flags );

    return rval;
}

ErrorCode ReadRTT::read_all_flags(const char* filename, std::vector<int> n_flags, std::string flag_id, rtt_flags& flags)
{
    std::string start_flag = flag_id + "_flags";
    std::string end_flag   = "end_" + start_flag + "\0";
    std::string line;                      // the current line being read
    std::ifstream input_file( filename );  // filestream for rttfile
    // file ok?
    if( !input_file.good() )
    {
        std::cout << "Problems reading file = " << filename << std::endl;
        return MB_FAILURE;
    }
    // if it works
    if( input_file.is_open() )
    {
        while( std::getline( input_file, line ) )
        {
            if( line.compare( start_flag ) == 0 )
            {
                while( std::getline( input_file, line ) )
                {
                    // Read all the side block until we find the end
                    if( line.compare( end_flag ) == 0 ) break;
                    std::vector< std::string > token = ReadRTT::split_string( line, ' ' );
                    if( token.size() != 2 )
                    {
                        std::cout << "Error reading side flags" << std::endl;
                        return MB_FAILURE;
                    }
                    int flag_key    = std::stoi( token[0] ) - 1;
                    std::string key = token[1];
                    for( int i = 0; i < n_flags[flag_key]; i++ )
                    {
                        std::getline( input_file, line );
                        flags[key].push_back( line );
                    }
                }
            }
        }
        input_file.close();
    } 
    return MB_SUCCESS;
}

/*
 * process the FACES flag from the side_flags section
 */
ErrorCode ReadRTT::side_process_faces( rtt_flags side_flags, std::vector< side >& side_data )
{
    if( side_flags.find( "FACES" ) != side_flags.end() )
    {
        for( size_t i = 0; i < side_flags["FACES"].size(); i++ )
        {
            side data = ReadRTT::get_side_data( side_flags["FACES"][i] );
            side_data.push_back( data );
        }
    }
    if( side_data.size() == 0 ) return MB_FAILURE;
    return MB_SUCCESS;
}

/*
 * reads the cell data from the filename pointed to
 */
ErrorCode ReadRTT::read_cell_flags( const char* filename, rtt_flags& cell_flags )
{
    ErrorCode rval = MB_FAILURE;
    // read all the cell data
    rval = read_all_flags( filename, dim_data.ncell_flags, "cell", cell_flags );
    return rval;
}

/*
 * process the standard flag from the cell_flags section
 */
ErrorCode ReadRTT::cell_process_flag( rtt_flags cell_flags, std::string key, std::vector< cell >& cell_data )
{
    if( cell_flags.find( key ) != cell_flags.end() )
    {
        for( size_t i = 0; i < cell_flags[key].size(); i++ )
        {
            cell data = ReadRTT::get_cell_data( cell_flags[key][i] );
            cell_data.push_back( data );
        }
        if( cell_data.size() == 0 ) return MB_FAILURE;
    }
    return MB_SUCCESS;
}


/*
 * Reads the node data fromt the filename pointed to
 */
ErrorCode ReadRTT::read_nodes( const char* filename, std::vector< node >& node_data )
{
    std::string line;                      // the current line being read
    std::ifstream input_file( filename );  // filestream for rttfile
    // file ok?
    if( !input_file.good() )
    {
        std::cout << "Problems reading file = " << filename << std::endl;
        return MB_FAILURE;
    }

    // if it works
    if( input_file.is_open() )
    {
        while( std::getline( input_file, line ) )
        {
            if( line.compare( "nodes\0" ) == 0 )
            {
                // read lines until find end nodes
                while( std::getline( input_file, line ) )
                {
                    if( line.compare( "end_nodes\0" ) == 0 ) break;
                    node data = ReadRTT::get_node_data( line );
                    node_data.push_back( data );
                }
            }
        }
        input_file.close();
    }
    if( node_data.size() == 0 ) return MB_FAILURE;
    return MB_SUCCESS;
}

/*
 * Reads the facet data fromt the filename pointed to
 */
ErrorCode ReadRTT::read_facets( const char* filename, std::vector< facet >& facet_data )
{
    std::string line;                      // the current line being read
    std::ifstream input_file( filename );  // filestream for rttfile
    // file ok?
    if( !input_file.good() )
    {
        std::cout << "Problems reading file = " << filename << std::endl;
        return MB_FAILURE;
    }

    // if it works
    if( input_file.is_open() )
    {
        while( std::getline( input_file, line ) )
        {
            if( line.compare( "sides\0" ) == 0 )
            {
                // read lines until find end nodes
                while( std::getline( input_file, line ) )
                {
                    if( line.compare( "end_sides\0" ) == 0 ) break;
                    facet data = ReadRTT::get_facet_data( line );
                    facet_data.push_back( data );
                }
            }
        }
        input_file.close();
    }
    if( facet_data.size() == 0 ) return MB_FAILURE;
    return MB_SUCCESS;
}

/*
 * Reads the facet data fromt the filename pointed to
 */
ErrorCode ReadRTT::read_tets( const char* filename, std::vector< tet >& tet_data )
{
    std::string line;                      // the current line being read
    std::ifstream input_file( filename );  // filestream for rttfile
    // file ok?
    if( !input_file.good() )
    {
        std::cout << "Problems reading file = " << filename << std::endl;
        return MB_FAILURE;
    }
    // if it works
    if( input_file.is_open() )
    {
        while( std::getline( input_file, line ) )
        {
            if( line.compare( "cells\0" ) == 0 )
            {
                // read lines until find end nodes
                while( std::getline( input_file, line ) )
                {
                    if( line.compare( "end_cells\0" ) == 0 ) break;
                    tet data = ReadRTT::get_tet_data( line );
                    tet_data.push_back( data );
                }
            }
        }
        input_file.close();
    }
    if( tet_data.size() == 0 ) return MB_FAILURE;
    return MB_SUCCESS;
}

/*
 * given the open file handle read until we find
 */
ErrorCode ReadRTT::get_header_data( std::ifstream& input_file )
{
    std::string line;
    while( std::getline( input_file, line ) )
    {

        // tokenize the line
        std::istringstream iss( line );
        std::vector< std::string > split_string;
        do
        {
            std::string sub_string;
            iss >> sub_string;
            split_string.push_back( sub_string );
        } while( iss );

        // if we find version
        if( line.find( "version" ) != std::string::npos )
        {
            if( split_string[1].find( "v" ) != std::string::npos &&
                split_string[0].find( "version" ) != std::string::npos )
            {
                header_data.version = split_string[1];
            }
        }
        else if( line.find( "title" ) != std::string::npos )
        {
            header_data.title = split_string[1];
        }
        else if( line.find( "date" ) != std::string::npos )
        {
            header_data.date = split_string[1];
        }
        else if( line.find( "contiguity" ) != std::string::npos )
        {
            header_data.contiguity = split_string[1];
        }
        else if( line.find( "end_header" ) != std::string::npos )
        {
            return MB_SUCCESS;
        }
    }
    // otherwise we never found the end_header keyword
    return MB_FAILURE;
}

ErrorCode ReadRTT::parse_dims( std::ifstream& input_file )
{
    if( !input_file.good() || !input_file.is_open() )
    {
        std::cout << "Problems reading file" << std::endl;
        return MB_FAILURE;
    }

    std::string line;
    std::vector< std::string > tokens;
    while( std::getline( input_file, line ) )
    {
        if( line == "" ) continue;
        if( line.find( "end_dims" ) != std::string::npos ) break;

        tokens = ReadRTT::split_string( line, ' ' );
        if( tokens[0] == "coor_units" )
        {
            dim_data.coor_units = tokens[1];
        }
        else if( tokens[0] == "prob_time_units" )
        {
            dim_data.prob_time_units = tokens[1];
        }
        else if( tokens[0] == "ncell_defs" )
        {
            dim_data.ncell_defs = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nnodes_max" )
        {
            dim_data.nnodes_max = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nsides_max" )
        {
            dim_data.nsides_max = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nnodes_sides_max" )
        {
            dim_data.nnodes_sides_max = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "ndim" )
        {
            dim_data.ndim = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "n_dim_topo" )
        {
            dim_data.n_dim_topo = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nnodes" )
        {
            dim_data.nnodes = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nnode_flag_types" )
        {
            dim_data.nnode_flag_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nnode_flags" )
        {
            for( size_t i = 1; i < tokens.size(); i++ )
            {
                dim_data.nnode_flags.push_back( std::atoi( tokens[i].c_str() ) );
            }
        }
        else if( tokens[0] == "nnode_data" )
        {
            dim_data.nnode_data = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nsides" )
        {
            dim_data.nsides = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nside_types" )
        {
            dim_data.nside_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "side_types" )
        {
            dim_data.side_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nside_flag_types" )
        {
            dim_data.nside_flag_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "nside_flags" )
        {
            for( size_t i = 1; i < tokens.size(); i++ )
            {
                dim_data.nside_flags.push_back( std::atoi( tokens[i].c_str() ) );
            }
        }
        else if( tokens[0] == "nside_data" )
        {
            dim_data.nside_data = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "ncells" )
        {
            dim_data.ncells = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "ncell_types" )
        {
            dim_data.ncell_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "cell_types" )
        {
            dim_data.cell_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "ncell_flag_types" )
        {
            dim_data.ncell_flag_types = std::atoi( tokens[1].c_str() );
        }
        else if( tokens[0] == "ncell_flags" )
        {
            for( size_t i = 1; i < tokens.size(); i++ )
            {
                dim_data.ncell_flags.push_back( std::atoi( tokens[i].c_str() ) );
            }
        }
        else if( tokens[0] == "ncell_data" )
        {
            dim_data.ncell_data = std::atoi( tokens[1].c_str() );
        }
    }
    // Check that the data is valid and has the expected number of entries
    dim_data.validate();

    return MB_SUCCESS;
}

/*
 * parsing the cell deifinitions
 */
ErrorCode ReadRTT::read_cell_defs( std::ifstream& input_file )
{
    if( !input_file.good() || !input_file.is_open() )
    {
        std::cout << "Problems reading file" << std::endl;
        return MB_FAILURE;
    }

    std::string line;
    std::vector< std::string > tokens;
    while( std::getline( input_file, line ) )
    {
        if( line == "" ) continue;
        if( line.find( "end_cell_defs" ) != std::string::npos ) break;
        cell_def new_cell_def;
        // Tokenize the line
        tokens            = ReadRTT::split_string( line, ' ' );
        new_cell_def.id   = std::atoi( tokens[0].c_str() );
        new_cell_def.name = tokens[1];
        // Getting the number of nodes and sides
        std::getline( input_file, line );
        // Tokenize the line
        tokens              = ReadRTT::split_string( line, ' ' );
        new_cell_def.nnodes = std::atoi( tokens[0].c_str() );
        new_cell_def.nsides = std::atoi( tokens[1].c_str() );
        // Side type index
        std::getline( input_file, line );
        // Tokenize the line
        tokens = ReadRTT::split_string( line, ' ' );
        for( int i = 0; i < new_cell_def.nsides; i++ )
        {
            int side_type = std::atoi( tokens[i].c_str() );
            // Ensure side types exists
            if( cell_def_data.find( side_type ) == cell_def_data.end() )
            {
                std::cout << "Error: side type " << side_type << " not found in cell definitions." << std::endl;
                return MB_FAILURE;
            }
            new_cell_def.side_type.push_back( side_type );
        }
        // Read the nodes per side
        for( int i = 0; i < new_cell_def.nsides; i++ )
        {
            std::vector< int > side_nodes;
            std::getline( input_file, line );
            // Tokenize the line
            tokens = ReadRTT::split_string( line, ' ' );
            for( int j = 0; j < new_cell_def.nnodes; j++ )
            {
                side_nodes.push_back( std::atoi( tokens[i].c_str() ) );
            }
            new_cell_def.sides_nodes.push_back( side_nodes );
        }

        cell_def_data.insert( std::make_pair( new_cell_def.id, new_cell_def ) );
    }

    return MB_SUCCESS;
}

/*
 * given the string sidedata, get the id number, senses and names of the sides
 */
ReadRTT::side ReadRTT::get_side_data( std::string sidedata )
{
    side new_side;
    std::vector< std::string > tokens;
    tokens = ReadRTT::split_string( sidedata, ' ' );

    // set the side id
    if( tokens.size() != 2 )
    {
        MB_SET_ERR_RET_VAL( "Error, too many tokens found from side_data", new_side );
    }
    // create the new side
    new_side.id = std::atoi( tokens[0].c_str() );

    std::vector< std::string > cell_names = ReadRTT::split_string( tokens[1], '/' );
    // get the boundary
    boundary new_bnd = ReadRTT::split_name( cell_names[0] );
    // set the surface sense and name
    new_side.senses[0] = new_bnd.sense;
    new_side.names[0]  = new_bnd.name;
    //
    if( cell_names.size() > 1 )
    {
        boundary bnd       = ReadRTT::split_name( cell_names[1] );
        new_side.senses[1] = bnd.sense;
        new_side.names[1]  = bnd.name;
    }
    else
    {
        new_side.senses[1] = 0;
        new_side.names[1]  = "\0";
    }

    return new_side;
}

/*
 * given the string celldata, get the id number and name of each cell
 */
ReadRTT::cell ReadRTT::get_cell_data( std::string celldata )
{
    cell new_cell;
    std::vector< std::string > tokens;
    tokens = ReadRTT::split_string( celldata, ' ' );

    // set the side id
    if( tokens.size() != 2 )
    {
        MB_SET_ERR_RET_VAL( "Error, too many tokens found from cell_data", new_cell );
    }
    // create the new side
    new_cell.id   = std::atoi( tokens[0].c_str() );
    new_cell.name = tokens[1];

    return new_cell;
}

/*
 * given the string nodedata, get the id number and coordinates of the node
 */
ReadRTT::node ReadRTT::get_node_data( std::string nodedata )
{
    node new_node;
    std::vector< std::string > tokens;
    tokens = ReadRTT::split_string( nodedata, ' ' );

    // set the side id
    if( tokens.size() != 5 )
    {
        MB_SET_ERR_RET_VAL( "Error, too many tokens found from get_node_data", new_node );
    }
    new_node.id = std::atoi( tokens[0].c_str() );
    new_node.x  = std::atof( tokens[1].c_str() );
    new_node.y  = std::atof( tokens[2].c_str() );
    new_node.z  = std::atof( tokens[3].c_str() );
    return new_node;
}

/*
 * given the string nodedata, get the id number, connectivity and sense data
 */
ReadRTT::facet ReadRTT::get_facet_data( std::string facetdata )
{
    facet new_facet;
    std::vector< std::string > tokens;
    tokens = ReadRTT::split_string( facetdata, ' ' );

    // ensure we have the correct number of tokens
    int base_token_size = 0;
    int idx_offset      = 0;
    // branch on the rtt version number
    if( header_data.version == "v1.0.0" )
    {
        base_token_size = 4;
    }
    else if( header_data.version == "v1.0.1" )
    {
        base_token_size = 5;
        idx_offset      = 1;
    }
    else
    {
        MB_SET_ERR_RET_VAL( "Error, version number not understood", new_facet );
    }

    if( (int)tokens.size() != base_token_size + dim_data.nside_flag_types )
    {
        std::cout << facetdata << std::endl;
        std::cout << header_data.version << " " << (int)tokens.size() << " " << base_token_size << " "
                  << dim_data.nside_flag_types << std::endl;
        MB_SET_ERR_RET_VAL( "Error, too many tokens found from get_facet_data", new_facet );
        exit( 1 );
    }

    // set the side id
    new_facet.id              = std::atoi( tokens[0].c_str() );
    new_facet.connectivity[0] = std::atoi( tokens[idx_offset + 1].c_str() );
    new_facet.connectivity[1] = std::atoi( tokens[idx_offset + 2].c_str() );
    new_facet.connectivity[2] = std::atoi( tokens[idx_offset + 3].c_str() );
    new_facet.side_id         = std::atoi( tokens[idx_offset + 4].c_str() );
    new_facet.surface_number  = std::atoi( tokens[idx_offset + 5].c_str() );

    return new_facet;
}

/*
 * given the string tetdata, get the id number, connectivity and mat num of the
 * tet
 */
ReadRTT::tet ReadRTT::get_tet_data( std::string tetdata )
{
    tet new_tet;
    std::vector< std::string > tokens;
    tokens = ReadRTT::split_string( tetdata, ' ' );

    // ensure we have the correct number of tokens
    int base_token_size = 0;
    int n_nodes         = 0;
    // branch on the rtt version number
    if( header_data.version == "v1.0.0" )
    {
        base_token_size = 1;
        // for old format we assume 4 nodes (Tet)
        n_nodes = 4;
    }
    else if( header_data.version == "v1.0.1" )
    {
        base_token_size = 2;

        // get the cell type
        new_tet.type_id = std::atoi( tokens[1].c_str() );
        // check that the cell type exists
        if( cell_def_data.find( new_tet.type_id ) == cell_def_data.end() )
        {
            std::cout << "Error: cell type " << new_tet.type_id << " not found in cell definitions." << std::endl;
            MB_SET_ERR_RET_VAL( "Error, cell type not found in cell definitions", new_tet );
        }
        n_nodes     = cell_def_data[new_tet.type_id].nnodes;
        int n_sides = cell_def_data[new_tet.type_id].nsides;
        if( n_nodes != 4 || n_sides != 4 )
        {
            std::cout << "Error: cell type " << new_tet.type_id << " does not have 4 nodes and 4 sides." << std::endl;
            MB_SET_ERR_RET_VAL( "Error, cell is not tetrahedricon", new_tet );
        }
    }
    else
    {
        MB_SET_ERR_RET_VAL( "Error, version number not understood", new_tet );
    }
    // set the side id
    new_tet.id = std::atoi( tokens[0].c_str() );

    if( (int)tokens.size() != base_token_size + n_nodes + dim_data.ncell_flag_types )
    {
        MB_SET_ERR_RET_VAL( "Error, unexpected number of tokens found from get_tet_data", new_tet );
    }
    for( int i = 0; i < n_nodes; i++ )
    {
        new_tet.connectivity[i] = std::atoi( tokens[i + base_token_size].c_str() );
    }
    // set the material number
    new_tet.material_number = std::atoi( tokens[base_token_size + n_nodes].c_str() );

    return new_tet;
}

/*
 * splits string into sense and name, to later facilitate the building
 * of sense data, strips off the tailing @ if it exists
 */
ReadRTT::boundary ReadRTT::split_name( std::string atilla_cellname )
{
    boundary new_boundary;
    // default initialisation
    new_boundary.sense = 0;
    new_boundary.name  = "\0";
    // +ve sense
    if( atilla_cellname.find( "+" ) != std::string::npos )
    {
        new_boundary.sense = 1;
        // look for the @# we do not want it
        std::size_t found = atilla_cellname.find( "@" );
        if( found != std::string::npos )
            new_boundary.name = atilla_cellname.substr( 3, found );
        else
            new_boundary.name = atilla_cellname.substr( 3, atilla_cellname.length() );
    }
    else if( atilla_cellname.find( "-" ) != std::string::npos )
    {
        // negative sense
        new_boundary.sense = -1;
        new_boundary.name  = atilla_cellname.substr( 3, atilla_cellname.length() );
    }
    return new_boundary;
}

/*
 * splits a string in a vector of strings split by spaces
 */
std::vector< std::string > ReadRTT::split_string( std::string string_to_split, char split_char )
{
    std::istringstream ss( string_to_split );
    std::vector< std::string > tokens;
    while( !ss.eof() )
    {
        std::string x;                      // here's a nice, empty string
        std::getline( ss, x, split_char );  // try to read the next field into it
        tokens.push_back( x );
    }

    // remove empty tokens
    std::vector< std::string >::iterator it;
    for( it = tokens.begin(); it != tokens.end(); )
    {
        std::string string = *it;
        if( string.compare( "\0" ) == 0 )
            it = tokens.erase( it );
        else
            ++it;
    }
    return tokens;
}

/*
 * Generate the parent-child links bwetween the cell and surface meshsets
 */
void ReadRTT::generate_parent_child_links( int num_ents[4],
                                           std::vector< EntityHandle > entity_map[4],
                                           std::vector< side > side_data,
                                           std::vector< cell > cell_data )
{
    ErrorCode rval;  // return value
    // loop over the number of surfaces
    for( int i = 0; i < num_ents[2]; i++ )
    {
        // get the surface handle
        EntityHandle surf_handle = entity_map[2][i];
        // there are volumes that share this face
        for( unsigned int shared = 0; shared <= 1; shared++ )
        {
            std::string parent_name = side_data[i].names[shared];
            // find the @ sign
            unsigned pos = parent_name.find( "@" );
            parent_name  = parent_name.substr( 0, pos );

            // loop over tets looking for matching name
            for( int j = 0; j < num_ents[3]; j++ )
            {
                // if match found
                if( cell_data[j].name.compare( parent_name ) == 0 )
                {
                    EntityHandle cell_handle = entity_map[3][j];
                    // parent
                    rval = MBI->add_parent_child( cell_handle, surf_handle );
                    if( rval != MB_SUCCESS )
                    {
                        std::cerr << "Failed to add parent child relationship" << std::endl;
                    }
                }
            }
        }
    }
    return;
}

/*
 * sets the sense of the surfaces wrt to volumes using geom topo tool
 */
void ReadRTT::set_surface_senses( int num_ents[4],
                                  std::vector< EntityHandle > entity_map[4],
                                  std::vector< side > side_data,
                                  std::vector< cell > cell_data )
{

    ErrorCode rval;  // return value
    // loop over the number of surfaces
    for( int i = 0; i < num_ents[2]; i++ )
    {
        EntityHandle surf_handle = entity_map[2][i];
        // there are 2 volumes that share this face
        for( unsigned int shared = 0; shared <= 1; shared++ )
        {
            std::string parent_name = side_data[i].names[shared];
            unsigned pos            = parent_name.find( "@" );
            parent_name             = parent_name.substr( 0, pos );
            // loop over tets looking for matching name
            for( int j = 0; j < num_ents[3]; j++ )
            {
                // if match found
                if( cell_data[j].name.compare( parent_name ) == 0 )
                {
                    EntityHandle cell_handle = entity_map[3][j];
                    // in rtt mesh +represents the inside and -represents outside
                    // in moab reverse is outside and forward is inside
                    if( side_data[i].senses[shared] == 1 )
                        rval = myGeomTool->set_sense( surf_handle, cell_handle, SENSE_FORWARD );
                    else if( side_data[i].senses[shared] == -1 )
                        rval = myGeomTool->set_sense( surf_handle, cell_handle, SENSE_REVERSE );
                    else
                        rval = myGeomTool->set_sense( surf_handle, 0, SENSE_REVERSE );

                    if( rval != MB_SUCCESS )
                    {
                        std::cerr << "Failed to set sense appropriately" << std::endl;
                    }
                }
            }
        }
    }
    return;
}

/*
 * Add all entities that are to be part of the graveyard
 */
ErrorCode ReadRTT::setup_group_data( std::vector< EntityHandle > entity_map[4] )
{
    ErrorCode rval;  // error codes
    EntityHandle handle;
    handle = create_group( "graveyard_comp", 1 );

    // add any volume to group graveyard, it is ignored by dag
    EntityHandle vol_handle = entity_map[3][0];
    rval                    = MBI->add_entities( handle, &vol_handle, 1 );
    return rval;
}

/*
 * create a new group of with the name group name and id
 */
EntityHandle ReadRTT::create_group( std::string group_name, int id )
{
    ErrorCode rval;
    // category tags
    const char geom_categories[][CATEGORY_TAG_SIZE] = { "Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0" };

    EntityHandle handle;
    rval = MBI->create_meshset( MESHSET_SET, handle );
    if( MB_SUCCESS != rval ) return rval;

    rval = MBI->tag_set_data( name_tag, &handle, 1, group_name.c_str() );
    if( MB_SUCCESS != rval ) return MB_FAILURE;

    rval = MBI->tag_set_data( id_tag, &handle, 1, &id );
    if( MB_SUCCESS != rval ) return MB_FAILURE;

    rval = MBI->tag_set_data( category_tag, &handle, 1, &geom_categories[4] );
    if( MB_SUCCESS != rval ) return MB_FAILURE;

    return handle;
}

}  // namespace moab
