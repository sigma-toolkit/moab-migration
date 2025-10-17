/** @example GeomSetHierarchy.cpp
 * This example demonstrates geometric topology hierarchy traversal.
 * It shows how to load a mesh file with geometric topology information,
 * traverse the geometric hierarchy from regions down to vertices,
 * access parent-child relationships between geometric entities,
 * and determine geometric sense (forward/reverse) between entities.
 * This tool is useful for understanding the geometric structure
 * of CAD-based meshes and their hierarchical relationships.
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBCN.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include <iostream>
#include <cassert>

const char* ent_names[] = { "Vertex", "Edge", "Face", "Region" };

int main( int argc, char** argv )
{
    if( 1 == argc )
    {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 0;
    }

    // instantiate & load a file
    moab::Interface* mb = new moab::Core();
    MB_CHK_SET_ERR( mb->load_file( argv[1] ), "Failed to load file" );

    // get the geometric topology tag handle
    moab::Tag geom_tag, gid_tag;
    MB_CHK_SET_ERR( mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geom_tag ),
                    "Failed to get geometric dimension tag" );
    gid_tag = mb->globalId_tag();
    assert( NULL != gid_tag );

    // traverse the model, from dimension 3 downward
    moab::Range psets, chsets;
    std::vector< moab::EntityHandle > sense_ents;
    std::vector< int > senses;
    int dim, pgid, chgid;
    void* dim_ptr = &dim;
    int sense;

    moab::GeomTopoTool gt( mb, true );

    for( dim = 3; dim >= 0; dim-- )
    {
        // get parents at this dimension
        chsets.clear();
        MB_CHK_SET_ERR( mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, &dim_ptr, 1, chsets, 1, false ),
                        "Failed to get entities by type and tag" );

        // for each child, get parents and do something with them
        moab::Range::iterator ch_it, p_it;
        for( ch_it = chsets.begin(); ch_it != chsets.end(); ++ch_it )
        {
            // get the children and put in child set list
            psets.clear();
            MB_CHK_SET_ERR( mb->get_parent_meshsets( *ch_it, psets ), "Failed to get parent meshsets" );

            MB_CHK_SET_ERR( mb->tag_get_data( gid_tag, &( *ch_it ), 1, &chgid ), "Failed to get tag data" );

            // print # parents
            std::cout << ent_names[dim] << " " << chgid << " has " << psets.size() << " parents." << std::endl;

            if( 2 == dim )
            {
                for( p_it = psets.begin(); p_it != psets.end(); ++p_it )
                {
                    MB_CHK_SET_ERR( mb->tag_get_data( gid_tag, &( *p_it ), 1, &pgid ), "Failed to get tag data" );
                    rval = gt.get_sense( *ch_it, *p_it, sense );
                    if( moab::MB_SUCCESS != rval ) continue;
                    std::cout << ent_names[dim + 1] << " " << pgid << ", " << ent_names[dim] << " " << chgid
                              << " sense is: ";
                    if( 1 == sense )
                        std::cout << "FORWARD" << std::endl;
                    else
                        std::cout << "REVERSE" << std::endl;
                }
            }
            else if( 1 == dim )
            {
                sense_ents.clear();
                senses.clear();
                rval = gt.get_senses( *ch_it, sense_ents, senses );
                if( moab::MB_SUCCESS != rval ) continue;
                for( unsigned int i = 0; i < sense_ents.size(); i++ )
                {
                    MB_CHK_SET_ERR( mb->tag_get_data( gid_tag, &sense_ents[i], 1, &pgid ), "Failed to get tag data" );
                    std::cout << ent_names[dim + 1] << " " << pgid << ", " << ent_names[dim] << " " << chgid
                              << " sense is: ";
                    if( -1 == senses[i] )
                        std::cout << "REVERSED" << std::endl;
                    else if( 0 == senses[i] )
                        std::cout << "BOTH" << std::endl;
                    else if( 1 == senses[i] )
                        std::cout << "FORWARD" << std::endl;
                    else
                        std::cout << "(invalid)" << std::endl;
                }
            }
        }
    }

    delete mb;

    return 0;
}
