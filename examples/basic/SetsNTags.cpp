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
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <string>

// Using declarations for cleaner code
using moab::Core;
using moab::ErrorCode;
using moab::Range;
using moab::Tag;

namespace
{
// Default mesh file path
const char* const DEFAULT_MESH_FILE = MESH_DIR "/hex01.vtk";

// Tag names for conventional tags from MBTagConventions.hpp
const std::vector< const char* > TAG_NAMES = { MATERIAL_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME };

// Print usage information
void print_usage( const char* program_name )
{
    std::cout << "Usage: " << program_name << " [meshfile]\n"
              << "  meshfile - Path to the mesh file (default: " << DEFAULT_MESH_FILE << ")\n";
}
}  // namespace

int main( int argc, char** argv )
{
    try
    {
        // Parse command line arguments
        if( argc > 1 && ( std::string( argv[1] ) == "-h" || std::string( argv[1] ) == "--help" ) )
        {
            print_usage( argv[0] );
            return 0;
        }

        const std::string mesh_file = ( argc > 1 ) ? argv[1] : std::string( DEFAULT_MESH_FILE );

        // Create MOAB instance with smart pointer for automatic cleanup
        auto moab = std::make_unique< Core >();
        if( !moab )
        {
            std::cerr << "Error: Failed to create MOAB instance\n";
            return 1;
        }

        // Load the mesh file
        MB_CHK_SET_ERR( moab->load_file( mesh_file.c_str() ), "Failed to load mesh file: " << mesh_file );

        Range sets, set_entities;
        Tag tag_handle;

        // Process each tag type
        for( const auto& tag_name : TAG_NAMES )
        {
            // Get the tag handle for this tag name
            MB_CHK_SET_ERR( moab->tag_get_handle( tag_name, 1, moab::MB_TYPE_INTEGER, tag_handle ),
                            "Failed to get tag handle for: " << tag_name );

            // Get all sets with this tag
            sets.clear();
            MB_CHK_SET_ERR( moab->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &tag_handle, nullptr, 1, sets ),
                            "Failed to get entities for tag: " << tag_name );

            std::cout << "\nFound " << sets.size() << " sets with tag: " << tag_name << "\n";

            // Process each set
            for( const auto& set_handle : sets )
            {
                // Get the set ID
                int set_id;
                MB_CHK_SET_ERR( moab->tag_get_data( tag_handle, &set_handle, 1, &set_id ),
                                "Failed to get tag data for set" );

                // Get all entities in the set (recursively)
                set_entities.clear();
                MB_CHK_SET_ERR( moab->get_entities_by_handle( set_handle, set_entities, true ),
                                "Failed to get entities in set" );

                std::cout << "  " << tag_name << " " << set_id << " contains " << set_entities.size() << " entities\n";

                // Print entity information
                set_entities.print( "   " );
                set_entities.clear();
            }
        }

        return 0;
    }
    catch( const std::exception& e )
    {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
    catch( ... )
    {
        std::cerr << "\nError: Unknown exception occurred\n";
        return 1;
    }
}
