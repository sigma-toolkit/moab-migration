/** @example HelloMOAB.cpp
 * @brief A simple example demonstrating basic MOAB functionality.
 *
 * This example demonstrates the fundamental operations in MOAB:
 * - Creating a MOAB instance
 * - Loading a mesh from a file
 * - Querying entities by type and dimension
 * - Basic error handling with MB_CHK_ERR
 *
 * The program reads a mesh file (default: 3k-tri-sphere.vtk) and reports
 * the number of vertices, edges, faces, and elements in the mesh.
 *
 * To run: ./HelloMOAB [meshfile]
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include <iostream>
#include <memory>
#include <string>

// Using declarations for cleaner code
using moab::Core;
using moab::ErrorCode;
using moab::MBEDGE;
using moab::MBVERTEX;
using moab::Range;

// Constants
namespace
{
const char* const DEFAULT_MESH_FILE = "3k-tri-sphere.vtk";
}  // namespace

int main( int argc, char** argv )
{
    try
    {
        // Create MOAB instance with smart pointer for automatic cleanup
        auto moab = std::make_unique< Core >();
        if( !moab )
        {
            std::cerr << "Error: Failed to create MOAB instance\n";
            return 1;
        }

        // Parse input filename
        const std::string mesh_file =
            ( argc > 1 ) ? argv[1] : std::string( MESH_DIR ) + "/" + std::string( DEFAULT_MESH_FILE );

        // Load the mesh from file with error checking
        MB_CHK_SET_ERR( moab->load_mesh( mesh_file.c_str() ), "Error: Could not load mesh file: " << mesh_file );

        // Get entities by type/dimension
        Range vertices, edges, faces, elements;
        MB_CHK_ERR( moab->get_entities_by_type( 0, MBVERTEX, vertices ) );
        MB_CHK_ERR( moab->get_entities_by_type( 0, MBEDGE, edges ) );
        MB_CHK_ERR( moab->get_entities_by_dimension( 0, 2, faces ) );
        MB_CHK_ERR( moab->get_entities_by_dimension( 0, 3, elements ) );

        // Output the number of entities with formatted output
        const auto print_entity_count = []( const std::string& name, const auto& range ) {
            std::cout << "Number of " << name << ": " << range.size() << "\n";
        };

        print_entity_count( "vertices", vertices );
        print_entity_count( "edges   ", edges );
        print_entity_count( "faces   ", faces );
        print_entity_count( "elements", elements );

        return 0;
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    catch( ... )
    {
        std::cerr << "Error: Unknown exception occurred\n";
        return 1;
    }
}
