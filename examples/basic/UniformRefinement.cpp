/** @example UniformRefinement.cpp
 * @brief Example demonstrating uniform mesh refinement using MOAB's AHF data structures.
 *
 * This example shows how to:
 * - Load a mesh file
 * - Use the NestedRefine class for uniform mesh refinement
 * - Control refinement levels and degrees
 * - Write the refined mesh hierarchy to a file
 *
 * @note If no refinement degrees are specified, two levels are used by default.
 *
 * @code
 * Usage: ./UniformRefinement [filename] [level1_degree] [level2_degree] ...
 * Example: ./UniformRefinement input.h5m 2 3 2  # Three levels with degrees 2, 3, 2
 * @endcode
 */

#include "moab/Core.hpp"
#include "moab/NestedRefine.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>

// Using declarations for cleaner code
using moab::Core;
using moab::EntityHandle;
using moab::ErrorCode;
using moab::NestedRefine;

namespace
{
constexpr const char* DEFAULT_OUTPUT_FILE = "mesh_hierarchy.h5m";
constexpr int DEFAULT_REFINEMENT_LEVELS   = 2;

void print_usage( const char* program_name )
{
    std::cout << "Usage: " << program_name << " filename [level1_degree level2_degree ...]\n"
              << "  filename - Path to the input mesh file\n"
              << "  levelN_degree - (Optional) Degree of refinement for each level\n"
              << "                  If not specified, uses " << DEFAULT_REFINEMENT_LEVELS
              << " levels with default degrees\n";
}

std::vector< int > parse_refinement_levels( int argc, char* argv[] )
{
    std::vector< int > level_degrees;

    if( argc > 2 )
    {
        level_degrees.reserve( argc - 2 );
        for( int i = 2; i < argc; ++i )
        {
            try
            {
                int degree = std::stoi( argv[i] );
                if( degree < 1 )
                {
                    throw std::invalid_argument( "Refinement degree must be at least 1" );
                }
                level_degrees.push_back( degree );
            }
            catch( const std::exception& e )
            {
                throw std::runtime_error( "Invalid refinement degree: " + std::string( e.what() ) );
            }
        }
    }
    else
    {
        // Default: two levels with default degrees
        level_degrees.resize( DEFAULT_REFINEMENT_LEVELS );
    }

    return level_degrees;
}
}  // namespace

int main( int argc, char* argv[] )
{
    try
    {
        // Parse command line arguments
        if( argc < 2 || std::string( argv[1] ) == "--help" || std::string( argv[1] ) == "-h" )
        {
            print_usage( argv[0] );
            return argc < 2 ? 1 : 0;
        }

        const std::string input_file = argv[1];
        auto level_degrees           = parse_refinement_levels( argc, argv );
        const int num_levels         = static_cast< int >( level_degrees.size() );

        std::cout << "Input file: " << input_file << "\n"
                  << "Refinement levels: " << num_levels << "\n"
                  << "Level degrees: ";

        for( const auto& degree : level_degrees )
        {
            std::cout << degree << " ";
        }
        std::cout << "\n";

        // Initialize MOAB and load the mesh
        auto moab = std::make_unique< Core >();
        if( !moab )
        {
            throw std::runtime_error( "Failed to create MOAB instance" );
        }

        // Load the input mesh file
        MB_CHK_SET_ERR_CONT( moab->load_file( input_file.c_str() ), "Failed to load mesh file: " << input_file );

        // Set up the refinement
        NestedRefine refiner( moab.get() );

        std::vector< EntityHandle > mesh_sets;
        std::cout << "Starting mesh hierarchy generation for " << num_levels << " levels...\n";
        // Perform the refinement
        MB_CHK_SET_ERR( refiner.generate_mesh_hierarchy( num_levels, level_degrees.data(), mesh_sets ),
                        "Failed to generate mesh hierarchy" );
        std::cout << "Successfully generated mesh hierarchy\n";

        // Write the refined mesh to file
        MB_CHK_SET_ERR( moab->write_file( DEFAULT_OUTPUT_FILE ),
                        "Failed to write output file: " << DEFAULT_OUTPUT_FILE );

        std::cout << "Output written to: " << DEFAULT_OUTPUT_FILE << "\n";
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
