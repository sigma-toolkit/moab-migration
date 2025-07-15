/** @example UniformRefinement.cpp
 * This example demonstrates uniform mesh refinement using MOAB's AHF data structures.
 * It shows how to load a mesh file, use the NestedRefine class to perform uniform mesh refinement,
 * control the number of refinement levels and degree at each level,
 * and write the refined mesh hierarchy to a file.
 *
 * Usage:
 *   ./UniformRefinement [filename] [level1_degree] [level2_degree] ...
 *
 * If no refinement degrees are specified, two levels are used by default.
 */

#include <iostream>
#include "moab/Core.hpp"
#include "moab/NestedRefine.hpp"

using namespace moab;

int main( int argc, char* argv[] )
{
    Core mb;
    Interface* mbImpl = &mb;
    ErrorCode error;

    if( argc == 1 )
    {
        std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
        return 1;
    }

    const char* filename = argv[1];
    error                = mbImpl->load_file( filename );
    MB_CHK_ERR( error );

    NestedRefine uref( &mb );

    // Usage: The level degrees array controls number of refinemetns and
    // the degree of refinement at each level.
    // Example: int level_degrees[4] = {2,3,2,3};
    std::vector< int > level_degrees;
    if( argc > 2 )
    {
        level_degrees.resize( argc - 2 );
        for( int i = 2; i < argc; ++i )
            level_degrees[i - 2] = atoi( argv[i] );
    }
    else
    {
        level_degrees.resize( 2 );
    }
    int num_levels = static_cast< int >( level_degrees.size() );

    std::cout << "Starting hierarchy generation" << std::endl;
    std::vector< EntityHandle > set;
    error = uref.generate_mesh_hierarchy( num_levels, level_degrees.data(), set );
    MB_CHK_ERR( error );
    std::cout << "Finished hierarchy generation" << std::endl;

    std::stringstream file;
    file << "mesh_hierarchy.h5m";
    error = mbImpl->write_file( file.str().c_str() );
    MB_CHK_ERR( error );
    return 0;
}
