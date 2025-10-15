/** @example ReadWriteTest.cpp
 * @brief Example demonstrating parallel mesh reading and writing with performance timing.
 *
 * This example shows how to:
 * - Read a mesh file in parallel with specific options
 * - Write a mesh file in parallel with specific options
 * - Measure and report read/write performance times
 * - Handle parallel mesh operations with shared entities
 * - Use different partition methods for parallel processing
 *
 * The example is designed for stress testing of MOAB's parallel
 * reader/writer capabilities and provides timing information
 * for performance analysis.
 *
 * @note This example requires MOAB to be compiled with MPI support.
 *
 * @par Usage:
 * @code
 * mpiexec -np <num_procs> ReadWriteTest [input] [output] -O <read_opts> -o <write_opts>
 * @endcode
 *
 * @par Example:
 * @code
 * mpiexec -np 4 ReadWriteTest input.nc output.h5m \
 *   -O "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=T,U" \
 *   -o "PARALLEL=WRITE_PART;VARIABLE=T,U"
 * @endcode
 */

#include "moab/Core.hpp"
#include <iostream>
#include <chrono>
#include <string>
#include <memory>
#include <cstdlib>

// Only include MPI-specific headers if MOAB was built with MPI support
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include <mpi.h>
#endif

// Using declarations for cleaner code
using moab::Core;
using moab::EntityHandle;
using moab::ErrorCode;
using moab::ParallelComm;

namespace
{
// Default values for command line arguments
struct Config
{
    std::string input_file;
    std::string output_file;
    std::string read_opts;
    std::string write_opts;
    bool use_defaults = false;
};

// Print usage information
void print_usage( const char* program_name )
{
    std::cout << "Usage: mpiexec -np <num_procs> " << program_name << " [input_file] [output_file] [options]\n"
              << "Options:\n"
              << "  -O <read_opts>  Options for reading the input file\n"
              << "  -o <write_opts> Options for writing the output file\n"
              << "\nExample:\n"
              << "  mpiexec -np 4 " << program_name << " input.nc output.h5m \\\n"
              << "    -O \"PARALLEL=READ_PART;PARTITION_METHOD=SQIJ\" \\n"
              << "    -o \"PARALLEL=WRITE_PART\"\n";
}

// Parse command line arguments
Config parse_arguments( int argc, char** argv )
{
    Config config;

    // Set default values
#ifdef MOAB_HAVE_NETCDF
    config.input_file   = std::string( MESH_DIR ) + "/io/fv3x46x72.t.3.nc";
    config.output_file  = "ReadWriteTestOut.h5m";
    config.read_opts    = "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=T,U";
    config.write_opts   = "PARALLEL=WRITE_PART";
    config.use_defaults = true;
#else
    if( argc < 3 )
    {
        print_usage( argv[0] );
        std::exit( 1 );
    }
#endif

    // Parse command line arguments if provided
    if( argc >= 3 )
    {
        config.input_file   = argv[1];
        config.output_file  = argv[2];
        config.use_defaults = false;
    }

    // Parse options
    for( int i = 3; i < argc; ++i )
    {
        std::string arg = argv[i];
        if( arg == "-O" && i + 1 < argc )
        {
            config.read_opts = argv[++i];
        }
        else if( arg == "-o" && i + 1 < argc )
        {
            config.write_opts = argv[++i];
        }
        else if( arg == "-h" || arg == "--help" )
        {
            print_usage( argv[0] );
            std::exit( 0 );
        }
    }

    return config;
}

// Print configuration information
void print_config( const Config& config, int rank, int nprocs )
{
    if( rank == 0 )
    {
        std::cout << "\n=== MOAB Parallel Read/Write Test ===\n"
                  << "Input file:  " << config.input_file << "\n"
                  << "Output file: " << config.output_file << "\n"
                  << "Read options: " << ( config.read_opts.empty() ? "(default)" : config.read_opts ) << "\n"
                  << "Write options: " << ( config.write_opts.empty() ? "(default)" : config.write_opts ) << "\n"
                  << "Running on " << nprocs << " MPI processes\n"
                  << std::endl;
    }
}
}  // namespace

int main( int argc, char** argv )
{
    // Initialize MPI
    int mpi_initialized = 0;
    MPI_Initialized( &mpi_initialized );
    if( !mpi_initialized )
    {
        MPI_Init( &argc, &argv );
    }

    int rank   = 0;
    int nprocs = 1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

    try
    {
        // Parse command line arguments
        Config config = parse_arguments( argc, argv );
        if( rank == 0 )
        {
            print_config( config, rank, nprocs );
        }

        // Create MOAB instance with smart pointer for automatic cleanup
        auto moab = std::make_unique< Core >();
        if( !moab )
        {
            std::cerr << "Error: Failed to create MOAB instance" << std::endl;
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }

        // Create parallel communication interface
        auto pcomm = std::make_unique< ParallelComm >( moab.get(), MPI_COMM_WORLD );

        // Create a mesh set to store the loaded mesh
        EntityHandle mesh_set;
        MB_CHK_ERR( moab->create_meshset( moab::MESHSET_SET, mesh_set ) );

        // Time the read operation
        auto read_start = std::chrono::high_resolution_clock::now();

        if( rank == 0 )
        {
            std::cout << "Reading mesh file..." << std::endl;
        }

        // Read the mesh file
        MB_CHK_SET_ERR( moab->load_file( config.input_file.c_str(), &mesh_set, config.read_opts.c_str() ),
                        "Failed to read input file: " << config.input_file );

        auto read_end                                = std::chrono::high_resolution_clock::now();
        std::chrono::duration< double > read_elapsed = read_end - read_start;

        if( rank == 0 )
        {
            std::cout << "Read completed in " << read_elapsed.count() << " seconds" << std::endl;
            std::cout << "Writing mesh file..." << std::endl;
        }

        // Time the write operation
        auto write_start = std::chrono::high_resolution_clock::now();

        // Write the mesh file
        MB_CHK_SET_ERR( moab->write_file( config.output_file.c_str(), nullptr, config.write_opts.c_str(), &mesh_set,
                                          1 ),
                        "Failed to write output file: " << config.output_file );

        auto write_end                                = std::chrono::high_resolution_clock::now();
        std::chrono::duration< double > write_elapsed = write_end - write_start;

        if( rank == 0 )
        {
            std::cout << "Write completed in " << write_elapsed.count() << " seconds" << std::endl;
            std::cout << "\n=== Test completed successfully ===" << std::endl;
        }

        // Clean up
        pcomm.reset();
        moab.reset();

        // Only finalize MPI if we initialized it
        if( !mpi_initialized )
        {
            MPI_Finalize();
        }

        return 0;
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }
    catch( ... )
    {
        std::cerr << "Unknown error occurred on rank " << rank << std::endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }
}
