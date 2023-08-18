/** @example ExchangeHalos.cpp
 * \brief Example program that shows the use case for performing tag data exchange
 * between parallel processors in order to sync data on shared entities. The reduction
 * operation on tag data is also shown where the user can perform any of the actions supported
 * by MPI_Op on data residing on shared entities. \n
 *
 * <b>This example </b>:
 *    -# Initialize MPI and instantiate MOAB
 *    -# Get user options: Input mesh file name, tag name (default: USERTAG), tag value
 * (default: 1.0)
 *    -# Create the root and partition sets
 *    -# Instantiate ParallelComm and read the mesh file in parallel using appropriate options
 *    -# Create two tags: USERTAG_EXC (exchange) and USERTAG_RED (reduction)
 *    -# Set tag data and exchange shared entity information between processors
 *      -# Get entities in all dimensions and set local (current rank, dimension) dependent data for
 *     exchange tag (USERTAG_EXC)
 *      -# Perform exchange of tag data so that data on shared entities are synced via
 * ParallelCommunicator.
 *    -#  Set tag data and reduce shared entity information between processors using MPI_SUM
 *      -#  Get higher dimensional entities in the current partition and set local (current rank)
 *     dependent data for reduce tag (USERTAG_EXC)
 *      -#  Perform the reduction operation (MPI_SUM) on shared entities via ParallelCommunicator.
 *    -#  Destroy the MOAB instance and finalize MPI
 *
 * <b>To run:</b> \n mpiexec -n 2 ./ExchangeHalos <mesh_file> <tag_name> <tag_value> \n
 * <b>Example:</b> \n mpiexec -n 2 ./ExchangeHalos ../MeshFiles/unittest/64bricks_1khex.h5m
 * USERTAG 100 \n
 *
 */

// MOAB includes
#include "moab/Core.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"

#ifndef MOAB_HAVE_MPI
#error "Please build MOAB with MPI..."
#endif

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

// C++ includes
#include <iostream>
#include <string>
#include <sstream>

using namespace moab;
using namespace std;

#define dbgprinti( MSG )                    \
    do                                      \
    {                                       \
        if( !proc_id ) cout << MSG << endl; \
    } while( false )

#define dbgprint( MSG )                             \
    do                                              \
    {                                               \
        if( !context.proc_id ) cout << MSG << endl; \
    } while( false )

#define dbgprintall( MSG )                                      \
    do                                                          \
    {                                                           \
        cout << "[" << context.proc_id << "]: " << MSG << endl; \
    } while( false )

#define runchk( CODE, MSG )         \
    do                              \
    {                               \
        ErrorCode err = CODE;       \
        MB_CHK_SET_ERR( err, MSG ); \
    } while( false )

#define runchk0( CODE, MSG )        \
    do                              \
    {                               \
        ErrorCode err = CODE;       \
        if( err ) dbgprinti( MSG ); \
        MB_CHK_ERR_CONT( err );     \
    } while( false )

struct RuntimeContext
{
  public:
    int dimension{ 2 };           /// dimension of the problem
    std::string input_filename;   /// input file name (nc format)
    std::string output_filename;  /// output file name (h5m format)
    int ghost_layers{ 2 };        /// number of ghost layers
    std::string scalar_tagname;   /// scalar tag name
    std::string vector_tagname;   /// vector tag name
    int vector_length{ 2 };       /// length of the vector tag components
    int num_max_exchange{ 10 };   /// total number of exchange iterations
    bool debug_output{ false };   /// write debug output information?
    int proc_id;                  /// process identifier
    int num_procs;                /// total number of processes
    double last_counter{};        /// last time counter between push/pop timer

    // MOAB objects
    Interface* moab_interface{};
    ParallelComm* parallel_communicator;
    EntityHandle fileset{}, partnset{};

    /// @brief Constructor: allocate MOAB interface and communicator, and initialize
    /// other data members with some default values
    RuntimeContext( MPI_Comm comm = MPI_COMM_WORLD )
        : input_filename( string( MESH_DIR ) + string( "/io/mpasx1.642.t.2.nc" ) ),
          output_filename( "exchangeHalos_output.h5m" ), scalar_tagname( "h_s" ), vector_tagname( "ke" )
    {
        // Create the moab instance
        moab_interface = new( std::nothrow ) Core;
        if( NULL == moab_interface ) exit( 1 );

        // Create sets for the mesh and partition.  Then pass these to the load_file functions to populate the mesh.
        runchk0( moab_interface->create_meshset( MESHSET_SET, fileset ), "Creating root set failed" );
        runchk0( moab_interface->create_meshset( MESHSET_SET, partnset ), "Creating partition set failed" );

        // Create the parallel communicator object with the partition handle associated with MOAB
        parallel_communicator = ParallelComm::get_pcomm( moab_interface, partnset, &comm );

        timer = new moab::CpuTimer();

        proc_id   = parallel_communicator->rank();
        num_procs = parallel_communicator->size();
    }

    /// @brief Destructor: deallocate MOAB interface and communicator
    ~RuntimeContext()
    {
        delete timer;
        delete parallel_communicator;
        delete moab_interface;
    }

    /// @brief Parse the runtime command line options
    /// @param argc - number of command line arguments
    /// @param argv - command line arguments as string list
    void ParseCLOptions( int argc, char* argv[] )
    {
        ProgOptions opts;
        // Input mesh
        opts.addOpt< std::string >( "input", "Input mesh filename to load in parallel", &input_filename );
        // Output mesh
        opts.addOpt< std::string >(
            "output", "Output mesh filename for verification (default=exchangeHalos_output.h5m)", &output_filename );
        // Dimension of the input mesh
        opts.addOpt< int >( "dimension", "Input mesh dimension (default = 2)", &dimension );
        // Scalar and Vector tag names
        opts.addOpt< std::string >( "stag", "Scalar tag name to exchange with neighboring tasks (default=h_s)",
                                    &scalar_tagname );
        opts.addOpt< std::string >( "vtag", "Vector tag name to exchange with neighboring tasks (default=ke)",
                                    &vector_tagname );
        // Vector tag length
        opts.addOpt< int >( "vtaglength", "Size of vector components per each entity (default=2)", &vector_length );
        // Number of halo (ghost) regions
        opts.addOpt< int >( "nghosts", "Number of ghost layers (halos) to exchange (default=2)", &ghost_layers );
        // Number of times to perform the halo exchange for timing
        opts.addOpt< int >( "nexchanges", "Number of ghost-halo exchange iterations to perform (default=10)",
                            &num_max_exchange );
        opts.addOpt< void >( "debug", "Should we write output file? (default=false)", &debug_output );

        opts.parseCommandLine( argc, argv );
    }

    void timer_push( std::string operation )
    {
        timer_ops = timer->time_since_birth();
        opName    = operation;
    }

    void timer_pop( int nruns = 1 )
    {
        double locElapsed = timer->time_since_birth() - timer_ops;
        double avgElapsed = 0;
        double maxElapsed = 0;
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, parallel_communicator->comm() );
        MPI_Reduce( &locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, parallel_communicator->comm() );
        if( proc_id == 0 )
        {
            avgElapsed /= num_procs;
            if( nruns > 1 )
                std::cout << "[LOG] Time taken to " << opName.c_str() << ", averaged over " << nruns
                          << " runs : max = " << maxElapsed / nruns << ", avg = " << avgElapsed / nruns << "\n";
            else
                std::cout << "[LOG] Time taken to " << opName.c_str() << " : max = " << maxElapsed
                          << ", avg = " << avgElapsed << "\n";

            last_counter = maxElapsed / nruns;
        }
        opName.clear();
    }

    void load_file() const;

    inline double last_elapsed() const
    {
        return last_counter;
    }

    moab::ErrorCode create_sv_tags( Tag& tagScalar, Tag& tagVector, Range& entities ) const;

  private:
    /// @brief Compute the centroids of elements in 2D lat/lon space
    /// @param ents
    /// @return centroids (as lat/lon)
    std::vector< double > compute_centroids( const Range& ents ) const;

    moab::CpuTimer* timer;
    double timer_ops{};
    std::string opName;
};

moab::ErrorCode RuntimeContext::create_sv_tags( Tag& tagScalar, Tag& tagVector, Range& entities ) const
{
    // Get element (centroid) coordinates so that we can evaluate some arbitrary data
    std::vector< double > entCoords = compute_centroids( entities );  // [entities * [lon, lat]]

    dbgprinti( "> Getting scalar tag handle " << scalar_tagname << "..." );
    double defSTagValue = -1.0;
    bool createdTScalar = false;
    // Create the exchange tag: default name = USERTAG_EXC
    runchk( moab_interface->tag_get_handle( scalar_tagname.c_str(), 1, MB_TYPE_DOUBLE, tagScalar,
                                            MB_TAG_CREAT | MB_TAG_DENSE, &defSTagValue, &createdTScalar ),
            "Retrieving scalar tag handle failed" );

    if( createdTScalar )
    {
        std::vector< double > tagValues( entities.size(), -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = index * 2;

            // double value =
            //     ( 2.0 + cos( entCoords[offset] ) * cos( entCoords[offset] ) * cos( 2.0 * entCoords[offset + 1] ) );
            double value =
                ( 2.0 + std::pow( sin( 2.0 * entCoords[offset + 1] ), 16.0 ) * cos( 16.0 * entCoords[offset] ) );

            index++;
            return value;
        } );
        // Set local scalar tag data for exchange
        runchk( moab_interface->tag_set_data( tagScalar, entities, tagValues.data() ),
                "Setting scalar tag data failed" );
    }

    dbgprinti( "> Getting vector tag handle " << vector_tagname << "..." );
    std::vector< double > defVTagValue( vector_length, -1.0 );
    bool createdTVector = false;
    // Create the exchange tag: default name = USERTAG_RED
    runchk( moab_interface->tag_get_handle( vector_tagname.c_str(), vector_length, MB_TYPE_DOUBLE, tagVector,
                                            MB_TAG_CREAT | MB_TAG_DENSE, defVTagValue.data(), &createdTVector ),
            "Retrieving vector tag handle failed" );

    if( createdTVector )
    {
        const int veclength = vector_length;
        std::vector< double > tagValues( entities.size() * veclength, -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = ( index / veclength ) * 2;

            double value =
                ( 2.0 + cos( entCoords[offset] ) * cos( entCoords[offset] ) * cos( 2.0 * entCoords[offset + 1] ) ) *
                ( index % veclength + 1.0 );  // assign some scalar multiple value for different vector components

            index++;
            return value;
        } );
        // Set local tag data for exchange
        runchk( moab_interface->tag_set_data( tagVector, entities, tagValues.data() ),
                "Setting vector tag data failed" );
    }

    return moab::MB_SUCCESS;
}

void RuntimeContext::load_file() const
{
    /// Parallel Read options:
    ///   PARALLEL = type {READ_PART}
    ///   PARTITION = PARALLEL_PARTITION : Partition as you read
    ///   PARALLEL_RESOLVE_SHARED_ENTS : Communicate to all processors to get the shared adjacencies
    ///   consistently in parallel PARALLEL_GHOSTS : a.b.c
    ///                   : a = 3 - highest dimension of entities
    ///                   : b = 0 -
    ///                   : c = 1 - number of layers
    string read_options        = "DEBUG_IO=0;";
    std::string::size_type idx = input_filename.rfind( '.' );
    std::string extension      = "";
    if( num_procs > 1 && idx != std::string::npos )
    {
        extension = input_filename.substr( idx + 1 );
        if( !extension.compare( "nc" ) )
            read_options += "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;"
                            "PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=;";  // NO_EDGES;NO_MIXED_ELEMENTS;RCBZOLTAN, TRIVIAL
        else if( !extension.compare( "h5m" ) )
            read_options += "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                            "PARALLEL_RESOLVE_SHARED_ENTS;";
        else
            read_options += "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;"
                            "PARALLEL_RESOLVE_SHARED_ENTS;";
    }

    // Load the file from disk with given options
    runchk0( moab_interface->load_file( input_filename.c_str(), &fileset, read_options.c_str() ),
             "MOAB::load_file failed" );
}

std::vector< double > RuntimeContext::compute_centroids( const Range& ents ) const
{
    double node[3];
    std::vector< double > eCentroids( ents.size() * 2 );  // [lon, lat]
    for( size_t ients = 0, offset = 0; ients < ents.size(); ++ients, offset += 2 )
    {
        const EntityHandle entity = ents[ients];
        runchk0( moab_interface->get_coords( &entity, 1, node ), "Getting entity coordinates failed" );

        // scale by magnitude so that mesh is on unit sphere
        double magnitude = std::sqrt( node[0] * node[0] + node[1] * node[1] + node[2] * node[2] );
        node[0] /= magnitude;
        node[1] /= magnitude;
        node[2] /= magnitude;

        // compute the spherical transformation onto unit sphere
        eCentroids[offset] = atan2( node[1], node[0] );
        if( eCentroids[offset] < 0.0 ) eCentroids[offset] += 2.0 * M_PI;
        eCentroids[offset + 1] = asin( node[2] );
    }
    return eCentroids;
}

//
// Start of main test program
//
int main( int argc, char** argv )
{
    // Initialize MPI first
    MPI_Init( &argc, &argv );

    {
        RuntimeContext context;
        dbgprint( "********** Exchange halos example **********\n" );

        // Get the input options
        context.ParseCLOptions( argc, argv );

        /////////////////////////////////////////////////////////////////////////
        // Print out the input parameters in use
        dbgprint( " -- Input Parameters -- " );
        dbgprint( "    Number of Processes  = " << context.num_procs );
        dbgprint( "    Input mesh           = " << context.input_filename );
        dbgprint( "    Ghost Layers         = " << context.ghost_layers );
        dbgprint( "    Scalar Tag name      = " << context.scalar_tagname );
        dbgprint( "    Vector Tag name      = " << context.vector_tagname );
        dbgprint( "    Vector Tag length    = " << context.vector_length << endl );
        /////////////////////////////////////////////////////////////////////////

        // Timer storage for all phases
        double elapsed_times[4];

        context.timer_push( "Read input file" );
        // Load the file from disk with given options
        context.load_file();
        context.timer_pop();
        elapsed_times[0] = context.last_elapsed();

        dbgprint( "- " );

        context.timer_push( "Setup ghost layers" );
        // Ensure that all processes understand about multi-shared vertices and entities
        // in case some adjacent parts are only m layers thick (where m < context.ghost_layers)
        // runchk( context.parallel_communicator->correct_thin_ghost_layers(), "Thin layer correction failed" );

        // Exchange ghost cells
        int ghost_dimension  = context.dimension;
        int bridge_dimension = context.dimension - 1;
        // Let us now get all ghost layers from adjacent parts
        runchk( context.parallel_communicator->exchange_ghost_cells(
                    ghost_dimension, bridge_dimension, context.ghost_layers, 0, true /* store_remote_handles */,
                    true /* wait_all */, &context.fileset ),
                "Exchange ghost cells failed" );  // true to store remote handles

        // Mesh is now loaded and ghost cells are available on each task.
        // Ensure to augment the ghost cells with essential tag data such as MATERIAL_SET etc if we need them
        // runchk( context.parallel_communicator->augment_default_sets_with_ghosts( fileset ), "Ghost cell data augment failed");
        context.timer_pop();
        elapsed_times[1] = context.last_elapsed();

        Range dimEnts;
        // Get all entities of dimension = dim
        runchk( context.moab_interface->get_entities_by_dimension( context.fileset, context.dimension, dimEnts ),
                "Getting 2D entities failed" );
        // Get only owned entities! The ghosted/shared entities will get their data when we exchange
        runchk( context.parallel_communicator->filter_pstatus( dimEnts, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                "Filtering pstatus failed" );

        // Aggregate the total number of elements in the mesh
        auto numEntities = dimEnts.size();
        // dbgprintall( " number of " << context.dimension << "D elements in local mesh = " << numEntities );
        int numTotalEntities = 0;
        MPI_Reduce( &numEntities, &numTotalEntities, 1, MPI_INT, MPI_SUM, 0,
                    context.parallel_communicator->proc_config().proc_comm() );
        dbgprint( "Total number of " << context.dimension << "D elements in the mesh = " << numTotalEntities );

        // Create two tag handles: Exchange and Reduction operations
        Tag tagScalar = nullptr;
        Tag tagVector = nullptr;
        runchk( context.create_sv_tags( tagScalar, tagVector, dimEnts ), "Unable to create scalar and vector tags" );

        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "- Writing to file *before* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_pre.h5m", "H5M", "DEBUG_IO=0;" ),
                    "Writing to disk failed" );
        }

        // Perform exchange tag data
        dbgprint( "> Exchanging tags between processors " );
        context.timer_push( "Exchange scalar tag data" );
        for( auto irun = 0; irun < context.num_max_exchange; ++irun )
        {
            // Exchange scalar tags between processors
            runchk( context.parallel_communicator->exchange_tags( tagScalar, dimEnts ),
                    "Exchanging scalar tag between processors failed" );
        }
        context.timer_pop( context.num_max_exchange );
        elapsed_times[2] = context.last_elapsed();

        context.timer_push( "Exchange vector tag data" );
        for( auto irun = 0; irun < context.num_max_exchange; ++irun )
        {
            // Exchange vector tags between processors
            runchk( context.parallel_communicator->exchange_tags( tagVector, dimEnts ),
                    "Exchanging vector tag between processors failed" );
        }
        context.timer_pop( context.num_max_exchange );
        elapsed_times[3] = context.last_elapsed();

        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "- Writing to file *after* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_post.h5m", "H5M", "DEBUG_IO=0;" ),
                    "Writing to disk failed" );
        }

        if( context.debug_output )
        {
            dbgprint( "> Writing out the final mesh and data in MOAB h5m format. File = " << context.output_filename );
            string write_options = ( context.num_procs > 1 ? "PARALLEL=WRITE_PART;DEBUG_IO=0;" : "" );
            // Write out to output file to visualize reduction/exchange of tag data
            runchk( context.moab_interface->write_file( context.output_filename.c_str(), "H5M", write_options.c_str() ),
                    "File write failed" );
        }

        dbgprint( "> Consolidated: [" << context.num_procs << ", " << context.ghost_layers << ", " << elapsed_times[0]
                                      << ", " << elapsed_times[1] << ", " << elapsed_times[2] << ", "
                                      << elapsed_times[3] << "]," );

        dbgprint( "\n********** ExchangeHalos Example DONE! **********" );
    }
    // Done, cleanup
    MPI_Finalize();

    return 0;
}
