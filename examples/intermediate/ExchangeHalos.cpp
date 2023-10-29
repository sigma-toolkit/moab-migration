/** @example ExchangeHalos Driver
 * \brief Example program that shows the use case for performing tag data exchange
 * between parallel processors in order to sync data on shared entities.
 *
 * <b>This example </b>:
 *    -# Initialize MPI and instantiates MOAB
 *    -# Gets user options: Input mesh file name, vector tag length, ghost layer size etc
 *    -# Create the root and partition sets
 *    -# Instantiate ParallelComm and read the mesh file in parallel using appropriate options
 *    -# Create the required number of ghost layers as requested by the user (default = 3)
 *    -# Get 2D MPAS polygonal entities in the mesh and filter to get only the "owned" entities
 *    -# Create two tags: scalar_variable (single data/cell) and vector_variable (multiple data/cell)
 *    -# Set tag data using analytical functions for both scalar and vector fields on owned entities
 *    -# Exchange shared entity information and tags between processors
 *      -# If debugging is turned on, store mesh file and tag on root process (will not contain data on shared entities)
 *      -# Perform exchange of scalar tag data on shared entities
 *      -# Perform exchange of vector tag data on shared entities
 *      -# If debugging is turned on, store mesh file and tag on root process (will now contain data on *all* entities)
 *    -#  Destroy the MOAB instance and finalize MPI
 *
 * <b>To run: </b>
 *      mpiexec -n np ./ExchangeHalos --input <mpas_mesh_file> --nghosts <ghostlayers> --vtaglength <vector component
 * size> \
 *                    --nexchanges <number of exchange runs>
 * <b>Example:</b>
 *      mpiexec -n 4 ./ExchangeHalos --input $MOAB_DIR/MeshFiles/unittest/io/mpasx1.642.t.2.nc --nghosts 3 --vtaglength 100
 *
 * NOTE: --debug option can be added to write out extra files in h5m format to visualize some output (written from root
 * task only)
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

using namespace moab;
using namespace std;
#define dbgprint( MSG )                                           \
    do                                                            \
    {                                                             \
        if( context.proc_id == 0 ) std::cout << MSG << std::endl; \
    } while( false )

#define runchk( CODE, MSG )         \
    do                              \
    {                               \
        moab::ErrorCode err = CODE; \
        MB_CHK_SET_ERR( err, MSG ); \
    } while( false )

#define runchk_cont( CODE, MSG )                               \
    do                                                         \
    {                                                          \
        moab::ErrorCode err = CODE;                            \
        MB_CHK_ERR_CONT( err );                                \
        if( err ) std::cout << "Error:: " << MSG << std::endl; \
    } while( false )

// Run time context structure

/// @brief The RunttimeContext is an example specific class to store
/// the run specific input data, MOAB datastructures used during the run
/// and provides other utility functions to profile operations etc
struct RuntimeContext
{
  public:
    int dimension{ 2 };           /// dimension of the problem
    std::string input_filename;   /// input file name (nc format)
    std::string output_filename;  /// output file name (h5m format)
    int ghost_layers{ 3 };        /// number of ghost layers
    std::string scalar_tagname;   /// scalar tag name
    std::string vector_tagname;   /// vector tag name
    int vector_length{ 3 };       /// length of the vector tag components
    int num_max_exchange{ 10 };   /// total number of exchange iterations
    bool debug_output{ false };   /// write debug output information?
    int proc_id{ 1 };             /// process identifier
    int num_procs{ 1 };           /// total number of processes
    double last_counter{ 0.0 };   /// last time counter between push/pop timer

    // MOAB objects
    moab::Interface* moab_interface{ nullptr };
    moab::ParallelComm* parallel_communicator{ nullptr };
    moab::EntityHandle fileset{ 0 }, partnset{ 0 };

    /// @brief Constructor: allocate MOAB interface and communicator, and initialize
    /// other data members with some default values
    explicit RuntimeContext( MPI_Comm comm = MPI_COMM_WORLD );

    /// @brief Destructor: deallocate MOAB interface and communicator
    ~RuntimeContext();

    /// @brief Parse the runtime command line options
    /// @param argc - number of command line arguments
    /// @param argv - command line arguments as string list
    void ParseCLOptions( int argc, char* argv[] );

    /// @brief Measure and start the timer to profile a task
    /// @param operation String name of the task being measured
    inline void timer_push( std::string operation );

    /// @brief Stop the timer and store the elapsed duration
    /// @param nruns Optional argument used to average the measured time
    void timer_pop( const int nruns = 1 );

    /// @brief Return the last elapsed time
    /// @return last_counter from timer_pop was called
    inline double last_elapsed() const;

    /// @brief Load a MOAB supported file (h5m or nc format) from disk
    ///        representing an MPAS mesh
    /// @param load_ghosts Optional boolean to specify whether to load ghosts
    ///                    when reading the file (only relevant for h5m)
    /// @return Error code if any (else MB_SUCCESS)
    moab::ErrorCode load_file( bool load_ghosts = false );

    /// @brief Create scalar and vector tags in the MOAB mesh instance
    /// @param tagScalar Tag reference to the scalar field
    /// @param tagVector Tag reference to the vector field
    /// @param entities Entities on which both the scalar and vector fields are defined
    /// @return Error code if any (else MB_SUCCESS)
    moab::ErrorCode create_sv_tags( moab::Tag& tagScalar, moab::Tag& tagVector, moab::Range& entities ) const;

    /// @brief Evaluate some closed-form Spherical Harmonic functions with an optional multiplier term
    /// @param lon Longitude in lat-lon space
    /// @param lat Latitude  in lat-lon space
    /// @param type Function type
    /// @param multiplier Optional multiplier to scale value (default=1.0)
    /// @return value of the evaluated function
    inline double evaluate_function( double lon, double lat, int type = 1, double multiplier = 1.0 ) const
    {
        switch( type )
        {
            case 1:
                return ( 2.0 + std::pow( sin( 2.0 * lat ), 16.0 ) * cos( 16.0 * lon ) ) * multiplier;
            default:
                return ( 2.0 + cos( lon ) * cos( lon ) * cos( 2.0 * lat ) ) * multiplier;
        }
    }

  private:
    /// @brief Compute the centroids of elements in 2D lat/lon space
    /// @param entities Entities to compute centroids
    /// @return Vector of centroids (as lat/lon)
    std::vector< double > compute_centroids( const moab::Range& entities ) const;

    moab::CpuTimer mTimer;
    double mTimerOps{ 0.0 };
    std::string mOpName;
};

//
// Start of main test program
//
int main( int argc, char** argv )
{
    // Initialize MPI first
    MPI_Init( &argc, &argv );

    {
        // Create our context for this example run
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

        // Read the input file specified by user, in parallel, using appropriate options
        // Supports reading partitioned h5m files and MPAS nc files directly with online Zoltan partitioning
        context.timer_push( "Read input file" );
        {
            // Load the file from disk with given options
            runchk( context.load_file( false ), "MOAB::load_file failed for filename: " << context.input_filename );
        }
        context.timer_pop();
        elapsed_times[0] = context.last_elapsed();

        // Let the actual measurements begin...
        dbgprint( "\n- Starting execution -\n" );

        // We need to set up the ghost layers requested by the user. First correct for thin layers and then
        // call `exchange_ghost_cells` to prepare the mesh for use with halo regions
        context.timer_push( "Setup ghost layers" );
        {
            // Loop over the number of ghost layers needed and ask MOAB for layers 1 at a time
            for( int ighost = 0; ighost < context.ghost_layers; ++ighost )
            {
                // Exchange ghost cells
                int ghost_dimension  = context.dimension;
                int bridge_dimension = context.dimension - 1;
                // Let us now get all ghost layers from adjacent parts
                runchk( context.parallel_communicator->exchange_ghost_cells(
                            ghost_dimension, bridge_dimension, ( ighost + 1 ), 0, true /* store_remote_handles */,
                            true /* wait_all */, &context.fileset ),
                        "Exchange ghost cells failed" );  // true to store remote handles

                // Ensure that all processes understand about multi-shared vertices and entities
                // in case some adjacent parts are only m layers thick (where m < context.ghost_layers)
                if( ighost < context.ghost_layers - 1 )
                    runchk( context.parallel_communicator->correct_thin_ghost_layers(),
                            "Thin layer correction failed" );
            }
        }
        context.timer_pop();
        elapsed_times[1] = context.last_elapsed();

        // Get the 2D MPAS elements and filter it so that we have only owned elements
        Range dimEnts;
        {
            // Get all entities of dimension = dim
            runchk( context.moab_interface->get_entities_by_dimension( context.fileset, context.dimension, dimEnts ),
                    "Getting 2D entities failed" );
            // Get only owned entities! The ghosted/shared entities will get their data when we exchange
            // So let us filter entities based on the status: NOT x NOT_OWNED = OWNED status :-)
            runchk( context.parallel_communicator->filter_pstatus( dimEnts, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                    "Filtering pstatus failed" );

            // Aggregate the total number of elements in the mesh
            auto numEntities     = dimEnts.size();
            int numTotalEntities = 0;
            MPI_Reduce( &numEntities, &numTotalEntities, 1, MPI_INT, MPI_SUM, 0,
                        context.parallel_communicator->proc_config().proc_comm() );

            // We expect the total number of elements to be constant, immaterial of number of processes.
            // If not, we have a bug!
            dbgprint( "Total number of " << context.dimension << "D elements in the mesh = " << numTotalEntities );
        }

        Tag tagScalar = nullptr;
        Tag tagVector = nullptr;
        // Create two tag handles: scalar_variable and vector_variable
        // Set these tags with appropriate closed form functional data
        // based on element centroid information
        runchk( context.create_sv_tags( tagScalar, tagVector, dimEnts ), "Unable to create scalar and vector tags" );

        // let us write out the local mesh before tag_exchange is called
        // we expect to see data only on the owned entities - and ghosted entities should have default values
        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "> Writing to file *before* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_pre.h5m", "H5M", "" ),
                    "Writing to disk failed" );
        }

        // Perform exchange of tag data between neighboring tasks
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

        // let us write out the local mesh after tag_exchange is called
        // we expect to see real data on both owned and ghost entities in halo regions (non-default values)
        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "> Writing to file *after* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_post.h5m", "H5M", "" ),
                    "Writing to disk failed" );
        }

        // Write out the final mesh with the tag data and mesh -- just for verification
        if( context.debug_output )
        {
            dbgprint( "> Writing out the final mesh and data in MOAB h5m format. File = " << context.output_filename );
            string write_options = ( context.num_procs > 1 ? "PARALLEL=WRITE_PART;DEBUG_IO=0;" : "" );
            // Write out to output file to visualize reduction/exchange of tag data
            runchk( context.moab_interface->write_file( context.output_filename.c_str(), "H5M", write_options.c_str() ),
                    "File write failed" );
        }

        // Consolidated timing results: the data is listed as follows
        // [ntasks,  nghosts,  load_mesh(I/O),  exchange_ghost_cells(setup), exchange_tags(scalar),
        // exchange_tags(vector)]
        dbgprint( "\n> Consolidated: [" << context.num_procs << ", " << context.ghost_layers << ", " << elapsed_times[0]
                                        << ", " << elapsed_times[1] << ", " << elapsed_times[2] << ", "
                                        << elapsed_times[3] << "]," );

        // execution finished
        dbgprint( "\n********** ExchangeHalos Example DONE! **********" );
    }
    // Done, cleanup
    MPI_Finalize();

    return 0;
}

/// Implementation details ///

///
RuntimeContext::RuntimeContext( MPI_Comm comm )
    : input_filename( std::string( MESH_DIR ) + std::string( "/io/mpasx1.642.t.2.nc" ) ),
      output_filename( "exchangeHalos_output.h5m" ), scalar_tagname( "scalar_variable" ),
      vector_tagname( "vector_variable" )
{
    // Create the moab instance
    moab_interface = new( std::nothrow ) moab::Core;
    if( NULL == moab_interface ) exit( 1 );

    // Create sets for the mesh and partition.  Then pass these to the load_file functions to populate the mesh.
    runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, fileset ), "Creating root set failed" );
    runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, partnset ), "Creating partition set failed" );

    // Create the parallel communicator object with the partition handle associated with MOAB
    parallel_communicator = moab::ParallelComm::get_pcomm( moab_interface, partnset, &comm );

    proc_id   = parallel_communicator->rank();
    num_procs = parallel_communicator->size();
}

RuntimeContext::~RuntimeContext()
{
    delete parallel_communicator;
    delete moab_interface;
}

void RuntimeContext::ParseCLOptions( int argc, char* argv[] )
{
    ProgOptions opts;
    // Input mesh
    opts.addOpt< std::string >( "input", "Input mesh filename to load in parallel. Default=data/default_mesh_holes.h5m",
                                &input_filename );
    // Output mesh
    opts.addOpt< void >( "debug", "Should we write output file? Default=false", &debug_output );
    opts.addOpt< std::string >( "output",
                                "Output mesh filename for verification (use --debug). Default=exchangeHalos_output.h5m",
                                &output_filename );
    // Dimension of the input mesh
    // Vector tag length
    opts.addOpt< int >( "vtaglength", "Size of vector components per each entity. Default=3", &vector_length );
    // Number of halo (ghost) regions
    opts.addOpt< int >( "nghosts", "Number of ghost layers (halos) to exchange. Default=3", &ghost_layers );
    // Number of times to perform the halo exchange for timing
    opts.addOpt< int >( "nexchanges", "Number of ghost-halo exchange iterations to perform. Default=10",
                        &num_max_exchange );

    opts.parseCommandLine( argc, argv );
}

void RuntimeContext::timer_push( std::string operation )
{
    mTimerOps = mTimer.time_since_birth();
    mOpName   = operation;
}

void RuntimeContext::timer_pop( const int nruns )
{
    double locElapsed = mTimer.time_since_birth() - mTimerOps;
    double avgElapsed = 0;
    double maxElapsed = 0;
    MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, parallel_communicator->comm() );
    MPI_Reduce( &locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, parallel_communicator->comm() );
    if( proc_id == 0 )
    {
        avgElapsed /= num_procs;
        if( nruns > 1 )
            std::cout << "[LOG] Time taken to " << mOpName.c_str() << ", averaged over " << nruns
                      << " runs : max = " << maxElapsed / nruns << ", avg = " << avgElapsed / nruns << "\n";
        else
            std::cout << "[LOG] Time taken to " << mOpName.c_str() << " : max = " << maxElapsed
                      << ", avg = " << avgElapsed << "\n";

        last_counter = maxElapsed / nruns;
    }
    mOpName.clear();
}

double RuntimeContext::last_elapsed() const
{
    return last_counter;
}

moab::ErrorCode RuntimeContext::create_sv_tags( moab::Tag& tagScalar,
                                                moab::Tag& tagVector,
                                                moab::Range& entities ) const
{
    // Get element (centroid) coordinates so that we can evaluate some arbitrary data
    std::vector< double > entCoords = compute_centroids( entities );  // [entities * [lon, lat]]

    if( proc_id == 0 ) std::cout << "> Getting scalar tag handle " << scalar_tagname << "..." << std::endl;
    double defSTagValue = -1.0;
    bool createdTScalar = false;
    // Get or create the scalar tag: default name = "scalar_variable"
    // Type: double, Components: 1, Layout: Dense (all entities potentially), Default: -1.0
    runchk( moab_interface->tag_get_handle( scalar_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, tagScalar,
                                            moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &defSTagValue, &createdTScalar ),
            "Retrieving scalar tag handle failed" );

    // we expect to create a new tag -- fail if Tag already exists since we do not want to overwrite data
    assert( createdTScalar );
    // set the data for scalar tag with an analytical Spherical Harmonic function
    {
        std::vector< double > tagValues( entities.size(), -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = index++ * 2;
            return evaluate_function( entCoords[offset], entCoords[offset + 1] );
        } );
        // Set local scalar tag data for exchange
        runchk( moab_interface->tag_set_data( tagScalar, entities, tagValues.data() ),
                "Setting scalar tag data failed" );
    }

    if( proc_id == 0 ) std::cout << "> Getting vector tag handle " << vector_tagname << "..." << std::endl;
    std::vector< double > defVTagValue( vector_length, -1.0 );
    bool createdTVector = false;
    // Get or create the scalar tag: default name = "vector_variable"
    // Type: double, Components: vector_length, Layout: Dense (all entities potentially), Default: [-1.0,..]
    runchk( moab_interface->tag_get_handle( vector_tagname.c_str(), vector_length, moab::MB_TYPE_DOUBLE, tagVector,
                                            moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, defVTagValue.data(),
                                            &createdTVector ),
            "Retrieving vector tag handle failed" );

    // we expect to create a new tag -- fail if Tag already exists since we do not want to overwrite data
    assert( createdTVector );
    // set the data for vector tag with an analytical Spherical Harmonic function
    // with an optional scaling for each component; just to make it look different :-)
    {
        const int veclength = vector_length;
        std::vector< double > tagValues( entities.size() * veclength, -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = ( index++ / veclength ) * 2;
            return this->evaluate_function( entCoords[offset], entCoords[offset + 1], 2, ( index % veclength + 1.0 ) );
        } );
        // Set local tag data for exchange
        runchk( moab_interface->tag_set_data( tagVector, entities, tagValues.data() ),
                "Setting vector tag data failed" );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode RuntimeContext::load_file( bool load_ghosts )
{
    /// Parallel Read options:
    ///   PARALLEL = type {READ_PART} : Read on all tasks
    ///   PARTITION_METHOD = RCBZOLTAN : Use Zoltan partitioner to compute an online partition and redistribute on the
    ///   fly PARTITION = PARALLEL_PARTITION : Partition as you read based on part information stored in h5m file
    ///   PARALLEL_RESOLVE_SHARED_ENTS : Communicate to all processors to get the shared adjacencies
    ///   consistently in parallel
    ///   PARALLEL_GHOSTS : a.b.c
    ///                   : a = 2 - highest dimension of entities (2D in this case)
    ///                   : b = 1 - dimension of entities to calculate adjacencies (vertex=0, edges=1)
    ///                   : c = 3 - number of ghost layers needed (3 in this case)
    std::string read_options   = "DEBUG_IO=0;";
    std::string::size_type idx = input_filename.rfind( '.' );
    std::string extension      = "";
    if( num_procs > 1 && idx != std::string::npos )
    {
        extension = input_filename.substr( idx + 1 );
        if( !extension.compare( "nc" ) )
            // PARTITION_METHOD= [RCBZOLTAN, TRIVIAL]
            read_options += "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;"
                            "PARALLEL_RESOLVE_SHARED_ENTS;NO_EDGES;NO_MIXED_ELEMENTS;VARIABLE=;";
        else if( !extension.compare( "h5m" ) )
            read_options +=
                "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                "PARALLEL_RESOLVE_SHARED_ENTS;" +
                ( load_ghosts ? "PARALLEL_THIN_GHOST_LAYER;PARALLEL_GHOSTS=2.1." + std::to_string( ghost_layers ) + ";"
                              : "" );
        else
        {
            std::cout << "Error unsupported file type (only h5m and nc) for this example: " << input_filename
                      << std::endl;
            return moab::MB_UNSUPPORTED_OPERATION;
        }
    }

    // Load the file from disk with given read options in parallel and associate all entities to fileset
    return moab_interface->load_file( input_filename.c_str(), &fileset, read_options.c_str() );
}

std::vector< double > RuntimeContext::compute_centroids( const moab::Range& entities ) const
{
    double node[3];
    std::vector< double > eCentroids( entities.size() * 2 );  // [lon, lat]
    size_t offset = 0;
    for( auto entity : entities )
    {
        // Get the element coordinates (centroid) on the real mesh
        runchk_cont( moab_interface->get_coords( &entity, 1, node ), "Getting entity coordinates failed" );

        // scale by magnitude so that element is on unit sphere
        double magnitude = std::sqrt( node[0] * node[0] + node[1] * node[1] + node[2] * node[2] );
        node[0] /= magnitude;
        node[1] /= magnitude;
        node[2] /= magnitude;

        // compute the spherical transformation onto unit sphere
        eCentroids[offset] = atan2( node[1], node[0] );
        if( eCentroids[offset] < 0.0 ) eCentroids[offset] += 2.0 * M_PI;
        eCentroids[offset + 1] = asin( node[2] );

        offset += 2;  // increment the offset
    }
    // return centroid list for elements
    return eCentroids;
}
