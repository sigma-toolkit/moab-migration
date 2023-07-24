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

#define dbgprint( MSG )                  \
    do                                   \
    {                                    \
        if( !rank ) cerr << MSG << endl; \
    } while( false )

#define dbgprintall( MSG )                           \
    do                                               \
    {                                                \
        cerr << "[" << rank << "]: " << MSG << endl; \
    } while( false )

struct RuntimeContext
{
  public:
    int dimension;
    std::string input_filename;
    int ghost_layers;
    std::string scalar_tagname;
    std::string vector_tagname;
    int vector_length;

    RuntimeContext()
        : dimension( 2 ), input_filename( string( MESH_DIR ) + string( "/io/mpasx1.642.t.2.nc" ) ), ghost_layers( 1 ),
          scalar_tagname( "h_s" ), vector_tagname( "ke" ), vector_length( 1 )
    {
    }
};

// Function to parse input parameters
ErrorCode parse_options( int argc, char** argv, RuntimeContext& context )
{
    ProgOptions opts;
    // Input mesh
    opts.addOpt< std::string >( "input", "Input mesh filename to load in parallel", &context.input_filename );
    // Dimension of the input mesh
    opts.addOpt< int >( "dimension", "Input mesh dimension (default = 2)", &context.dimension );
    // Tag names
    opts.addOpt< std::string >( "stag", "Scalar tag name to exchange with neighboring tasks", &context.scalar_tagname );
    opts.addOpt< std::string >( "vtag", "Vector tag name to exchange with neighboring tasks", &context.vector_tagname );
    opts.addOpt< int >( "vtaglength", "Size of vector components per each entity", &context.vector_length );
    // Ghost layers
    opts.addOpt< int >( "nghosts", "Number of ghost layers (halos) to exchange", &context.ghost_layers );

    opts.parseCommandLine( argc, argv );

    return MB_SUCCESS;
}

//
// Start of main test program
//
int main( int argc, char** argv )
{
    ErrorCode err;
    RuntimeContext context;

    // Initialize MPI first
    MPI_Init( &argc, &argv );

    // Create the moab instance
    Interface* mbi = new( std::nothrow ) Core;
    if( NULL == mbi ) exit( 1 );

    // Create sets for the mesh and partition.  Then pass these to the load_file functions to populate the mesh.
    EntityHandle fileset, partnset;
    err = mbi->create_meshset( MESHSET_SET, fileset );MB_CHK_SET_ERR( err, "Creating root set failed" );
    err = mbi->create_meshset( MESHSET_SET, partnset );MB_CHK_SET_ERR( err, "Creating partition set failed" );

    // Create the parallel communicator object with the partition handle associated with MOAB
    MPI_Comm comm                       = MPI_COMM_WORLD;
    ParallelComm* parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );

    const int rank = parallel_communicator->rank();
    const int size = parallel_communicator->size();

    dbgprint( "********** Exchange halos example **********\n" );

    // Get the input options
    err = parse_options( argc, argv, context );MB_CHK_SET_ERR( err, "Parsing command-line options failed" );

    /////////////////////////////////////////////////////////////////////////
    // Print out the input parameters in use
    dbgprint( " -- Input Parameters -- " );
    dbgprint( "    Input mesh        = " << context.input_filename );
    dbgprint( "    Ghost Layers      = " << context.ghost_layers );
    dbgprint( "    Scalar Tag name   = " << context.scalar_tagname );
    dbgprint( "    Vector Tag name   = " << context.vector_tagname );
    dbgprint( "    Vector Tag length = " << context.vector_length << endl );
    /////////////////////////////////////////////////////////////////////////

    /// Parallel Read options:
    ///   PARALLEL = type {READ_PART}
    ///   PARTITION = PARALLEL_PARTITION : Partition as you read
    ///   PARALLEL_RESOLVE_SHARED_ENTS : Communicate to all processors to get the shared adjacencies
    ///   consistently in parallel PARALLEL_GHOSTS : a.b.c
    ///                   : a = 3 - highest dimension of entities
    ///                   : b = 0 -
    ///                   : c = 1 - number of layers
    ///   PARALLEL_COMM = index
    // string read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;"
    //                       "PARTITION_DISTRIBUTE;PARALLEL_COMM=0";
    // string read_options = "PARALLEL=READ_PART;PARTITION=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;"
    //                       "PARTITION_DISTRIBUTE;PARALLEL_COMM=0";
    // string read_options = ( size > 1 ? ";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;DEBUG_IO=0;NO_EDGES;" : "" );
    string read_options = ( size > 1 ? ";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;"
                                       "NO_EDGES;NO_MIXED_ELEMENTS;DEBUG_IO=0;"
                                     : "" );
    // Load the file from disk with given options
    err = mbi->load_file( context.input_filename.c_str(), &fileset, read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file failed" );

    dbgprint( "- Writing to file " );
    err = mbi->write_file( "exchangeHalos_output_tmp.h5m", "H5M", "PARALLEL=WRITE_PART;DEBUG_IO=0;", &fileset, 1 );MB_CHK_ERR( err );
    dbgprint( "- " );

    // Ensure that all processes understand about multi-shared vertices and entities
    err = parallel_communicator->correct_thin_ghost_layers();MB_CHK_ERR( err );

    // Exchange ghost cells
    int ghost_dim = 2, bridge_dimension = 0, additional_entities = ghost_dim;
    // Let us get one layer at a time to avoid issues with thin partitions
    for( auto igh = 0; igh < context.ghost_layers; ++igh )
        err = parallel_communicator->exchange_ghost_cells( context.dimension, bridge_dimension, 1,
                                                           additional_entities, true /* store_remote_handles */,
                                                           true /* wait_all */ );MB_CHK_ERR( err );  // true to store remote handles

    // Ensure to augment the ghost cells with essential tag data
    // err = parallel_communicator->augment_default_sets_with_ghosts( fileset );MB_CHK_ERR( err );

    Range dimEnts;
    // Get all entities of dimension = dim
    err = mbi->get_entities_by_dimension( fileset, context.dimension, dimEnts );MB_CHK_ERR( err );
    err = parallel_communicator->filter_pstatus( dimEnts, PSTATUS_NOT_OWNED, PSTATUS_NOT );MB_CHK_ERR( err );

    // Aggregate the total number of elements in the mesh
    auto numEntities = dimEnts.size();
    int numTotalEntities;
    MPI_Reduce( &numEntities, &numTotalEntities, 1, MPI_INT, MPI_SUM, 0,
                parallel_communicator->proc_config().proc_comm() );
    dbgprint( "Total number of " << context.dimension << "D elements in the mesh = " << numTotalEntities );

    // Get element (centroid) coordinates so that we can evaluate some arbitrary data
    std::vector< double > entCoords( dimEnts.size() * 2 );  // [lon, lat]
    for( size_t ients = 0, offset = 0; ients < dimEnts.size(); ++ients, offset += 2 )
    {
        const EntityHandle entity = dimEnts[ients];
        double node[3];
        err = mbi->get_coords( &entity, 1, node );MB_CHK_ERR( err );

        // scale by magnitude so that mesh is on unit sphere
        double magnitude = std::sqrt( node[0] * node[0] + node[1] * node[1] + node[2] * node[2] );
        node[0] /= magnitude;
        node[1] /= magnitude;
        node[2] /= magnitude;

        // compute the spherical transformation onto unit sphere
        entCoords[offset] = atan2( node[1], node[0] );
        if( entCoords[offset] < 0.0 ) entCoords[offset] += 2.0 * M_PI;
        entCoords[offset + 1] = asin( node[2] );
    }

    // Create two tag handles: Exchange and Reduction operations
    Tag tagScalar, tagVector;
    bool createdTScalar, createdTVector;
    {
        dbgprint( "> Getting scalar tag handle " << context.scalar_tagname << "..." );
        double defSTagValue = -1.0;
        // Create the exchange tag: default name = USERTAG_EXC
        err = mbi->tag_get_handle( context.scalar_tagname.c_str(), 1, MB_TYPE_DOUBLE, tagScalar,
                                   MB_TAG_CREAT | MB_TAG_DENSE, &defSTagValue, &createdTScalar );MB_CHK_SET_ERR( err, "Retrieving scalar tag handle failed" );

        if( createdTScalar )
        {
            std::vector< double > tagValues( dimEnts.size(), -1.0 );
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
            err = mbi->tag_set_data( tagScalar, dimEnts, tagValues.data() );MB_CHK_SET_ERR( err, "Setting scalar tag data failed" );
        }

        dbgprint( "> Getting vector tag handle " << context.vector_tagname << "..." );
        std::vector< double > defVTagValue( context.vector_length, -1.0 );
        // Create the exchange tag: default name = USERTAG_RED
        err = mbi->tag_get_handle( context.vector_tagname.c_str(), context.vector_length, MB_TYPE_DOUBLE, tagVector,
                                   MB_TAG_CREAT | MB_TAG_DENSE, defVTagValue.data(), &createdTVector );MB_CHK_SET_ERR( err, "Retrieving vector tag handle failed" );

        if( createdTVector )
        {
            const int veclength = context.vector_length;
            std::vector< double > tagValues( dimEnts.size() * veclength, -1.0 );
            std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
                static int index = 0;
                const int offset = ( index / veclength ) * 2;

                double value =
                    ( 2.0 + cos( entCoords[offset] ) * cos( entCoords[offset] ) * cos( 2.0 * entCoords[offset + 1] ) ) *
                    ( index % veclength + 1.0 );  // assign some scalar multiple value for different vector components

                // if (index%veclength == 0) printf("Veclength = %d, offset = %d, value = %f\n", veclength, offset, value);
                index++;
                return value;
            } );
            // Set local tag data for exchange
            err = mbi->tag_set_data( tagVector, dimEnts, tagValues.data() );MB_CHK_SET_ERR( err, "Setting vector tag data failed" );
        }
    }

    // Perform exchange tag data
    dbgprint( "> Exchanging tags between processors " );
    // if( false )
    {
        // Exchange tags between processors
        err = parallel_communicator->exchange_tags( tagScalar, dimEnts );MB_CHK_SET_ERR( err, "Exchanging scalar tag between processors failed" );
        err = parallel_communicator->exchange_tags( tagVector, dimEnts );MB_CHK_SET_ERR( err, "Exchanging vector tag between processors failed" );
    }

    dbgprint( "> Writing out the final mesh and data in MOAB h5m format. File = exchangeHalos_output.h5m." );
    string write_options = ( size > 1 ? "PARALLEL=WRITE_PART;DEBUG_IO=0;" : "" );
    // string write_options = "PARALLEL=WRITE_PART;DEBUG_IO=2;";
    // Write out to output file to visualize reduction/exchange of tag data
    // err = mbi->write_file( "exchangeHalos_output.h5m", "H5M", "PARALLEL=WRITE_PART" );MB_CHK_ERR( err );
    // err = mbi->write_file( "exchangeHalos_output.h5m", "H5M", write_options.c_str() );MB_CHK_ERR( err );
    err = mbi->write_file( "exchangeHalos_output.h5m", "H5M", write_options.c_str(), &fileset, 1 );MB_CHK_ERR( err );

    // Done, cleanup
    delete parallel_communicator;
    delete mbi;

    dbgprint( "\n********** ExchangeHalos Example DONE! **********" );

    MPI_Finalize();
    return 0;
}
