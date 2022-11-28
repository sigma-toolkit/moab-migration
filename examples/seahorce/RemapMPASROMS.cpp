/** @example ReduceExchangeTags.cpp
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
 * <b>To run:</b> \n mpiexec -n 2 ./ReduceExchangeTags <mesh_file> <tag_name> <tag_value> \n
 * <b>Example:</b> \n mpiexec -n 2 ./ReduceExchangeTags ../MeshFiles/unittest/64bricks_1khex.h5m
 * USERTAG 100 \n
 *
 */

#include <iostream>
#include <string>
#include <sstream>
#include "moab/Core.hpp"

#ifndef MOAB_HAVE_MPI
#error " compile with MPI and HDF5 for this example to work \n";
#endif

#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"

// Remapping related includes
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"

using namespace moab;
using namespace std;

// Error routines for use with MPI API
#define MPICHKERR( CODE, MSG )       \
    do                               \
    {                                \
        if( 0 != ( CODE ) )          \
        {                            \
            cerr << ( MSG ) << endl; \
            MPI_Finalize();          \
        }                            \
    } while( false )

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

ErrorCode ScaleCoords( Interface* mb, Range& nodes, double R, bool is_cartesian = true );

//
// Start of main test program
//
int main( int argc, char** argv )
{
    ErrorCode err;
    int ierr, rank;
    string mpas_filename, roms_filename, output_filename;
    MPI_Comm comm = MPI_COMM_WORLD;
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
    //                       "PARTITION_DISTRIBUTE";  // ;PARALLEL_GHOSTS=3.0.1
    string read_options  = "";
    string write_options = "";  // "PARALLEL=WRITE_PART"
    const double radius  = 1.0;

    {
        ProgOptions opts;

        // set default values
        mpas_filename   = "mpas_grid.nc";
        roms_filename   = "roms_grid.h5m";
        output_filename = "output_mpas_roms_map2d.nc";

        opts.addOpt< std::string >( "mpas", "MPAS filename with 2D mesh and 3D dataset", &mpas_filename );
        opts.addOpt< std::string >( "roms", "ROMS filename with 2D mesh", &roms_filename );
        opts.addOpt< std::string >( "out", "Output filename for ROMS 2D mesh and remapped dataset", &output_filename );

        opts.parseCommandLine( argc, argv );
    }

    // Print usage if not enough arguments
    if( argc < 1 )
    {
        cerr << "Usage: ";
        cerr << argv[0] << " --mpas file_name --roms file_name --out file_name" << endl;
        ierr = MPI_Finalize();
        MPICHKERR( ierr, "MPI_Finalize failed; Aborting" );

        return 1;
    }

    // Initialize MPI first
    ierr = MPI_Init( &argc, &argv );MPICHKERR( ierr, "MPI_Init failed" );

    ierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );MPICHKERR( ierr, "MPI_Comm_rank failed" );

    dbgprint( "********** Remap MPAS-to-ROMS **********\n" );

    // Create the moab instance
    Interface* mbi = new( std::nothrow ) Core;
    if( NULL == mbi ) return 1;

    // Print out the input parameters
    dbgprint( " Input Parameters - " );
    dbgprint( "   Filenames:: " );
    dbgprint( "       MPAS: " << mpas_filename );
    dbgprint( "       ROMS: " << roms_filename );
    dbgprint( "     Output: " << output_filename << endl );

    // Create root sets for each mesh.  Then pass these
    // to the load_file functions to be populated.
    // EntityHandle mpasset, romsset;
    // err = mbi->create_meshset( MESHSET_SET, mpasset );MB_CHK_SET_ERR( err, "Creating root set failed" );
    // err = mbi->create_meshset( MESHSET_SET, romsset );MB_CHK_SET_ERR( err, "Creating root set failed" );

    EntityHandle partnset;
    err = mbi->create_meshset( MESHSET_SET, partnset );MB_CHK_SET_ERR( err, "Creating partition set failed" );
    // Create the parallel communicator object with the partition handle associated with MOAB
    ParallelComm* parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );

#ifdef MOAB_HAVE_MPI
    moab::TempestRemapper remapper( mbi, parallel_communicator );
#else
    moab::TempestRemapper remapper( mbi );
#endif
    remapper.meshValidate     = true;
    remapper.constructEdgeMap = false;
    remapper.initialize();

    EntityHandle& mpasset = remapper.GetMeshSet( moab::Remapper::SourceMesh );
    EntityHandle& romsset = remapper.GetMeshSet( moab::Remapper::TargetMesh );

    // Load the MPAS file from disk with given options
    {
        dbgprint( "Reading MPAS file from disk" );
        err = mbi->load_file( mpas_filename.c_str(), &mpasset, read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for MPAS mesh failed" );

        // Get all entities in the database
        moab::Range mpas_verts, mpas_elems;
        err = mbi->get_entities_by_dimension( mpasset, 0, mpas_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( mpasset, 2, mpas_elems );MB_CHK_ERR( err );
        dbgprint( "MPAS mesh contains " << mpas_verts.size() << " vertices and " << mpas_elems.size() << " elements");

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, mpas_verts, radius, true );MB_CHK_ERR( err );
        err = mbi->write_file( "mpas_modified_2d.h5m", "H5M", write_options.c_str(), &mpasset, 1 );MB_CHK_ERR( err );
        err = remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( err );
    }

    // Load the ROMS file from disk with given options
    {
        dbgprint( "Reading ROMS file from disk" );
        err = mbi->load_file( roms_filename.c_str(), &romsset, read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS mesh failed" );

        // Get all entities in the database
        moab::Range roms_verts, roms_elems;
        err = mbi->get_entities_by_dimension( romsset, 0, roms_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( romsset, 2, roms_elems );MB_CHK_ERR( err );
        dbgprint( "ROMS mesh contains " << roms_verts.size() << " vertices and " << roms_elems.size() << " elements");

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, roms_verts, radius, false );MB_CHK_ERR( err );
        err = mbi->write_file( "roms_modified_2d.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );
        err = remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh );MB_CHK_ERR( err );
    }

    // exit(1);

    dbgprint( "Constructing covering set for intersection" );
    const double epsrel = ReferenceTolerance;
    const double boxeps = 1e-6;
    err = remapper.ConstructCoveringSet( epsrel, 1.0, 1.0, boxeps, false );MB_CHK_ERR( err );

    // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
    dbgprint( "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes" );
    err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );

    EntityHandle& intersection_set = remapper.GetMeshSet( moab::Remapper::OverlapMesh );

    // print some diagnostic checks to see if the overlap grid resolved the input meshes correctly
    {
        moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::GaussQuadrature );

        double local_areas[3];  // Array for Initial area, and through Method 1 and Method 2
        // local_areas[0] = area_on_sphere_lHuiller ( mbi, runCtx->meshsets[1], radius );
        local_areas[0] = areaAdaptor.area_on_sphere( mbi, mpasset, radius );
        local_areas[1] = areaAdaptor.area_on_sphere( mbi, romsset, radius );
        local_areas[2] = areaAdaptor.area_on_sphere( mbi, intersection_set, radius );

        dbgprint( "initial area: source mesh = " << local_areas[0] << ", target mesh = " << local_areas[1]
                                                 << ", overlap mesh = " << local_areas[2] );
        dbgprint( "relative error w.r.t source = " << fabs( local_areas[0] - local_areas[2] ) / local_areas[0]
                                                   << ", and target = "
                                                   << fabs( local_areas[1] - local_areas[2] ) / local_areas[1] );
    }

    // Write out to output file to visualize reduction/exchange of tag data
    dbgprint( "Writing intersection mesh... " );
    err = mbi->write_file( "mesh_intersection.h5m", "H5M", write_options.c_str(), &intersection_set, 1 );MB_CHK_ERR( err );

    // Create two tag handles: Exchange and Reduction operations
    // dbgprint( "-Computing weights " << "..." );
    // Tag tagReduce, tagExchange;
    // {
    //     stringstream sstr;
    //     // Create the exchange tag: default name = USERTAG_EXC
    //     sstr << tagName << "_EXC";
    //     err = mbi->tag_get_handle( sstr.str().c_str(), 1, MB_TYPE_INTEGER, tagExchange, MB_TAG_CREAT | MB_TAG_DENSE,
    //                                &tagValue );MB_CHK_SET_ERR( err, "Retrieving tag handles failed" );

    //     // Create the exchange tag: default name = USERTAG_RED
    //     sstr.str( "" );
    //     sstr << tagName << "_RED";
    //     err = mbi->tag_get_handle( sstr.str().c_str(), 1, MB_TYPE_DOUBLE, tagReduce, MB_TAG_CREAT | MB_TAG_DENSE,
    //                                &tagValue );MB_CHK_SET_ERR( err, "Retrieving tag handles failed" );
    // }

    // // Perform exchange tag data
    // dbgprint( "-Exchanging tags between processors " );
    // {
    //     Range partEnts, dimEnts;
    //     for( int dim = 0; dim <= 3; dim++ )
    //     {
    //         // Get all entities of dimension = dim
    //         err = mbi->get_entities_by_dimension( rootset, dim, dimEnts, false );MB_CHK_ERR( err );

    //         vector< int > tagValues( dimEnts.size(), static_cast< int >( tagValue ) * ( rank + 1 ) * ( dim + 1 ) );
    //         // Set local tag data for exchange
    //         err = mbi->tag_set_data( tagExchange, dimEnts, &tagValues[0] );MB_CHK_SET_ERR( err, "Setting local tag data failed during exchange phase" );
    //         // Merge entities into parent set
    //         partEnts.merge( dimEnts );
    //     }

    //     // Exchange tags between processors
    //     err = parallel_communicator->exchange_tags( tagExchange, partEnts );MB_CHK_SET_ERR( err, "Exchanging tags between processors failed" );
    // }

    // Write out to output file to visualize reduction/exchange of tag data
    // err = mbi->write_file( output_filename.c_str(), "H5M", write_options.c_str() );MB_CHK_ERR( err );

    // Done, cleanup
    delete mbi;

    dbgprint( "\n********** Remap MPAS-to-ROMS DONE! **********" );

    MPI_Finalize();
    return 0;
}

static void spherical_to_cart( double lat, double lon, double R, double res[3] )
{
    lat *= 3.14159265358979323846 / 180;
    lon *= 3.14159265358979323846 / 180;
    res[0] = R * cos( lat ) * cos( lon );  // x coordinate
    res[1] = R * cos( lat ) * sin( lon );  // y
    res[2] = R * sin( lat );               // z
    // res[0] = R * sin( lat ) * cos( lon );  // x coordinate
    // res[1] = R * sin( lat ) * sin( lon );  // y
    // res[2] = R * cos( lat );               // z
}

ErrorCode ScaleCoords( Interface* mb, Range& nodes, double R, bool is_cartesian )
{
    ErrorCode rval;
    int rank = 0;
    double posi[3], posf[3];

    // one by one, get the node and project it on the sphere, with a radius given
    // the center of the sphere is at 0,0,0
    for( Range::iterator nit = nodes.begin(); nit != nodes.end(); ++nit )
    {
        EntityHandle nd = *nit;

        if( !is_cartesian )
        {
            rval = mb->get_coords( &nd, 1, posi );MB_CHK_ERR( rval );
            spherical_to_cart( posi[0], posi[1], R, posf );
            // dbgprint( nd << " lat=" << posi[0] << ", lon=" << posi[1] << "; Cartesian = [" << posf[0] << ", " << posf[1] << ", " << posf[2] << "]" );
        }
        else
        {
            rval = mb->get_coords( &nd, 1, posf );MB_CHK_ERR( rval );
        }

        double len = std::sqrt( posf[0] * posf[0] + posf[1] * posf[1] + posf[2] * posf[2] );
        if( len < 1e-12 )
        {
            dbgprint( nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2] << ": Failed with length == 0.");
            return MB_FAILURE;
        }

        // rescale to radius
        posf[0] *= R / len;
        posf[1] *= R / len;
        posf[2] *= R / len;
        rval    = mb->set_coords( &nd, 1, posf );MB_CHK_ERR( rval );
    }
    return MB_SUCCESS;
}
