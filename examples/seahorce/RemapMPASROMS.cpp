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

#include "mba.hpp"
#include "nanoflann.hpp"

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

moab::ErrorCode shepard_interpolate( int dimension,
                                     std::vector< double >& xyzd,
                                     std::vector< double >& fd,
                                     double power,
                                     std::vector< double >& xyzi,
                                     std::vector< double >& fi );

moab::ErrorCode modified_shepard_interpolate( int dimension,
                                     std::vector< double >& xyzd,
                                     std::vector< double >& fd,
                                     double power,
                                     std::vector< double >& xyzi,
                                     std::vector< double >& fi );

moab::ErrorCode compute_mba( std::vector< double >& xyzd,
                     std::vector< double >& fd,
                     std::vector< double >& xyzi,
                     std::vector< double >& fi );

moab::ErrorCode ComputeFieldProjectionMBA( moab::Interface* mbi,
                                           Mesh& meshOverlap,
                                           std::string varProject,
                                           moab::Range& srcelems,
                                           moab::Range& dstelems,
                                           bool normalize );

template < typename T >
struct PointCloud
{
    using coord_t = T;  //!< The type of each coordinate

    const static int projection    = 2;
    const static int in_dimension  = 3;
    const static int out_dimension = 2;
    const std::vector< T >& xyz;
    std::vector< T > xyz_T;
    const size_t count;

    PointCloud( const std::vector< T >& pxyz ) : xyz( pxyz ), count( pxyz.size() / in_dimension )
    {
        xyz_T.resize( out_dimension * count );

        init();
    }

    void init()
    {
        double cd[3];
        unsigned offset = 0;
        for( size_t i = 0; i < count; i++ )
        {
            cd[0] = xyz[offset];
            cd[1] = xyz[offset + 1];
            cd[2] = xyz[offset + 2];
            IntxUtils::transform_coordinates( cd, projection );
            xyz_T[i * out_dimension]     = cd[0];
            xyz_T[i * out_dimension + 1] = cd[1];
        }
    }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return count;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt( const size_t idx, const size_t dim ) const
    {
        return ( dim ? xyz_T[idx * out_dimension + 1] : xyz_T[idx * out_dimension] );
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template < class BBOX >
    bool kdtree_get_bbox( BBOX& /* bb */ ) const
    {
        return false;
    }
};

template < typename T >
struct PC3D
{
    using coord_t = T;  //!< The type of each coordinate

    const static int dimension  = 3;
    const std::vector< T >& xyz;
    const size_t count;

    PC3D( const std::vector< T >& pxyz ) : xyz( pxyz ), count( pxyz.size() / dimension )
    {
    }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return count;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt( const size_t idx, const size_t dim ) const
    {
        return xyz[idx * dimension + dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template < class BBOX >
    bool kdtree_get_bbox( BBOX& /* bb */ ) const
    {
        return false;
    }
};

//
// Start of main test program
//
int main( int argc, char** argv )
{
    ErrorCode err;
    int ierr, rank, size;
    string mpas_filename, roms_filename, output_filename;
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
    const double radius         = 1.0;
    bool ensureMonotonicity     = false;
    bool computeWeights         = false;
    bool computeShepards        = false;
    bool normalize              = false;
    std::string strMethod       = "";
    const char* varProject      = "bottomDepth";
    const double shepard_power  = 2;
    const bool useTranspose     = false;
    const double mpas_zh[60]    = { 10,      20,      30,      40,      50,      60,      70,      80,      90,
                                    100,     110,     120,     130,     140,     150,     160,     170.197, 180.761,
                                    191.821, 203.499, 215.923, 229.233, 243.584, 259.156, 276.152, 294.815, 315.424,
                                    338.312, 363.875, 392.58,  424.989, 461.767, 503.707, 551.749, 606.997, 670.729,
                                    744.398, 829.607, 928.043, 1041.37, 1171.04, 1318.09, 1482.9,  1664.99, 1863.01,
                                    2074.87, 2298.04, 2529.9,  2768.1,  3010.67, 3256.14, 3503.45, 3751.89, 4001.01,
                                    4250.53, 4500.26, 4750.12, 5000.05, 5250.01, 5499.99 };

    {
        ProgOptions opts;

        // set default values
        if( useTranspose )
        {
            roms_filename = "mpas_grid.h5m";
            mpas_filename = "roms_grid.h5m";
        }
        else
        {
            mpas_filename = "mpas_grid.h5m";
            roms_filename = "roms_grid.h5m";
        }

        output_filename = "output_mpas_roms_map2d.nc";

        opts.addOpt< std::string >( "mpas", "MPAS filename with 2D mesh and 3D dataset", &mpas_filename );
        opts.addOpt< std::string >( "roms", "ROMS filename with 2D mesh", &roms_filename );
        opts.addOpt< std::string >( "out", "Output filename for ROMS 2D mesh and remapped dataset", &output_filename );
        opts.addOpt< std::string >( "method", "Additional computational method arguments (invdist, bilin, intbilin)",
                                    &strMethod );
        opts.addOpt< void >( "shepard", "Use Shepard's inverse-distance weighting to compute projection",
                             &computeShepards );
        opts.addOpt< void >( "mono", "Ensure monotonicity in the weight generation", &ensureMonotonicity );
        opts.addOpt< void >( "normalize", "Re-normalize interpolant to preserve field integral", &normalize );
        opts.addOpt< void >( "weights,w", "Compute and output the weights", &computeWeights );

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

    MPI_Comm comm = MPI_COMM_WORLD;

    ierr = MPI_Comm_rank( comm, &rank );MPICHKERR( ierr, "MPI_Comm_rank failed" );
    ierr = MPI_Comm_size( comm, &size );MPICHKERR( ierr, "MPI_Comm_size failed" );

    // string mpas_read_options = size > 1 ? "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;NO_EDGES;VARIABLE=bottomDepth;"
    //     : "NO_EDGES;VARIABLE=bottomDepth;";
    string mpas_read_options = size > 1 ? "" : "";
    string roms_read_options = size > 1 ? "" : "";
    string write_options     = size > 1 ? "PARALLEL=WRITE_PART" : "";

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

    // construct the remapper
    moab::TempestRemapper remapper( mbi, parallel_communicator );
    remapper.meshValidate     = true;
    remapper.constructEdgeMap = false;
    remapper.initialize();

    EntityHandle& mpasset           = remapper.GetMeshSet( moab::Remapper::SourceMesh );
    EntityHandle& romsset           = remapper.GetMeshSet( moab::Remapper::TargetMesh );
    EntityHandle& mpas_covering_set = remapper.GetMeshSet( moab::Remapper::CoveringMesh );

    // Load the MPAS file from disk with given options
    {
        dbgprint( "Reading MPAS file from disk" );
        err = mbi->load_file( mpas_filename.c_str(), &mpasset, mpas_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for MPAS mesh failed" );

        // Tag depthtag;
        // err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, depthtag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
        // if( err != MB_SUCCESS ) dbgprint( "Error: " << err << "; Failed to get bottomDepth tag handle" );
        // MB_CHK_SET_ERR( err, "MPAS bottomDepth tag failed" );

        // Get all entities in the database
        moab::Range mpas_verts, mpas_elems;
        err = mbi->get_entities_by_dimension( mpasset, 0, mpas_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( mpasset, 2, mpas_elems );MB_CHK_ERR( err );
        dbgprint( "MPAS mesh contains " << mpas_verts.size() << " vertices and " << mpas_elems.size() << " elements");

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, mpas_verts, radius, !useTranspose );MB_CHK_ERR( err );
        err = mbi->write_file( "mpas_modified_2d.h5m", "H5M", write_options.c_str(), &mpasset, 1 );MB_CHK_ERR( err );
        err = remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( err );
    }

    // Load the ROMS file from disk with given options
    {
        dbgprint( "Reading ROMS file from disk" );
        err = mbi->load_file( roms_filename.c_str(), &romsset, roms_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS mesh failed" );

        // Get all entities in the database
        moab::Range roms_verts, roms_elems;
        err = mbi->get_entities_by_dimension( romsset, 0, roms_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( romsset, 2, roms_elems );MB_CHK_ERR( err );
        dbgprint( "ROMS mesh contains " << roms_verts.size() << " vertices and " << roms_elems.size() << " elements");

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, roms_verts, radius, useTranspose );MB_CHK_ERR( err );
        err = mbi->write_file( "roms_modified_2d.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );
        err = remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh );MB_CHK_ERR( err );
    }

    const double epsrel = ReferenceTolerance;
    const double boxeps = 1e-6;

    // first create the covering set
    err = mbi->create_meshset( moab::MESHSET_SET, mpas_covering_set );MB_CHK_SET_ERR( err, "Can't create new set" );
    if (false)
    {
        moab::Intx2MeshOnSphere mbintx( mbi );
        mbintx.set_error_tolerance( epsrel );
        mbintx.set_radius_source_mesh( radius );
        mbintx.set_radius_destination_mesh( radius );
        mbintx.set_box_error( boxeps );

        // moab::Range local_verts;
        // err = mbi->get_entities_by_dimension( romsset, 0, local_verts );MB_CHK_ERR( err );
        // err = mbintx.build_processor_euler_boxes( romsset, local_verts );MB_CHK_ERR( err );
        dbgprint( "Constructing covering set for intersection" );
        err = mbintx.construct_covering_set( mpasset, mpas_covering_set );MB_CHK_ERR( err );
    }
    else
    {
        // err = remapper.ConstructCoveringSet( epsrel, radius, radius, boxeps, false );MB_CHK_ERR( err );
        // mpas_covering_set = remapper.GetMeshSet( moab::Remapper::CoveringMesh );

        // Cull the MPAS set so that we don't have a global mesh

        // construct a kd-tree index:
        using KdTree =
            nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                 PC3D< double >, 3 /* dim */
                                                 >;

        moab::Range mpas_elems;
        err = mbi->get_entities_by_dimension( mpasset, 2, mpas_elems );MB_CHK_ERR( err );
        std::vector< double > mpas_xyz( mpas_elems.size() * 3 );
        err = mbi->get_coords( mpas_elems, mpas_xyz.data() );MB_CHK_ERR( err );

        moab::Range roms_elems;
        err = mbi->get_entities_by_dimension( romsset, 2, roms_elems );MB_CHK_ERR( err );

        double query_pt[3];  // dimension
        PC3D< double > cloud( mpas_xyz );
        KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

        moab::Range culled_elems;
        const size_t num_results =
            static_cast< size_t >( ( 2 * shepard_power + 1 ) * ( 2 * shepard_power + 1 ) - shepard_power );
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        nanoflann::KNNResultSet< double > resultSet( num_results );
        for( size_t i = 0; i < roms_elems.size(); i++ )
        {
            const moab::EntityHandle ehandle = roms_elems[i];
            err               = mbi->get_coords( &ehandle, 1, &query_pt[0] );MB_CHK_ERR( err );

            // Do a KNN search
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            for( size_t j = 0; j < num_results; ++j )
                culled_elems.insert( mpas_elems[srcindx[j]] );
        }
        err = mbi->add_entities( mpas_covering_set, culled_elems );MB_CHK_ERR( err );
        dbgprint( "Culled MPAS mesh contains " << culled_elems.size() << " elements" );

        err = remapper.ConvertMeshToTempest( moab::Remapper::CoveringMesh );MB_CHK_ERR( err );

        err = mbi->write_file( "mpas_covering_2d.h5m", "H5M", write_options.c_str(), &mpas_covering_set, 1 );MB_CHK_ERR( err );
    }

    if( computeShepards )
    {
        // call Shepard's interpolant to compute data
        moab::Range mpas_elems;
        // moab::Range& mpas_elems = remapper.GetMeshEntities( Remapper::SourceMesh );
        err = mbi->get_entities_by_dimension( mpas_covering_set, 2, mpas_elems, true );MB_CHK_ERR( err );

        moab::Range& roms_elems = remapper.GetMeshEntities( Remapper::TargetMesh );

        // compute the actual shepard's interpolation
        // dbgprint( "Computing the Shepard's interpolant now" );
        // err = modified_shepard_interpolate( 3, mpas_xyz, mpas_tdata, shepard_power , roms_xyz, roms_tdata );MB_CHK_ERR( err );

        Mesh meshOverlap;
        if( normalize )
        {
            Mesh& meshInput  = *remapper.GetMesh( moab::Remapper::CoveringMesh );
            Mesh& meshOutput = *remapper.GetMesh( moab::Remapper::TargetMesh );
            meshInput.ConstructEdgeMap();
            meshOutput.ConstructEdgeMap();

            // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
            std::cout << "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes\n";
            // err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );
            bool concaveMeshA = false, concaveMeshB = false, allowNoOverlap = true, verbose = false;
            int err = GenerateOverlapWithMeshes( meshInput, meshOutput, meshOverlap, "" /*outFilename*/, "Netcdf4", "exact",
                                                concaveMeshA, concaveMeshB, allowNoOverlap, verbose );
            if( err )
            {
                MB_CHK_SET_ERR( MB_FAILURE, "TempestRemap: Can't compute the intersection of meshes on the sphere" );
            }
        }

        // Now let us compute the mba hierarchy for each field
        err = ComputeFieldProjectionMBA( mbi, meshOverlap, "bottomDepth", mpas_elems, roms_elems, true /* bool normalize */ );MB_CHK_ERR( err );
        err = ComputeFieldProjectionMBA( mbi, meshOverlap, "salinity", mpas_elems, roms_elems, true /* bool normalize */ );MB_CHK_ERR( err );
        err = ComputeFieldProjectionMBA( mbi, meshOverlap, "temperature", mpas_elems, roms_elems, true /* bool normalize */ );MB_CHK_ERR( err );
    }
    else
    {
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

        // compute the mapping weights
        if( computeWeights )
        {
            dbgprint( "\nSetup computation of weights" );
            // Call to generate the remapping weights with the tempest meshes
            moab::TempestOnlineMap weightMap( &remapper );

            GenerateOfflineMapAlgorithmOptions mapOptions;
            mapOptions.nPin             = 1;
            mapOptions.nPout            = 1;
            mapOptions.fSourceConcave   = false;
            mapOptions.fTargetConcave   = false;
            mapOptions.strMethod        = strMethod; // invdist, bilin
            mapOptions.fMonotone        = ensureMonotonicity;
            mapOptions.fNoCorrectAreas  = false;
            mapOptions.fNoCheck         = true;
            mapOptions.strOutputMapFile = output_filename; // ask TR to write it out
            mapOptions.strOutputFormat  = "Netcdf4";

            dbgprint( "Compute weights with TempestRemap" );
            err = weightMap.GenerateRemappingWeights( "fv",         // std::string strInputType
                                                    "fv",         // std::string strOutputType,
                                                    mapOptions,   // const GenerateOfflineMapAlgorithmOptions& options
                                                    "GLOBAL_ID",  // const std::string& source_tag_name
                                                    "GLOBAL_ID"   // const std::string& target_tag_name
            );MB_CHK_ERR( err );

            // check the generated weights and output information
            {
                const double dNormalTolerance = 1.0E-8;
                const double dStrictTolerance = 1.0E-12;
                weightMap.CheckMap( true, true, ensureMonotonicity, dNormalTolerance, dStrictTolerance );
            }

            {
                // Write the map to disk
                dbgprint( "\nWrite the weights to " << output_filename );

                typedef std::map< std::string, std::string > AttributeMap;
                typedef AttributeMap::value_type AttributePair;

                AttributeMap mapAttributes;

                mapAttributes.insert( AttributePair( "domain_a", mpas_filename ) );
                mapAttributes.insert( AttributePair( "domain_b", roms_filename ) );
                mapAttributes.insert( AttributePair( "grid_file_src", mpas_filename ) );
                mapAttributes.insert( AttributePair( "grid_file_dst", roms_filename ) );
                mapAttributes.insert( AttributePair( "grid_file_ovr", "mesh_intersection.h5m" ) );
                mapAttributes.insert(
                    AttributePair( "concave_src", ( mapOptions.fSourceConcave ) ? ( "true" ) : ( "false" ) ) );
                mapAttributes.insert(
                    AttributePair( "concave_dst", ( mapOptions.fTargetConcave ) ? ( "true" ) : ( "false" ) ) );
                if( mapOptions.strSourceMeta != "" )
                {
                    mapAttributes.insert( AttributePair( "meta_src", mapOptions.strSourceMeta ) );
                }
                if( mapOptions.strTargetMeta != "" )
                {
                    mapAttributes.insert( AttributePair( "meta_dst", mapOptions.strTargetMeta ) );
                }
                mapAttributes.insert( AttributePair( "type_src", "fv" ) );
                mapAttributes.insert( AttributePair( "type_dst", "fv" ) );
                mapAttributes.insert( AttributePair( "np_src", std::to_string( (long long)mapOptions.nPin ) ) );
                mapAttributes.insert( AttributePair( "np_dst", std::to_string( (long long)mapOptions.nPout ) ) );
                mapAttributes.insert( AttributePair( "mono", ( mapOptions.fMonotone ) ? ( "true" ) : ( "false" ) ) );
                mapAttributes.insert( AttributePair( "nobubble", "false" ) );
                mapAttributes.insert(
                    AttributePair( "nocorrectareas", "false" ) );
                mapAttributes.insert( AttributePair( "noconserve", "false" ) );
                mapAttributes.insert( AttributePair( "sparse_constraints", "false" ) );
                mapAttributes.insert( AttributePair( "method", mapOptions.strMethod ) );
                mapAttributes.insert( AttributePair( "version", "RemapMPASROMS v0.1" ) );

                weightMap.Write( mapOptions.strOutputMapFile, mapAttributes, NcFile::Netcdf4Classic );

                // // Write the map file to disk in parallel using either HDF5 or SCRIP interface
                // err = weightMap.WriteParallelMap( output_filename.c_str() );MB_CHK_ERR( err );
            }

            // Now apply the map to compute bottom depth in ROMS mesh
            if( true )
            {
                // MPAS variable to project = "bottomDepth"
                moab::Tag dtag;
                err = mbi->tag_get_handle( varProject, 1, moab::MB_TYPE_DOUBLE, dtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

                // Now let us apply the weights onto the vector and project onto target mesh
                err = weightMap.ApplyWeights( dtag, dtag, useTranspose );MB_CHK_ERR( err );
            }
        }

        remapper.clear();
    }

    if (useTranspose)
    {
        err = mbi->write_file( "roms_2d_projected.h5m", "H5M", write_options.c_str(), &mpasset, 1 );MB_CHK_ERR( err );
    }
    else
    {
        err = mbi->write_file( "roms_2d_projected.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );
    }

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
            spherical_to_cart( posi[1], posi[0], R, posf );
            // dbgprint( nd << " lat=" << posi[1] << ", lon=" << posi[0] << "; Cartesian = [" << posf[0] << ", " << posf[1] << ", " << posf[2] << "]" );
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

moab::ErrorCode shepard_interpolate( int dimension,
                             std::vector< double >& xyzd,
                             std::vector< double >& fd,
                             double power,
                             std::vector< double >& xyzi,
                             std::vector< double >& fi )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    size_t nd  = xyzd.size() / dimension;
    size_t ni  = xyzi.size() / dimension;
    double ind = 1.0 / static_cast< double >( nd );

    for( size_t i = 0; i < ni; i++ )
    {
        std::vector< double > w( nd, 0.0 );
        int z;
        const int ioffset = i * dimension;
        if( power < 1 )
        {
            for( size_t j = 0; j < nd; j++ )
                w[j] = ind;
        }
        else
        {
            z = -1;
            for( size_t j = 0; j < nd; j++ )
            {
                double t = 0.0;
                const int joffset = j * dimension;
                for( int i2 = 0; i2 < dimension; i2++ )
                {
                    t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                }
                w[j] = std::sqrt( t );
                if( w[j] < 1e-12 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < nd; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
            }
            else
            {
                double s = 0.0;
                for( size_t j = 0; j < nd; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], power );
                    s += w[j];
                }

                for( size_t j = 0; j < nd; j++ )
                    w[j] /= s;
            }
        }

        fi[i] = 0.0;
        for( size_t k = 0; k < nd; k++ )
        {
            fi[i] += w[k] * fd[k];
        }
        w.clear();
    }

    return moab::MB_SUCCESS;
}



moab::ErrorCode compute_mba( std::vector< double >& xyzd,
                             std::vector< double >& fd,
                             std::vector< double >& xyzi,
                             std::vector< double >& fi )
{
    const size_t nd = fd.size();
    const size_t ni = fi.size();

    // Bounding box containing the data points.
    mba::point< 3 > lo = { -1, -1, -1 };
    mba::point< 3 > hi = { 1, 1, 1 };

    // Initial grid size.
    const size_t init_grid_size = static_cast< size_t >( std::max( 10.0, std::sqrt( nd ) / 8 ) );
    // mba::index< 3 > grid        = { init_grid_size, init_grid_size, 2 };
    mba::index< 3 > grid        = { 80, 80, 2 };

    std::vector< mba::point< 3 > > coords(nd);
    size_t offset = 0;
    for( size_t k = 0; k < nd; k++, offset += 3 )
        coords[k] = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };

    // Algorithm setup.
    mba::MBA< 3 > interp( lo, hi, grid, coords, fd, 5 /*levels*/, 1e-14 /*tolerance*/, 0.5 /*min_fill*/ );
    // mba::linear_approximation< 3 > interp( coords.begin(), coords.end(), fd.begin() );

    // Get interpolated value at arbitrary location.
    offset = 0;
    for( size_t k = 0; k < ni; k++, offset+=3 )
        fi[k] = interp( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );

    return moab::MB_SUCCESS;
}

moab::ErrorCode modified_shepard_interpolate( int dimension,
                                              std::vector< double >& xyzd,
                                              std::vector< double >& fd,
                                              double power,
                                              std::vector< double >& xyzi,
                                              std::vector< double >& fi )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7302837/
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    power = 2;
    // First call the regular inverse-distance weighting method
    shepard_interpolate( dimension, xyzd, fd, power, xyzi, fi );

    power = 4;

    size_t nd  = xyzd.size() / dimension;
    size_t ni  = xyzi.size() / dimension;
    std::vector< double > w( nd, 0.0 );

    double fisecnum = 0.0;
    for( size_t k = 0; k < nd; k++ )
        fisecnum += fd[k];

    assert(power >= 1.0);

    for( size_t i = 0; i < ni; i++ )
    {
        const int ioffset = i * dimension;
        int z;
        double s = 0.0, is = 0.0;
        {
            z = -1;
            for( size_t j = 0; j < nd; j++ )
            {
                double t          = 0.0;
                const int joffset = j * dimension;
                for( int i2 = 0; i2 < dimension; i2++ )
                {
                    t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                }
                w[j] = std::sqrt( t );
                if( w[j] < 1e-6 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < nd; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
                s = is = 1.0;

            }
            else
            {
                for( size_t j = 0; j < nd; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], power );
                    s += w[j];
                    is += 1.0/w[j];
                }

                for( size_t j = 0; j < nd; j++ )
                    w[j] /= s;
            }
        }

        fi[i] += nd * ( fisecnum - nd * fi[i] ) / ( nd * nd - s * is );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode modified_shepard_interpolate2( int dimension,
                                              std::vector< double >& xyzd,
                                              std::vector< double >& fd,
                                              double power,
                                              std::vector< double >& xyzi,
                                              std::vector< double >& fi )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7302837/
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    // size_t nd  = xyzd.size() / dimension;
    size_t ni  = xyzi.size() / dimension;

    // construct a kd-tree index:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PointCloud< double > >,
                                                        PointCloud< double >, 3 /* dim */
                                                        >;

    PointCloud< double > cloud( xyzd );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    double query_pt[3];  // dimension

    const size_t num_results = static_cast< size_t >( ( 2 * power + 1 ) * ( 2 * power + 1 ) - power );
    double ind               = ( 1.0 / num_results );
    std::vector< size_t > srcindx( num_results );
    std::vector< double > srcdist( num_results );
    nanoflann::KNNResultSet< double > resultSet( num_results );
    for( size_t i = 0; i < ni; i++ )
    {
        const int ioffset = i * dimension;
        {
            for( int dd = 0; dd < dimension; dd++ )
                query_pt[dd] = xyzi[ioffset + dd];
            // IntxUtils::transform_coordinates( query_pt, PointCloud< double >::projection );

            // Do a KNN search
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            // tree.knnSearch( &query_pt[0], num_results, &srcindx[0], &srcdist[0] );

            // std::cout << "knnSearch(nn=" << num_results << "): \n";
            for( size_t ll = 0; ll < num_results; ++ll )
            {
                double dist = 0.0;
                for( int dd = 0; dd < dimension; dd++ )
                    dist += ( xyzi[ioffset + dd] - xyzd[srcindx[ll] * dimension + dd] ) *
                            ( xyzi[ioffset + dd] - xyzd[srcindx[ll] * dimension + dd] );
                dist        = std::sqrt( dist );
                std::cout << i << "\tret_index=" << srcindx[ll] << " out_dist_sqr=" << srcdist[ll] << ", " << dist << std::endl;
                srcdist[ll] = dist;
            }
        }

        std::vector< double > w( num_results, 0.0 );
        int z;
        if( power < 1 )
        {
            for( size_t j = 0; j < num_results; j++ )
                w[j] = ind;
        }
        else
        {
            z = -1;
            for( size_t j = 0; j < num_results; j++ )
            {
                // double t          = 0.0;
                // const int joffset = j * dimension;
                // for( int i2 = 0; i2 < dimension; i2++ )
                // {
                //     t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                // }
                w[j] = srcdist[j];// std::sqrt( t );
                if( w[j] < 1e-12 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < num_results; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
            }
            else
            {
                double s = 0.0;
                for( size_t j = 0; j < num_results; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], power );
                    s += w[j];
                }

                for( size_t j = 0; j < num_results; j++ )
                    w[j] /= s;
            }
        }

        fi[i] = 0.0;
        for( size_t k = 0; k < num_results; k++ )
        {
            fi[i] += w[k] * fd[k];
        }
        w.clear();
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ComputeFieldProjectionMBA( moab::Interface* mbi,
                                           Mesh& meshOverlap,
                                           std::string varProject,
                                           moab::Range& srcelems,
                                           moab::Range& dstelems,
                                           bool normalize )
{
    moab::ErrorCode err;
    moab::Tag dtag;
    err = mbi->tag_get_handle( varProject.c_str(), 1, moab::MB_TYPE_DOUBLE, dtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

    // get the source data from tag
    std::vector< double > src_tdata( srcelems.size() ), dst_tdata( dstelems.size() );
    err = mbi->tag_get_data( dtag, srcelems, src_tdata.data() );MB_CHK_ERR( err );

    // get the coordinates of the elements
    std::vector< double > src_xyz( srcelems.size() * 3 ), dst_xyz( dstelems.size() * 3 );
    err = mbi->get_coords( srcelems, src_xyz.data() );MB_CHK_ERR( err );
    err = mbi->get_coords( dstelems, dst_xyz.data() );MB_CHK_ERR( err );

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0, normFactor = 1.0;
    if( normalize )
    {
        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += src_tdata[iSourceFace] * meshOverlap.vecFaceArea[j];
        }
    }

    std::cout << "Computing the MBA interpolant now for field " + varProject + "\n";
    err = compute_mba( src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );

    if( normalize )
    {
        // Loop through all overlap-target faces and compute integral
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iTargetFace = meshOverlap.vecTargetFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iTargetFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralOut += dst_tdata[iTargetFace] * meshOverlap.vecFaceArea[j];
        }

        normFactor = dTotalFieldIntegralIn / dTotalFieldIntegralOut;
        for( size_t ind = 0; ind < dst_tdata.size(); ind++ )
            dst_tdata[ind] *= normFactor;
    }

    // now set the data on ROMS instance of MOAB tag
    std::cout << "Setting tag data to destination mesh\n";
    err = mbi->tag_set_data( dtag, dstelems, dst_tdata.data() );MB_CHK_ERR( err );

    return moab::MB_SUCCESS;
}
