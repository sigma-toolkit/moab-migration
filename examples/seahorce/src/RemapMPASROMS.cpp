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
#include "moab/MOABConfig.h"
#undef MOAB_HAVE_MPI

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#include "moab/Remapping/MBA.hpp"
#include "moab/nanoflann.hpp"

#ifdef MOAB_HAVE_MPI
// #error " compile with MPI and HDF5 for this example to work \n";
#include "moab/ParallelComm.hpp"
#endif

#include "moab/ProgOptions.hpp"

// Remapping related includes
#include "FiniteVolumeTools.h"
#include "moab/IntxMesh/IntxUtils.hpp"
// #include "moab/Remapping/TempestRemapper.hpp"

// #include "moab/Remapping/mlinterp.hpp"

#include "pchip.hpp"

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

#define dbgprint( MSG )      \
    do                       \
    {                        \
        cout << MSG << endl; \
    } while( false )

#define dbgprintall( MSG )                           \
    do                                               \
    {                                                \
        cout << "[" << rank << "]: " << MSG << endl; \
    } while( false )

struct RemappingContext
{
    Mesh meshInput;
    Mesh meshOutput;
    Mesh meshOverlap;
    OfflineMap weightMap;
};

ErrorCode ScaleCoords( Interface* mb,
                       std::vector< moab::EntityHandle >& nodes,
                       double R,
                       bool is_cartesian = true,
                       bool is_threed    = false );

moab::ErrorCode ShepardInterpolator( int dimension,
                                     std::vector< double >& xyzd,
                                     std::vector< double >& fd,
                                     std::vector< double >& xyzi,
                                     std::vector< double >& fi );

moab::ErrorCode ModifiedShepardInterpolator( int dimension,
                                             std::vector< double >& xyzd,
                                             std::vector< double >& fd,
                                             std::vector< double >& xyzi,
                                             std::vector< double >& fi,
                                             int order );

moab::ErrorCode ComputeMBAInterpolant( std::vector< double >& xyzd,
                                       std::vector< double >& fd,
                                       std::vector< double >& xyzi,
                                       std::vector< double >& fi,
                                       bool is_threed,
                                       int order );

moab::ErrorCode ComputeNNInterpolant( const std::vector< double >& src_xyz,
                                      const std::vector< double >& src_tdata,
                                      const std::vector< double >& dst_xyz,
                                      std::vector< double >& dst_tdata );

moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& xyzd,
                                            std::vector< double >& fd,
                                            std::vector< double >& xyzi,
                                            std::vector< double >& fi );

moab::ErrorCode ExtrudeMPASPolygonsToPolyhedra( Interface* mb,
                                                std::vector< double >& layer_thickness,
                                                moab::EntityHandle& poly2dset,
                                                moab::EntityHandle& outputset );

moab::ErrorCode ExtrudePolygonsToPolyhedra( Interface* mb,
                                            std::vector< double >& layer_thickness,
                                            moab::EntityHandle& poly2dset,
                                            moab::EntityHandle& outputset,
                                            const bool is_mpas = false );

moab::ErrorCode ComputeFieldProjections( moab::Interface* mbi,
                                         RemappingContext& context,
                                         std::string varProjectSrc,
                                         std::string varProjectDst,
                                         std::vector< moab::EntityHandle >& srcelems,
                                         std::vector< moab::EntityHandle >& dstelems,
                                         bool is_three_dimensional,
                                         bool is_three2x1_dimensional,
                                         bool normalize                                = true,
                                         const double constantoffset                   = 0.0,
                                         const std::string strMethod                   = "mba",
                                         int order                                     = 3,
                                         std::vector< moab::EntityHandle >* src3delems = nullptr,
                                         std::vector< moab::EntityHandle >* dst3delems = nullptr );

// For CAAS
double ApplyCAASLimiting( OfflineMap& mapOperator,
                          Mesh& meshInput,
                          Mesh& meshOverlap,
                          const int nPin,
                          DataArray1D< double >& dataInDouble,
                          DataArray1D< double >& dataOutDouble,
                          bool useCAASLocal );

ErrorCode CloneToTRMesh( moab::Interface* m_interface, Mesh& mesh, EntityHandle mesh_set );

// 3D settings
constexpr int mpas_zreflevels = 60;
constexpr int mpas_zlevels    = 60;
constexpr int roms_zlevels    = 60;
constexpr int nvars           = 2;
int src_zlayers = 0, dst_zlayers = 0;

// tag name data
const char* mpas_twod_tagnames[nvars]       = { "salinity", "temperature" };
const char* mpas_threed_cum_tagnames[nvars] = { "salinity_3d", "temperature_3d" };
const char* mpas_threed_tagnames[nvars]     = { "Salinity3d", "Temperature3d" };
const char* roms_twod_tagnames[nvars]       = { "Salinity2DROMS", "Temperature2DROMS" };
const char* roms_threed_tagnames[nvars]     = { "Salinity3dROMS", "Temperature3dROMS" };

moab::ErrorCode ComputeTempestRemapWeights( moab::Interface* mbi,
                                            RemappingContext& context,
                                            EntityHandle src_set,
                                            EntityHandle tgt_set,
                                            std::string strMethod,
                                            bool ensureMonotonicity );

template < typename T >
struct PC3D
{
    using coord_t = T;  //!< The type of each coordinate

    const static int dimension = 3;
    const std::vector< T >& xyz;
    const size_t count;

    PC3D( const std::vector< T >& pxyz ) : xyz( pxyz ), count( pxyz.size() / dimension ) {}

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

    RemappingContext context;
    int rank, size;
    string mpas_filename, roms_filename, roms_3d_filename;
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
    const double radius          = 1.0;
    int dimension                = 2;
    bool ensureMonotonicity      = false;
    bool computeTR               = false;
    bool computeShepard          = false;
    bool computeMBA              = false;
    bool nearestNeighbor         = false;
    bool threetwooneD            = false;
    bool normalize               = false;
    bool use_3dprojection        = false;
    bool generateExtrusions      = false;
    bool oneDfirst               = false;
    std::string strMethod        = "";
    std::string bathymetryMethod = "bilin";
    int bathymetryOrder          = 3;
    int fieldOrder               = 3;

    constexpr bool useCAAS = true;

    {
        ProgOptions opts;

        // set default values
        mpas_filename    = "mpas_grid.h5m";
        roms_filename    = "roms_grid.h5m";
        roms_3d_filename = "roms_3d_grid.h5m";

        // Input and output meshes
        opts.addOpt< std::string >( "mpas", "MPAS filename with 2D mesh and 3D dataset", &mpas_filename );
        opts.addOpt< std::string >( "roms", "ROMS filename with 2D mesh", &roms_filename );

        // Problem setup
        opts.addOpt< int >( "dimension", "Compute 2D surface or 3D volumetric coupling (default=2)", &dimension );
        opts.addOpt< void >( "setup", "Compute full mesh extrusions needed for coupling in 3D", &generateExtrusions );
        opts.addOpt< std::string >( "method",
                                    "Additional computational method arguments (fv, invdist, bilin, intbilin, "
                                    "delaunay, shepard, mba). default=MBA",
                                    &strMethod );
        opts.addOpt< void >( "mono", "Ensure monotonicity in the weight generation (only for TR-FV)",
                             &ensureMonotonicity );
        opts.addOpt< void >( "321D", "Compute three-dimensional projections using a 2Dx1D approach", &threetwooneD );
        opts.addOpt< void >( "1D2D", "Use 1Dx2D as opposed to 2Dx1D for 321D projection", &oneDfirst );
        ;
        opts.addOpt< void >(
            "normalize",
            "Re-normalize interpolant to preserve global field integral (only 2D and requires mesh intersection)",
            &normalize );
        opts.addOpt< int >( "bathymetryOrder",
                            "Specify order of Bathymetry reconstruction. \n"
                            "\tTR: method='' -> FV order, method='invdist,bilin,intbilin' -> order 2, \n"
                            "\tShepard: shepard_power=order\n"
                            "\tMBA: order=1 -> bilinear, else order 3\n"
                            "(default=MBA3)",
                            &bathymetryOrder );
        opts.addOpt< int >( "fieldOrder",
                            "Specify order for Temperature and Salinity field projection. \n"
                            "\tTR: method='' -> FV order, method='invdist,bilin,intbilin' -> order 2, \n"
                            "\tShepard: shepard_power=order\n"
                            "\tMBA: order=1 -> bilinear, else order 3\n"
                            "(default=MBA1)",
                            &fieldOrder );

        opts.parseCommandLine( argc, argv );

        if( strMethod == "shepard" )
        {
            computeShepard = true;
        }
        else if( strMethod == "mba" )
        {
            computeMBA = true;
            // if( order == 1 ) strMethod = "mba:linear";
            // else strMethod = "mba:cubic";
        }
        else if( strMethod == "nn" )
        {
            nearestNeighbor = true;
        }
        else
        {
            computeTR = true;
            if( strMethod == "fv" ) strMethod = "";  // no sub-method necessary
        }

        if( threetwooneD )
        {
            strMethod = "bilin";
            computeTR = true;
        }

        if( dimension == 3 ) use_3dprojection = true;
        if( !computeMBA && !computeShepard && !computeTR ) computeMBA = true;
    }

    //
    src_zlayers = use_3dprojection ? mpas_zlevels : 1;
    dst_zlayers = use_3dprojection ? roms_zlevels : 1;

    // only MBA is right now tested with 3D projections?
    if( use_3dprojection ) computeMBA = true;

    // Print usage if not enough arguments
    if( argc < 1 )
    {
        cerr << "Usage: ";
        cerr << argv[0] << " --mpas file_name --roms file_name --out file_name" << endl;
#ifdef MOAB_HAVE_MPI
        ierr = MPI_Finalize();
        MPICHKERR( ierr, "MPI_Finalize failed; Aborting" );
#endif

        return 1;
    }

#ifdef MOAB_HAVE_MPI
    // Initialize MPI first
    ierr = MPI_Init( &argc, &argv );
    MPICHKERR( ierr, "MPI_Init failed" );

    MPI_Comm comm = MPI_COMM_WORLD;

    ierr = MPI_Comm_rank( comm, &rank );
    MPICHKERR( ierr, "MPI_Comm_rank failed" );
    ierr = MPI_Comm_size( comm, &size );
    MPICHKERR( ierr, "MPI_Comm_size failed" );
#else
    rank = 0;
    size = 1;
#endif
    string mpas_read_options = size > 1 ? "" : "";
    string roms_read_options = size > 1 ? "" : "";
    string write_options     = size > 1 ? "PARALLEL=WRITE_PART" : "";

    if( !rank ) dbgprint( "********** Remap MPAS-to-ROMS **********\n" );

    // Create the moab instance
    Interface* mbi = new( std::nothrow ) Core;
    if( NULL == mbi ) return 1;

    // Print out the input parameters
    if( !rank )
    {
        dbgprint( " -- Runtime Parameters -- " );
        dbgprint( "   MPAS mesh file: " << mpas_filename );
        dbgprint( "   ROMS mesh file: " << roms_filename );

        dbgprint( "        Dimension: " << dimension );
        dbgprint( "        Algorithm: " << ( nearestNeighbor ? "Nearest Neighbor mapping"
                                             : computeTR     ? "TempestRemap Conservative mapping"
                                             : computeMBA    ? "Multilevel B-spline Approximation"
                                                             : "Shepard interpolant" ) );
        if( computeTR ) dbgprint( "           Method: " << strMethod );
        dbgprint( " Bathymetry Order: " << bathymetryOrder );
        dbgprint( "      Field Order: " << fieldOrder );
        dbgprint( endl );
    }

    // construct the remapper
    // #ifdef MOAB_HAVE_MPI
    //     EntityHandle partnset;
    //     err = mbi->create_meshset( MESHSET_SET, partnset );MB_CHK_SET_ERR( err, "Creating partition set failed" );
    //     // Create the parallel communicator object with the partition handle associated with MOAB
    //     ParallelComm* parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );
    //     moab::TempestRemapper remapper( mbi, parallel_communicator );
    // #else
    //     moab::TempestRemapper remapper( mbi );
    // #endif
    //     remapper.meshValidate     = false;
    //     remapper.constructEdgeMap = false;
    //     remapper.initialize();

    EntityHandle mpasset, mpas_covering_set, romsset;
    err = mbi->create_meshset( moab::MESHSET_SET, mpasset );MB_CHK_SET_ERR( err, "Can't create new set" );
    err = mbi->create_meshset( moab::MESHSET_SET, mpas_covering_set );MB_CHK_SET_ERR( err, "Can't create new set" );
    err = mbi->create_meshset( moab::MESHSET_SET, romsset );MB_CHK_SET_ERR( err, "Can't create new set" );

    // Load the MPAS file from disk with given options
    std::vector< moab::EntityHandle > mpas_verts, mpas_elems;
    std::vector< moab::EntityHandle > mpas3d_verts, mpas3d_elems;
    {
        dbgprint( "Reading MPAS file from disk" );
        err = mbi->load_file( mpas_filename.c_str(), &mpasset, mpas_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for MPAS mesh failed" );

        // Tag depthtag;
        // err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, depthtag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
        // if( err != MB_SUCCESS ) dbgprint( "Error: " << err << "; Failed to get bottomDepth tag handle" );
        // MB_CHK_SET_ERR( err, "MPAS bottomDepth tag failed" );

        // Get all entities in the database
        err = mbi->get_entities_by_dimension( mpasset, 0, mpas_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( mpasset, 2, mpas_elems );MB_CHK_ERR( err );
        dbgprint( "MPAS mesh contains " << mpas_verts.size() << " vertices and " << mpas_elems.size() << " elements" );

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, mpas_verts, radius, true, false );MB_CHK_ERR( err );
        err = mbi->write_file( "mpas_modified_2d.h5m", "H5M", write_options.c_str(), &mpasset, 1 );MB_CHK_ERR( err );
    }

    // Load the ROMS file from disk with given options
    std::vector< moab::EntityHandle > roms_verts, roms_elems;
    std::vector< moab::EntityHandle > roms3d_verts, roms3d_elems;
    {
        dbgprint( "Reading ROMS file from disk" );
        err = mbi->load_file( roms_filename.c_str(), &romsset, roms_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS mesh failed" );

        // Get all entities in the database
        err = mbi->get_entities_by_dimension( romsset, 0, roms_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( romsset, 2, roms_elems );MB_CHK_ERR( err );
        dbgprint( "ROMS mesh contains " << roms_verts.size() << " vertices and " << roms_elems.size() << " elements" );

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, roms_verts, radius, false, false );MB_CHK_ERR( err );
        err = mbi->write_file( "roms_modified_2d.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );
    }

    // Cull the MPAS set so that we don't have a global mesh
    {
        const int nring_neighborhood = 1;
        // construct a kd-tree index:
        using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                            PC3D< double >, 3 /* dim */
                                                            >;

        moab::Range orig_mpas_elems;
        err = mbi->get_entities_by_dimension( mpasset, 2, orig_mpas_elems );MB_CHK_ERR( err );
        std::vector< double > mpas_xyz( orig_mpas_elems.size() * 3 );
        err = mbi->get_coords( orig_mpas_elems, mpas_xyz.data() );MB_CHK_ERR( err );

        double query_pt[3];  // dimension
        PC3D< double > cloud( mpas_xyz );
        KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

        const size_t num_results = static_cast< size_t >(
            ( 2 * nring_neighborhood + 1 ) * ( 2 * nring_neighborhood + 1 ) - nring_neighborhood );
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        nanoflann::KNNResultSet< double > resultSet( num_results );
        mpas_elems.clear();
        mpas_verts.clear();
        moab::Range lelems, lverts;
        for( size_t i = 0; i < roms_elems.size(); i++ )
        {
            const moab::EntityHandle ehandle = roms_elems[i];
            err                              = mbi->get_coords( &ehandle, 1, &query_pt[0] );MB_CHK_ERR( err );

            // Do a KNN search
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            for( size_t j = 0; j < num_results; ++j )
                lelems.insert( orig_mpas_elems[srcindx[j]] );
        }
        err = mbi->add_entities( mpas_covering_set, lelems );MB_CHK_ERR( err );

        err = mbi->write_file( "mpas_covering_2d.h5m", "H5M", write_options.c_str(), &mpas_covering_set, 1 );MB_CHK_ERR( err );

        err = mbi->get_connectivity( lelems, lverts, true );MB_CHK_ERR( err );
        err = mbi->add_entities( mpas_covering_set, lverts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( mpas_covering_set, 2, mpas_elems );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( mpas_covering_set, 0, mpas_verts );MB_CHK_ERR( err );
        dbgprint( "Culled MPAS mesh contains " << mpas_elems.size() << " elements and " << mpas_verts.size()
                                               << " vertices." );
    }

    if( normalize || computeTR || true )
    {
        // call compute 2D map
        err = ComputeTempestRemapWeights( mbi, context, mpas_covering_set, romsset, bathymetryMethod,
                                          ensureMonotonicity );MB_CHK_SET_ERR( err, "Canot compute 2D remapping weights" );
    }

    constexpr bool ProjectMPASBathymetryToROMS = true;

    // let us perform 3D extrusions as needed
    EntityHandle root_set  = 0;
    EntityHandle mpasset3d = 0, romsset3d = 0;

    if( use_3dprojection || threetwooneD )
    {
        if( generateExtrusions )
        {
            constexpr bool useConstantRefAxialThickness = true;
            std::vector< double > zmh_xyz3d, zrh_xyz3d, ref_zmh_z1d;
            zmh_xyz3d.resize( mpas_elems.size() * src_zlayers );
            zrh_xyz3d.resize( roms_elems.size() * dst_zlayers );
            ref_zmh_z1d.resize( mpas_zreflevels );

            if( useConstantRefAxialThickness )
            {
                moab::Tag mztag;
                err = mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag,
                                           moab::MB_TAG_SPARSE );MB_CHK_SET_ERR( err, "Can't get tag handle: refBottomDepth" );

                err = mbi->tag_get_data( mztag, &root_set, 1, ref_zmh_z1d.data() );MB_CHK_SET_ERR( err, "Can't get refBottomDepth data" );
                // constexpr double mpas_ref_levels[] = {
                //     10,      20,      30,      40,      50,      60,      70,      80,      90,      100,
                //     110,     120,     130,     140,     150,     160,     170.197, 180.761, 191.821, 203.499,
                //     215.923, 229.233, 243.584, 259.156, 276.152, 294.815, 315.424, 338.312, 363.875, 392.58,
                //     424.989, 461.767, 503.707, 551.749, 606.997, 670.729, 744.398, 829.607, 928.043, 1041.37,
                //     1171.04, 1318.09, 1482.9,  1664.99, 1863.01, 2074.87, 2298.04, 2529.9,  2768.1,  3010.67,
                //     3256.14, 3503.45, 3751.89, 4001.01, 4250.53, 4500.26, 4750.12, 5000.05, 5250.01, 5499.99 };
                for( size_t i = 0; i < mpas_elems.size(); ++i )
                {
                    const int offset = i * src_zlayers;
                    // zmh_xyz3d[offset] = 0;
                    // for( int j = 1; j <= src_zlayers; ++j )
                    //     zmh_xyz3d[offset + j] = ref_zmh_z1d[j - 1] + zmh_xyz3d[offset + j - 1];
                    // for( int j = 0; j < src_zlayers; ++j )
                    //     zmh_xyz3d[offset + j] = ref_zmh_z1d[j];

                    zmh_xyz3d[offset] = ref_zmh_z1d[0];
                    for( int j = 1; j < src_zlayers; ++j )
                        zmh_xyz3d[offset + j] = ref_zmh_z1d[j] - ref_zmh_z1d[j - 1];
                }
            }
            else
            {
                moab::Tag mhtag;
                err = mbi->tag_get_handle( "layerThickness_3d", src_zlayers, moab::MB_TYPE_DOUBLE, mhtag,
                                           moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get tag handle: layerThickness_3d" );
                err = mbi->tag_get_data( mhtag, mpas_elems.data(), mpas_elems.size(), zmh_xyz3d.data() );MB_CHK_SET_ERR( err, "Can't get layerThickness_3d data" );
            }

            dbgprint( "\nExtruding MPAS polygonal mesh ..." );
            // err = ExtrudeMPASPolygonsToPolyhedra( mbi, zmh_xyz3d, mpas_covering_set, mpasset3d );MB_CHK_SET_ERR( err, "Can't extrude MPAS polygons" );
            err = ExtrudePolygonsToPolyhedra( mbi, zmh_xyz3d, mpas_covering_set, mpasset3d, true );MB_CHK_SET_ERR( err, "Can't extrude MPAS polygons" );

            if( ProjectMPASBathymetryToROMS )
            {
                // Project the bottom Bathymetry data from MPAS to ROMS so that we can impose it.
                err = ComputeFieldProjections( mbi, context, "bottomDepth", "bottomDepth", mpas_elems, roms_elems,
                                               false /*use_3dprojection*/, false /* 2Dx1D */,
                                               false /* bool normalize */, 2000.0, bathymetryMethod, bathymetryOrder );MB_CHK_SET_ERR( err, "Can't create new set" );

                std::vector< double > zrh_xyz2d( roms_elems.size() );
                moab::Tag rhtag;
                err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
                err = mbi->tag_get_data( rhtag, roms_elems.data(), roms_elems.size(), zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

                // err = mbi->write_file( "roms_3d_2dsurface.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );

                for( size_t i = 0; i < roms_elems.size(); ++i )
                {
                    const double pbathymetry = zrh_xyz2d[i];
                    const double delz        = pbathymetry / dst_zlayers;
                    const int offset         = i * dst_zlayers;
                    for( int j = 0; j < dst_zlayers; ++j )
                    {
                        zrh_xyz3d[offset + j] = delz;
                        // printf( "Thickness value for ROMS element %zu, %d = %f, %f\n", i, dst_zlayers, zrh_xyz2d[i], delz );
                    }
                }

                dbgprint( "\nExtruding ROMS structured quad mesh..." );
                err = ExtrudePolygonsToPolyhedra( mbi, zrh_xyz3d, romsset, romsset3d, false );MB_CHK_SET_ERR( err, "Can't extrude ROMS polygons" );
            }
            else
            {
                err = mbi->create_meshset( moab::MESHSET_SET, romsset3d );MB_CHK_SET_ERR( err, "Can't create new set" );

                dbgprint( "Reading ROMS 3D mesh file from disk" );
                err = mbi->load_file( roms_3d_filename.c_str(), &romsset3d, roms_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS 3D mesh failed" );
            }

            err = mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems );MB_CHK_ERR( err );

            // std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size() << " elements.\n";
            err = mbi->write_file( "mpas_full_3d.h5m", "H5M", write_options.c_str(), &mpasset3d, 1 );MB_CHK_ERR( err );

            err = mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems );MB_CHK_ERR( err );

            // std::cout << "3D ROMS: " << roms3d_verts.size() << " vertices and " << roms3d_elems.size() << " elements.\n";
            // Rescale the radius of both to compute the intersection
            // err = ScaleCoords( mbi, roms3d_verts, radius, false, true );MB_CHK_ERR( err );
            err = mbi->write_file( "roms_full_3d.h5m", "H5M", write_options.c_str(), &romsset3d, 1 );MB_CHK_ERR( err );
        }
        else
        {
            err = mbi->create_meshset( moab::MESHSET_SET, mpasset3d );MB_CHK_SET_ERR( err, "Can't create new set" );
            err = mbi->create_meshset( moab::MESHSET_SET, romsset3d );MB_CHK_SET_ERR( err, "Can't create new set" );

            err = mbi->load_file( "mpas_full_3d.h5m", &mpasset3d, mpas_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for MPAS 3D mesh failed" );
            err = mbi->load_file( "roms_full_3d.h5m", &romsset3d, roms_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS 3D mesh failed" );

            err = mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems );MB_CHK_ERR( err );

            std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size()
                      << " elements.\n";
            std::cout << "3D ROMS: " << roms3d_verts.size() << " vertices and " << roms3d_elems.size()
                      << " elements.\n";

            // // get MPAS and ROMS height factors for elements
            // {
            //     moab::Tag mztag;
            //     err = mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag,
            //                                moab::MB_TAG_SPARSE );MB_CHK_SET_ERR( err, "Can't get tag handle: refBottomDepth" );

            //     err = mbi->tag_get_data( mztag, &root_set, 1, ref_zmh_z1d.data() );MB_CHK_SET_ERR( err, "Can't get refBottomDepth data" );
            //     for( size_t i = 0; i < mpas_elems.size(); ++i )
            //     {
            //         const int offset = i * src_zlayers;
            //         zmh_xyz3d[offset] = ref_zmh_z1d[0];
            //         for( int j = 1; j < src_zlayers; ++j )
            //             zmh_xyz3d[offset + j] = ref_zmh_z1d[j] - ref_zmh_z1d[j - 1];
            //     }

            //     // Now ROMS
            //     std::vector< double > zrh_xyz2d( roms_elems.size() );
            //     moab::Tag rhtag;
            //     err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
            //     err = mbi->tag_get_data( rhtag, roms_elems.data(), roms_elems.size(), zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

            //     for( size_t i = 0; i < roms_elems.size(); ++i )
            //     {
            //         const double pbathymetry = zrh_xyz2d[i];
            //         const double delz        = pbathymetry / dst_zlayers;
            //         const int offset         = i * dst_zlayers;
            //         for( int j = 0; j < dst_zlayers; ++j )
            //         {
            //             zrh_xyz3d[offset + j] = delz;
            //         }
            //     }
            // }
        }
    }
    else
    {
        // Project the bottom Bathymetry data from MPAS to ROMS so that we can impose it.
        err = ComputeFieldProjections( mbi, context, "bottomDepth", "bottomDepth", mpas_elems, roms_elems,
                                       false /*use_3dprojection*/, false /* 2Dx1D */, false /* bool normalize */,
                                       2000.0, bathymetryMethod, bathymetryOrder );MB_CHK_SET_ERR( err, "Can't create new set" );
    }

    bool twoDfirst = !oneDfirst;
#define VERTICAL_INTERPOLATION

    if( threetwooneD )
    {

        // Now let us compute the mba hierarchy for each field
        // err =
        //     ComputeFieldProjections( mbi, context,
        //                              mpas_threed_tagnames[0],
        //                              roms_threed_tagnames[0], mpas_elems,
        //                              roms_elems, true, false /* 2Dx1D */, normalize /* bool normalize */,
        //                              35.0, "nn", fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
        // err =
        //     ComputeFieldProjections( mbi, context,
        //                              mpas_threed_tagnames[1],
        //                              roms_threed_tagnames[1], mpas_elems,
        //                              roms_elems, true, false /* 2Dx1D */, normalize /* bool normalize */,
        //                              6, "nn", fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );

        const size_t mpassize = mpas_elems.size();
        const size_t romssize = roms_elems.size();

        // context.weightMap.SetEnforcementBounds( "Lp", &context.meshInput, &context.meshOverlap, nullptr, nullptr, 1 );

        // get the handle to the weight matrix
        const SparseMatrix< double >& weights = context.weightMap.GetSparseMatrix();

        double defaultvalue = -1.0;
        moab::Tag mpas_soltags_3d[2], roms_soltags_elem[2];
        err = mbi->tag_get_handle( mpas_threed_cum_tagnames[0], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
                                   mpas_soltags_3d[0], moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't create salinity tag" );
        err = mbi->tag_get_handle( mpas_threed_cum_tagnames[1], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
                                   mpas_soltags_3d[1], moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't create temperature tag" );
        err = mbi->tag_get_handle( roms_threed_tagnames[0], 1, moab::MB_TYPE_DOUBLE, roms_soltags_elem[0],
                                   moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defaultvalue );MB_CHK_SET_ERR( err, "Can't create salinity tag on ROMS3D" );
        err = mbi->tag_get_handle( roms_threed_tagnames[1], 1, moab::MB_TYPE_DOUBLE, roms_soltags_elem[1],
                                   moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defaultvalue );MB_CHK_SET_ERR( err, "Can't create temperature tag on ROMS3D" );

        std::vector< double > src_salinity_data( mpassize * mpas_zreflevels ),
            src_temperature_data( mpassize * mpas_zreflevels );
        // get the source data from tag
        err = mbi->tag_get_data( mpas_soltags_3d[0], mpas_elems.data(), mpassize, src_salinity_data.data() );MB_CHK_ERR( err );
        err = mbi->tag_get_data( mpas_soltags_3d[1], mpas_elems.data(), mpassize, src_temperature_data.data() );MB_CHK_ERR( err );

        std::vector< double > tgt_salinity_data( romssize * roms_zlevels ),
            tgt_temperature_data( romssize * roms_zlevels );

        std::vector< double > zmh_ztmp( mpas_zreflevels );
        moab::Tag mztag;
        err =
            mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag, moab::MB_TAG_SPARSE );MB_CHK_SET_ERR( err, "Can't get tag handle: refBottomDepth" );

        err = mbi->tag_get_data( mztag, &root_set, 1, zmh_ztmp.data() );MB_CHK_SET_ERR( err, "Can't get refBottomDepth data" );

        moab::Tag rhtag;
        err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
        std::vector< double > zrh_xyz2d( romssize );
        err = mbi->tag_get_data( rhtag, roms_elems.data(), romssize, zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

        std::vector< double > zmh_z( mpas_zreflevels );
        zmh_z[0] = zmh_ztmp[0] / 2;
        for( size_t j = 1; j < mpas_zlevels; ++j )
        {
            zmh_z[j] = zmh_ztmp[j - 1] + 0.5 * ( zmh_ztmp[j] - zmh_ztmp[j - 1] );
        }

        if( twoDfirst )
        {
            // First compute the projections in 2D for each MPAS layer.
            // This will give the MPAS projected data on to ROMS mesh - on MPAS axial levels
            // NOTE: implicit assumption is that ROMS levels >= MPAS levels. Should fix how
            // the projection is invoked by perhaps skipping the tag_set_data

            std::vector< double > tgtsrc_salinity_data( romssize * mpas_zreflevels ),
                tgtsrc_temperature_data( romssize * mpas_zreflevels );

            DataArray1D< double > dataInDoubleS( mpassize ), dataInDoubleT( mpassize );
            unsigned offsetm = 0, offsetr = 0;
            for( int ii = 0; ii < mpas_zlevels; ii++ )
            {
                dbgprint( "Computing projection for MPAS level: " << ii );

                DataArray1D< double > dataOutDoubleS( romssize, false ), dataOutDoubleT( romssize, false );

                std::vector< EntityHandle > mpas_slice( mpas3d_elems.begin() + offsetm,
                                                        mpas3d_elems.begin() + offsetm + mpassize );
                std::vector< EntityHandle > roms_slice( roms3d_elems.begin() + offsetr,
                                                        roms3d_elems.begin() + offsetr + romssize );

                for( size_t j = 0; j < mpassize; ++j )
                {
                    const size_t offset = j * mpas_zlevels;
                    dataInDoubleS[j]    = src_salinity_data[ii + offset];
                    dataInDoubleT[j]    = src_temperature_data[ii + offset];
                }

                // dataInDoubleS.AttachToData( src_salinity_data.data() + offsetm );
                dataOutDoubleS.AttachToData( tgtsrc_salinity_data.data() + offsetr );

                // Compute the projection for the salinity field
                weights.Apply( dataInDoubleS, dataOutDoubleS );

                if( useCAAS )
                    ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap, 1 /*nPin*/,
                                       dataInDoubleS, dataOutDoubleS, false /*useCAASLocal*/ );

#ifndef VERTICAL_INTERPOLATION
                err = mbi->tag_set_data( roms_soltags_elem[0], roms_slice.data(), roms_slice.size(), dataOutDoubleS );MB_CHK_SET_ERR( err, "Can't set salinity tag data" );
#endif
                // dataInDoubleT.AttachToData( src_temperature_data.data() + offsetm );
                dataOutDoubleT.AttachToData( tgtsrc_temperature_data.data() + offsetr );

                // Compute the projection for the temperature field
                weights.Apply( dataInDoubleT, dataOutDoubleT );

                if( useCAAS )
                    ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap, 1 /*nPin*/,
                                       dataInDoubleT, dataOutDoubleT, false /*useCAASLocal*/ );

#ifndef VERTICAL_INTERPOLATION
                err = mbi->tag_set_data( roms_soltags_elem[1], roms_slice.data(), roms_slice.size(), dataOutDoubleT );MB_CHK_SET_ERR( err, "Can't set temperature tag data" );
#endif
                // err = ComputeFieldProjections( mbi, context, mpas_threed_tagnames[0], roms_threed_tagnames[0], mpas_slice,
                //                                roms_slice, false, true /* 2Dx1D */, normalize /* bool normalize */, 35.0,
                //                                strMethod, fieldOrder );MB_CHK_SET_ERR( err, "Can't project salinity data" );
                // err = ComputeFieldProjections( mbi, context, mpas_threed_tagnames[1], roms_threed_tagnames[1], mpas_slice,
                //                                roms_slice, false, true /* 2Dx1D */, normalize /* bool normalize */, 5.0,
                //                                strMethod, fieldOrder );MB_CHK_SET_ERR( err, "Can't project temperature data" );

                offsetm += mpassize;
                offsetr += romssize;
            }

#ifdef VERTICAL_INTERPOLATION
            // First MPAS z-levels

            // Now ROMS

            // Now project each data layer in axial direction.
            // NOTE: Embarassingly parallel
            std::vector< double > roms_zvalsS( mpas_zlevels ), roms_zvalsT( mpas_zlevels );

            for( size_t i = 0; i < romssize; ++i )
            {
                const double pbathymetry = zrh_xyz2d[i];
                const double delz        = pbathymetry / dst_zlayers;

                for( size_t j = 0; j < mpas_zlevels; ++j )
                {
                    const size_t offset = j * romssize;
                    roms_zvalsS[j]      = tgtsrc_salinity_data[i + offset];
                    roms_zvalsT[j]      = tgtsrc_temperature_data[i + offset];
                }

                for( int j = 0; j < dst_zlayers; ++j )
                {
                    const size_t offset   = j * romssize;
                    double roms_zlocation = ( j + 0.5 ) * delz;
                    // Project data from zmh_z to zrh_z
                    tgt_salinity_data[i + offset]    = pchipInterpolate( zmh_z, roms_zvalsS, roms_zlocation );
                    tgt_temperature_data[i + offset] = pchipInterpolate( zmh_z, roms_zvalsT, roms_zlocation );
                    // tgt_salinity_data[i + offset]    = linearInterpolate( zmh_z, roms_zvalsS, roms_zlocation );
                    // tgt_temperature_data[i + offset] = linearInterpolate( zmh_z, roms_zvalsT, roms_zlocation );
                }
            }
#endif
        }
        else  // Vertical first and horizontal next
        {

            std::vector< double > srctgt_salinity_data( mpassize * roms_zlevels ),
                srctgt_temperature_data( mpassize * roms_zlevels );

#ifdef VERTICAL_INTERPOLATION
            // Now project each data layer in axial direction.
            // NOTE: Embarassingly parallel
            std::vector< double > mpasroms_zvalsS( mpas_zlevels ),
                mpasroms_zvalsT( mpas_zlevels );  // Data on MPAS 2D but ROMS 1D vertical

            for( size_t i = 0; i < mpassize; ++i )
            {
                const size_t offset = i * mpas_zlevels;
                for( size_t j = 0; j < mpas_zlevels; ++j )
                {
                    mpasroms_zvalsS[j] = src_salinity_data[j + offset];
                    mpasroms_zvalsT[j] = src_temperature_data[j + offset];
                    // roms_zvalsS[mpas_zlevels - 1 - j] = tgtsrc_salinity_data[i + offset];
                    // roms_zvalsT[mpas_zlevels - 1 - j] = tgtsrc_temperature_data[i + offset];
                }

                const size_t offsetr = i * roms_zlevels;
                for( int j = 0; j < roms_zlevels; ++j )
                {
                    const double pbathymetry = zrh_xyz2d[i];
                    const double delz        = pbathymetry / dst_zlayers;

                    double roms_zlocation                = ( j + 0.5 ) * delz;
                    srctgt_salinity_data[j + offsetr]    = pchipInterpolate( zmh_z, mpasroms_zvalsS, roms_zlocation );
                    srctgt_temperature_data[j + offsetr] = pchipInterpolate( zmh_z, mpasroms_zvalsT, roms_zlocation );
                    // tgt_salinity_data[i + offset]    = linearInterpolate( zmh_z, mpasroms_zvalsT, roms_zlocation );
                    // tgt_temperature_data[i + offset] = linearInterpolate( zmh_z, mpasroms_zvalsT, roms_zlocation );
                }
            }
#endif

            DataArray1D< double > dataInDoubleS( mpassize ), dataInDoubleT( mpassize );
            unsigned offsetr = 0;
            for( int ii = 0; ii < roms_zlevels; ii++ )
            {
                dbgprint( "Computing projection for ROMS level: " << ii );

                DataArray1D< double > dataOutDoubleS( romssize, false ), dataOutDoubleT( romssize, false );

                // std::vector< EntityHandle > mpas_slice( mpas3d_elems.begin() + offsetr,
                //                                         mpas3d_elems.begin() + offsetr + mpassize );
                std::vector< EntityHandle > roms_slice( roms3d_elems.begin() + offsetr,
                                                        roms3d_elems.begin() + offsetr + romssize );

                for( size_t j = 0; j < mpassize; ++j )
                {
                    const size_t offset = j * roms_zlevels;
                    dataInDoubleS[j]    = srctgt_salinity_data[ii + offset];
                    dataInDoubleT[j]    = srctgt_temperature_data[ii + offset];
                }

                // dataInDoubleS.AttachToData( src_salinity_data.data() + offsetr );
                dataOutDoubleS.AttachToData( tgt_salinity_data.data() + offsetr );

                // Compute the projection for the salinity field
                weights.Apply( dataInDoubleS, dataOutDoubleS );

                if( useCAAS )
                    ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap, 1 /*nPin*/,
                                       dataInDoubleS, dataOutDoubleS, true /*useCAASLocal*/ );

#ifndef VERTICAL_INTERPOLATION
                err = mbi->tag_set_data( roms_soltags_elem[0], roms_slice.data(), roms_slice.size(), dataOutDoubleS );MB_CHK_SET_ERR( err, "Can't set salinity tag data" );
#endif
                // dataInDoubleT.AttachToData( src_temperature_data.data() + offsetr );
                dataOutDoubleT.AttachToData( tgt_temperature_data.data() + offsetr );

                // Compute the projection for the temperature field
                weights.Apply( dataInDoubleT, dataOutDoubleT );

                if( useCAAS )
                    ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap, 1 /*nPin*/,
                                       dataInDoubleT, dataOutDoubleT, true /*useCAASLocal*/ );

#ifndef VERTICAL_INTERPOLATION
                err = mbi->tag_set_data( roms_soltags_elem[1], roms_slice.data(), roms_slice.size(), dataOutDoubleT );MB_CHK_SET_ERR( err, "Can't set temperature tag data" );
#endif
                // err = ComputeFieldProjections( mbi, context, mpas_threed_tagnames[0], roms_threed_tagnames[0], mpas_slice,
                //                                roms_slice, false, true /* 2Dx1D */, normalize /* bool normalize */, 35.0,
                //                                strMethod, fieldOrder );MB_CHK_SET_ERR( err, "Can't project salinity data" );
                // err = ComputeFieldProjections( mbi, context, mpas_threed_tagnames[1], roms_threed_tagnames[1], mpas_slice,
                //                                roms_slice, false, true /* 2Dx1D */, normalize /* bool normalize */, 5.0,
                //                                strMethod, fieldOrder );MB_CHK_SET_ERR( err, "Can't project temperature data" );

                offsetr += romssize;
            }
        }

#ifdef VERTICAL_INTERPOLATION
        err = mbi->tag_set_data( roms_soltags_elem[0], roms3d_elems.data(), roms3d_elems.size(),
                                 tgt_salinity_data.data() );MB_CHK_SET_ERR( err, "Can't set salinity tag data" );
        err = mbi->tag_set_data( roms_soltags_elem[1], roms3d_elems.data(), roms3d_elems.size(),
                                 tgt_temperature_data.data() );MB_CHK_SET_ERR( err, "Can't set temperature tag data" );
#endif
    }
    else if( use_3dprojection || computeMBA )
    {
        // Now let us compute the mba hierarchy for each field
        err =
            ComputeFieldProjections( mbi, context,
                                     ( use_3dprojection ? mpas_threed_tagnames[0] : mpas_twod_tagnames[0] ),
                                     ( use_3dprojection ? roms_threed_tagnames[0] : roms_twod_tagnames[0] ), mpas_elems,
                                     roms_elems, use_3dprojection, false /* 2Dx1D */, normalize /* bool normalize */,
                                     35.0, strMethod, fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
        err =
            ComputeFieldProjections( mbi, context,
                                     ( use_3dprojection ? mpas_threed_tagnames[1] : mpas_twod_tagnames[1] ),
                                     ( use_3dprojection ? roms_threed_tagnames[1] : roms_twod_tagnames[1] ), mpas_elems,
                                     roms_elems, use_3dprojection, false /* 2Dx1D */, normalize /* bool normalize */,
                                     8.5, strMethod, fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
    }
    else if( computeShepard )
    {
        // Now let us compute the Shepard's interpolant to compute data for each field
        err =
            ComputeFieldProjections( mbi, context,
                                     ( use_3dprojection ? mpas_threed_tagnames[0] : mpas_twod_tagnames[0] ),
                                     ( use_3dprojection ? roms_threed_tagnames[0] : roms_twod_tagnames[0] ), mpas_elems,
                                     roms_elems, use_3dprojection, false /* 2Dx1D */, normalize /* bool normalize */,
                                     35.0, strMethod, fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
        err =
            ComputeFieldProjections( mbi, context,
                                     ( use_3dprojection ? mpas_threed_tagnames[1] : mpas_twod_tagnames[1] ),
                                     ( use_3dprojection ? roms_threed_tagnames[1] : roms_twod_tagnames[1] ), mpas_elems,
                                     roms_elems, use_3dprojection, false /* 2Dx1D */, normalize /* bool normalize */,
                                     8.5, strMethod, fieldOrder, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
    }
    else
    {
        // Now apply the map to compute the field projections on the ROMS mesh
        // Now let us apply the weights onto the vector and project onto target mesh
        // err = weightMap.ApplyWeights( stag, stag, false );MB_CHK_ERR( err );
        // err = weightMap.ApplyWeights( ttag, ttag, false );MB_CHK_ERR( err );

        // get the handle to the weight matrix
        const SparseMatrix< double >& weights = context.weightMap.GetSparseMatrix();
        DataArray1D< double > dataInDouble( mpas_elems.size() );
        DataArray1D< double > dataOutDouble( roms_elems.size() );

        // assert( useConservativeSalinity );
        {
            constexpr double salinity_avg = 35.0;
            moab::Tag stag;
            err = mbi->tag_get_handle( "salinity", 1, moab::MB_TYPE_DOUBLE, stag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

            // Apply the map onto the salinity solution field
            err = mbi->tag_get_data( stag, mpas_elems.data(), mpas_elems.size(), dataInDouble );MB_CHK_ERR( err );
            for( size_t i = 0; i < dataInDouble.GetRows(); i++ )
                dataInDouble[i] -= salinity_avg;

            // Compute the projection for the salinity field
            weights.Apply( dataInDouble, dataOutDouble );

            // Scale data values
            for( size_t i = 0; i < dataOutDouble.GetRows(); i++ )
                dataOutDouble[i] += salinity_avg;

            // Set the salinity solution field on the ROMS mesh
            err = mbi->tag_set_data( stag, roms_elems.data(), roms_elems.size(), dataOutDouble );MB_CHK_ERR( err );
        }

        // Apply the map onto the temperature solution field
        // assert( useConservativeTemperature );
        {
            constexpr double temperature_avg = 8.5;
            moab::Tag ttag;
            err = mbi->tag_get_handle( "temperature", 1, moab::MB_TYPE_DOUBLE, ttag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

            // Apply the map onto the salinity solution field
            err = mbi->tag_get_data( ttag, mpas_elems.data(), mpas_elems.size(), dataInDouble );MB_CHK_ERR( err );
            for( size_t i = 0; i < dataInDouble.GetRows(); i++ )
                dataInDouble[i] -= temperature_avg;

            // Compute the projection for the salinity field
            weights.Apply( dataInDouble, dataOutDouble );

            // Scale data values
            for( size_t i = 0; i < dataOutDouble.GetRows(); i++ )
                dataOutDouble[i] += temperature_avg;

            // Set the temperature solution field on the ROMS mesh
            err = mbi->tag_set_data( ttag, roms_elems.data(), roms_elems.size(), dataOutDouble );MB_CHK_ERR( err );
        }

        // remapper.clear();
    }
    dbgprint( std::endl );

    const std::string mpas_output_file = ( use_3dprojection ? "mpas_3d_source.h5m" : "mpas_2d_source.h5m" );
    dbgprint( "Writing out the source mesh with fields to '" << mpas_output_file << "'" );
    err = mbi->write_file( mpas_output_file.c_str(), "H5M", write_options.c_str(),
                           ( use_3dprojection || threetwooneD ? &mpasset3d : &mpas_covering_set ), 1 );MB_CHK_ERR( err );

    const std::string roms_output_file = ( use_3dprojection ? "roms_3d_projected.h5m" : "roms_2d_projected.h5m" );
    dbgprint( "Writing out the target mesh with projected fields to '" << roms_output_file << "'" );
    err = mbi->write_file( roms_output_file.c_str(), "H5M", write_options.c_str(),
                           ( use_3dprojection || threetwooneD ? &romsset3d : &romsset ), 1 );MB_CHK_ERR( err );

    // Done, cleanup
    delete mbi;

    if( !rank ) dbgprint( "\n********** Remap MPAS-to-ROMS DONE! **********" );

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

ErrorCode ScaleCoords( Interface* mb,
                       std::vector< moab::EntityHandle >& nodes,
                       double R,
                       bool is_cartesian,
                       bool is_threed )
{
    ErrorCode rval;
    double posi[3], posf[3], len = 0;

    // one by one, get the node and project it on the sphere, with a radius given
    // the center of the sphere is at 0,0,0
    for( auto nit = nodes.begin(); nit != nodes.end(); ++nit )
    {
        EntityHandle nd = *nit;

        if( !is_cartesian )
        {
            rval = mb->get_coords( &nd, 1, posi );MB_CHK_ERR( rval );
            const double lat = posi[1] * 3.14159265358979323846 / 180;
            const double lon = posi[0] * 3.14159265358979323846 / 180;
            posf[0]          = R * cos( lat ) * cos( lon );  // x coordinate
            posf[1]          = R * cos( lat ) * sin( lon );  // y
            posf[2]          = R * sin( lat );               // z
            // spherical_to_cart( posi[1], posi[0], R, posf );
            // dbgprint( nd << " lat=" << posi[1] << ", lon=" << posi[0] << "; Cartesian = [" << posf[0] << ", " << posf[1] << ", " << posf[2] << "]" );
        }
        else
        {
            rval = mb->get_coords( &nd, 1, posf );MB_CHK_ERR( rval );
        }

        len = std::sqrt( posf[0] * posf[0] + posf[1] * posf[1] + posf[2] * posf[2] );
        if( len < 1e-12 )
        {
            dbgprint( nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2]
                         << ": Failed with length == 0." );
            return MB_FAILURE;
        }

        // rescale to radius
        posf[0] *= R / len;
        posf[1] *= R / len;
        if( is_threed && !is_cartesian )
            posf[2] = posi[2];
        else
            posf[2] *= R / len;
        // if (is_threed)
        //     dbgprint( nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2]  );
        rval = mb->set_coords( &nd, 1, posf );MB_CHK_ERR( rval );
    }
    return MB_SUCCESS;
}

moab::ErrorCode ComputeTempestRemapWeights( moab::Interface* mbi,
                                            RemappingContext& context,
                                            EntityHandle src_set,
                                            EntityHandle tgt_set,
                                            std::string strMethod,
                                            bool ensureMonotonicity )
{
    // err = remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( err );
    // err = remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh );MB_CHK_ERR( err );
    CloneToTRMesh( mbi, context.meshInput, src_set );
    CloneToTRMesh( mbi, context.meshOutput, tgt_set );

    context.meshInput.ConstructEdgeMap();
    context.meshOutput.ConstructEdgeMap();

    // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
    std::cout << "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes\n";
    // err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );
    bool concaveMeshA = false, concaveMeshB = false, allowNoOverlap = true, verbose = false;
    int ierr =
        GenerateOverlapWithMeshes( context.meshInput, context.meshOutput, context.meshOverlap, "" /*outFilename*/,
                                   "Netcdf4", "exact", concaveMeshA, concaveMeshB, allowNoOverlap, verbose );
    if( ierr )
    {
        MB_CHK_SET_ERR( MB_FAILURE, "TempestRemap: Can't compute the intersection of meshes on the sphere" );
    }

    dbgprint( "\nSetup computation of weights" );
    // Call to generate the remapping weights with the tempest meshes

    const std::string output_filename = "output_mpas_roms_map2d.nc";
    GenerateOfflineMapAlgorithmOptions mapOptions;
    mapOptions.nPin             = 1;
    mapOptions.nPout            = 1;
    mapOptions.fSourceConcave   = false;
    mapOptions.fTargetConcave   = false;
    mapOptions.strMethod        = strMethod;  // invdist, bilin, intbilin, delaunay
    mapOptions.fMonotone        = ensureMonotonicity;
    mapOptions.fNoCorrectAreas  = false;
    mapOptions.fNoCheck         = true;
    mapOptions.strOutputMapFile = output_filename;  // ask TR to write it out
    mapOptions.strOutputFormat  = "Netcdf4";

    dbgprint( "Compute weights with TempestRemap" );
    ierr = GenerateOfflineMapWithMeshes( context.meshInput,    // Mesh inputMesh
                                         context.meshOutput,   // Mesh outputMesh,
                                         context.meshOverlap,  // Mesh overlapMesh,
                                         "fv",                 // std::string inputDiscretization,
                                         "fv",                 // std::string outputDiscretization,
                                         mapOptions,           // const GenerateOfflineMapAlgorithmOptions& options
                                         context.weightMap );

    // check the generated weights and output information
    {
        const double dNormalTolerance = 1.0E-8;
        const double dStrictTolerance = 1.0E-12;
        context.weightMap.CheckMap( true, true, ensureMonotonicity, dNormalTolerance, dStrictTolerance );
    }

    // Write the map to disk
#ifdef WRITE_MAP_FILE
    {
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
        mapAttributes.insert( AttributePair( "nocorrectareas", "false" ) );
        mapAttributes.insert( AttributePair( "noconserve", "false" ) );
        mapAttributes.insert( AttributePair( "sparse_constraints", "false" ) );
        mapAttributes.insert( AttributePair( "method", mapOptions.strMethod ) );
        mapAttributes.insert( AttributePair( "version", "RemapMPASROMS v0.1" ) );

        dbgprint( "\nWrite the weights to " << output_filename );
        context.weightMap.Write( mapOptions.strOutputMapFile, mapAttributes, NcFile::Netcdf4Classic );

        // // Write the map file to disk in parallel using either HDF5 or SCRIP interface
        // err = weightMap.WriteParallelMap( output_filename.c_str() );MB_CHK_ERR( err );
    }
#endif
    return moab::MB_SUCCESS;
}

moab::ErrorCode ShepardInterpolator( int dimension,
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
                double t          = 0.0;
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

#define DelaunaySparseSerial   MOAB_FC_FUNC( delaunaysparses, DELAUNAYSPARSES )
#define DelaunaySparseParallel MOAB_FC_FUNC( delaunaysparsep, DELAUNAYSPARSEP )

// serial: optional arguments and compute interpolant values
extern "C" void DelaunaySparseSerial( int* d,
                                      int* n,
                                      double pts[],
                                      int* m,
                                      double q[],
                                      int simps[],
                                      double weights[],
                                      int ierr[],
                                      int* ir,
                                      double interp_in[],
                                      double interp_out[],
                                      double* eps,
                                      double* extrap,
                                      double rnorm[],
                                      int* ibudget,
                                      bool* chain,
                                      bool* exact );

// parallel: optional arguments and compute interpolant values
extern "C" void DelaunaySparseParallel( int* d,
                                        int* n,
                                        double pts[],
                                        int* m,
                                        double q[],
                                        int simps[],
                                        double weights[],
                                        int ierr[],
                                        int* ir,
                                        double interp_in[],
                                        double interp_out[],
                                        double* eps,
                                        double* extrap,
                                        double rnorm[],
                                        int* ibudget,
                                        bool* chain,
                                        bool* exact,
                                        int* pmode );

moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& /*xyzd*/,
                                            std::vector< double >& /*fd*/,
                                            std::vector< double >& /*xyzi*/,
                                            std::vector< double >& /*fi*/ )
{
    // int nd = fd.size();
    // int ni = fi.size();

    return moab::MB_SUCCESS;
}

moab::ErrorCode ComputeMBAInterpolant( std::vector< double >& xyzd,
                                       std::vector< double >& fd,
                                       std::vector< double >& xyzi,
                                       std::vector< double >& fi,
                                       bool is_threed,
                                       int order )
{
    const size_t nd = fd.size();
    const size_t ni = fi.size();

    if( false )
    {
        if( is_threed )
        {
            std::vector< double > xd( nd ), yd( nd ), zd( nd ), xi( ni ), yi( ni ), zi( ni );
#pragma omp parallel for
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                xd[k]               = xyzd[offset];
                yd[k]               = xyzd[offset + 1];
                zd[k]               = xyzd[offset + 2];
            }
#pragma omp parallel for
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                xi[k]               = xyzi[offset];
                yi[k]               = xyzi[offset + 1];
                zi[k]               = xyzi[offset + 2];
            }
            // mlinterp::interp( &nd, ni, fd, fi, xd, xi, yd, yi, zd, zi );

            return MB_SUCCESS;
        }
        else
        {
            std::vector< double > xd( nd ), yd( nd ), xi( ni ), yi( ni );
#pragma omp parallel for
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                xd[k]               = xyzd[offset];
                yd[k]               = xyzd[offset + 1];
            }
#pragma omp parallel for
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                xi[k]               = xyzi[offset];
                yi[k]               = xyzi[offset + 1];
            }
            // mlinterp::interp( &nd, ni, fd, fi, xd, xi, yd, yi );

            return MB_SUCCESS;
        }
    }

    constexpr bool use_recursive = false;
    std::vector< double > filocal;
    if( use_recursive ) filocal.resize( fi.size() );

    // Algorithm setup.
    if( order == 1 )
    {
        mba::linear_approximation< 3 >* interp;
        if( use_recursive )
        {
            ErrorCode err = ComputeNNInterpolant( xyzd, fd, xyzi, filocal );MB_CHK_ERR( err );

            std::vector< mba::point< 3 > > coords( ni );
#pragma omp parallel for
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                coords[k]           = mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] };
            }

            interp = new mba::linear_approximation< 3 >( coords.begin(), coords.end(), filocal.begin() );
        }
        else
        {
            std::vector< mba::point< 3 > > coords( nd );
#pragma omp parallel for
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                coords[k]           = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };
                // coords[k] = mba::point< 2 >{ xyzd[offset], xyzd[offset + 1] };
            }

            interp = new mba::linear_approximation< 3 >( coords.begin(), coords.end(), fd.begin() );
        }

        // Get interpolated value at arbitrary location.
#pragma omp parallel for
        for( size_t k = 0; k < ni; k++ )
        {
            const size_t offset = k * 3;
            fi[k]               = ( *interp )( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );
        }
        // fi[k] = interp( mba::point< 2 >{ xyzi[offset], xyzi[offset + 1] } );

        delete interp;
    }
    else
    {
        std::vector< mba::point< 3 > > coords( nd );
#pragma omp parallel for shared( coords )
        for( size_t k = 0; k < nd; k++ )
        {
            const size_t offset = k * 3;
            coords[k]           = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };
        }

        // construct a kd-tree index:
        using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                            PC3D< double >, 3 /* dim */
                                                            >;

        PC3D< double > cloud_src( xyzd ), cloud_tgt( xyzi );
        KdTree tree_src( 3 /*dim*/, cloud_src, { 5 /* max leaf */ } );
        KdTree tree_tgt( 3 /*dim*/, cloud_tgt, { 5 /* max leaf */ } );
        KdTree::BoundingBox bbox_src, bbox_tgt;
        tree_src.computeBoundingBox( bbox_src );
        tree_tgt.computeBoundingBox( bbox_tgt );

        printf( "Found bounding boxes: (%f, %f), (%f, %f), (%f, %f)\n", bbox_src[0].low, bbox_src[0].high,
                bbox_src[1].low, bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

        printf( "Found bounding boxes: (%f, %f), (%f, %f), (%f, %f)\n", bbox_tgt[0].low, bbox_tgt[0].high,
                bbox_tgt[1].low, bbox_tgt[1].high, bbox_tgt[2].low, bbox_tgt[2].high );

        int nlevels = 10;

        // Bounding box containing the data points.
        // mba::point< 3 > lo = { bbox_src[0].low * 0.9, bbox_src[1].low * 1.1, bbox_src[2].low * 1.1 };
        // mba::point< 3 > hi = { bbox_src[0].high, bbox_src[1].high, bbox_src[2].high };
        mba::point< 3 > lo = { -1.0, -1.0, -1E5 };
        mba::point< 3 > hi = { 1.0, 1.0, 1000.0 };

        // Initial grid size.
        // const size_t init_grid_size = static_cast< size_t >( std::max( 10.0, std::sqrt( nd ) / 8 ) );
        // mba::index< 3 > grid        = { init_grid_size, init_grid_size, 2 };
        mba::index< 3 > grid = { 13, 13, mpas_zlevels / 2 };

        // if( is_threed )
        // {
        //     grid[2] = mpas_zlevels;
        //     nlevels = 3;
        // }

        mba::MBA< 3 >* interp = nullptr;
        if( use_recursive )
        {
            ErrorCode err = ComputeNNInterpolant( xyzd, fd, xyzi, filocal );MB_CHK_ERR( err );

            const auto& tmpTree = tree_src;

            std::function< double( mba::point< 3 > ) > initFn = [&tmpTree, fd]( mba::point< 3 > query_pt ) {
                size_t srcindx;
                double srcdist;
                nanoflann::KNNResultSet< double > resultSet( 1 );

                // Do a KNN search
                resultSet.init( &srcindx, &srcdist );
                // printf( "Querying point: %f %f %f\n", query_pt[0], query_pt[1], query_pt[2] );
                tmpTree.findNeighbors( resultSet, query_pt.data() );

                return fd[srcindx];
            };

            interp = new mba::MBA< 3 >( lo, hi, grid, coords, fd, nlevels /*levels*/, 1e-8 /*tolerance*/,
                                        0.25 /*min_fill*/, initFn );
        }
        else
        {
            interp = new mba::MBA< 3 >( lo, hi, grid, coords, fd, nlevels /*levels*/, 1e-8 /*tolerance*/,
                                        0.25 /*min_fill*/ );
        }

        // Get interpolated value at arbitrary location.
        std::cout << "\nEvaluating the interpolant now...\n";
#pragma omp parallel for shared( fi, interp )
        for( size_t k = 0; k < ni; k++ )
        {
            const size_t offset = k * 3;
            fi[k]               = ( *interp )( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );
        }

        delete interp;
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ModifiedShepardInterpolator( int dimension,
                                             std::vector< double >& xyzd,
                                             std::vector< double >& fd,
                                             std::vector< double >& xyzi,
                                             std::vector< double >& fi,
                                             int order )
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
    const double shepard_power = order + 1;
    // First call the regular inverse-distance weighting method
    ShepardInterpolator( dimension, xyzd, fd, shepard_power, xyzi, fi );

    size_t nd = xyzd.size() / dimension;
    size_t ni = xyzi.size() / dimension;
    std::vector< double > w( nd, 0.0 );

    double fisecnum = 0.0;
    for( size_t k = 0; k < nd; k++ )
        fisecnum += fd[k];

    assert( shepard_power >= 1.0 );

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
                    w[j] = 1.0 / std::pow( w[j], shepard_power );
                    s += w[j];
                    is += 1.0 / w[j];
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
    const double shepard_power = 2;
    // size_t nd  = xyzd.size() / dimension;
    size_t ni = xyzi.size() / dimension;

    // construct a kd-tree index:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                        PC3D< double >, 3 /* dim */
                                                        >;

    PC3D< double > cloud( xyzd );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    double query_pt[3];  // dimension

    const size_t num_results =
        static_cast< size_t >( ( 2 * shepard_power + 1 ) * ( 2 * shepard_power + 1 ) - shepard_power );
    double ind = ( 1.0 / num_results );
    std::vector< size_t > srcindx( num_results );
    std::vector< double > srcdist( num_results );
    for( size_t i = 0; i < ni; i++ )
    {
        const int ioffset = i * dimension;
        {
            for( int dd = 0; dd < dimension; dd++ )
                query_pt[dd] = xyzi[ioffset + dd];
            // IntxUtils::transform_coordinates( query_pt, PointCloud< double >::projection );

            // Do a KNN search
            nanoflann::KNNResultSet< double > resultSet( num_results );
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
                dist = std::sqrt( dist );
                std::cout << i << "\tret_index=" << srcindx[ll] << " out_dist_sqr=" << srcdist[ll] << ", " << dist
                          << std::endl;
                srcdist[ll] = dist;
            }
        }

        std::vector< double > w( num_results, 0.0 );
        int z;
        if( shepard_power < 1 )
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
                w[j] = srcdist[j];  // std::sqrt( t );
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
                    w[j] = 1.0 / std::pow( w[j], shepard_power );
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

moab::ErrorCode ComputeFieldProjections( moab::Interface* mbi,
                                         RemappingContext& context,
                                         std::string varProjectSrc,
                                         std::string varProjectDst,
                                         std::vector< moab::EntityHandle >& srcelems,
                                         std::vector< moab::EntityHandle >& dstelems,
                                         bool is_three_dimensional,
                                         bool is_three2x1_dimensional,
                                         bool normalize,
                                         const double constantoffset,
                                         const std::string strMethod,
                                         int order,
                                         std::vector< moab::EntityHandle >* src3delems,
                                         std::vector< moab::EntityHandle >* dst3delems )
{
    moab::ErrorCode err;
    moab::Tag dmtag;
    if( is_three_dimensional ) assert( src3delems && dst3delems );
    const std::vector< moab::EntityHandle >& source_range = is_three_dimensional ? *src3delems : srcelems;
    const std::vector< moab::EntityHandle >& target_range = is_three_dimensional ? *dst3delems : dstelems;

    // err = mbi->tag_get_handle( varProject.c_str(), src_zlayers, moab::MB_TYPE_DOUBLE, dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
    err = mbi->tag_get_handle( varProjectSrc.c_str(), 1, moab::MB_TYPE_DOUBLE, dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

    // get the source data from tag
    std::vector< double > src_tdata( source_range.size() ), dst_tdata( target_range.size() );
    err = mbi->tag_get_data( dmtag, source_range.data(), source_range.size(), src_tdata.data() );MB_CHK_ERR( err );

    // get the coordinates of the elements
    std::vector< double > src_xyz( source_range.size() * 3 ), dst_xyz( target_range.size() * 3 );
    {
        err = mbi->get_coords( source_range.data(), source_range.size(), src_xyz.data() );MB_CHK_ERR( err );
        err = mbi->get_coords( target_range.data(), target_range.size(), dst_xyz.data() );MB_CHK_ERR( err );
    }

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0, normFactor = 1.0;
    if( normalize && ( !is_three_dimensional || is_three2x1_dimensional ) && context.meshOverlap.faces.size() )
    {
        for( size_t j = 0; j < src_tdata.size(); j++ )
            src_tdata[j] -= constantoffset;

        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < context.meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = context.meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += src_tdata[iSourceFace] * context.meshOverlap.vecFaceArea[j];
        }
    }

    if( !strMethod.compare( "" ) || !strMethod.compare( "bilin" ) || !strMethod.compare( "intbilin" ) ||
        !strMethod.compare( "delaunay" ) )
    {
        // get the handle to the weight matrix
        const SparseMatrix< double >& weights = context.weightMap.GetSparseMatrix();
        DataArray1D< double > dataInDouble( source_range.size(), false );
        dataInDouble.AttachToData( src_tdata.data() );
        DataArray1D< double > dataOutDouble( target_range.size(), false );
        dataOutDouble.AttachToData( dst_tdata.data() );

        // Compute the projection for the bottomDepth field
        weights.Apply( dataInDouble, dataOutDouble );
    }
    else if( !strMethod.compare( "nn" ) )
    {
        std::cout << "\nComputing Nearest-neighbor interpolant (order=1, degree=0) for field " << varProjectSrc
                  << std::endl;
        err = ComputeNNInterpolant( src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }
    else if( !strMethod.compare( "delaunay" ) )
    {
        std::cout << "\nComputing Delaunay Piecewise-Linear interpolant (order=2, degree=1) for field " << varProjectSrc
                  << std::endl;
        err = ComputeDelaunayInterpolant( src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }
    else if( !strMethod.compare( "mba" ) )
    {
        std::cout << "\nComputing MBA interpolant (order=4, degree=3) for field " << varProjectSrc << std::endl;
        err = ComputeMBAInterpolant( src_xyz, src_tdata, dst_xyz, dst_tdata, is_three_dimensional, order );MB_CHK_ERR( err );
    }
    else if( !strMethod.compare( "shepard" ) )
    {
        std::cout << "\nComputing Shepard interpolant (order=" << order << ") for field " << varProjectSrc << std::endl;
        err = ModifiedShepardInterpolator( 3, src_xyz, src_tdata, dst_xyz, dst_tdata, order );MB_CHK_ERR( err );
    }

    if( normalize && ( !is_three_dimensional || is_three2x1_dimensional ) && context.meshOverlap.faces.size() )
    {
        // Loop through all overlap-target faces and compute integral
        for( size_t j = 0; j < context.meshOverlap.faces.size(); j++ )
        {
            int iTargetFace = context.meshOverlap.vecTargetFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iTargetFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralOut += dst_tdata[iTargetFace] * context.meshOverlap.vecFaceArea[j];
        }

        normFactor = dTotalFieldIntegralIn / dTotalFieldIntegralOut;
        for( size_t ind = 0; ind < dst_tdata.size(); ind++ )
        {
            dst_tdata[ind] += constantoffset;
            dst_tdata[ind] *= normFactor;
        }
    }

    // now set the data on ROMS instance of MOAB tag
    moab::Tag drtag;
    err = mbi->tag_get_handle( varProjectDst.c_str(), 1, moab::MB_TYPE_DOUBLE, drtag,
                               moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( err );

    err = mbi->tag_set_data( drtag, target_range.data(), target_range.size(), dst_tdata.data() );MB_CHK_ERR( err );

    return moab::MB_SUCCESS;
}

ErrorCode CloneToTRMesh( moab::Interface* m_interface, Mesh& mesh, EntityHandle mesh_set )
{
    ErrorCode rval;
    moab::Range elems, verts;

    NodeVector& nodes = mesh.nodes;
    FaceVector& faces = mesh.faces;

    rval = m_interface->get_entities_by_dimension( mesh_set, 2, elems );MB_CHK_ERR( rval );

    // resize the number of elements in Tempest mesh
    faces.resize( elems.size() );

    // let us now get the vertices from all the elements
    rval = m_interface->get_connectivity( elems, verts );MB_CHK_ERR( rval );
    if( verts.size() == 0 )
    {
        rval = m_interface->get_entities_by_dimension( mesh_set, 0, verts );MB_CHK_ERR( rval );
    }
    // assert(verts.size() > 0); // If not, this may be an invalid mesh ! possible for unbalanced
    // loads

    std::map< EntityHandle, int > indxMap;
    {
        int j = 0;
        for( Range::iterator it = verts.begin(); it != verts.end(); it++ )
            indxMap[*it] = j++;
    }

    for( size_t iface = 0; iface < elems.size(); ++iface )
    {
        Face& face           = faces[iface];
        EntityHandle ehandle = elems[iface];

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity( ehandle, connectface, nnodesf );MB_CHK_ERR( rval );

        face.edges.resize( nnodesf );
        for( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            int indx = indxMap[connectface[iverts]];
            assert( indx >= 0 );
            face.SetNode( iverts, indx );
        }
    }

    size_t nnodes = verts.size();
    nodes.resize( nnodes );

    // Set the data for the vertices
    std::vector< double > coordx( nnodes ), coordy( nnodes ), coordz( nnodes );
    rval = m_interface->get_coords( verts, &coordx[0], &coordy[0], &coordz[0] );MB_CHK_ERR( rval );
    for( size_t inode = 0; inode < nnodes; ++inode )
    {
        Node& node = nodes[inode];
        node.x     = coordx[inode];
        node.y     = coordy[inode];
        node.z     = coordz[inode];
    }
    coordx.clear();
    coordy.clear();
    coordz.clear();

    mesh.RemoveZeroEdges();
    mesh.RemoveCoincidentNodes();

    // Generate reverse node array and edge map
    mesh.ConstructEdgeMap( false );
    // mesh.ConstructReverseNodeArray();

    // mesh.Validate();

    return MB_SUCCESS;
}

moab::ErrorCode ExtrudePolygonsToPolyhedra( Interface* mb,
                                            std::vector< double >& layer_thickness,
                                            moab::EntityHandle& poly2dset,
                                            moab::EntityHandle& outputset,
                                            const bool is_mpas )
{
    ErrorCode rval;
    if( outputset == 0 )
    {
        rval = mb->create_meshset( moab::MESHSET_SET, outputset );MB_CHK_SET_ERR( rval, "Can't create new set" );
    }

    // Get verts entities, by type
    Range verts, edges, faces;
    // rval = mb->get_entities_by_type( poly2dset, MBVERTEX, verts );MB_CHK_ERR( rval );
    rval = mb->get_entities_by_dimension( poly2dset, 0, verts );MB_CHK_ERR( rval );

    // Get faces, by dimension, so we stay generic to entity type
    rval = mb->get_entities_by_dimension( poly2dset, 2, faces );MB_CHK_ERR( rval );
    // std::cout << "Number of 2D entities: vertices = " << verts.size() << ", faces = " << faces.size() << endl;

    // add the initial faces to the first set
    if( is_mpas )
    {
        rval = mb->add_entities( outputset, verts );MB_CHK_ERR( rval );
        rval = mb->add_entities( outputset, faces );MB_CHK_ERR( rval );
    }

    // Create all edges
    rval = mb->get_adjacencies( faces, 1, true, edges, Interface::UNION );MB_CHK_ERR( rval );

    const size_t nverts = verts.size();
    const size_t nedges = edges.size();
    const size_t nfaces = faces.size();
    const int nlayers   = static_cast< int >( layer_thickness.size() / nfaces );
    const size_t nquads = nedges * nlayers;
    std::vector< double > coords( 3 * nverts );

    // output some information
    {
        dbgprint( " Input 2D " << ( is_mpas ? "MPAS" : "ROMS" ) << " Mesh details ::" );
        dbgprint( "\tNumber of Vertices = " << nverts );
        dbgprint( "\t          Edges    = " << edges.size() );
        dbgprint( "\t          Faces    = " << nfaces );
    }

    // get the vertex coordinates for the polygonal mesh
    rval = mb->get_coords( verts, &coords[0] );MB_CHK_ERR( rval );

    Tag gidTag = mb->globalId_tag();

    Tag parentTag;
    rval = mb->tag_get_handle( "ColumnParent", 1, moab::MB_TYPE_INTEGER, parentTag,
                               moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( rval );

    std::vector< int > gidData( nverts ), gidParentVertexData( nverts ), gidParentFaceData( nfaces );
    rval = mb->tag_get_data( gidTag, verts, gidParentVertexData.data() );MB_CHK_ERR( rval );
    rval = mb->tag_get_data( gidTag, faces, gidParentFaceData.data() );MB_CHK_ERR( rval );

    std::vector< int > layerchildren;
    // create first vertices
    Range* newVerts = new Range[nlayers + 1];
    newVerts[0]     = verts;  // just for convenience
    for( int ii = 0; ii < nlayers; ii++ )
    {
        for( size_t i = 0; i < nverts; i++ )
        {
            // coords[3 * i + 2] -= 0.5;
            // if( false )
            {
                Range eladjs;
                const EntityHandle vtx = verts[i];
                rval                   = mb->get_adjacencies( &vtx, 1, 2, false, eladjs, Interface::UNION );MB_CHK_ERR( rval );

                if( eladjs.size() )
                {
                    double thickness = 0.0;
                    double invweight = 0.0;
                    for( size_t k = 0; k < eladjs.size(); ++k )
                    {
                        int il = faces.index( eladjs[k] );
                        if( il < 0 ) continue;

                        invweight += 1.0;
                        thickness += layer_thickness[il * nlayers + ii];
                        if( thickness < 0 )
                        {
                            printf( "Thickness value for layer %d: element %zu  = %f\n", ii, k,
                                    layer_thickness[il * nlayers + ii] );
                            exit( 1 );
                        }
                    }
                    thickness /= invweight;

                    if( !is_mpas && thickness < 1e-10 )
                    {
                        printf( "Thickness value for layer %d: vertex %zu, adj = %zu = %f\n", ii, i, eladjs.size(),
                                thickness );
                        exit( 1 );
                    }

                    // Subtract or Add depending on the z-Direction to extrude
                    coords[3 * i + 2] -= thickness;
                }
                else
                    printf( "Vertex %zu has no adjacencies, coord = %f\n", i, coords[3 * i + 2] );
            }
        }

        rval = mb->create_vertices( &coords[0], nverts, newVerts[ii + 1] );MB_CHK_ERR( rval );
        std::iota( gidData.begin(), gidData.end(), nverts * ( 1 + ii ) );
        rval = mb->tag_set_data( gidTag, newVerts[ii + 1], gidData.data() );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( parentTag, newVerts[ii + 1], gidParentVertexData.data() );MB_CHK_ERR( rval );

        rval = mb->add_entities( outputset, newVerts[ii + 1] );MB_CHK_ERR( rval );
    }

    EntityHandle start_elem;
    std::vector< EntityHandle > allPolygons;
    if( is_mpas )
    {
        // for each edge, we will create nlayers quads
        ReadUtilIface* read_iface;
        rval = mb->query_interface( read_iface );MB_CHK_SET_ERR( rval, "Error in query_interface" );

        // Create quads
        EntityHandle* connect;
        rval = read_iface->get_element_connect( nquads, 4, MBQUAD, 0, start_elem, connect );MB_CHK_SET_ERR( rval, "Error in get_element_connect" );
        Range quads( start_elem, start_elem + nquads );

        // ---------------------------------------------------------------------------
        int indexConn = 0;
        for( size_t j = 0; j < nedges; j++ )
        {
            EntityHandle edge = edges[j];

            const EntityHandle* conn2 = NULL;
            int nnodes;
            rval = mb->get_connectivity( edge, conn2, nnodes );MB_CHK_ERR( rval );
            if( 2 != nnodes ) MB_CHK_ERR( MB_FAILURE );

            int i0 = verts.index( conn2[0] );
            int i1 = verts.index( conn2[1] );
            for( int ii = 0; ii < nlayers; ii++ )
            {
                connect[indexConn++] = newVerts[ii][i0];
                connect[indexConn++] = newVerts[ii][i1];
                connect[indexConn++] = newVerts[ii + 1][i1];
                connect[indexConn++] = newVerts[ii + 1][i0];
            }
        }

        // rval = mb->add_entities( outputset, quads );MB_CHK_ERR( rval );

        std::vector< int > allPolygonsGID( nquads );
        // allPolygonsGID.resize( nfaces * ( nlayers + 1 ) + nquads );

        // std::vector< int > allPolygonsGID;
        // GIDS for lateral quads will be at the end of the list
        // std::iota( allPolygonsGID.begin(), allPolygonsGID.end(), static_cast< int >( nfaces * ( nlayers + 1 ) ) + 1 );
        // TODO: this fails. Need to fix
        // rval = mb->tag_set_data( gidTag, quads, allPolygonsGID.data() );MB_CHK_ERR( rval );

        // next allocate for the x-y extruded faces
        allPolygons.resize( nfaces * ( nlayers + 1 ) );
        // rval = mb->tag_get_data( gidTag, faces, allPolygonsGID.data() );MB_CHK_ERR( rval );
        for( size_t i = 0; i < nfaces; i++ )
        {
            allPolygons[i] = faces[i];
        }
    }

    // vertices are parallel to the base vertices
    int gidElem                               = 1;
    int ipolygon                              = nfaces;
    int indexVerts[MAXEDGES]                  = { 0 };  // polygons with at most MAXEDGES edges
    EntityHandle polyhedronConn[MAXEDGES + 2] = { 0 };
    EntityHandle vertexConn[MAXEDGES * 2]     = { 0 };
    // edges will be used to determine the lateral faces of polyhedra (prisms)
    int indexEdges[MAXEDGES] = { 0 };  // index of edges in base polygon
    std::vector< int > minlevelFace, maxlevelFace;
    moab::Tag minlvlTag, maxlvlTag;
    moab::Tag mpas_soltags[nvars], mpas_soltags_new[nvars];
    std::vector< double > src_data( mpas_zreflevels * nvars );
    if( is_mpas )
    {
        rval = mb->tag_get_handle( "minLevelCell", 1, moab::MB_TYPE_INTEGER, minlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( rval );
        minlevelFace.resize( nfaces );

        rval = mb->tag_get_data( minlvlTag, faces, minlevelFace.data() );MB_CHK_ERR( rval );

        rval = mb->tag_get_handle( "maxLevelCell", 1, moab::MB_TYPE_INTEGER, maxlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( rval );
        maxlevelFace.resize( nfaces );
        rval = mb->tag_get_data( maxlvlTag, faces, maxlevelFace.data() );MB_CHK_ERR( rval );

        for( auto it = 0; it < nvars; ++it )
        {
            rval = mb->tag_get_handle( mpas_threed_cum_tagnames[it], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
                                       mpas_soltags[it], moab::MB_TAG_DENSE );MB_CHK_ERR( rval );

            rval = mb->tag_get_handle( mpas_threed_tagnames[it], 1, moab::MB_TYPE_DOUBLE, mpas_soltags_new[it],
                                       moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( rval );
        }
    }

    for( int ii = 0; ii < nlayers; ii++ )
        for( size_t j = 0; j < nfaces; j++ )
        {
            const EntityHandle polyg = faces[j];
            const EntityType etype   = mb->type_from_handle( polyg );
            const int polyGID        = gidParentFaceData[j];

            // printf( "Polygon %d has %d nodes\n", j, nnodes );

            const EntityHandle* connp = nullptr;
            int nnodes;
            rval = mb->get_connectivity( polyg, connp, nnodes );MB_CHK_ERR( rval );

            std::vector< int > vecents( nnodes + 2, polyGID );
            if( is_mpas )
            {
                const int orig_nodes = nnodes;

                // account for padded polygons
                while( connp[nnodes - 2] == connp[nnodes - 1] && nnodes > 3 )
                    nnodes--;

                std::copy( connp, connp + nnodes, vertexConn );
                // we had padded entities
                if( orig_nodes != nnodes )
                {
                    rval = mb->set_connectivity( polyg, vertexConn, nnodes );MB_CHK_ERR( rval );
                }

                for( int i = 0; i < nnodes; i++ )
                {
                    indexVerts[i]             = verts.index( connp[i] );
                    int i1                    = ( i + 1 ) % nnodes;
                    EntityHandle edgeVerts[2] = { connp[i], connp[i1] };
                    // get edge adjacent to these vertices
                    Range adjEdges;
                    rval = mb->get_adjacencies( edgeVerts, 2, 1, false, adjEdges );MB_CHK_ERR( rval );
                    if( adjEdges.size() < 1 ) MB_CHK_SET_ERR( MB_FAILURE, " did not find edge " );
                    indexEdges[i] = edges.index( adjEdges[0] );
                    if( indexEdges[i] < 0 ) MB_CHK_SET_ERR( MB_FAILURE, "did not find edge in range" );
                }

                for( auto it = 0; it < nvars; ++it )
                {
                    // get the source data from tag
                    rval = mb->tag_get_data( mpas_soltags[it], &polyg, 1, src_data.data() + it * src_zlayers );MB_CHK_ERR( rval );
                }
            }
            else
            {
                for( int i = 0; i < nnodes; i++ )
                {
                    // vertexConn[nnodes + i] = connp[i];
                    indexVerts[i] = verts.index( connp[i] );
                }
            }

            {
                // only add this extruded MPAS element if it is within the accepted layer mask
                if( is_mpas && ( ii + 1 < minlevelFace[j] || ii + 1 > maxlevelFace[j] ) ) continue;

                // create a polygon on each layer
                if( is_mpas )
                {
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[nnodes + i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1
                }
                else
                {
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[nnodes + i] = newVerts[ii][indexVerts[i]];  // vertices in layer ii+1
                }

                EntityHandle polyhedron;
                if( is_mpas )
                {
                    rval =
                        mb->create_element( etype, &vertexConn[nnodes], nnodes, allPolygons[nfaces * ( ii + 1 ) + j] );MB_CHK_ERR( rval );
                    // allPolygonsGID[nfaces * ( ii + 1 ) + j] = ipolygon++;

                    // now create a polyhedra with top, bottom and lateral swept faces
                    // first face is the bottom
                    polyhedronConn[0] = allPolygons[nfaces * ii + j];
                    // next add lateral quads, in order of edges, using the start_elem
                    // first layer of quads has EntityHandle from start_elem to start_elem+nedges-1
                    // second layer of quads has EntityHandle from start_elem + nedges to
                    // start_elem+2*nedges-1 ,etc
                    for( int i = 0; i < nnodes; i++ )
                    {
                        polyhedronConn[1 + i] = start_elem + ii + nlayers * indexEdges[i];
                    }
                    // second face is the top
                    polyhedronConn[1 + nnodes] = allPolygons[nfaces * ( ii + 1 ) + j];

                    // Create polyhedron
                    rval = mb->create_element( MBPOLYHEDRON, polyhedronConn, 2 + nnodes, polyhedron );MB_CHK_ERR( rval );

                    rval = mb->add_entities( outputset, polyhedronConn, nnodes + 2 );MB_CHK_ERR( rval );

                    rval = mb->tag_set_data( parentTag, polyhedronConn, nnodes + 2, vecents.data() );MB_CHK_ERR( rval );

                    for( auto it = 0; it < nvars; ++it )
                    {
                        // get the source data from tag
                        rval = mb->tag_set_data( mpas_soltags_new[it], &polyhedron, 1,
                                                 src_data.data() + it * src_zlayers + ii );MB_CHK_ERR( rval );
                    }

                    rval = mb->tag_set_data( gidTag, &allPolygons[nfaces * ( ii + 1 ) + j], 1, &ipolygon );MB_CHK_ERR( rval );
                    ipolygon++;
                }
                else
                {
                    EntityType extrudedType;
                    switch( etype )
                    {
                        case MBTRI:
                            extrudedType = MBPRISM;
                            break;
                        case MBQUAD:
                            extrudedType = MBHEX;
                            break;
                        default:
                            MB_CHK_SET_ERR( MB_FAILURE, "Unsupported standard 2D element type found for extrusion." );
                    }
                    rval = mb->create_element( extrudedType, vertexConn, 2 * nnodes, polyhedron );MB_CHK_ERR( rval );
                }

                rval = mb->add_entities( outputset, &polyhedron, 1 );MB_CHK_ERR( rval );
                rval = mb->tag_set_data( gidTag, &polyhedron, 1, &gidElem );MB_CHK_ERR( rval );
                gidElem++;

                rval = mb->tag_set_data( parentTag, &polyhedron, 1, &polyGID );MB_CHK_ERR( rval );
            }
        }

    // output some information
    {
        dbgprint( " Output 3D " << ( is_mpas ? "MPAS" : "ROMS" ) << " Mesh details ::" );
        dbgprint( "\tNumber of Vertices = " << nverts * ( nlayers + 1 ) );
        dbgprint( "\t          Elements = " << gidElem - 1 );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ComputeNNInterpolant( const std::vector< double >& src_xyz,
                                      const std::vector< double >& src_tdata,
                                      const std::vector< double >& dst_xyz,
                                      std::vector< double >& dst_tdata )
{
    constexpr double power       = 2.0;
    constexpr size_t num_results = 1;
    // construct a kd-tree index:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                        PC3D< double >, 3 /* dim */
                                                        >;

    double query_pt[3];  // dimension
    PC3D< double > cloud( src_xyz );
    KdTree tree( 3 /*dim*/, cloud, { 15 /* max leaf */ } );

    size_t offset = 0;
    for( size_t i = 0; i < dst_tdata.size(); i++, offset += 3 )
    {
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        nanoflann::KNNResultSet< double > resultSet( num_results );

        query_pt[0] = dst_xyz[offset];
        query_pt[1] = dst_xyz[offset + 1];
        query_pt[2] = dst_xyz[offset + 2];

        // Do a KNN search
        resultSet.init( srcindx.data(), srcdist.data() );
        tree.findNeighbors( resultSet, query_pt );

        double value = 0.0, weights = 0.0;
        for( size_t j = 0; j < num_results; ++j )
        {
            value += src_tdata[srcindx[j]] / std::pow( srcdist[j], power );
            weights += 1.0 / std::pow( srcdist[j], power );
        }

        // dst_tdata[i] = src_tdata[srcindx[0]];
        dst_tdata[i] = value / weights;
    }

    return moab::MB_SUCCESS;
}

// double moab::TempestOnlineMap::ApplyCAASLimiting( std::vector< double >& dataInDouble,
//                                                   std::vector< double >& dataOutDouble,
//                                                   bool useCAASLocal );
double ApplyCAASLimiting( OfflineMap& mapOperator,
                          Mesh& meshInput,
                          Mesh& meshOverlap,
                          const int nPin,
                          DataArray1D< double >& dataInDouble,
                          DataArray1D< double >& dataOutDouble,
                          bool useCAASLocal )
{
    const size_t nSourceCount                   = dataInDouble.GetRows();
    const size_t nTargetCount                   = dataOutDouble.GetRows();
    const DataArray1D< double >& m_dSourceAreas = mapOperator.GetSourceAreas();
    const DataArray1D< double >& m_dTargetAreas = mapOperator.GetTargetAreas();

    // Announce input mass
    double dSourceMass = 0.0;
    double dSourceMin  = dataInDouble[0];
    double dSourceMax  = dataInDouble[0];
    for( size_t i = 0; i < nSourceCount; i++ )
    {
        dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
        dSourceMax = fmax( dSourceMax, dataInDouble[i] );
        dSourceMin = fmin( dSourceMin, dataInDouble[i] );
    }

    // Apply the offline map to the data
    {
        DataArray1D< double > x( nTargetCount );
        DataArray1D< double > dataLowerBound( nTargetCount );
        DataArray1D< double > dataUpperBound( nTargetCount );

        double dMassDiff = dSourceMass;

        double dTargetMin = dataOutDouble[0];
        double dTargetMax = dataOutDouble[0];
        for( size_t i = 0; i < nTargetCount; i++ )
        {
            dMassDiff -= dataOutDouble[i] * m_dTargetAreas[i];
            dTargetMax = fmax( dTargetMax, dataOutDouble[i] );
            dTargetMin = fmin( dTargetMin, dataOutDouble[i] );
        }

        // Early exit if the values are monotone already.
        if( dTargetMax <= dSourceMax && dTargetMin <= dSourceMin ) return 0.0;

        if( useCAASLocal )
        {
            double dMinI;
            double dMaxI;

            int nTargetFaces = nTargetCount;
            std::vector< double > vecLocalUpperBound( nTargetCount );
            std::vector< double > vecLocalLowerBound( nTargetCount );

            std::vector< std::vector< int > > vecSourceOvTarget( nTargetFaces );
            for( size_t i = 0; i < meshOverlap.faces.size(); i++ )
            {

                int ixT = meshOverlap.vecTargetFaceIx[i];
                int ixS = meshOverlap.vecSourceFaceIx[i];
                vecSourceOvTarget[ixT].push_back( ixS );
            }

            //FV to FV
            {
                for( size_t i = 0; i < nTargetCount; i++ )
                {
                    if( !vecSourceOvTarget[i].size() ) continue;
                    dMaxI = dataInDouble[vecSourceOvTarget[i][0]];
                    dMinI = dataInDouble[vecSourceOvTarget[i][0]];

                    //Compute max over interstecting source faces

                    for( size_t j = 0; j < vecSourceOvTarget[i].size(); j++ )
                    {
                        int k = vecSourceOvTarget[i][j];
                        dMaxI = fmax( dMaxI, dataInDouble[k] );
                        dMinI = fmin( dMinI, dataInDouble[k] );
                    }

                    if( useCAASLocal )
                    {
                        double dMaxIAdj = dMaxI;
                        double dMinIAdj = dMinI;

                        AdjacentFaceVector vecAdjFaces;

                        GetAdjacentFaceVectorByEdge( meshInput, vecSourceOvTarget[i][0], ( nPin + 1 ) * ( nPin + 1 ),
                                                     vecAdjFaces );

                        //Compute max over neighboring faces
                        for( size_t j = 0; j < vecAdjFaces.size(); j++ )
                        {
                            int k = vecAdjFaces[j].first;

                            dMaxIAdj = fmax( dMaxIAdj, dataInDouble[k] );
                            dMinIAdj = fmin( dMinIAdj, dataInDouble[k] );
                        }

                        vecLocalLowerBound[i] = dMinIAdj;
                        vecLocalUpperBound[i] = dMaxIAdj;
                    }
                    else
                    {
                        vecLocalLowerBound[i] = dMinI;
                        vecLocalUpperBound[i] = dMaxI;
                    }
                }
            }

            for( size_t i = 0; i < dataLowerBound.GetRows(); i++ )
            {
                dataLowerBound[i] = vecLocalLowerBound[i] - dataOutDouble[i];
                dataUpperBound[i] = vecLocalUpperBound[i] - dataOutDouble[i];
            }

        }     // if( useCAASLocal )
        else  // useCAASGlobal
        {
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = dSourceMin - dataLowerBound[i];
                dataUpperBound[i] = dSourceMax - dataUpperBound[i];
            }
        }

        // Invoke CAAS application on the offline map
        mapOperator.CAAS( dataOutDouble, dataLowerBound, dataUpperBound, dMassDiff );
    }

    // Announce output mass
    double dTargetMass = 0.0;
    double dTargetMin  = dataOutDouble[0];
    double dTargetMax  = dataOutDouble[0];
    for( size_t i = 0; i < nTargetCount; i++ )
    {
        dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
        if( dataOutDouble[i] < dTargetMin )
        {
            dTargetMin = dataOutDouble[i];
        }
        if( dataOutDouble[i] > dTargetMax )
        {
            dTargetMax = dataOutDouble[i];
        }
    }

    return ( dTargetMass - dSourceMass );
}

double ApplyCAASLimiting_ABC( OfflineMap& mapOperator,
                              Mesh& meshInput,
                              Mesh& meshOverlap,
                              const int nPin,
                              DataArray1D< double >& dataInDouble,
                              DataArray1D< double >& dataOutDouble,
                              bool useCAAS,
                              bool useCAASLocal )
{

    const int nSourceCount                      = dataInDouble.GetRows();
    const int nTargetCount                      = dataOutDouble.GetRows();
    const DataArray1D< double >& m_dSourceAreas = mapOperator.GetSourceAreas();
    const DataArray1D< double >& m_dTargetAreas = mapOperator.GetTargetAreas();

    // Announce input mass
    double dSourceMass = 0.0;
    double dSourceMin  = dataInDouble[0];
    double dSourceMax  = dataInDouble[0];
    for( int i = 0; i < nSourceCount; i++ )
    {
        dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
        if( dataInDouble[i] < dSourceMin )
        {
            dSourceMin = dataInDouble[i];
        }
        if( dataInDouble[i] > dSourceMax )
        {
            dSourceMax = dataInDouble[i];
        }
    }

    // Apply the offline map to the data

    if( useCAASLocal || useCAAS )
    {
        DataArray1D< double > l = dataOutDouble;
        DataArray1D< double > u = dataOutDouble;
        DataArray1D< double > x( nTargetCount );
        double b = dSourceMass;

        for( size_t i = 0; i < l.GetRows(); i++ )
        {
            b -= dataOutDouble[i] * m_dTargetAreas[i];
        }

        if( useCAASLocal )
        {
            // int GLLSizeIn  = 0;  // FV
            // int GLLSizeOut = 0;  // FV
            // int pOut       = dataGLLNodesOut.GetSize( 0 );
            // int qOut       = dataGLLNodesOut.GetSize( 1 );
            // int pIn        = dataGLLNodesIn.GetSize( 0 );
            // int qIn        = dataGLLNodesIn.GetSize( 1 );
            // int pIn          = nPin;
            double f_maxI    = 0.0;
            double f_minI    = 0.0;
            int nTargetFaces = nTargetCount;

            std::vector< std::vector< int > > SourceOvTarget( nTargetFaces );

            for( size_t i = 0; i < meshOverlap.faces.size(); i++ )
            {
                int ixT = meshOverlap.vecTargetFaceIx[i];
                int ixS = meshOverlap.vecSourceFaceIx[i];
                SourceOvTarget[ixT].push_back( ixS );
            }

            std::vector< double > local_UB( nTargetCount );
            std::vector< double > local_LB( nTargetCount );

            for( int i = 0; i < nTargetCount; i++ )
            {
                AdjacentFaceVector vecAdjFaces;

                if( !SourceOvTarget[i].size() ) continue;

                GetAdjacentFaceVectorByEdge( meshInput, SourceOvTarget[i][0], ( nPin + 1 ) * ( nPin + 1 ),
                                             vecAdjFaces );
                f_maxI = dataInDouble[vecAdjFaces[0].first];
                f_minI = dataInDouble[vecAdjFaces[0].first];
                for( size_t j = 0; j < vecAdjFaces.size(); j++ )
                {
                    int k  = vecAdjFaces[j].first;
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                for( size_t j = 0; j < SourceOvTarget[i].size(); j++ )
                {

                    int k  = SourceOvTarget[i][j];
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                // f_minI=fmax(f_minI,0.0);

                local_UB[i] = f_maxI;
                local_LB[i] = f_minI;
            }

            double mt = 0.0;
            for( int i = 0; i < nTargetCount; i++ )
            {
                mt += m_dTargetAreas[i] * ( local_LB[i] - dataOutDouble[i] );
            }

            for( size_t i = 0; i < l.GetRows(); i++ )
            {
                l[i] = local_LB[i] - l[i];
                u[i] = local_UB[i] - u[i];
            }

            // Adjust mass of lower bound if greater than b
            double mL = 0.0;

            for( int i = 0; i < nTargetCount; i++ )
            {
                mL += m_dTargetAreas[i] * l[i];
            }
            if( mL > b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mL   = mL - m_dTargetAreas[i] * l[i] + m_dTargetAreas[i] * ( dSourceMin - dataOutDouble[i] );
                    l[i] = dSourceMin - dataOutDouble[i];
                    if( mL < b )
                    {
                        break;
                    }
                }
            }

            // Adjust mass of upper bound if less than b
            double mU = 0.0;

            if( mU < b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mU   = mU - m_dTargetAreas[i] * u[i] + m_dTargetAreas[i] * ( dSourceMax - dataOutDouble[i] );
                    u[i] = dSourceMax - dataOutDouble[i];
                    if( mU > b )
                    {
                        break;
                    }
                }
            }
        }  // if( useCAASLocal )
        else if( useCAAS )
        {
            for( size_t i = 0; i < l.GetRows(); i++ )
            {
                l[i] = dSourceMin - l[i];
                u[i] = dSourceMax - u[i];
            }
        }

        // Invoke CAAS application on the offline map
        mapOperator.CAAS( x, l, u, b );

        // Add correction
        for( size_t i = 0; i < l.GetRows(); i++ )
        {
            dataOutDouble[i] += x[i];
        }
    }

    // Announce output mass
    double dTargetMass = 0.0;
    double dTargetMin  = dataOutDouble[0];
    double dTargetMax  = dataOutDouble[0];
    for( int i = 0; i < nTargetCount; i++ )
    {
        dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
        if( dataOutDouble[i] < dTargetMin )
        {
            dTargetMin = dataOutDouble[i];
        }
        if( dataOutDouble[i] > dTargetMax )
        {
            dTargetMax = dataOutDouble[i];
        }
    }

    return ( dTargetMass - dSourceMass );
}