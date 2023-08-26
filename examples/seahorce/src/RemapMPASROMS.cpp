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
#include <numeric>  // std::iota

#include "moab/MOABConfig.h"
// #undef MOAB_HAVE_MPI

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/ProgOptions.hpp"

// Other includes
#include "RemapMPASROMS.hpp"
#include "ComputeNN.hpp"
#include "ComputeShepard.hpp"
#include "PCHIP.hpp"
#include "spline.h"
#include "ComputeTR.hpp"
#include "ComputeMBA.hpp"
#include "moab/Remapping/mlinterp.hpp"
#include "MeshUtilities.hpp"

using namespace moab;
using namespace std;

moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& xyzd,
                                            std::vector< double >& fd,
                                            std::vector< double >& xyzi,
                                            std::vector< double >& fi );

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

//
// Start of main test program
//
int main( int argc, char** argv )
{
    ErrorCode err;

    RemappingContext context;
    int ierr, rank, size;
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
    int bathymetryOrder          = 1;
    int fieldOrder               = 3;

    int src_zlayers = 0, dst_zlayers = 0;

    constexpr bool useCAAS = false;

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

        // if( threetwooneD )
        // {
        //     strMethod = "bilin";
        //     computeTR = true;
        // }

        if( dimension == 3 ) use_3dprojection = true;
        if( !computeMBA && !computeShepard && !computeTR ) computeMBA = true;
    }

    // set the number of z-layers
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

    if( normalize || computeTR )
    {
        // call to compute the 2D map and store to disk
        err = ComputeTempestRemapWeights( mbi, context, mpas_covering_set, romsset, bathymetryMethod,
                                          ensureMonotonicity );MB_CHK_SET_ERR( err, "Cannot compute 2D remapping weights" );
    }
    else
    {
        // load the computed 2D map files
        err = LoadTempestRemapWeights( mbi, context, mpas_covering_set, romsset, bathymetryMethod );MB_CHK_SET_ERR( err, "Cannot load 2D remapping weights" );
    }

    // let us perform 3D extrusions as needed
    EntityHandle root_set  = 0;
    EntityHandle mpasset3d = 0, romsset3d = 0;

    {
        // Initialize all important data

        // MPAS reference-z-levels
        {
            // constexpr double mpas_ref_levels[] = {
            //     10,      20,      30,      40,      50,      60,      70,      80,      90,      100,
            //     110,     120,     130,     140,     150,     160,     170.197, 180.761, 191.821, 203.499,
            //     215.923, 229.233, 243.584, 259.156, 276.152, 294.815, 315.424, 338.312, 363.875, 392.58,
            //     424.989, 461.767, 503.707, 551.749, 606.997, 670.729, 744.398, 829.607, 928.043, 1041.37,
            //     1171.04, 1318.09, 1482.9,  1664.99, 1863.01, 2074.87, 2298.04, 2529.9,  2768.1,  3010.67,
            //     3256.14, 3503.45, 3751.89, 4001.01, 4250.53, 4500.26, 4750.12, 5000.05, 5250.01, 5499.99 };
            moab::Tag mztag;
            err = mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag,
                                       moab::MB_TAG_SPARSE );MB_CHK_SET_ERR( err, "Can't get tag handle: refBottomDepth" );

            err = mbi->tag_get_data( mztag, &root_set, 1, context.mpas_zref_heights );MB_CHK_SET_ERR( err, "Can't get refBottomDepth data" );
        }

        // Project the bottom Bathymetry data from MPAS to ROMS so that we can impose it.
        err = ComputeFieldProjections( mbi, context, "bottomDepth", "bottomDepth", mpas_elems, roms_elems,
                                    false /*use_3dprojection*/, false /* 2Dx1D */, false /* bool normalize */, 2000.0,
                                    bathymetryMethod, bathymetryOrder );MB_CHK_SET_ERR( err, "Can't create new set" );
    }

    std::vector< double > zmh_xyz3d, zrh_xyz3d;
    if( use_3dprojection || threetwooneD )
    {
        zmh_xyz3d.resize( mpas_elems.size() * src_zlayers );
        zrh_xyz3d.resize( roms_elems.size() * dst_zlayers );

        constexpr bool useConstantRefAxialThickness = true;
        if( useConstantRefAxialThickness )
        {
            for( size_t i = 0; i < mpas_elems.size(); ++i )
            {
                const int offset = i * src_zlayers;
                // zmh_xyz3d[offset] = 0;
                // for( int j = 1; j <= src_zlayers; ++j )
                //     zmh_xyz3d[offset + j] = context.mpas_zref_heights[j - 1] + zmh_xyz3d[offset + j - 1];
                // for( int j = 0; j < src_zlayers; ++j )
                //     zmh_xyz3d[offset + j] = context.mpas_zref_heights[j];

                zmh_xyz3d[offset] = context.mpas_zref_heights[0];
                for( int j = 1; j < src_zlayers; ++j )
                    zmh_xyz3d[offset + j] = context.mpas_zref_heights[j] - context.mpas_zref_heights[j - 1];
            }
        }
        else
        {
            moab::Tag mhtag;
            err = mbi->tag_get_handle( "layerThickness_3d", src_zlayers, moab::MB_TYPE_DOUBLE, mhtag,
                                       moab::MB_TAG_DENSE );
            MB_CHK_SET_ERR( err, "Can't get tag handle: layerThickness_3d" );
            err = mbi->tag_get_data( mhtag, mpas_elems.data(), mpas_elems.size(), zmh_xyz3d.data() );
            MB_CHK_SET_ERR( err, "Can't get layerThickness_3d data" );
        }

        if( generateExtrusions )
        {
            dbgprint( "\nExtruding MPAS polygonal mesh ..." );
            err = ExtrudePolygonsToPolyhedra( mbi, zmh_xyz3d, mpas_covering_set, mpasset3d, true, src_zlayers );MB_CHK_SET_ERR( err, "Can't extrude MPAS polygons" );

            {
                std::vector< double > zrh_xyz2d( roms_elems.size() );
                moab::Tag rhtag;
                err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
                err = mbi->tag_get_data( rhtag, roms_elems.data(), roms_elems.size(), zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

                err = mbi->write_file( "roms_3d_2dsurface.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );

#define USE_STRETCHING_FUNCTION

#ifdef USE_STRETCHING_FUNCTION
                /// Call stretching functions
                // h = np.linspace( 10, 200, 10 );
                constexpr double zeta = 0.0;

                // return N evenly s-coordinate w-points.
                auto Sw = [&]() {
                    std::vector< double > s( roms_zlevels+1 );
                    double del = 1.0 / roms_zlevels;
                    s[0] = 0.0;
                    for( auto i = 0; i < roms_zlevels; ++i )
                        s[i + 1] = s[i] - del;
                    return s;
                };

                /// reference: https://github.com/seahorce-scidac/seahorce-notebooks/blob/main/ROMS_Vertical_Grid_examples.ipynb
                auto vstretching_1 = [&]( double s ) {
                    constexpr double theta_s = 5.0;
                    constexpr double theta_b = 0.5;
                    return ( 1 - theta_b ) * ( sinh( s * theta_s ) / sinh( theta_s ) ) +
                           theta_b * ( -0.5 + 0.5 * tanh( theta_s * ( s + 0.5 ) ) / tanh( 0.5 * theta_s ) );
                };
                auto vstretching_2 = [&]( double s ) {
                    constexpr double theta_s = 5.0;
                    // constexpr double theta_b = 0.5;
                    return ( 1 - cosh( theta_s * s ) ) / ( cosh( theta_s ) - 1.0 );
                };
                auto vstretching_4 = [&]( double s ) {
                    constexpr double theta_s = 5.0;
                    constexpr double theta_b = 0.5;
                    double C                 = ( 1.0 - cosh( theta_s * s ) ) / ( cosh( theta_s ) - 1.0 );
                    return ( exp( theta_b * C ) - 1.0 ) / ( 1.0 - exp( -theta_b ) );
                };

                auto vtransform_1 = [&]( double s, double h ) {
                    constexpr double hc = 3150;  // hc has to be less than or equal to min(h)
                    double vstretch_val = vstretching_1( s );
                    return hc * ( s - vstretch_val ) + vstretch_val * h;
                };

                auto vtransform_2 = [&]( double s, double h ) {
                    constexpr double hc = 3000;  // hc has to be less than or equal to min(h)
                    double vstretch_val = vstretching_4( s );
                    return ( hc * s + vstretch_val * h ) / (hc + h);
                };
#endif
                constexpr int transform_id = 2; // 1 or 2

                for( size_t i = 0; i < roms_elems.size(); ++i )
                {
                    const double pbathymetry = zrh_xyz2d[i];
                    const double delz        = pbathymetry / dst_zlayers;
                    const int offset         = i * dst_zlayers;
#ifdef USE_STRETCHING_FUNCTION
                    std::vector< double > s_w = Sw();
                    // auto z_w = zeta + ( zeta + h ) * vtransform_1( sw );
                    {
                        // h = delz * i
                        std::vector< double > z_w( roms_zlevels + 1 );
                        double h = 0.0;
                        z_w[0]   = zeta + ( zeta + h ) * vtransform_1( s_w[0], h );
                        // if( i < 10 )
                        //     printf( "Element = %d, Level = %d, s_w = %f, zheight = %f\n", i, 0, s_w[0], z_w[0] );
                        for( auto ilevel = 0; ilevel < roms_zlevels; ++ilevel )
                        {
                            h = delz * ( ilevel + 1 );
                            if( transform_id == 1 )
                                z_w[ilevel + 1] = (vtransform_1( s_w[ilevel + 1], h ) +
                                                  zeta * ( 1 + vtransform_1( s_w[ilevel + 1], h ) / h ));
                            else
                                z_w[ilevel + 1] = (zeta + ( zeta + h ) * vtransform_2( s_w[ilevel + 1], h ));
                            // zrh_xyz3d[offset + ilevel] = 0.5 * ( z_w[ilevel] + z_w[ilevel + 1] );
                            zrh_xyz3d[offset + ilevel] = -(z_w[ilevel + 1] - z_w[ilevel]);
                            // if (i < 10) printf( "Element = %d, Level = %d, s_w = %f, zheight = %f, default = %f\n", i, ilevel, s_w[ilevel+1], zrh_xyz3d[offset + ilevel], delz );
                        }

                        // std::generate( z_w.begin(), z_w.end(), []() { return zeta + ( zeta + h ) * vtransform_1( sw ) } );
                    }
#else
                    for( int j = 0; j < dst_zlayers; ++j )
                    {
                        zrh_xyz3d[offset + j] = delz;
                        // printf( "Thickness value for ROMS element %zu, %d = %f, %f\n", i, dst_zlayers, zrh_xyz2d[i], delz );
                    }
#endif
                }

                dbgprint( "\nExtruding ROMS structured quad mesh..." );
                err = ExtrudePolygonsToPolyhedra( mbi, zrh_xyz3d, romsset, romsset3d, false, dst_zlayers );MB_CHK_SET_ERR( err, "Can't extrude ROMS polygons" );
            }

            err = mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems );MB_CHK_ERR( err );

            // std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size() << " elements.\n";
            err = mbi->write_file( "mpas_full_3d.h5m", "H5M", write_options.c_str(), &mpasset3d, 1 );MB_CHK_ERR( err );

            err = mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems );MB_CHK_ERR( err );

            {
                moab::Tag rh3tag;
                err = mbi->tag_get_handle( "ROMSlayerThickness", 1, moab::MB_TYPE_DOUBLE, rh3tag,
                                           moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_SET_ERR( err, "Can't get tag handle: ROMSlayerThickness" );
                err = mbi->tag_set_data( rh3tag, roms3d_elems.data(), roms3d_elems.size(), zrh_xyz3d.data() );MB_CHK_SET_ERR( err, "Can't get layerThickness_3d data" );
            }

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

            moab::Tag rh3tag;
            err = mbi->tag_get_handle( "ROMSlayerThickness", 1, moab::MB_TYPE_DOUBLE, rh3tag,
                                       moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_SET_ERR( err, "Can't get tag handle: ROMSlayerThickness" );

            err = mbi->tag_get_data( rh3tag, roms3d_elems.data(), roms3d_elems.size(), zrh_xyz3d.data() );MB_CHK_SET_ERR( err, "Can't get layerThickness_3d data" );

            // // get MPAS and ROMS height factors for elements
            // {
            //     moab::Tag mztag;
            //     err = mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag,
            //                                moab::MB_TAG_SPARSE );MB_CHK_SET_ERR( err, "Can't get tag handle: refBottomDepth" );

            //     err = mbi->tag_get_data( mztag, &root_set, 1, context.mpas_zref_heights.data() );MB_CHK_SET_ERR( err, "Can't get refBottomDepth data" );
            //     for( size_t i = 0; i < mpas_elems.size(); ++i )
            //     {
            //         const int offset = i * src_zlayers;
            //         zmh_xyz3d[offset] = context.mpas_zref_heights[0];
            //         for( int j = 1; j < src_zlayers; ++j )
            //             zmh_xyz3d[offset + j] = context.mpas_zref_heights[j] - context.mpas_zref_heights[j - 1];
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

    bool twoDfirst = !oneDfirst;

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

        moab::Tag rhtag;
        err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
        std::vector< double > zrh_xyz2d( romssize );
        err = mbi->tag_get_data( rhtag, roms_elems.data(), romssize, zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

        std::vector< int > minlevelFace, maxlevelFace;
        moab::Tag minlvlTag, maxlvlTag;
        err = mbi->tag_get_handle( "minLevelCell", 1, moab::MB_TYPE_INTEGER, minlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        minlevelFace.resize( mpassize );
        err = mbi->tag_get_data( minlvlTag, mpas_elems.data(), mpassize, minlevelFace.data() );MB_CHK_ERR( err );

        err = mbi->tag_get_handle( "maxLevelCell", 1, moab::MB_TYPE_INTEGER, maxlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        maxlevelFace.resize( mpassize );
        err = mbi->tag_get_data( maxlvlTag, mpas_elems.data(), mpassize, maxlevelFace.data() );MB_CHK_ERR( err );

        std::vector< double > zmh_z( mpas_zreflevels ), zmh_z_rev( mpas_zreflevels );
        zmh_z_rev[0] = 0.5*context.mpas_zref_heights[0];
        // double totalz = -context.mpas_zref_heights[0];
        // printf( "zmh_z_rev[%d] = %f and element z-coord = %f \n", 0, zmh_z_rev[0], context.mpas_zref_heights[0] );
        for( size_t j = 1; j < mpas_zlevels; ++j )
        {
            zmh_z_rev[j] = 0.5 * ( context.mpas_zref_heights[j - 1] + context.mpas_zref_heights[j] );
            // zmh_z[j] = ( totalz - 0.5 * context.mpas_zref_heights[j] );
            // printf( "zmh_z[%d] = %f and element z-coord = %f \n", j, zmh_z_rev[j], -context.mpas_zref_heights[j] );
        }
        // for( int j = mpas_zlevels - 1; j >= 0; --j )
        for( size_t j = 0; j < mpas_zlevels; ++j )
        {
            // zmh_z[j] = zmh_z_rev[mpas_zlevels - 1 - j];
            zmh_z[j] = zmh_z_rev[j];
            // printf( "zmh_z_rev[%d] = %f and element z-coord = %f \n", j, zmh_z_rev[j], zmh_z[j] );
        }

        // for( size_t j = 0; j < mpas_zlevels; ++j )
        // {
        //     double coords[3];
        //     {
        //         mbi->get_coords( &mpas3d_elems[j * mpassize], 1, coords );
        //         printf( "zmh_z[%d] = %f and element z-coord = %f\n", j, zmh_z[j], coords[2] );
        //     }
        // }

        std::vector< double >
            tgt_salinity_data( romssize * roms_zlevels ),
             tgt_temperature_data( romssize * roms_zlevels );

        if( twoDfirst )
        {
            // First compute the projections in 2D for each MPAS layer.
            // This will give the MPAS projected data on to ROMS mesh - on MPAS axial levels
            // NOTE: implicit assumption is that ROMS levels >= MPAS levels. Should fix how
            // the projection is invoked by perhaps skipping the tag_set_data

            std::vector< double > tgtsrc_salinity_data( romssize * mpas_zreflevels ),
                tgtsrc_temperature_data( romssize * mpas_zreflevels );

#pragma omp parallel for shared( tgtsrc_salinity_data, tgtsrc_temperature_data, weights )
            for( int ii = 0; ii < mpas_zlevels; ii++ )
            {
                DataArray1D< double > dataInDoubleS( mpassize ), dataInDoubleT( mpassize );
                unsigned offsetr = romssize * ii;

                dbgprint( "Computing projection for MPAS level: " << ii );

                DataArray1D< double > dataOutDoubleS( romssize, false ), dataOutDoubleT( romssize, false );

#ifndef VERTICAL_INTERPOLATION
                // std::vector< EntityHandle > mpas_slice( mpas3d_elems.begin() + offsetm,
                //                                         mpas3d_elems.begin() + offsetm + mpassize );
                std::vector< EntityHandle > roms_slice( roms3d_elems.begin() + offsetr,
                                                        roms3d_elems.begin() + offsetr + romssize );
#endif
                for( size_t j = 0; j < mpassize; ++j )
                {
                    const int offset = j * mpas_zlevels;
                    dataInDoubleS[j]    = maxlevelFace[j] - 1 < ii ? src_salinity_data[maxlevelFace[j] - 1 + offset]
                                                                   : src_salinity_data[ii + offset];
                    dataInDoubleT[j]    = maxlevelFace[j] - 1 < ii ? src_temperature_data[maxlevelFace[j] - 1 + offset]
                                                                   : src_temperature_data[ii + offset];
                    // if( ii == 0 || ii == mpas_zlevels - 1 )
                    // {
                    //     std::cout << "Level " << ii << ": Salinity = " << dataInDoubleS[j]
                    //               << ", Temperature = " << dataInDoubleT[j] << std::endl;
                    // }
                }

                dataOutDoubleS.AttachToData( tgtsrc_salinity_data.data() + offsetr );

                // Compute the projection for the salinity field
                weights.Apply( dataInDoubleS, dataOutDoubleS );

                if( useCAAS )
                    ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap, 1 /*nPin*/,
                                       dataInDoubleS, dataOutDoubleS, false /*useCAASLocal*/ );

#ifndef VERTICAL_INTERPOLATION
                err = mbi->tag_set_data( roms_soltags_elem[0], roms_slice.data(), roms_slice.size(), dataOutDoubleS );MB_CHK_SET_ERR( err, "Can't set salinity tag data" );
#endif
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
            }

#ifdef VERTICAL_INTERPOLATION
            // First MPAS z-levels
            // Now ROMS

            // Now project each data layer in axial direction.
            // NOTE: Embarassingly parallel
#pragma omp parallel for shared( zrh_xyz2d, tgtsrc_salinity_data, tgtsrc_temperature_data, tgt_salinity_data, \
                                     tgt_temperature_data )
            for( size_t i = 0; i < romssize; ++i )
            {
                std::vector< double > roms_zvalsS( mpas_zreflevels ), roms_zvalsT( mpas_zreflevels );
                // const double pbathymetry = zrh_xyz2d[i];
                // const double delz        = pbathymetry / roms_zlevels;
                // if (i<10) printf( "--- ROMS delz = %f\n", delz );

                for( int imzl = 0; imzl < mpas_zreflevels; ++imzl )
                {
                    const int offset               = imzl * romssize;
                    // roms_zvalsS[mpas_zreflevels - 1 - imzl] = tgtsrc_salinity_data[i + offset];
                    // roms_zvalsT[mpas_zreflevels - 1 - imzl] = tgtsrc_temperature_data[i + offset];
                    roms_zvalsS[imzl] = tgtsrc_salinity_data[i + offset];
                    roms_zvalsT[imzl] = tgtsrc_temperature_data[i + offset];
                }

#ifndef VERTICAL_INTERPOLANT_LINEAR
                tk::spline splS( zmh_z, roms_zvalsS, tk::spline::cspline_hermite, true );
                tk::spline splT( zmh_z, roms_zvalsT, tk::spline::cspline_hermite, true );
#endif

                double roms_ztotal = 0.0;
                for( int irzl = 0; irzl < roms_zlevels; ++irzl )
                {
                    const int offset   = irzl * romssize;
                    // double roms_zlocation = ( irzl + 0.5 ) * delz;
                    double roms_zlocation = roms_ztotal + zrh_xyz3d[irzl]/2;
                    roms_ztotal += zrh_xyz3d[irzl];

                    // Project data from zmh_z to zrh_z
#ifdef VERTICAL_INTERPOLANT_LINEAR
                    mlinterp::interp( &mpas_zlevels, 1,                              // Number of points
                                      roms_zvalsS.data(), &tgt_salinity_data[i + offset],  // Output axis (y)
                                      zmh_z.data(), &roms_zlocation                        // Input axis (x)
                    );
                    mlinterp::interp( &mpas_zlevels, 1,                                 // Number of points
                                      roms_zvalsT.data(), &tgt_temperature_data[i + offset],  // Output axis (y)
                                      zmh_z.data(), &roms_zlocation                           // Input axis (x)
                    );
#else
                    tgt_salinity_data[i + offset]    = splS( roms_zlocation );
                    tgt_temperature_data[i + offset] = splT( roms_zlocation );
#endif

                    // tgt_salinity_data[i + offset]    = pchipInterpolate( zmh_z, roms_zvalsS, roms_zlocation );
                    // tgt_temperature_data[i + offset] = pchipInterpolate( zmh_z, roms_zvalsT, roms_zlocation );
                    // tgt_salinity_data[i + offset]    = linearInterpolate( zmh_z, roms_zvalsS, roms_zlocation );
                    // tgt_temperature_data[i + offset] = linearInterpolate( zmh_z, roms_zvalsT, roms_zlocation );
                    if( tgt_salinity_data[i + offset] < 34.86 || tgt_salinity_data[i + offset] > 35.32 )
                    {
                        printf( "---(%zu, %d) ROMS z = %f, offset = %d, value = %f, zmh_z: [%f, %f]\n", i, irzl,
                                roms_zlocation, offset, tgt_salinity_data[i + offset], zmh_z[mpas_zreflevels - 1],
                                zmh_z[mpas_zreflevels-2] );
                    }
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
                for( int imzl = 0; imzl < mpas_zlevels; ++imzl )
                {

                    // mpasroms_zvalsS[mpas_zlevels - 1 - imzl] = imzl + 1 > maxlevelFace[i]
                    //                                             ? src_salinity_data[maxlevelFace[i] - 1 + offset]
                    //                                             : src_salinity_data[imzl + offset];
                    // mpasroms_zvalsT[mpas_zlevels - 1 - imzl] = imzl + 1 > maxlevelFace[i]
                    //                                             ? src_temperature_data[maxlevelFace[i] - 1 + offset]
                    //                                             : src_temperature_data[imzl + offset];

                    mpasroms_zvalsS[imzl] = imzl + 1 > maxlevelFace[i]
                                                                ? src_salinity_data[maxlevelFace[i] - 1 + offset]
                                                                : src_salinity_data[imzl + offset];
                    mpasroms_zvalsT[imzl] = imzl + 1 > maxlevelFace[i]
                                                                ? src_temperature_data[maxlevelFace[i] - 1 + offset]
                                                                : src_temperature_data[imzl + offset];
                }

#ifndef VERTICAL_INTERPOLANT_LINEAR
                tk::spline splS( zmh_z, mpasroms_zvalsS, tk::spline::cspline_hermite, true );
                tk::spline splT( zmh_z, mpasroms_zvalsT, tk::spline::cspline_hermite, true );
#endif

                const size_t offsetr = i * roms_zlevels;
                double roms_ztotal   = 0.0;
                for( size_t irzl = 0; irzl < roms_zlevels; ++irzl )
                {
                    // const double pbathymetry = zrh_xyz2d[i];
                    // const double delz        = pbathymetry / dst_zlayers;

                    // double roms_zlocation = ( irzl + 0.5 ) * delz;
                    double roms_zlocation = roms_ztotal + zrh_xyz3d[irzl] / 2;
                    roms_ztotal += zrh_xyz3d[irzl];
                    // Project data from zmh_z to zrh_z
#ifdef VERTICAL_INTERPOLANT_LINEAR
                    mlinterp::interp( &mpas_zlevels, 1,                                    // Number of points
                                      roms_zvalsS.data(), &srctgt_salinity_data[irzl + offsetr],  // Output axis (y)
                                      zmh_z.data(), &roms_zlocation                        // Input axis (x)
                    );
                    mlinterp::interp( &mpas_zlevels, 1,                                       // Number of points
                                      roms_zvalsT.data(), &srctgt_temperature_data[irzl + offsetr],  // Output axis (y)
                                      zmh_z.data(), &roms_zlocation                           // Input axis (x)
                    );
#else
                    srctgt_salinity_data[irzl + offsetr]    = splS( roms_zlocation );
                    srctgt_temperature_data[irzl + offsetr] = splT( roms_zlocation );
#endif
                    // srctgt_salinity_data[irzl + offsetr]    = pchipInterpolate( zmh_z, mpasroms_zvalsS, roms_zlocation );
                    // srctgt_temperature_data[irzl + offsetr] = pchipInterpolate( zmh_z, mpasroms_zvalsT, roms_zlocation );
                }
            }
#endif

            for( size_t ii = 0; ii < roms_zlevels; ii++ )
            {
                DataArray1D< double > dataInDoubleS( mpassize ), dataInDoubleT( mpassize );
                dbgprint( "Computing projection for ROMS level: " << ii );

                DataArray1D< double > dataOutDoubleS( romssize, false ), dataOutDoubleT( romssize, false );

                const size_t offsetr = ii * romssize;
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

moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& /*xyzd*/,
                                            std::vector< double >& /*fd*/,
                                            std::vector< double >& /*xyzi*/,
                                            std::vector< double >& /*fi*/ )
{
    // int nd = fd.size();
    // int ni = fi.size();

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
        // constexpr int nlevels         = 7;
        // std::array< size_t, 3 > grid = { 15, 15, mpas_zlevels };
        constexpr int nlevels        = 9;
        std::array< size_t, 3 > grid = { 4, 4, mpas_zlevels/4 };
        std::array< double, 6 > bbox = { -1.0, -1.0, -1E6, 1.0, 1.0, 1000 };
        std::cout << "\nComputing MBA interpolant (order=4, degree=3) for field " << varProjectSrc << std::endl;
        err = ComputeMBAInterpolant( src_xyz, src_tdata, dst_xyz, dst_tdata, is_three_dimensional, order, grid, bbox,
                                     nlevels );MB_CHK_ERR( err );
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
