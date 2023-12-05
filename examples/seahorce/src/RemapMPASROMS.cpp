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

#include "moab/ReadUtilIface.hpp"

// Other includes
#include "RemapMPASROMS.hpp"
#include "ComputeNN.hpp"
#include "ComputeShepard.hpp"
#include "ComputeDelaunay.hpp"
#include "HermiteCubicCurve.hpp"
#include "ComputeTR.hpp"
#include "ComputeMBA.hpp"
// #include "ComputeRBF.hpp"

#ifdef USE_GMLS
#include "ComputeGMLS.hpp"
#endif
#include "moab/Remapping/mlinterp.hpp"
#include "MeshUtilities.hpp"

using namespace moab;
// using namespace std;

// Utility macros
#define dbgprint( MSG )                                           \
    do                                                            \
    {                                                             \
        if( context.proc_id == 0 ) std::cout << MSG << std::endl; \
    } while( false )

//
// Start of main test program
//
int main( int argc, char** argv )
{
    // constexpr double radius = 6371220.0;
    constexpr double radius = 1.0;

    // Initialize MPI first
    MPI_Init( &argc, &argv );

    {
        // Let us create the runtime context with default params
        RuntimeContext context;

        int src_zlayers = 0;
        int dst_zlayers = 0;
        const int rank  = context.proc_id;
        const int size  = context.num_procs;

        // get the moab instance
        Interface* mbi = context.moab_interface;

        // write options for h5m files
        const std::string write_options = size > 1 ? "PARALLEL=WRITE_PART" : "";

        // Get the command-line options for the run
        context.ParseCLOptions( argc, argv );

        // set the number of z-layers
        src_zlayers = context.use_3dprojection ? mpas_zlevels : 1;
        dst_zlayers = context.use_3dprojection ? roms_zlevels : 1;

        // print details about runtime parameters
        context.describe();

        // Load the MPAS file from disk with given options
        std::vector< moab::EntityHandle > mpas_verts, mpas_elems;
        std::vector< moab::EntityHandle > mpas3d_verts, mpas3d_elems;
        {
            dbgprint( "Reading MPAS file from disk" );
            context.timer_push( "Load ROMS 2D mesh file" );
            runchk( mbi->load_file( context.mpas_filename.c_str(), &context.mpasset ),
                    "MOAB::load_file for MPAS mesh failed" );
            // Get all entities in the database
            runchk( mbi->get_entities_by_dimension( context.mpasset, 0, mpas_verts ) );
            runchk( mbi->get_entities_by_dimension( context.mpasset, 2, mpas_elems ) );
            // Rescale the radius to unit sphere to compute the intersection
            runchk( ScaleCoords( mbi, mpas_verts, radius, true ) );
            context.timer_pop();

            dbgprint( "MPAS mesh contains " << mpas_verts.size() << " vertices and " << mpas_elems.size()
                                            << " elements" );

#ifdef VERBOSE_OUTPUT
            runchk( mbi->write_file( "mpas_modified_2d.h5m", "H5M", write_options.c_str(), &context.mpasset, 1 ) );
#endif
        }

        // Load the ROMS file from disk with given options
        std::vector< moab::EntityHandle > roms_verts, roms_elems;
        std::vector< moab::EntityHandle > roms3d_verts, roms3d_elems;
        {
            dbgprint( "Reading ROMS file from disk" );
            context.timer_push( "Load ROMS 2D mesh file" );
            runchk( mbi->load_file( context.roms_filename.c_str(), &context.romsset ),
                    "MOAB::load_file for ROMS mesh failed" );
            // Get all entities in the database
            runchk( mbi->get_entities_by_dimension( context.romsset, 0, roms_verts ) );
            runchk( mbi->get_entities_by_dimension( context.romsset, 2, roms_elems ) );
            // Rescale the radius to unit sphere to compute the intersection
            runchk( ScaleCoords( mbi, roms_verts, radius, false ) );
            context.timer_pop();

            dbgprint( "ROMS mesh contains " << roms_verts.size() << " vertices and " << roms_elems.size()
                                            << " elements" );

#ifdef VERBOSE_OUTPUT
            runchk( mbi->write_file( "roms_modified_2d.h5m", "H5M", write_options.c_str(), &context.romsset, 1 ) );
#endif
        }

        // Cull the MPAS set so that we don't have a global mesh
        {
            context.timer_push( "Cull MPAS surface mesh: covering region" );
            const int nring_neighborhood = 2;
            // construct a kd-tree index:
            using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                                PC3D< double >, 3 /* dim */
                                                                >;

            moab::Range orig_mpas_elems;
            runchk( mbi->get_entities_by_dimension( context.mpasset, 2, orig_mpas_elems ) );
            std::vector< double > mpas_xyz( orig_mpas_elems.size() * 3 );
            runchk( mbi->get_coords( orig_mpas_elems, mpas_xyz.data() ) );

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
                runchk( mbi->get_coords( &ehandle, 1, &query_pt[0] ) );

                // Do a KNN search
                resultSet.init( srcindx.data(), srcdist.data() );
                tree.findNeighbors( resultSet, query_pt );

                for( size_t j = 0; j < num_results; ++j )
                    lelems.insert( orig_mpas_elems[srcindx[j]] );
            }
            runchk( mbi->add_entities( context.mpas_covering_set, lelems ) );
            context.timer_pop();

#ifdef VERBOSE_OUTPUT
            runchk( mbi->write_file( "mpas_covering_2d.h5m", "H5M", write_options.c_str(), &context.mpas_covering_set,
                                     1 ) );
#endif

            runchk( mbi->get_connectivity( lelems, lverts, true ) );
            runchk( mbi->add_entities( context.mpas_covering_set, lverts ) );
            runchk( mbi->get_entities_by_dimension( context.mpas_covering_set, 2, mpas_elems ) );
            runchk( mbi->get_entities_by_dimension( context.mpas_covering_set, 0, mpas_verts ) );
            dbgprint( "Culled MPAS mesh contains " << mpas_elems.size() << " elements and " << mpas_verts.size()
                                                   << " vertices." );
        }

        if( context.computeTRMaps )
        {
            context.timer_push( "Compute TempestRemap weights for method: " +
                                RuntimeContext::GetMethod( context.field_methods["Bathymetry"].first ) );
            // call to compute the 2D map and store to disk
            runchk( ComputeTempestRemapWeights( context, context.mpas_covering_set, context.romsset ),
                    "Cannot compute 2D remapping weights" );
            context.timer_pop();

            // context.timer_push( "Load TempestRemap weights for method: " + context.bathymetryMethod );
            // // load the computed 2D map files
            // runchk( LoadTempestRemapWeights( context, context.mpas_covering_set, context.romsset,
            //                                  context.bathymetryMethod ),
            //         "Cannot load 2D remapping weights" );
            // context.timer_pop();
        }

        // let us perform 3D extrusions as needed
        EntityHandle root_set  = 0;
        EntityHandle mpasset3d = 0, romsset3d = 0;


        {
            // Initialize all important data

            // Project the bottom Bathymetry and SeaSurfaceHeight data from MPAS to ROMS so that we can impose it.
            for (int iv=0; iv < nstandardvars; ++iv)
            {
                context.timer_push( "Project " + std::string(mpas_twod_standardtagnames[iv]) + " field" );
                runchk( context.ComputeFieldProjections( 2, mpas_twod_standardtagnames[iv], roms_twod_standardtagnames[iv],
                                                        mpas_elems, roms_elems ),
                        "Can't project " + std::string(mpas_twod_standardtagnames[iv]) + " field" );
                context.timer_pop();
            }

            // Next project the forcing functions from MPAS to ROMS so that we can impose it as boundary conditions
            for (int iv=0; iv < nforcingvars; ++iv)
            {
                context.timer_push( "Project " + std::string(mpas_twod_forcingtagnames[iv]) + " field" );
                runchk( context.ComputeFieldProjections( 2, mpas_twod_forcingtagnames[iv], mpas_twod_forcingtagnames[iv],
                                                        mpas_elems, roms_elems ),
                        "Can't project " + std::string(mpas_twod_forcingtagnames[iv]) + " field" );
                context.timer_pop();
            }

        }

        std::vector< double > zmh_xyz3d, zrh_xyz3d;
        if( context.use_3dprojection || context.threetwooneD )
        {
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
                runchk( mbi->tag_get_handle( "refBottomDepth", mpas_zreflevels, moab::MB_TYPE_DOUBLE, mztag,
                                             moab::MB_TAG_SPARSE ),
                        "Can't get tag handle: refBottomDepth" );

                runchk( mbi->tag_get_data( mztag, &root_set, 1, context.mpas_zref_heights ),
                        "Can't get refBottomDepth data" );
                for( int j = 0; j < src_zlayers; ++j )
                {
                    context.mpas_zref_heights[j] /= axial_scaling;
                }
            }

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
                        zmh_xyz3d[offset + j] = ( context.mpas_zref_heights[j] - context.mpas_zref_heights[j - 1] );
                }
            }
            else
            {
                moab::Tag mhtag;
                runchk( mbi->tag_get_handle( "layerThickness_3d", src_zlayers, moab::MB_TYPE_DOUBLE, mhtag,
                                             moab::MB_TAG_DENSE ),
                        "Can't get tag handle: layerThickness_3d" );
                runchk( mbi->tag_get_data( mhtag, mpas_elems.data(), mpas_elems.size(), zmh_xyz3d.data() ),
                        "Can't get layerThickness_3d data" );
            }

            if( context.generateExtrusions )
            {
                dbgprint( "\nExtruding MPAS polygonal mesh ..." );
                context.timer_push( "Extrude MPAS 2D polygonal mesh to 3D polyhedral mesh" );
                runchk( ExtrudePolygonsToPolyhedra( context, zmh_xyz3d, context.mpas_covering_set, mpasset3d, true,
                                                    src_zlayers ),
                        "Can't extrude MPAS polygons" );
                context.timer_pop();

                {
                    std::vector< double > zrh_xyz2d( roms_elems.size() );
                    moab::Tag rhtag;
                    runchk( mbi->tag_get_handle( "Bathymetry", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE ),
                            "Can't get bottomDepth tag" );
                    runchk( mbi->tag_get_data( rhtag, roms_elems.data(), roms_elems.size(), zrh_xyz2d.data() ),
                            "Can't get bottomDepth tag data" );

#ifdef VERBOSE_OUTPUT
                    runchk(
                        mbi->write_file( "roms_3d_2dsurface.h5m", "H5M", write_options.c_str(), &context.romsset, 1 ) );
#endif

#ifdef USE_STRETCHING_FUNCTION
                    /// Call stretching functions
                    // h = np.linspace( 10, 200, 10 );
                    constexpr double zeta = 0.0;

                    // return N evenly s-coordinate w-points.
                    auto Sw = [&]() {
                        std::vector< double > s( roms_zlevels + 1 );
                        double del = 1.0 / roms_zlevels;
                        s[0]       = 0.0;
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
                        constexpr double hc = 100 / axial_scaling;  // hc has to be less than or equal to min(h)
                        double vstretch_val = vstretching_1( s );
                        return hc * ( s - vstretch_val ) + vstretch_val * h;
                    };

                    auto vtransform_2 = [&]( double s, double h ) {
                        constexpr double hc = 100 / axial_scaling;  // hc has to be less than or equal to min(h)
                        double vstretch_val = vstretching_4( s );
                        return ( hc * s + vstretch_val * h ) / ( hc + h );
                    };
#endif
                    constexpr int transform_id = 2;  // 1 or 2

                    for( size_t i = 0; i < roms_elems.size(); ++i )
                    {
                        const double pbathymetry = zrh_xyz2d[i] / axial_scaling;
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
                                    z_w[ilevel + 1] = ( vtransform_1( s_w[ilevel + 1], h ) +
                                                        zeta * ( 1 + vtransform_1( s_w[ilevel + 1], h ) / h ) );
                                else
                                    z_w[ilevel + 1] = ( zeta + ( zeta + h ) * vtransform_2( s_w[ilevel + 1], h ) );
                                // zrh_xyz3d[offset + ilevel] = 0.5 * ( z_w[ilevel] + z_w[ilevel + 1] );
                                zrh_xyz3d[offset + ilevel] = -( z_w[ilevel + 1] - z_w[ilevel] );
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
                    context.timer_push( "Extrude ROMS 2D-QUAD mesh to 3D-HEX mesh with adaptive z-layers" );
                    runchk( ExtrudePolygonsToPolyhedra( context, zrh_xyz3d, context.romsset, romsset3d, false,
                                                        dst_zlayers ),
                            "Can't extrude ROMS polygons" );
                    context.timer_pop();
                }

                runchk( mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts ) );
                runchk( mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems ) );

                // std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size() << " elements.\n";
                runchk( mbi->write_file( "mpas_full_3d.h5m", "H5M", write_options.c_str(), &mpasset3d, 1 ) );

                runchk( mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts ) );
                runchk( mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems ) );

                {
                    moab::Tag rh3tag;
                    runchk( mbi->tag_get_handle( "ROMSlayerThickness", 1, moab::MB_TYPE_DOUBLE, rh3tag,
                                                 moab::MB_TAG_DENSE | moab::MB_TAG_CREAT ),
                            "Can't get tag handle: ROMSlayerThickness" );
                    runchk( mbi->tag_set_data( rh3tag, roms3d_elems.data(), roms3d_elems.size(), zrh_xyz3d.data() ),
                            "Can't get layerThickness_3d data" );
                }

                // std::cout << "3D ROMS: " << roms3d_verts.size() << " vertices and " << roms3d_elems.size() << " elements.\n";
                // Rescale the radius of both to compute the intersection
                // runchk( ScaleCoords( mbi, roms3d_verts, radius, false ) );
                runchk( mbi->write_file( "roms_full_3d.h5m", "H5M", write_options.c_str(), &romsset3d, 1 ) );
            }
            else
            {
                runchk( mbi->create_meshset( moab::MESHSET_SET, mpasset3d ), "Can't create new set" );
                runchk( mbi->create_meshset( moab::MESHSET_SET, romsset3d ), "Can't create new set" );

                runchk( mbi->load_file( "mpas_full_3d.h5m", &mpasset3d ), "MOAB::load_file for MPAS 3D mesh failed" );
                runchk( mbi->load_file( "roms_full_3d.h5m", &romsset3d ), "MOAB::load_file for ROMS 3D mesh failed" );

                runchk( mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts ) );
                runchk( mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems ) );
                runchk( mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts ) );
                runchk( mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems ) );

                std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size()
                          << " elements.\n";
                std::cout << "3D ROMS: " << roms3d_verts.size() << " vertices and " << roms3d_elems.size()
                          << " elements.\n";

                moab::Tag rh3tag;
                runchk( mbi->tag_get_handle( "ROMSlayerThickness", 1, moab::MB_TYPE_DOUBLE, rh3tag,
                                             moab::MB_TAG_DENSE | moab::MB_TAG_CREAT ),
                        "Can't get tag handle: ROMSlayerThickness" );

                runchk( mbi->tag_get_data( rh3tag, roms3d_elems.data(), roms3d_elems.size(), zrh_xyz3d.data() ),
                        "Can't get layerThickness_3d data" );
            }
        }

        bool twoDfirst = !context.oneDfirst;

        if( context.threetwooneD )
        {
            const size_t mpassize = mpas_elems.size();
            const size_t romssize = roms_elems.size();

            // context.weightMap.SetEnforcementBounds( "Lp", &context.meshInput, &context.meshOverlap, nullptr, nullptr, 1 );

            // get the handle to the weight matrix
            const SparseMatrix< double >& weights = context.weightMap.GetSparseMatrix();

            double defaultvalue = -1.0;
            moab::Tag mpas_soltags_3d[nvars], roms_soltags_elem[nvars];
            for (int iv=0; iv < nvars; ++iv)
            {
                runchk( mbi->tag_get_handle( mpas_tagnames[iv], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
                                             mpas_soltags_3d[iv], moab::MB_TAG_DENSE ),
                        "Can't get MPAS " << mpas_tagnames[iv] << " tag" );

                runchk( mbi->tag_get_handle( roms_tagnames[iv], 1, moab::MB_TYPE_DOUBLE, roms_soltags_elem[iv],
                                             moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defaultvalue ),
                        "Can't create ROMS " << roms_tagnames[iv] << " tag" );
            }
            // runchk( mbi->tag_get_handle( mpas_tagnames[0], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
            //                              mpas_soltags_3d[0], moab::MB_TAG_DENSE ),
            //         "Can't create salinity tag" );
            // runchk( mbi->tag_get_handle( mpas_tagnames[1], mpas_zreflevels, moab::MB_TYPE_DOUBLE,
            //                              mpas_soltags_3d[1], moab::MB_TAG_DENSE ),
            //         "Can't create temperature tag" );
            // runchk( mbi->tag_get_handle( roms_threed_tagnames[0], 1, moab::MB_TYPE_DOUBLE, roms_soltags_elem[0],
            //                              moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defaultvalue ),
            //         "Can't create salinity tag on ROMS3D" );
            // runchk( mbi->tag_get_handle( roms_threed_tagnames[1], 1, moab::MB_TYPE_DOUBLE, roms_soltags_elem[1],
            //                              moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defaultvalue ),
            //         "Can't create temperature tag on ROMS3D" );

            const int svaroffset = mpassize * mpas_zreflevels;
            std::vector< double > src_data( nvars * svaroffset );
            // get the source data from tags
            for( int iv = 0; iv < nvars; ++iv )
                runchk( mbi->tag_get_data( mpas_soltags_3d[iv], mpas_elems.data(), mpassize,
                                           src_data.data() + svaroffset * iv ) );

            moab::Tag rhtag;
            runchk( mbi->tag_get_handle( "Bathymetry", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE ),
                    "Can't get Bathymetry tag handle" );
            std::vector< double > zrh_xyz2d( romssize );
            runchk( mbi->tag_get_data( rhtag, roms_elems.data(), romssize, zrh_xyz2d.data() ),
                    "Can't get Bathymetry tag data" );

            std::vector< int > minlevelFace, maxlevelFace;
            moab::Tag minlvlTag, maxlvlTag;
            runchk( mbi->tag_get_handle( "minLevelCell", 1, moab::MB_TYPE_INTEGER, minlvlTag, moab::MB_TAG_DENSE ) );
            minlevelFace.resize( mpassize );
            runchk( mbi->tag_get_data( minlvlTag, mpas_elems.data(), mpassize, minlevelFace.data() ) );

            runchk( mbi->tag_get_handle( "maxLevelCell", 1, moab::MB_TYPE_INTEGER, maxlvlTag, moab::MB_TAG_DENSE ) );
            maxlevelFace.resize( mpassize );
            runchk( mbi->tag_get_data( maxlvlTag, mpas_elems.data(), mpassize, maxlevelFace.data() ) );

            std::vector< double > zmh_z( mpas_zreflevels );
            zmh_z[0] = 0.5 * context.mpas_zref_heights[0];
            for( size_t j = 1; j < mpas_zlevels; ++j )
            {
                zmh_z[j] = 0.5 * ( context.mpas_zref_heights[j - 1] + context.mpas_zref_heights[j] );
            }

            const int tvaroffset = romssize * roms_zlevels;
            std::vector< double > tgt_data( nvars * tvaroffset );
            if( twoDfirst )
            {
                // First compute the projections in 2D for each MPAS layer.
                // This will give the MPAS projected data on to ROMS mesh - on MPAS axial levels
                // NOTE: implicit assumption is that ROMS levels >= MPAS levels. Should fix how
                // the projection is invoked by perhaps skipping the tag_set_data
                const int tsvaroffset = romssize * mpas_zreflevels;
                std::vector< double > tgtsrc_data( nvars * tsvaroffset );

                context.timer_push( "Compute 3D projection: 2Dx1D algorithm" );
#pragma omp parallel for shared( tgtsrc_data, weights )
                for( int ii = 0; ii < mpas_zlevels; ii++ )
                {
                    std::cout << "Computing projection for MPAS level: " + std::to_string( ii ) + "\n";

                    for( int iv = 0; iv < nvars; ++iv )
                    {
                        DataArray1D< double > dataInDouble( mpassize );
                        DataArray1D< double > dataOutDouble( romssize, false );
                        unsigned offsetr = romssize * ii;

                        for( size_t j = 0; j < mpassize; ++j )
                        {
                            const int offset = svaroffset * iv + j * mpas_zlevels;
                            dataInDouble[j] = maxlevelFace[j] - 1 < ii || src_data[ii + offset] < 1e-8
                                                   ? src_data[maxlevelFace[j] - 1 + offset]
                                                   : src_data[ii + offset];
                        }

                        dataOutDouble.AttachToData( tgtsrc_data.data() + offsetr + tsvaroffset * iv );

                        // Compute the projection for the salinity field
                        weights.Apply( dataInDouble, dataOutDouble );

                        if( context.useCAAS )
                            ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap,
                                               context.field_methods["Bathymetry"].second /*nPin*/, dataInDouble,
                                               dataOutDouble, true /*useCAASLocal*/ );
                    }
                }

#ifdef VERTICAL_INTERPOLATION
                // We already have the 2D projection on all MPAS z-level to ROMS (on MPAS z-levels)
                // Now perform the vertical interpolation from MPAS-z-levels to ROMS-z-levels

                // Now project each data layer in axial direction.
                // NOTE: Embarassingly parallel
#pragma omp parallel for shared( zrh_xyz3d, tgtsrc_data, tgt_data )
                for( size_t i = 0; i < romssize; ++i )
                {
                    std::vector< double > roms_zvals( mpas_zreflevels );
                    for( int iv = 0; iv < nvars; ++iv )
                    {
                        // const double pbathymetry = zrh_xyz2d[i];
                        // const double delz        = pbathymetry / roms_zlevels;
                        // if (i<10) printf( "--- ROMS delz = %f\n", delz );

#ifndef VERTICAL_INTERPOLANT_LINEAR
                        // tk::spline splS( zmh_z, roms_zvalsS, tk::spline::cspline_hermite, true );
                        // tk::spline splT( zmh_z, roms_zvalsT, tk::spline::cspline_hermite, true );
                        HermiteCubicCurve< double > splHC;
#endif
                        for( int imzl = 0; imzl < mpas_zreflevels; ++imzl )
                        {
                            const int offset = imzl * romssize;
                            // roms_zvalsS[mpas_zreflevels - 1 - imzl] = tgtsrc_salinity_data[i + offset];
                            // roms_zvalsT[mpas_zreflevels - 1 - imzl] = tgtsrc_temperature_data[i + offset];

                            roms_zvals[imzl] = tgtsrc_data[tsvaroffset * iv + i + offset];
#ifndef VERTICAL_INTERPOLANT_LINEAR
                            splHC.add( zmh_z[imzl], roms_zvals[imzl] );
#endif
                        }
#ifndef VERTICAL_INTERPOLANT_LINEAR
                        splHC.finish();
#endif
                        double rcoords[3];
                        auto entROMS = roms_elems[i];
                        mbi->get_coords( &entROMS, 1, rcoords );

                        double roms_ztotal = rcoords[2];

                        for( int irzl = 0; irzl < roms_zlevels; ++irzl )
                        {
                            const double* zrh_xyz3d_loc = zrh_xyz3d.data() + i * roms_zlevels;
                            const int offset            = irzl * romssize;
                            // double roms_zlocation = ( irzl + 0.5 ) * delz;
                            double roms_zlocation = roms_ztotal + zrh_xyz3d_loc[irzl] / 2;
                            roms_ztotal += zrh_xyz3d_loc[irzl];

                            // Project data from zmh_z to zrh_z
#ifdef VERTICAL_INTERPOLANT_LINEAR
                            mlinterp::interp( &mpas_zlevels, 1,  // Number of points
                                              roms_zvals.data(),
                                              &tgt_data[iv * tvaroffset + i + offset],  // Output axis (y)
                                              zmh_z.data(), &roms_zlocation             // Input axis (x)
                            );
#else
                            // tgt_data[iv * tvaroffset + i + offset] = splS( roms_zlocation );
                            tgt_data[iv * tvaroffset + i + offset] = splHC.at( roms_zlocation );
#endif
                        }
                    }
                }
#endif
                context.timer_pop();
            }
//             else  // Vertical first and horizontal next
//             {
//                 std::vector< double > srctgt_salinity_data( mpassize * roms_zlevels ),
//                     srctgt_temperature_data( mpassize * roms_zlevels );

// #ifdef VERTICAL_INTERPOLATION
//                 // Now project each data layer in axial direction.
//                 // NOTE: Embarassingly parallel
//                 context.timer_push( "Compute 3D projection: 1Dx2D algorithm" );
// #pragma omp parallel for shared( zrh_xyz3d, src_salinity_data, src_temperature_data, tgt_salinity_data, \
//                                      tgt_temperature_data, srctgt_salinity_data, srctgt_temperature_data )
//                 for( size_t i = 0; i < mpassize; ++i )
//                 {
//                     std::vector< double > mpasroms_zvalsS( mpas_zlevels ),
//                         mpasroms_zvalsT( mpas_zlevels );  // Data on MPAS 2D but ROMS 1D vertical
//                     const size_t offset = i * mpas_zlevels;
//                     for( int imzl = 0; imzl < mpas_zlevels; ++imzl )
//                     {
//                         mpasroms_zvalsS[imzl] = imzl + 1 > maxlevelFace[i] || src_salinity_data[imzl + offset] < 1e-8
//                                                     ? src_salinity_data[maxlevelFace[i] - 1 + offset]
//                                                     : src_salinity_data[imzl + offset];
//                         mpasroms_zvalsT[imzl] = imzl + 1 > maxlevelFace[i] || src_temperature_data[imzl + offset] < 1e-8
//                                                     ? src_temperature_data[maxlevelFace[i] - 1 + offset]
//                                                     : src_temperature_data[imzl + offset];
//                     }

// #ifndef VERTICAL_INTERPOLANT_LINEAR
//                     // tk::spline splS( zmh_z, mpasroms_zvalsS, tk::spline::cspline_hermite, true );
//                     // tk::spline splT( zmh_z, mpasroms_zvalsT, tk::spline::cspline_hermite, true );
//                     HermiteCubicCurve< double > splHC_S;
//                     HermiteCubicCurve< double > splHC_T;
//                     for( int imzl = 0; imzl < mpas_zlevels; imzl++ )
//                     {
//                         splHC_S.add( zmh_z[imzl], mpasroms_zvalsS[imzl] );
//                         splHC_T.add( zmh_z[imzl], mpasroms_zvalsT[imzl] );
//                     }
//                     splHC_S.finish();
//                     splHC_T.finish();
// #endif

//                     double rcoords[3] = { 0, 0, 0 };

//                     // auto entROMS = roms_elems[i];
//                     // mbi->get_coords( &entROMS, 1, rcoords );
//                     // for( int irzl = 0; irzl < roms_zlevels; ++irzl )
//                     // {
//                     //     const double* zrh_xyz3d_loc = zrh_xyz3d.data() + i * roms_zlevels;
//                     //     const int offset            = irzl * romssize;
//                     //     // double roms_zlocation = ( irzl + 0.5 ) * delz;
//                     //     double roms_zlocation = roms_ztotal + zrh_xyz3d_loc[irzl] / 2;
//                     //     roms_ztotal += zrh_xyz3d_loc[irzl];
//                     // }

//                     const size_t offsetr = i * roms_zlevels;
//                     double roms_ztotal   = rcoords[2];
//                     for( size_t irzl = 0; irzl < roms_zlevels; ++irzl )
//                     {
//                         // const double pbathymetry = zrh_xyz2d[i];
//                         // const double delz        = pbathymetry / dst_zlayers;

//                         // double roms_zlocation = ( irzl + 0.5 ) * delz;
//                         double roms_zlocation = roms_ztotal + zrh_xyz3d[irzl] / 2;
//                         roms_ztotal += zrh_xyz3d[irzl];
//                         // Project data from zmh_z to zrh_z
// #ifdef VERTICAL_INTERPOLANT_LINEAR
//                         mlinterp::interp( &mpas_zlevels, 1,  // Number of points
//                                           roms_zvalsS.data(), &srctgt_salinity_data[irzl + offsetr],  // Output axis (y)
//                                           zmh_z.data(), &roms_zlocation                               // Input axis (x)
//                         );
//                         mlinterp::interp( &mpas_zlevels, 1,  // Number of points
//                                           roms_zvalsT.data(),
//                                           &srctgt_temperature_data[irzl + offsetr],  // Output axis (y)
//                                           zmh_z.data(), &roms_zlocation              // Input axis (x)
//                         );
// #else
//                         // srctgt_salinity_data[irzl + offsetr]    = splS( roms_zlocation );
//                         // srctgt_temperature_data[irzl + offsetr] = splT( roms_zlocation );
//                         srctgt_salinity_data[irzl + offsetr]    = splHC_S.at( roms_zlocation );
//                         srctgt_temperature_data[irzl + offsetr] = splHC_T.at( roms_zlocation );
// #endif
//                         // srctgt_salinity_data[irzl + offsetr]    = pchipInterpolate( zmh_z, mpasroms_zvalsS, roms_zlocation );
//                         // srctgt_temperature_data[irzl + offsetr] = pchipInterpolate( zmh_z, mpasroms_zvalsT, roms_zlocation );
//                     }
//                 }
// #endif

// #pragma omp parallel for shared( srctgt_salinity_data, srctgt_temperature_data, tgt_salinity_data, \
//                                      tgt_temperature_data, weights )
//                 for( size_t ii = 0; ii < roms_zlevels; ii++ )
//                 {
//                     DataArray1D< double > dataInDoubleS( mpassize ), dataInDoubleT( mpassize );
//                     std::cout << "Computing projection for ROMS level: " + std::to_string( ii ) + "\n";

//                     DataArray1D< double > dataOutDoubleS( romssize, false ), dataOutDoubleT( romssize, false );

//                     const size_t offsetr = ii * romssize;
//                     // std::vector< EntityHandle > mpas_slice( mpas3d_elems.begin() + offsetr,
//                     //                                         mpas3d_elems.begin() + offsetr + mpassize );
//                     std::vector< EntityHandle > roms_slice( roms3d_elems.begin() + offsetr,
//                                                             roms3d_elems.begin() + offsetr + romssize );

//                     for( size_t j = 0; j < mpassize; ++j )
//                     {
//                         const size_t offset = j * roms_zlevels;
//                         dataInDoubleS[j]    = srctgt_salinity_data[ii + offset];
//                         dataInDoubleT[j]    = srctgt_temperature_data[ii + offset];
//                     }

//                     // dataInDoubleS.AttachToData( src_salinity_data.data() + offsetr );
//                     dataOutDoubleS.AttachToData( tgt_salinity_data.data() + offsetr );

//                     // Compute the projection for the salinity field
//                     weights.Apply( dataInDoubleS, dataOutDoubleS );

//                     if( context.useCAAS )
//                         ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap,
//                                            context.field_methods["Bathymetry"].second /*nPin*/, dataInDoubleS,
//                                            dataOutDoubleS, true /*useCAASLocal*/ );

// // #ifndef VERTICAL_INTERPOLATION
// //                     runchk( mbi->tag_set_data( roms_soltags_elem[0], roms_slice.data(), roms_slice.size(),
// //                                                dataOutDoubleS ),
// //                             "Can't set salinity tag data" );
// // #endif
//                     // dataInDoubleT.AttachToData( src_temperature_data.data() + offsetr );
//                     dataOutDoubleT.AttachToData( tgt_temperature_data.data() + offsetr );

//                     // Compute the projection for the temperature field
//                     weights.Apply( dataInDoubleT, dataOutDoubleT );

//                     if( context.useCAAS )
//                         ApplyCAASLimiting( context.weightMap, context.meshInput, context.meshOverlap,
//                                            context.field_methods["Bathymetry"].second /*nPin*/, dataInDoubleT,
//                                            dataOutDoubleT, true /*useCAASLocal*/ );

// // #ifndef VERTICAL_INTERPOLATION
// //                     runchk( mbi->tag_set_data( roms_soltags_elem[1], roms_slice.data(), roms_slice.size(),
// //                                                dataOutDoubleT ),
// //                             "Can't set temperature tag data" );
// // #endif
//                 }
//                 context.timer_pop();
//             }

#ifdef VERTICAL_INTERPOLATION
            for( int iv = 0; iv < nvars; ++iv )
                runchk( mbi->tag_set_data( roms_soltags_elem[iv], roms3d_elems.data(), roms3d_elems.size(),
                                           tgt_data.data() + iv * tvaroffset ),
                        "Can't set tag data" );
#endif
        }
        else  //if( context.use_3dprojection || context.computeMBAInterpolant )
        {
            if( context.use_3dprojection && ( context.field_methods["Salinity"].first == DelaunayInterpolant ||
                                              context.field_methods["Temperature"].first == DelaunayInterpolant ) )
            {
                // get the coordinates of the elements
                std::vector< double > src_xyz( mpas3d_elems.size() * 3 );
                runchk( mbi->get_coords( mpas3d_elems.data(), mpas3d_elems.size(), src_xyz.data() ) );
                // setup the necessary data-structures for computing the Delaunay interpolant
                runchk( SetupDelaunayInterpolant( context, src_xyz ) );
            }

            for( int iv = 0; iv < nvars; ++iv )
            {
                if( context.use_3dprojection )
                {
                    context.timer_push( "Compute 3D projection: " +
                                        RuntimeContext::GetMethod( context.field_methods[roms_tagnames[iv]].first ) +
                                        " algorithm" );
                    // Now let us compute the 3D field projections for each field
                    runchk( context.ComputeFieldProjections( 3, mpas_ele_tagnames[iv], roms_tagnames[iv], mpas3d_elems,
                                                             roms3d_elems ) );
                }
                else
                {
                    context.timer_push( "Compute 2D projection: " +
                                        RuntimeContext::GetMethod( context.field_methods[roms_tagnames[iv]].first ) +
                                        " algorithm" );
                    // Now let us compute the 2D field projections for each field
                    runchk( context.ComputeFieldProjections( 2, mpas_tagnames[iv], roms_tagnames[iv], mpas_elems,
                                                             roms_elems ) );
                }
                context.timer_pop();
            }
        }
        dbgprint( std::endl );

        const std::string mpas_output_file = ( context.use_3dprojection ? "mpas_3d_source.h5m" : "mpas_2d_source.h5m" );
        dbgprint( "Writing out the source mesh with fields to '" << mpas_output_file << "'" );
        runchk( mbi->write_file(
            mpas_output_file.c_str(), "H5M", write_options.c_str(),
            ( context.use_3dprojection || context.threetwooneD ? &mpasset3d : &context.mpas_covering_set ), 1 ) );

        const std::string roms_output_file =
            ( context.use_3dprojection ? "roms_3d_projected.h5m" : "roms_2d_projected.h5m" );
        dbgprint( "Writing out the target mesh with projected fields to '" << roms_output_file << "'" );
        runchk( mbi->write_file( roms_output_file.c_str(), "H5M", write_options.c_str(),
                                 ( context.use_3dprojection || context.threetwooneD ? &romsset3d : &context.romsset ),
                                 1 ) );

        if( !rank ) dbgprint( "\n********** Remap MPAS-to-ROMS DONE! **********" );
    }

    // Done, cleanup
    MPI_Finalize();

    return 0;
}

moab::ErrorCode RuntimeContext::ComputeFieldProjections( int dimension,
                                                         std::string varProjectSrc,
                                                         std::string varProjectDst,
                                                         const std::vector< moab::EntityHandle >& source_range,
                                                         const std::vector< moab::EntityHandle >& target_range,
                                                         const double constantoffset )
{
    moab::ErrorCode err;
    moab::Tag dmtag;

    moab::Interface* mbi = this->moab_interface;

    // determine if this is either one of the standard 2D fields or forcing functions
    bool standard_fields = false;
    for (int iv = 0; iv < nstandardvars; ++iv)
        if ( !varProjectDst.compare( roms_twod_standardtagnames[iv] ) ) standard_fields = true;
    for (int iv = 0; iv < nforcingvars; ++iv)
        if ( !varProjectDst.compare( mpas_twod_forcingtagnames[iv] ) ) standard_fields = true;

    const bool is_three_dimensional = (dimension == 3);
    const bool normalize            = this->normalize;
    const RemappingMethod rmethod   = this->field_methods[varProjectDst].first;
    const int order                 = this->field_methods[varProjectDst].second;

    // get the source data from tag
    std::vector< double > src_tdata( source_range.size() ), dst_tdata( target_range.size() );
    if( standard_fields )
    {
        err = mbi->tag_get_handle( varProjectSrc.c_str(), 1, moab::MB_TYPE_DOUBLE, dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        err = mbi->tag_get_data( dmtag, source_range.data(), source_range.size(), src_tdata.data() );MB_CHK_ERR( err );
    }
    else
    {
        err = mbi->tag_get_handle( varProjectSrc.c_str(), mpas_zreflevels, moab::MB_TYPE_DOUBLE, dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        std::vector< double > src_tdata_layers( source_range.size() * mpas_zreflevels );
        err = mbi->tag_get_data( dmtag, source_range.data(), source_range.size(), src_tdata_layers.data() );MB_CHK_ERR( err );
        for( size_t il = 0; il < source_range.size(); ++il )
            src_tdata[il] = src_tdata_layers[il * mpas_zreflevels];
    }

    // get the coordinates of the elements
    std::vector< double > src_xyz( source_range.size() * 3 ), dst_xyz( target_range.size() * 3 );
    if( !( rmethod == TempestRemapFV || rmethod == TempestRemapBilinear || rmethod == TempestRemapInvDist ||
           rmethod == TempestRemapDelaunay || rmethod == TempestRemapIntegratedBilinear ) )
    {
        err = mbi->get_coords( source_range.data(), source_range.size(), src_xyz.data() );MB_CHK_ERR( err );
        err = mbi->get_coords( target_range.data(), target_range.size(), dst_xyz.data() );MB_CHK_ERR( err );
    }

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0, normFactor = 1.0;
    if( normalize && ( !is_three_dimensional ) && this->meshOverlap.faces.size() )
    {
        for( size_t j = 0; j < src_tdata.size(); j++ )
            src_tdata[j] -= constantoffset;

        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < this->meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = this->meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += src_tdata[iSourceFace] * this->meshOverlap.vecFaceArea[j];
        }
    }

    if( rmethod == TempestRemapFV || rmethod == TempestRemapBilinear || rmethod == TempestRemapInvDist ||
        rmethod == TempestRemapDelaunay || rmethod == TempestRemapIntegratedBilinear )
    {
        // get the handle to the weight matrix
        const SparseMatrix< double >& weights = this->weightMap.GetSparseMatrix();
        DataArray1D< double > dataInDouble( source_range.size(), false );
        dataInDouble.AttachToData( src_tdata.data() );
        DataArray1D< double > dataOutDouble( target_range.size(), false );
        dataOutDouble.AttachToData( dst_tdata.data() );

        // Compute the projection for the Bathymetry field
        weights.Apply( dataInDouble, dataOutDouble );
    }
    else if( rmethod == NearestNeighborInterpolant )
    {
        std::cout << "\nComputing Nearest-neighbor interpolant (order=1, degree=0) for field " << varProjectSrc
                  << std::endl;
        err = ComputeNNInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }
    else if( rmethod == DelaunayInterpolant )
    {
        std::cout << "\nComputing Delaunay Piecewise-Linear interpolant (order=2, degree=1) for field " << varProjectSrc
                  << std::endl;
        err = ComputeDelaunayInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }
    else if( rmethod == MultilevelBsplineApproximation )
    {
        if( this->dimension == 3 )
        {
            // constexpr int nlevels         = 7;
            // std::array< size_t, 3 > grid = { 16, 16, mpas_zlevels };
            // constexpr int nlevels        = 9;
            // std::array< size_t, 3 > grid = { 4, 4, mpas_zlevels / 4 };
            constexpr int nlevels        = 10;
            std::array< size_t, 3 > grid = { 2, 2, mpas_zlevels / 4 };
            std::array< double, 6 > bbox = { -1.0, -1.0, -1E6, 1.0, 1.0, 5E2 };
            // std::array< double, 6 > bbox = { -1.0, -1.0, -1.0, 1.0, 1.0, 1 };
            std::cout << "\nComputing MBA interpolant (order=4, degree=3) for field " << varProjectSrc << std::endl;
            err = ComputeMBAInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata, is_three_dimensional, order,
                                         grid, bbox, nlevels );MB_CHK_ERR( err );
        }
        else
        {
            constexpr int nlevels        = 10;
            std::array< size_t, 3 > grid = { 5, 5, 5 };
            std::array< double, 6 > bbox = { -1.0, -1.0, -1E2, 1.0, 1.0, 1E2 };
            std::cout << "\nComputing MBA interpolant (order=4, degree=3) for field " << varProjectSrc << std::endl;
            err = ComputeMBAInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata, is_three_dimensional, order,
                                         grid, bbox, nlevels );MB_CHK_ERR( err );
        }
    }
    else if( rmethod == ShepardInterpolant )
    {
        std::cout << "\nComputing Shepard interpolant (order=" << order << ") for field " << varProjectSrc << std::endl;
#ifdef MOAB_HAVE_ALGLIB
        err = ShepardInterpolatorAlgLib( this->dimension, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
#endif
        // err = modified_shepard_interpolate2( this->dimension, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
        err = ComputeHierarchicalShepardInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata, order );MB_CHK_ERR( err );
        // err = RBFAlgLib( *this, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
        //  err = ComputeGMLSInterpolant( *this, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }

    if( normalize && ( !is_three_dimensional ) && this->meshOverlap.faces.size() )
    {
        // Loop through all overlap-target faces and compute integral
        for( size_t j = 0; j < this->meshOverlap.faces.size(); j++ )
        {
            int iTargetFace = this->meshOverlap.vecTargetFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iTargetFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralOut += dst_tdata[iTargetFace] * this->meshOverlap.vecFaceArea[j];
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
