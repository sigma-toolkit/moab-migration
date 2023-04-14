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
#include "moab/ReadUtilIface.hpp"

#include "moab/Remapping/MBA.hpp"
#include "moab/nanoflann.hpp"

#ifndef MOAB_HAVE_MPI
#error " compile with MPI and HDF5 for this example to work \n";
#endif

#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"

// Remapping related includes
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/Remapping/TempestRemapper.hpp"

#include "spline.h"

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

moab::ErrorCode ShepardInterpolator( int dimension,
                                     std::vector< double >& xyzd,
                                     std::vector< double >& fd,
                                     std::vector< double >& xyzi,
                                     std::vector< double >& fi );

moab::ErrorCode ModifiedShepardInterpolator( int dimension,
                                             std::vector< double >& xyzd,
                                             std::vector< double >& fd,
                                             std::vector< double >& xyzi,
                                             std::vector< double >& fi );

moab::ErrorCode ComputeMBAInterpolant( std::vector< double >& xyzd,
                                       std::vector< double >& fd,
                                       std::vector< double >& xyzi,
                                       std::vector< double >& fi );

moab::ErrorCode ExtrudePolygonsToPolyhedra( Interface* mb,
                                            std::vector< double >& layer_thickness,
                                            moab::EntityHandle& poly2dset,
                                            moab::EntityHandle& outputset,
                                            bool is_mpas );

moab::ErrorCode ComputeFieldProjections( moab::Interface* mbi,
                                         Mesh& meshOverlap,
                                         std::string varProjectSrc,
                                         std::string varProjectDst,
                                         moab::Range& srcelems,
                                         moab::Range& dstelems,
                                         bool is_three_dimensional,
                                         bool normalize              = true,
                                         const double constantoffset = 0.0,
                                         bool useMBA                 = true,
                                         moab::Range* src3delems     = nullptr,
                                         moab::Range* dst3delems     = nullptr );

ErrorCode CloneToTRMesh( moab::Interface* m_interface, Mesh& mesh, EntityHandle mesh_set );

// 3D settings
constexpr int mpas_zlevels = 60;
constexpr int roms_zlevels = 5;
int src_zlayers = 0, dst_zlayers = 0;

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
    const double radius     = 1.0;
    bool ensureMonotonicity = false;
    bool computeWeights     = false;
    bool computeShepard     = false;
    bool computeMBA         = false;

    bool normalize          = false;
    bool use_3dprojection   = false;
    bool generateExtrusions = false;
    std::string strMethod   = "";
    const bool useTranspose = false;

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
                             &computeShepard );
        opts.addOpt< void >( "mba", "Use Shepard's inverse-distance weighting to compute projection", &computeMBA );
        opts.addOpt< void >( "mono", "Ensure monotonicity in the weight generation", &ensureMonotonicity );
        opts.addOpt< void >( "normalize", "Re-normalize interpolant to preserve field integral", &normalize );
        opts.addOpt< void >( "weights,w", "Compute and output the weights", &computeWeights );
        opts.addOpt< void >( "threed", "Compute full 3D temperature and salinity projection", &use_3dprojection );
        opts.addOpt< void >( "threed_setup", "Compute full 3D temperature and salinity projection",
                             &generateExtrusions );

        opts.parseCommandLine( argc, argv );
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
        ierr = MPI_Finalize();
        MPICHKERR( ierr, "MPI_Finalize failed; Aborting" );

        return 1;
    }

    // Initialize MPI first
    ierr = MPI_Init( &argc, &argv );
    MPICHKERR( ierr, "MPI_Init failed" );

    MPI_Comm comm = MPI_COMM_WORLD;

    ierr = MPI_Comm_rank( comm, &rank );
    MPICHKERR( ierr, "MPI_Comm_rank failed" );
    ierr = MPI_Comm_size( comm, &size );
    MPICHKERR( ierr, "MPI_Comm_size failed" );

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

    EntityHandle partnset;
    err = mbi->create_meshset( MESHSET_SET, partnset );MB_CHK_SET_ERR( err, "Creating partition set failed" );
    // Create the parallel communicator object with the partition handle associated with MOAB
    ParallelComm* parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );

    // construct the remapper
    moab::TempestRemapper remapper( mbi, parallel_communicator );
    remapper.meshValidate     = false;
    remapper.constructEdgeMap = false;
    remapper.initialize();

    EntityHandle mpasset, mpas_covering_set, romsset;
    err = mbi->create_meshset( moab::MESHSET_SET, mpasset );MB_CHK_SET_ERR( err, "Can't create new set" );
    err = mbi->create_meshset( moab::MESHSET_SET, mpas_covering_set );MB_CHK_SET_ERR( err, "Can't create new set" );
    err = mbi->create_meshset( moab::MESHSET_SET, romsset );MB_CHK_SET_ERR( err, "Can't create new set" );

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
        dbgprint( "MPAS mesh contains " << mpas_verts.size() << " vertices and " << mpas_elems.size() << " elements" );

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, mpas_verts, radius, !useTranspose );MB_CHK_ERR( err );
        err = mbi->write_file( "mpas_modified_2d.h5m", "H5M", write_options.c_str(), &mpasset, 1 );MB_CHK_ERR( err );
    }

    // Load the ROMS file from disk with given options
    moab::Range roms_verts, roms_elems, roms3d_verts, roms3d_elems;
    {
        dbgprint( "Reading ROMS file from disk" );
        err = mbi->load_file( roms_filename.c_str(), &romsset, roms_read_options.c_str() );MB_CHK_SET_ERR( err, "MOAB::load_file for ROMS mesh failed" );

        // Get all entities in the database
        err = mbi->get_entities_by_dimension( romsset, 0, roms_verts );MB_CHK_ERR( err );
        err = mbi->get_entities_by_dimension( romsset, 2, roms_elems );MB_CHK_ERR( err );
        dbgprint( "ROMS mesh contains " << roms_verts.size() << " vertices and " << roms_elems.size() << " elements" );

        // Rescale the radius of both to compute the intersection
        err = ScaleCoords( mbi, roms_verts, radius, useTranspose );MB_CHK_ERR( err );
        err = mbi->write_file( "roms_modified_2d.h5m", "H5M", write_options.c_str(), &romsset, 1 );MB_CHK_ERR( err );
        err = remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh );MB_CHK_ERR( err );
    }

    // Cull the MPAS set so that we don't have a global mesh
    moab::Range mpas_verts, mpas_elems, mpas3d_verts, mpas3d_elems;
    {
        const int nring_neighborhood = 3;
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
        for( size_t i = 0; i < roms_elems.size(); i++ )
        {
            const moab::EntityHandle ehandle = roms_elems[i];
            err                              = mbi->get_coords( &ehandle, 1, &query_pt[0] );MB_CHK_ERR( err );

            // Do a KNN search
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            for( size_t j = 0; j < num_results; ++j )
                mpas_elems.insert( orig_mpas_elems[srcindx[j]] );
        }
        err = mbi->add_entities( mpas_covering_set, mpas_elems );MB_CHK_ERR( err );

        err = mbi->write_file( "mpas_covering_2d.h5m", "H5M", write_options.c_str(), &mpas_covering_set, 1 );MB_CHK_ERR( err );

        err = mbi->get_connectivity( mpas_elems, mpas_verts, true );MB_CHK_ERR( err );
        err = mbi->add_entities( mpas_covering_set, mpas_verts );MB_CHK_ERR( err );
        dbgprint( "Culled MPAS mesh contains " << mpas_elems.size() << " elements and " << mpas_verts.size()
                                               << " vertices." );

        err = remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( err );
    }

    Mesh meshInput, meshOutput, meshOverlap;
    if( normalize || !( computeMBA || computeShepard ) )
    {
        CloneToTRMesh( mbi, meshInput, mpas_covering_set );
        CloneToTRMesh( mbi, meshOutput, romsset );

        meshInput.ConstructEdgeMap();
        meshOutput.ConstructEdgeMap();

        // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
        std::cout << "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes\n";
        // err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );
        bool concaveMeshA = false, concaveMeshB = false, allowNoOverlap = true, verbose = false;
        int ierr = GenerateOverlapWithMeshes( meshInput, meshOutput, meshOverlap, "" /*outFilename*/, "Netcdf4",
                                              "exact", concaveMeshA, concaveMeshB, allowNoOverlap, verbose );
        if( ierr )
        {
            MB_CHK_SET_ERR( MB_FAILURE, "TempestRemap: Can't compute the intersection of meshes on the sphere" );
        }
    }

    // let us perform 3D extrusions as needed
    EntityHandle mpasset3d = 0, romsset3d = 0;
    if( use_3dprojection )
    {
        if( generateExtrusions )
        {
            printf( "Getting MPAS height tag now\n" );
            moab::Tag mhtag;
            err = mbi->tag_get_handle( "layerThickness_3d", src_zlayers, moab::MB_TYPE_DOUBLE, mhtag,
                                       moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get tag handle: layerThickness_3d" );
            std::vector< double > zmh_xyz3d( mpas_elems.size() * src_zlayers ),
                zrh_xyz3d( roms_elems.size() * dst_zlayers );
            err = mbi->tag_get_data( mhtag, mpas_elems, zmh_xyz3d.data() );MB_CHK_SET_ERR( err, "Can't get layerThickness_3d data" );

            printf( "Extruding MPAS polyhedra now\n" );
            err = ExtrudePolygonsToPolyhedra( mbi, zmh_xyz3d, mpas_covering_set, mpasset3d, true );MB_CHK_SET_ERR( err, "Can't extrude MPAS polygons" );

            // Project the bottom Bathymetry data from MPAS to ROMS so that we can impose it.
            err = ComputeFieldProjections( mbi, meshOverlap, "bottomDepth", "bottomDepth", mpas_elems, roms_elems,
                                           false /*use_3dprojection*/, false /* bool normalize */, 2000.0, computeMBA );MB_CHK_SET_ERR( err, "Can't create new set" );

            std::vector< double > zrh_xyz2d( roms_elems.size() );
            moab::Tag rhtag;
            err = mbi->tag_get_handle( "bottomDepth", 1, moab::MB_TYPE_DOUBLE, rhtag, moab::MB_TAG_DENSE );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag" );
            err = mbi->tag_get_data( rhtag, roms_elems, zrh_xyz2d.data() );MB_CHK_SET_ERR( err, "Can't get bottomDepth tag data" );

            for( size_t i = 0; i < roms_elems.size(); ++i )
            {
                const double pbathymetry = zrh_xyz2d[i];
                const double delz        = pbathymetry / dst_zlayers;
                const int offset         = i * dst_zlayers;
                for( int j = 0; j < dst_zlayers; ++j )
                    zrh_xyz3d[offset + j] = delz;
                // printf( "Thickness value for ROMS element %d, %d = %f, %f\n", i, dst_zlayers, zrh_xyz2d[i] , delz );
            }

            printf( "Extruding ROMS polyhedra now\n" );
            err = ExtrudePolygonsToPolyhedra( mbi, zrh_xyz3d, romsset, romsset3d, false );MB_CHK_SET_ERR( err, "Can't extrude ROMS polygons" );

            err = mbi->get_entities_by_dimension( mpasset3d, 0, mpas3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( mpasset3d, 3, mpas3d_elems );MB_CHK_ERR( err );

            std::cout << "3D MPAS: " << mpas3d_verts.size() << " vertices and " << mpas3d_elems.size() << " elements.\n";
            err = mbi->write_file( "mpas_full_3d.h5m", "H5M", write_options.c_str(), &mpasset3d, 1 );MB_CHK_ERR( err );

            err = mbi->get_entities_by_dimension( romsset3d, 0, roms3d_verts );MB_CHK_ERR( err );
            err = mbi->get_entities_by_dimension( romsset3d, 3, roms3d_elems );MB_CHK_ERR( err );

            std::cout << "3D ROMS: " << roms3d_verts.size() << " vertices and " << roms3d_elems.size() << " elements.\n";
            err = mbi->write_file( "roms_full_3d.h5m", "H5M", write_options.c_str(), &romsset3d, 1 );MB_CHK_ERR( err );
        }
        else
        {
            mpas_verts.clear();
            roms_verts.clear();
            mpas_elems.clear();
            roms_elems.clear();

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
        }
    }

    if( use_3dprojection || computeMBA )
    {
        // Now let us compute the mba hierarchy for each field
        err = ComputeFieldProjections( mbi, meshOverlap, ( use_3dprojection ? "salinity_3d" : "salinity" ),
                                       ( use_3dprojection ? "salinity_3d_roms" : "salinity_roms" ), mpas_elems,
                                       roms_elems, use_3dprojection, true /* bool normalize */, 35.0, computeMBA, &mpas3d_elems, &roms3d_elems );MB_CHK_ERR( err );
        err = ComputeFieldProjections( mbi, meshOverlap, ( use_3dprojection ? "temperature_3d" : "temperature" ),
                                       ( use_3dprojection ? "temperature_3d_roms" : "temperature_roms" ), mpas_elems,
                                       roms_elems, use_3dprojection, true /* bool normalize */, 8.5, computeMBA,
                                       &mpas3d_elems, &roms3d_elems );
        MB_CHK_ERR( err );
    }
    else if( computeShepard )
    {
        // call Shepard's interpolant to compute data
        // compute the actual shepard's interpolation
        dbgprint( "Computing the Shepard's interpolant now" );

        // Now let us compute the mba hierarchy for each field
        err = ComputeFieldProjections( mbi, meshOverlap, ( use_3dprojection ? "salinity_3d" : "salinity" ),
                                       ( use_3dprojection ? "salinity_3d_roms" : "salinity_roms" ), mpas_elems,
                                       roms_elems, use_3dprojection, true /* bool normalize */, 35.0, computeMBA,
                                       &mpas3d_elems, &roms3d_elems );
        MB_CHK_ERR( err );
        err = ComputeFieldProjections( mbi, meshOverlap, ( use_3dprojection ? "temperature_3d" : "temperature" ),
                                       ( use_3dprojection ? "temperature_3d_roms" : "temperature_roms" ), mpas_elems,
                                       roms_elems, use_3dprojection, true /* bool normalize */, 8.5, computeMBA,
                                       &mpas3d_elems, &roms3d_elems );
        MB_CHK_ERR( err );
    }
    else
    {
        // compute the conservative 2D mapping weights
        {
            dbgprint( "\nSetup computation of weights" );
            // Call to generate the remapping weights with the tempest meshes
            OfflineMap weightMap;

            GenerateOfflineMapAlgorithmOptions mapOptions;
            mapOptions.nPin             = 1;
            mapOptions.nPout            = 1;
            mapOptions.fSourceConcave   = false;
            mapOptions.fTargetConcave   = false;
            mapOptions.strMethod        = strMethod;  // invdist, bilin
            mapOptions.fMonotone        = ensureMonotonicity;
            mapOptions.fNoCorrectAreas  = false;
            mapOptions.fNoCheck         = true;
            mapOptions.strOutputMapFile = output_filename;  // ask TR to write it out
            mapOptions.strOutputFormat  = "Netcdf4";

            dbgprint( "Compute weights with TempestRemap" );
            ierr = GenerateOfflineMapWithMeshes( meshInput,    // Mesh inputMesh
                                                 meshOutput,   // Mesh outputMesh,
                                                 meshOverlap,  // Mesh overlapMesh,
                                                 "fv",         // std::string inputDiscretization,
                                                 "fv",         // std::string outputDiscretization,
                                                 mapOptions,   // const GenerateOfflineMapAlgorithmOptions& options
                                                 weightMap );MB_CHK_ERR( err );

            // check the generated weights and output information
            {
                const double dNormalTolerance = 1.0E-8;
                const double dStrictTolerance = 1.0E-12;
                weightMap.CheckMap( true, true, ensureMonotonicity, dNormalTolerance, dStrictTolerance );
            }

            // Write the map to disk
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
                weightMap.Write( mapOptions.strOutputMapFile, mapAttributes, NcFile::Netcdf4Classic );

                // // Write the map file to disk in parallel using either HDF5 or SCRIP interface
                // err = weightMap.WriteParallelMap( output_filename.c_str() );MB_CHK_ERR( err );
            }

            // Now apply the map to compute the field projections on the ROMS mesh
            {
                // Now let us apply the weights onto the vector and project onto target mesh
                // err = weightMap.ApplyWeights( stag, stag, useTranspose );MB_CHK_ERR( err );
                // err = weightMap.ApplyWeights( ttag, ttag, useTranspose );MB_CHK_ERR( err );

                // get the handle to the weight matrix
                const SparseMatrix< double >& weights = weightMap.GetSparseMatrix();
                DataArray1D< double > dataInDouble( mpas_elems.size() );
                DataArray1D< double > dataOutDouble( roms_elems.size() );

                // assert( useConservativeSalinity );
                {
                    constexpr double salinity_avg = 35.0;
                    moab::Tag stag;
                    err = mbi->tag_get_handle( "salinity", 1, moab::MB_TYPE_DOUBLE, stag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

                    // Apply the map onto the salinity solution field
                    err = mbi->tag_get_data( stag, mpas_elems, dataInDouble );MB_CHK_ERR( err );
                    for( size_t i = 0; i < dataInDouble.GetRows(); i++ )
                        dataInDouble[i] -= salinity_avg;

                    // Compute the projection for the salinity field
                    weights.Apply( dataInDouble, dataOutDouble );

                    // Scale data values
                    for( size_t i = 0; i < dataOutDouble.GetRows(); i++ )
                        dataOutDouble[i] += salinity_avg;

                    // Set the salinity solution field on the ROMS mesh
                    err = mbi->tag_set_data( stag, roms_elems, dataOutDouble );MB_CHK_ERR( err );
                }

                // Apply the map onto the temperature solution field
                // assert( useConservativeTemperature );
                {
                    constexpr double temperature_avg = 8.5;
                    moab::Tag ttag;
                    err = mbi->tag_get_handle( "temperature", 1, moab::MB_TYPE_DOUBLE, ttag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

                    // Apply the map onto the salinity solution field
                    err = mbi->tag_get_data( ttag, mpas_elems, dataInDouble );MB_CHK_ERR( err );
                    for( size_t i = 0; i < dataInDouble.GetRows(); i++ )
                        dataInDouble[i] -= temperature_avg;

                    // Compute the projection for the salinity field
                    weights.Apply( dataInDouble, dataOutDouble );

                    // Scale data values
                    for( size_t i = 0; i < dataOutDouble.GetRows(); i++ )
                        dataOutDouble[i] += temperature_avg;

                    // Set the temperature solution field on the ROMS mesh
                    err = mbi->tag_set_data( ttag, roms_elems, dataOutDouble );MB_CHK_ERR( err );
                }
            }
        }

        remapper.clear();
    }

    dbgprint( "Writing out the target mesh with projected fields to 'roms_2d_projected.h5m'" );
    if( useTranspose )
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
            dbgprint( nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2]
                         << ": Failed with length == 0." );
            return MB_FAILURE;
        }

        // rescale to radius
        posf[0] *= R / len;
        posf[1] *= R / len;
        posf[2] *= R / len;
        rval = mb->set_coords( &nd, 1, posf );MB_CHK_ERR( rval );
    }
    return MB_SUCCESS;
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

moab::ErrorCode ComputeMBAInterpolant( std::vector< double >& xyzd,
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
    // const size_t init_grid_size = static_cast< size_t >( std::max( 10.0, std::sqrt( nd ) / 8 ) );
    // mba::index< 3 > grid        = { init_grid_size, init_grid_size, 2 };
    mba::index< 3 > grid = { 100, 100, 100 };

    std::vector< mba::point< 3 > > coords( nd );
    size_t offset = 0;
    for( size_t k = 0; k < nd; k++, offset += 3 )
        coords[k] = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };

    // Algorithm setup.
    mba::MBA< 3 > interp( lo, hi, grid, coords, fd, 8 /*levels*/, 1e-14 /*tolerance*/, 0.5 /*min_fill*/ );
    // mba::linear_approximation< 3 > interp( coords.begin(), coords.end(), fd.begin() );

    // Get interpolated value at arbitrary location.
    offset = 0;
    for( size_t k = 0; k < ni; k++, offset += 3 )
        fi[k] = interp( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );

    return moab::MB_SUCCESS;
}

moab::ErrorCode ModifiedShepardInterpolator( int dimension,
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
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PointCloud< double > >,
                                                        PointCloud< double >, 3 /* dim */
                                                        >;

    PointCloud< double > cloud( xyzd );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    double query_pt[3];  // dimension

    const size_t num_results =
        static_cast< size_t >( ( 2 * shepard_power + 1 ) * ( 2 * shepard_power + 1 ) - shepard_power );
    double ind = ( 1.0 / num_results );
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
                                         Mesh& meshOverlap,
                                         std::string varProjectSrc,
                                         std::string varProjectDst,
                                         moab::Range& srcelems,
                                         moab::Range& dstelems,
                                         bool is_three_dimensional,
                                         bool normalize,
                                         const double constantoffset,
                                         bool useMBA,
                                         moab::Range* src3delems,
                                         moab::Range* dst3delems )
{
    moab::ErrorCode err;
    moab::Tag dmtag;
    // err = mbi->tag_get_handle( varProject.c_str(), src_zlayers, moab::MB_TYPE_DOUBLE, dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
    err = mbi->tag_get_handle( varProjectSrc.c_str(), ( is_three_dimensional ? src_zlayers : 1 ), moab::MB_TYPE_DOUBLE,
                               dmtag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );

    // get the source data from tag
    std::vector< double > src_tdata( srcelems.size() * ( is_three_dimensional ? src_zlayers : 1 ) ),
        dst_tdata( dstelems.size() * ( is_three_dimensional ? dst_zlayers : 1 ) );
    err = mbi->tag_get_data( dmtag, srcelems, src_tdata.data() );MB_CHK_ERR( err );

    // get the coordinates of the elements
    std::vector< double > src_xyz( srcelems.size() * 3 ), dst_xyz( dstelems.size() * 3 );
    {
        err = mbi->get_coords( srcelems, src_xyz.data() );MB_CHK_ERR( err );
        err = mbi->get_coords( dstelems, dst_xyz.data() );MB_CHK_ERR( err );
    }

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0, normFactor = 1.0;
    if( normalize && !is_three_dimensional )
    {
        for( size_t j = 0; j < src_tdata.size(); j++ )
            src_tdata[j] -= constantoffset;

        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += src_tdata[iSourceFace] * meshOverlap.vecFaceArea[j];
        }
    }

    if( useMBA )
    {
        std::cout << "\nComputing the MBA interpolant now for field " << varProjectSrc << std::endl;
        err = ComputeMBAInterpolant( src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }
    else
    {
        std::cout << "\nComputing the Shepard(power=2) interpolant now for field " << varProjectSrc << std::endl;
        err = ModifiedShepardInterpolator( 3, src_xyz, src_tdata, dst_xyz, dst_tdata );MB_CHK_ERR( err );
    }

    if( normalize && !is_three_dimensional )
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
        {
            dst_tdata[ind] += constantoffset;
            dst_tdata[ind] *= normFactor;
        }
    }

    // now set the data on ROMS instance of MOAB tag
    std::cout << "Setting tag data to destination mesh\n";
    moab::Tag drtag;
    err = mbi->tag_get_handle( varProjectDst.c_str(), ( is_three_dimensional ? dst_zlayers : 1 ), moab::MB_TYPE_DOUBLE,
                               drtag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( err );
    err = mbi->tag_set_data( drtag, dstelems, dst_tdata.data() );MB_CHK_ERR( err );

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

    for( unsigned iface = 0; iface < elems.size(); ++iface )
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

    unsigned nnodes = verts.size();
    nodes.resize( nnodes );

    // Set the data for the vertices
    std::vector< double > coordx( nnodes ), coordy( nnodes ), coordz( nnodes );
    rval = m_interface->get_coords( verts, &coordx[0], &coordy[0], &coordz[0] );MB_CHK_ERR( rval );
    for( unsigned inode = 0; inode < nnodes; ++inode )
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
                                            bool is_mpas )
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
    std::cout << "Number of 2D entities: vertices = " << verts.size() << ", faces = " << faces.size() << endl;

    // add the initial faces to the first set
    rval = mb->add_entities( outputset, verts );MB_CHK_ERR( rval );
    // rval = mb->add_entities( outputset, faces );MB_CHK_ERR( rval );

    // Create all edges
    rval = mb->get_adjacencies( faces, 1, true, edges, Interface::UNION );MB_CHK_ERR( rval );

    if( true )
    {
        std::cout << "Input Polygonal Mesh details::" << std::endl;
        std::cout << "\tNumber of Vertices = " << verts.size() << std::endl;
        std::cout << "\t          Edges    = " << edges.size() << std::endl;
        std::cout << "\t          Faces    = " << faces.size() << std::endl;
    }

    const size_t nverts  = verts.size();
    const size_t nedges  = edges.size();
    const size_t nfaces  = faces.size();
    const size_t nlayers = layer_thickness.size() / nfaces;
    const size_t nquads  = nedges * nlayers;
    std::vector< double > coords( 3 * nverts );

    // get the vertex coordinates for the polygonal mesh
    rval = mb->get_coords( verts, &coords[0] );MB_CHK_ERR( rval );

    Tag gidTag = mb->globalId_tag();
    std::vector< int > gidData;
    std::vector< EntityHandle > vecents( verts.size() );
    std::copy( verts.begin(), verts.end(), std::back_inserter( vecents ) );

    Tag parentTag;
    rval = mb->tag_get_handle( "ColumnParent", 1, moab::MB_TYPE_HANDLE, parentTag,
                               moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( rval );
    std::vector< int > layerchildren;
    // create first vertices
    Range* newVerts = new Range[nlayers + 1];
    newVerts[0]     = verts;  // just for convenience
    for( size_t ii = 0; ii < nlayers; ii++ )
    {
        // for( size_t i = 0; i < nverts; i++ )
        //     coords[3 * i + 2] -= 0.5;  // constant offset for layers (testing)

        for( size_t i = 0; i < nverts; i++ )
        {
            Range eladjs;
            const EntityHandle vtx = verts[i];
            rval                   = mb->get_adjacencies( &vtx, 1, 2, true, eladjs, Interface::UNION );MB_CHK_ERR( rval );

            if( eladjs.size() )
            {
                double thickness = 0;
                for( size_t k = 0; k < eladjs.size(); ++k )
                {
                    int il = faces.index( eladjs[k] );
                    thickness += layer_thickness[il * nlayers + ii];
                }
                thickness /= eladjs.size();

                // printf( "Thickness value for layer %d: vertex %d = %f\n", ii, i, thickness );

                // Subtract or Add depending on the z-Direction to extrude
                coords[3 * i + 2] -= thickness;
            }
            else
                printf( "Vertex %zu has no adjacencies, coord = %f\n", i, coords[3 * i + 2] );
        }

        rval = mb->create_vertices( &coords[0], nverts, newVerts[ii + 1] );MB_CHK_ERR( rval );
        gidData.resize( nverts );
        std::iota( gidData.begin(), gidData.end(), nverts + ii * nverts );
        rval = mb->tag_set_data( gidTag, newVerts[ii + 1], gidData.data() );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( parentTag, newVerts[ii + 1], vecents.data() );MB_CHK_ERR( rval );

        rval = mb->add_entities( outputset, newVerts[ii + 1] );MB_CHK_ERR( rval );
    }

    EntityHandle start_elem, *connect;
    std::vector< EntityHandle > allPolygons;
    std::vector< int > allPolygonsGID;
    if( is_mpas )
    {
        // for each edge, we will create nlayers quads
        ReadUtilIface* read_iface;
        rval = mb->query_interface( read_iface );MB_CHK_SET_ERR( rval, "Error in query_interface" );

        // Create quads
        rval = read_iface->get_element_connect( nquads, 4, MBQUAD, 0, start_elem, connect );MB_CHK_SET_ERR( rval, "Error in get_element_connect" );
        Range quads( start_elem, start_elem + nquads );

        // ---------------------------------------------------------------------------
        int indexConn = 0;
        for( unsigned j = 0; j < nedges; j++ )
        {
            EntityHandle edge = edges[j];

            const EntityHandle* conn2 = NULL;
            int nnodes;
            rval = mb->get_connectivity( edge, conn2, nnodes );MB_CHK_ERR( rval );
            if( 2 != nnodes ) MB_CHK_ERR( MB_FAILURE );

            int i0 = verts.index( conn2[0] );
            int i1 = verts.index( conn2[1] );
            for( size_t ii = 0; ii < nlayers; ii++ )
            {
                connect[indexConn++] = newVerts[ii][i0];
                connect[indexConn++] = newVerts[ii][i1];
                connect[indexConn++] = newVerts[ii + 1][i1];
                connect[indexConn++] = newVerts[ii + 1][i0];
            }
        }

        // rval = mb->add_entities( outputset, quads );MB_CHK_ERR( rval );

        // std::vector< int > allPolygonsGID( nfaces * ( nlayers + 1 ) + nquads );
        allPolygonsGID.resize( nquads );
        // GIDS for lateral quads will be at the end of the list
        std::iota( allPolygonsGID.begin(), allPolygonsGID.end(), static_cast< int >( nfaces * ( nlayers + 1 ) ) + 1 );
        // TODO: this fails. Need to fix
        // rval = mb->tag_set_data( gidTag, quads, allPolygonsGID.data() );MB_CHK_ERR( rval );

        // next allocate for the x-y extruded faces
        allPolygons.resize( nfaces * ( nlayers + 1 ) );
        allPolygonsGID.resize( nfaces * ( nlayers + 1 ) );
        rval = mb->tag_get_data( gidTag, faces, allPolygonsGID.data() );MB_CHK_ERR( rval );
        for( unsigned int i = 0; i < nfaces; i++ )
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

    for( unsigned int j = 0; j < nfaces; j++ )
    {
        const EntityHandle polyg = faces[j];
        const EntityType etype   = mb->type_from_handle( polyg );

        // printf( "Polygon %d has %d nodes\n", j, nnodes );

        const EntityHandle* connp = nullptr;
        int nnodes;
        rval = mb->get_connectivity( polyg, connp, nnodes );MB_CHK_ERR( rval );

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

            vecents.resize( nnodes + 2, polyg );
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
        }
        else
        {
            std::copy( connp, connp + nnodes, vertexConn );
            for( int i = 0; i < nnodes; i++ )
                indexVerts[i] = verts.index( connp[i] );
        }

        for( unsigned ii = 0; ii < nlayers; ii++ )
        {
            // create a polygon on each layer
            for( int i = 0; i < nnodes; i++ )
                vertexConn[nnodes + i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1

            EntityHandle polyhedron;
            if( is_mpas )
            {
                rval = mb->create_element( etype, &vertexConn[nnodes], nnodes, allPolygons[nfaces * ( ii + 1 ) + j] );MB_CHK_ERR( rval );
                allPolygonsGID[nfaces * ( ii + 1 ) + j] = ipolygon++;

                // now create a polyhedra with top, bottom and lateral swept faces
                // first face is the bottom
                polyhedronConn[0] = allPolygons[nfaces * ii + j];
                // second face is the top
                polyhedronConn[1] = allPolygons[nfaces * ( ii + 1 ) + j];
                // next add lateral quads, in order of edges, using the start_elem
                // first layer of quads has EntityHandle from start_elem to start_elem+nedges-1
                // second layer of quads has EntityHandle from start_elem + nedges to
                // start_elem+2*nedges-1 ,etc
                for( int i = 0; i < nnodes; i++ )
                {
                    polyhedronConn[2 + i] = start_elem + ii + nlayers * indexEdges[i];
                }

                rval = mb->add_entities( outputset, polyhedronConn, nnodes + 2 );MB_CHK_ERR( rval );

                // Create polyhedron
                rval = mb->create_element( MBPOLYHEDRON, polyhedronConn, 2 + nnodes, polyhedron );MB_CHK_ERR( rval );

                rval = mb->tag_set_data( parentTag, polyhedronConn, nnodes + 2, vecents.data() );MB_CHK_ERR( rval );
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

            rval = mb->tag_set_data( parentTag, &polyhedron, 1, &polyg );MB_CHK_ERR( rval );
        }
    }

    // add subsequent layers of vertices
    if( is_mpas )
    {
        // rval = mb->add_entities( outputset, allPolygons.data(), allPolygons.size() );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( gidTag, allPolygons.data(), allPolygons.size(), allPolygonsGID.data() );MB_CHK_ERR( rval );
    }
    // add to the second neumann set the last layer of polygons
    // std::vector<EntityHandle> allPolygons(nfaces*(nlayers+1)) are all polygons
    // the last layer starts at allPolygons[nfaces*nlayers], size nfaces
    // rval = mb->add_entities( outputset, &allPolygons[nfaces * nlayers], nfaces );MB_CHK_ERR( rval );

    // give global ids to polyhedra, they will be used for side sets;
    // the order will be the same as the order in the file, by ranges of polyhedra
    // Range polyhs;
    // rval = mb->get_entities_by_type( outputset, MBPOLYHEDRON, polyhs );MB_CHK_ERR( rval );
    // std::vector< int > gidVals( polyhs.size() );
    // std::iota( gidVals.begin(), gidVals.end(), 1 );
    // rval = mb->tag_set_data( gidTag, polyhs, gidVals.data() );MB_CHK_ERR( rval );

    // rval = mb->write_file( std::string( prefix + "_polygon_3d.h5m" ).c_str(), "H5M", "DEBUG_IO=5;KEEP", &outputset,
    //                        1 );
    // rval = mb->write_file( std::string( prefix + "_polygon_3d.vtk" ).c_str(), "VTK", "", &outputset, 1 );MB_CHK_ERR( rval );

    return moab::MB_SUCCESS;
}
