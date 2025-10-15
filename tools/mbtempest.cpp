/*
 * Usage: MOAB-Tempest tool
 *
 * Generate a Cubed-Sphere mesh: ./mbtempest -t 0 -res 25 -f cubed_sphere_mesh.exo
 * Generate a RLL mesh: ./mbtempest -t 1 -res 25 -f rll_mesh.exo
 * Generate a Icosahedral-Sphere mesh: ./mbtempest -t 2 -res 25 <-dual> -f icosa_mesh.exo
 *
 * Now you can compute the intersections between the meshes too!
 *
 * Generate the overlap mesh: ./mbtempest -t 3 -l cubed_sphere_mesh.exo -l rll_mesh.exo -f
 * overlap_mesh.exo
 *
 */

// standard C++ includes
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

// MOAB includes
#include "moab/Core.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "DebugOutput.hpp"

//#ifndef MOAB_HAVE_MPI
//    #error mbtempest tool requires MPI configuration
//#endif

#ifdef MOAB_HAVE_MPI
// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

struct ToolContext
{
    moab::Core* mbcore;
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm;
#endif
    const int proc_id, n_procs;                  // MPI process id and number of processes
    moab::DebugOutput outputFormatter;           // Output formatter
    int blockSize;                               // Resolution of the mesh
    int nlayers;                                 // Number of ghost layers in the source mesh
    std::vector< std::string > inFilenames;      // Input mesh filenames
    std::vector< Mesh* > meshes;                 // Mesh references in TempestRemap format
    std::vector< moab::EntityHandle > meshsets;  // Meshsets in MOAB format
    std::vector< int > disc_orders;              // Discretization orders
    std::vector< std::string > disc_methods;     // Discretization method names (fv, cgll, dgll, etc.)
    std::vector< std::string > doftag_names;     // DoF tag names (GLOBAL_ID, etc.)
    std::string fvMethod;          // FV method name (invdist, delaunay, bilin, intbilin, intbilingb, none)
    std::string outFilename;       // Output filename
    std::string intxFilename;      // Output intersection mesh filename
    std::string baselineFile;      // Output baseline file (for verification)
    std::string variableToVerify;  // Variable name for verification
    moab::TempestRemapper::TempestMeshType
        meshType;  // Mesh type (CS, RLL, ICO, OVERLAP_FILES, OVERLAP_MEMORY, OVERLAP_MOAB)
    bool skip_io;  // Skip I/O operations; strictly to measure the performance of the intersection/map/SpMV algorithms
    bool computeDual;        // Compute the dual of the mesh
    bool computeWeights;     // Compute the projection weights
    bool verifyWeights;      // Verify the projection weights
    bool enforceConvexity;   // Enforce convexity of the input meshes
    int ensureMonotonicity;  // Ensure monotonicity in the weight generation process
    bool rrmGrids;      // Flag specifying that we are dealing with regionally refined grids (possibly nonoverlapping)
    bool kdtreeSearch;  // Use Kd-tree based search for computing mesh intersections instead of advancing front
    bool fCheck;        // Check the generated map for conservation and consistency
    bool fVolumetric;   // Apply a volumetric projection to compute the weights
    bool useGnomonicProjection;                     // Use Gnomonic plane projections to compute coverage mesh
    moab::TempestOnlineMap::CAASType cassType;      // CAAS filter type
    GenerateOfflineMapAlgorithmOptions mapOptions;  // Offline map options
    bool print_diagnostics;                         // Print diagnostics
    double boxeps;                                  // Box error tolerance
    double epsrel;                                  // Relative error tolerance

#ifdef MOAB_HAVE_MPI
    ToolContext( moab::Core* icore, moab::ParallelComm* p_pcomm )
        : mbcore( icore ), pcomm( p_pcomm ), proc_id( pcomm->rank() ), n_procs( pcomm->size() ),
          outputFormatter( std::cout, pcomm->rank(), 0 ),
#else
    ToolContext( moab::Core* icore )
        : mbcore( icore ), proc_id( 0 ), n_procs( 1 ), outputFormatter( std::cout, 0, 0 ),
#endif
          blockSize( 5 ), nlayers( 0 ), fvMethod( "none" ), outFilename( "outputFile.nc" ), intxFilename( "" ),
          baselineFile( "" ), variableToVerify( "" ), meshType( moab::TempestRemapper::DEFAULT ), skip_io( false ),
          computeDual( false ), computeWeights( false ), verifyWeights( false ), enforceConvexity( false ),
          ensureMonotonicity( 0 ), rrmGrids( false ), kdtreeSearch( true ), fCheck( false ), fVolumetric( false ),
          useGnomonicProjection( false ), cassType( moab::TempestOnlineMap::CAAS_NONE ), print_diagnostics( false ),
          boxeps( 1e-7 ),               // Box error tolerance default value
          epsrel( ReferenceTolerance )  // ReferenceTolerance is defined in Defines.h in TempestRemap
    {
        inFilenames.resize( 2 );
        doftag_names.resize( 2 );
        timer = new moab::CpuTimer();

        outputFormatter.set_prefix( "[MBTempest]: " );
    }

    ~ToolContext()
    {
        // for (unsigned i=0; i < meshes.size(); ++i) delete meshes[i];
        meshes.clear();
        inFilenames.clear();
        disc_orders.clear();
        disc_methods.clear();
        doftag_names.clear();
        outFilename.clear();
        intxFilename.clear();
        baselineFile.clear();
        variableToVerify.clear();
        meshsets.clear();
        delete timer;
    }

    void timer_push( std::string operation )
    {
        timer_ops = timer->time_since_birth();
        opName    = operation;
    }

    void timer_pop()
    {
        double locElapsed = timer->time_since_birth() - timer_ops, avgElapsed = 0, maxElapsed = 0;
#ifdef MOAB_HAVE_MPI
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, pcomm->comm() );
        MPI_Reduce( &locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, pcomm->comm() );
#else
        maxElapsed = locElapsed;
        avgElapsed = locElapsed;
#endif
        if( !proc_id )
        {
            avgElapsed /= n_procs;
            std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed
                      << ", avg = " << avgElapsed << "\n";
        }
        // std::cout << "\n[LOG" << proc_id << "] Time taken to " << opName << " = " <<
        // timer->time_since_birth() - timer_ops << std::endl;
        opName.clear();
    }

    void ParseCLOptions( int argc, char* argv[] )
    {
        ProgOptions opts;
        int imeshType                  = 0;
        std::string expectedFName      = "output.exo";
        std::string expectedMethod     = "fv";
        std::string expectedFVMethod   = "none";
        std::string expectedDofTagName = "GLOBAL_ID";
        int expectedOrder              = 1;
        int useCAAS                    = 0;
        int nlayer_input               = 0;

        if( !proc_id )
        {
            std::cout << "Command line options provided to mbtempest:\n  ";
            for( int iarg = 0; iarg < argc; ++iarg )
                std::cout << argv[iarg] << " ";
            std::cout << std::endl << std::endl;
        }

        opts.addOpt< int >( "type,t",
                            "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP_FILES=3, "
                            "OVERLAP_MEMORY=4, OVERLAP_MOAB=5])",
                            &imeshType );

        opts.addOpt< int >( "res,r", "Resolution of the mesh (default=5)", &blockSize );

        opts.addOpt< void >( "dual,d", "Output the dual of the mesh (relevant only for ICO mesh type)", &computeDual );

        opts.addOpt< std::string >( "file,f", "Output computed mesh or remapping weights to specified filename",
                                    &outFilename );
        opts.addOpt< std::string >(
            "load,l", "Input mesh filenames for source and target meshes. (relevant only when computing weights)",
            &expectedFName );

        opts.addOpt< void >( "advfront,a",
                             "Use the advancing front intersection instead of the Kd-tree based algorithm "
                             "to compute mesh intersections. (relevant only when computing weights)" );

        opts.addOpt< std::string >( "intx,i", "Output TempestRemap intersection mesh filename", &intxFilename );

        opts.addOpt< void >( "weights,w",
                             "Compute and output the weights using the overlap mesh (generally "
                             "relevant only for OVERLAP mesh)",
                             &computeWeights );

        opts.addOpt< void >(
            "verbose,v", "Print verbose diagnostic messages during intersection and map computation (default=false)",
            &print_diagnostics );

        opts.addOpt< std::string >( "method,m", "Discretization method for the source and target solution fields",
                                    &expectedMethod );

        opts.addOpt< int >( "order,o", "Discretization orders for the source and target solution fields",
                            &expectedOrder );

        opts.addOpt< std::string >( "global_id,g",
                                    "Tag name that contains the global DoF IDs for source and target solution fields",
                                    &expectedDofTagName );

        opts.addOpt< std::string >( "fvmethod",
                                    "Sub-type method for FV-FV projections (invdist, delaunay, bilin, "
                                    "intbilin, intbilingb, none. Default: none)",
                                    &expectedFVMethod );

        opts.addOpt< void >( "noconserve",
                             "Do not apply conservation to the resultant weights (relevant only "
                             "when computing weights)",
                             &mapOptions.fNoConservation );

        opts.addOpt< void >( "volumetric",
                             "Apply a volumetric projection to compute the weights (relevant only "
                             "when computing weights)",
                             &fVolumetric );

        opts.addOpt< void >( "skip_output", "For performance studies, skip all I/O operations.", &skip_io );

        opts.addOpt< void >( "gnomonic", "Use Gnomonic plane projections to compute coverage mesh.",
                             &useGnomonicProjection );

        opts.addOpt< void >( "enforce_convexity", "check convexity of input meshes to compute mesh intersections",
                             &enforceConvexity );

        opts.addOpt< void >( "nobubble", "do not use bubble on interior of spectral element nodes",
                             &mapOptions.fNoBubble );

        opts.addOpt< void >(
            "sparseconstraints",
            "Use sparse solver for constraints when we have high-valence (typical with high-res RLL mesh)",
            &mapOptions.fSparseConstraints );

        opts.addOpt< void >( "rrmgrids",
                             "At least one of the meshes is a regionally refined grid (relevant to "
                             "accelerate intersection computation)",
                             &rrmGrids );

        opts.addOpt< void >( "checkmap", "Check the generated map for conservation and consistency", &fCheck );

        opts.addOpt< void >( "verify",
                             "Verify the accuracy of the maps by projecting analytical functions "
                             "from source to target "
                             "grid by applying the maps",
                             &verifyWeights );

        opts.addOpt< std::string >( "var",
                                    "Tag name of the variable to use in the verification study "
                                    "(error metrics for user defined variables may not be available)",
                                    &variableToVerify );

        opts.addOpt< int >( "monotonicity", "Ensure monotonicity in the weight generation. Options=[0,1,2,3]",
                            &ensureMonotonicity );

        opts.addOpt< int >( "ghost", "Number of ghost layers in coverage mesh (default=1)", &nlayer_input );

        opts.addOpt< double >( "boxeps", "The tolerance for boxes (default=1e-7)", &boxeps );

        opts.addOpt< int >( "caas", "apply CAAS nonlinear filter after linear map application", &useCAAS );

        opts.addOpt< std::string >( "baseline", "Output baseline file", &baselineFile );

        opts.parseCommandLine( argc, argv );

        // By default - use Kd-tree based search; if user asks for advancing front, disable Kd-tree
        // algorithm
        kdtreeSearch = opts.numOptSet( "advfront,a" ) == 0;

        switch( imeshType )
        {
            case 0:
                meshType = moab::TempestRemapper::CS;
                break;

            case 1:
                meshType = moab::TempestRemapper::RLL;
                break;

            case 2:
                meshType = moab::TempestRemapper::ICO;
                break;

            case 3:
                meshType = moab::TempestRemapper::OVERLAP_FILES;
                break;

            case 4:
                meshType = moab::TempestRemapper::OVERLAP_MEMORY;
                break;

            case 5:
                meshType = moab::TempestRemapper::OVERLAP_MOAB;
                break;

            default:
                meshType = moab::TempestRemapper::DEFAULT;
                break;
        }

        // decipher whether we want to use CAAS filter when applying the projection
        switch( useCAAS )
        {
            case moab::TempestOnlineMap::CAAS_GLOBAL:
                cassType = moab::TempestOnlineMap::CAAS_GLOBAL;
                break;
            case moab::TempestOnlineMap::CAAS_LOCAL:
                cassType = moab::TempestOnlineMap::CAAS_LOCAL;
                break;
            case moab::TempestOnlineMap::CAAS_LOCAL_ADJACENT:
                cassType = moab::TempestOnlineMap::CAAS_LOCAL_ADJACENT;
                break;
            case moab::TempestOnlineMap::CAAS_QLT:
                cassType = moab::TempestOnlineMap::CAAS_QLT;
                break;
            default:
                cassType = moab::TempestOnlineMap::CAAS_NONE;
                break;
        }

        if( meshType > moab::TempestRemapper::ICO )  // compute overlap mesh and maps possibly
        {
            opts.getOptAllArgs( "load,l", inFilenames );
            opts.getOptAllArgs( "order,o", disc_orders );
            opts.getOptAllArgs( "method,m", disc_methods );
            opts.getOptAllArgs( "global_id,i", doftag_names );

            assert( inFilenames.size() == 2 );
            assert( ensureMonotonicity >= 0 && ensureMonotonicity <= 3 );

            // get discretization order parameters
            if( disc_orders.size() == 0 ) disc_orders.resize( 2, 1 );
            if( disc_orders.size() == 1 ) disc_orders.push_back( disc_orders[0] );

            // get discretization method parameters
            if( disc_methods.size() == 0 ) disc_methods.resize( 2, "fv" );
            if( disc_methods.size() == 1 ) disc_methods.push_back( disc_methods[0] );

            // get DoF tagname parameters
            if( doftag_names.size() == 0 ) doftag_names.resize( 2, "GLOBAL_ID" );
            if( doftag_names.size() == 1 ) doftag_names.push_back( doftag_names[0] );

            assert( disc_orders.size() == 2 );
            assert( disc_methods.size() == 2 );
            assert( doftag_names.size() == 2 );

            // for computing maps and overlaps, set discretization orders
            mapOptions.nPin           = disc_orders[0];
            mapOptions.nPout          = disc_orders[1];
            mapOptions.fSourceConcave = false;
            mapOptions.fTargetConcave = false;

            mapOptions.strMethod = "";

            if( expectedFVMethod != "none" )
            {
                mapOptions.strMethod += expectedFVMethod + ";";
                fvMethod = expectedFVMethod;

                // These FV projection methods are non-conservative; specify it explicitly
                mapOptions.fNoConservation = true;
            }
            switch( ensureMonotonicity )
            {
                case 0:
                    mapOptions.fMonotone = false;
                    break;
                case 3:
                    mapOptions.strMethod += "mono3;";
                    break;
                case 2:
                    mapOptions.strMethod += "mono2;";
                    break;
                default:
                    mapOptions.fMonotone = true;
            }
            mapOptions.fNoCorrectAreas = false;
            mapOptions.fNoCheck        = !fCheck;

            //assert( fVolumetric && fInverseDistanceMap == false );  // both options cannot be active
            if( fVolumetric ) mapOptions.strMethod += "volumetric;";

            // For global meshes, this default should work out of the box.
            if( !fvMethod.compare( "bilin" ) )
                nlayers = 3;
            else
                nlayers = ( mapOptions.nPin > 1 ? mapOptions.nPin + 1 : 0 );
            if( nlayer_input ) nlayers = std::max( nlayer_input, nlayers );
        }

        // clear temporary string name
        expectedFName.clear();

        mapOptions.strOutputMapFile = outFilename;
        mapOptions.strOutputFormat  = "Netcdf4";
    }

  private:
    moab::CpuTimer* timer;
    double timer_ops;
    std::string opName;
};

// Forward declare some methods
static moab::ErrorCode CreateTempestMesh( ToolContext&, moab::TempestRemapper& remapper, Mesh* );
inline double sample_slow_harmonic( double dLon, double dLat );
inline double sample_fast_harmonic( double dLon, double dLat );
inline double sample_constant( double dLon, double dLat );
inline double sample_stationary_vortex( double dLon, double dLat );

std::string get_file_read_options( ToolContext& ctx, std::string filename )
{
    // if running in serial, blank option may suffice in general
    std::string opts = "";
    if( ctx.n_procs > 1 )
    {
        size_t lastindex      = filename.find_last_of( "." );
        std::string extension = filename.substr( lastindex + 1, filename.size() );
        if( extension == "h5m" ) return "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        //    "PARALLEL_GHOSTS=2.0.3;PARALLEL_THIN_GHOST_LAYER;";
        else if( extension == "nc" )
        {
            // set default set of options
            opts = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;";

            if( ctx.proc_id == 0 )
            {
                {
                    NcFile ncFile( filename.c_str(), NcFile::ReadOnly );
                    if( !ncFile.is_valid() )
                        _EXCEPTION1( "Unable to open grid file \"%s\" for reading", filename.c_str() );

                    // Check for dimension names "grid_size", "grid_rank" and "grid_corners"
                    int iSCRIPFormat = 0, iMPASFormat = 0, iESMFFormat = 0;
                    for( int i = 0; i < ncFile.num_dims(); i++ )
                    {
                        std::string strDimName = ncFile.get_dim( i )->name();
                        if( strDimName == "grid_size" || strDimName == "grid_corners" || strDimName == "grid_rank" )
                            iSCRIPFormat++;
                        if( strDimName == "nodeCount" || strDimName == "elementCount" ||
                            strDimName == "maxNodePElement" )
                            iESMFFormat++;
                        if( strDimName == "nCells" || strDimName == "nEdges" || strDimName == "nVertices" ||
                            strDimName == "vertexDegree" )
                            iMPASFormat++;
                    }
                    ncFile.close();

                    if( iESMFFormat == 3 )  // Input from a NetCDF ESMF file
                    {
                        opts = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=;";
                    }
                    if( iSCRIPFormat == 3 )  // Input from a NetCDF SCRIP file
                    {
                        opts = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;";
                    }
                    if( iMPASFormat == 4 )  // Input from a NetCDF SCRIP file
                    {
                        opts = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;"
                               "PARALLEL_RESOLVE_SHARED_ENTS;NO_EDGES;NO_MIXED_ELEMENTS;VARIABLE=;";
                    }
                }
            }

#ifdef MOAB_HAVE_MPI
            int line_size = opts.size();
            MPI_Bcast( &line_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
            if( ctx.proc_id != 0 ) opts.resize( line_size );
            MPI_Bcast( const_cast< char* >( opts.data() ), line_size, MPI_CHAR, 0, MPI_COMM_WORLD );
#endif
        }
        else
            return "PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;";
    }
    return opts;
}
//#define MOAB_DBG
int main( int argc, char* argv[] )
{
    NcError error( NcError::verbose_nonfatal );
    std::stringstream sstr;
    std::string historyStr;

    int proc_id = 0, nprocs = 1;
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
#endif

    moab::Core* mbCore = new( std::nothrow ) moab::Core;

    if( nullptr == mbCore )
    {
        return 1;
    }

    // Build the history string
    for( int ia = 0; ia < argc; ++ia )
        historyStr += std::string( argv[ia] ) + " ";

    ToolContext* runCtx;
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm = new moab::ParallelComm( mbCore, MPI_COMM_WORLD, 0 );

    runCtx                   = new ToolContext( mbCore, pcomm );
    const char* writeOptions = ( nprocs > 1 ? "PARALLEL=WRITE_PART" : "" );
#else
    runCtx                   = new ToolContext( mbCore );
    const char* writeOptions = "";
#endif
    runCtx->ParseCLOptions( argc, argv );

    const double radius_src  = 1.0 /*2.0*acos(-1.0)*/;
    const double radius_dest = 1.0 /*2.0*acos(-1.0)*/;

    moab::DebugOutput& outputFormatter = runCtx->outputFormatter;

#ifdef MOAB_HAVE_MPI
    moab::TempestRemapper remapper( mbCore, pcomm );
#else
    moab::TempestRemapper remapper( mbCore );
#endif
    remapper.meshValidate     = true;
    remapper.constructEdgeMap = true;
    remapper.initialize();

    // Default area_method = lHuiller; Options: Girard, lHuiller, GaussQuadrature (if TR is available)
    moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::GaussQuadrature );

    Mesh* tempest_mesh = new Mesh();
    MB_CHK_ERR( CreateTempestMesh( *runCtx, remapper, tempest_mesh ) );

    if( runCtx->meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Compute intersections with MOAB
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert( runCtx->meshes.size() == 3 );

#ifdef MOAB_HAVE_MPI
        MB_CHK_ERR( pcomm->check_all_shared_handles() );
#endif

        // Load the meshes and validate
        MB_CHK_ERR( remapper.ConvertTempestMesh( moab::Remapper::SourceMesh ) );
        MB_CHK_ERR( remapper.ConvertTempestMesh( moab::Remapper::TargetMesh ) );
        MB_CHK_ERR( remapper.ConvertTempestMesh( moab::Remapper::OverlapMesh ) );
        if( !runCtx->skip_io )
        {
            MB_CHK_ERR( mbCore->write_mesh( "tempest_intersection.h5m", &runCtx->meshsets[2], 1 ) );
        }

        // print verbosely about the problem setting
        size_t velist[6], gvelist[6];
        {
            moab::Range rintxverts, rintxelems;
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 0, rintxverts ) );
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 2, rintxelems ) );
            velist[0] = rintxverts.size();
            velist[1] = rintxelems.size();

            moab::Range bintxverts, bintxelems;
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 0, bintxverts ) );
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, bintxelems ) );
            velist[2] = bintxverts.size();
            velist[3] = bintxelems.size();
        }

        moab::EntityHandle intxset;  // == remapper.GetMeshSet(moab::Remapper::OverlapMesh);

        // Compute intersections with MOAB
        {
            // Create the intersection on the sphere object
            runCtx->timer_push( "setup the intersector" );

            moab::Intx2MeshOnSphere* mbintx = new moab::Intx2MeshOnSphere( mbCore );
            mbintx->set_error_tolerance( runCtx->epsrel );
            mbintx->set_box_error( runCtx->boxeps );
            mbintx->set_radius_source_mesh( radius_src );
            mbintx->set_radius_destination_mesh( radius_dest );
#ifdef MOAB_HAVE_MPI
            mbintx->set_parallel_comm( pcomm );
#endif
            MB_CHK_ERR( mbintx->FindMaxEdges( runCtx->meshsets[0], runCtx->meshsets[1] ) );

#ifdef MOAB_HAVE_MPI
            moab::Range local_verts;
            MB_CHK_ERR( mbintx->build_processor_euler_boxes( runCtx->meshsets[1], local_verts ) );

            runCtx->timer_pop();

            moab::EntityHandle covering_set;
            runCtx->timer_push( "communicate the mesh" );
            // TODO:: needs clarification
            // we compute just intersection here, no need for extra ghost layers anyway
            // ghost layers are needed in coverage for bilinear map, which does not actually need intersection
            // this will be fixed in the future, bilinear map needs just coverage, not intersection
            // so I am not passing the ghost layer here, even though there is an option in runCtx for a ghost layer
            MB_CHK_ERR( mbintx->construct_covering_set(
                runCtx->meshsets[0], covering_set ) );  // lots of communication if mesh is distributed very differently
            runCtx->timer_pop();

            // print verbosely about the problem setting
            {
                moab::Range cintxverts, cintxelems;
                MB_CHK_ERR( mbCore->get_entities_by_dimension( covering_set, 0, cintxverts ) );
                MB_CHK_ERR( mbCore->get_entities_by_dimension( covering_set, 2, cintxelems ) );
                velist[4] = cintxverts.size();
                velist[5] = cintxelems.size();
            }

            MPI_Reduce( velist, gvelist, 6, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD );

#else
            moab::EntityHandle covering_set = runCtx->meshsets[0];
            for( int i = 0; i < 6; i++ )
                gvelist[i] = velist[i];
#endif

            if( !proc_id )
            {
                outputFormatter.printf( 0, "The source set contains %lu vertices and %lu elements \n", gvelist[0],
                                        gvelist[0] );
                outputFormatter.printf( 0, "The covering set contains %lu vertices and %lu elements \n", gvelist[2],
                                        gvelist[2] );
                outputFormatter.printf( 0, "The target set contains %lu vertices and %lu elements \n", gvelist[1],
                                        gvelist[1] );
            }

            // Now let's invoke the MOAB intersection algorithm in parallel with a
            // source and target mesh set representing two different decompositions
            runCtx->timer_push( "compute intersections with MOAB" );
            MB_CHK_SET_ERR( mbCore->create_meshset( moab::MESHSET_SET, intxset ), "Can't create new set" );
            MB_CHK_SET_ERR( mbintx->intersect_meshes( covering_set, runCtx->meshsets[1], intxset ),
                            "Can't compute the intersection of meshes on the sphere" );
            runCtx->timer_pop();

            // free the memory
            delete mbintx;
        }

        {
            moab::Range intxelems, intxverts;
            MB_CHK_ERR( mbCore->get_entities_by_dimension( intxset, 2, intxelems ) );
            MB_CHK_ERR( mbCore->get_entities_by_dimension( intxset, 0, intxverts, true ) );
            outputFormatter.printf( 0, "The intersection set contains %lu elements and %lu vertices \n",
                                    intxelems.size(), intxverts.size() );

            moab::IntxAreaUtils areaAdaptorHuiller( moab::IntxAreaUtils::lHuiller );  // lHuiller, GaussQuadrature
            double initial_sarea =
                areaAdaptorHuiller.area_on_sphere( mbCore, runCtx->meshsets[0],
                                                   radius_src );  // use the target to compute the initial area
            double initial_tarea =
                areaAdaptorHuiller.area_on_sphere( mbCore, runCtx->meshsets[1],
                                                   radius_dest );  // use the target to compute the initial area
            double intx_area = areaAdaptorHuiller.area_on_sphere( mbCore, intxset, radius_src );

            outputFormatter.printf( 0, "mesh areas: source = %12.10f, target = %12.10f, intersection = %12.10f \n",
                                    initial_sarea, initial_tarea, intx_area );
            outputFormatter.printf( 0, "relative error w.r.t source = %12.10e, target = %12.10e \n",
                                    fabs( intx_area - initial_sarea ) / initial_sarea,
                                    fabs( intx_area - initial_tarea ) / initial_tarea );
        }

        // Write out our computed intersection file
        if( !runCtx->skip_io )
        {
            MB_CHK_ERR( mbCore->write_mesh( "moab_intersection.h5m", &intxset, 1 ) );
        }

        if( runCtx->computeWeights )
        {
            runCtx->timer_push( "compute weights with the Tempest meshes" );
            // Call to generate an offline map with the tempest meshes
            OfflineMap weightMap;
            int err = GenerateOfflineMapWithMeshes( *runCtx->meshes[0], *runCtx->meshes[1], *runCtx->meshes[2],
                                                    runCtx->disc_methods[0],  // std::string strInputType
                                                    runCtx->disc_methods[1],  // std::string strOutputType,
                                                    runCtx->mapOptions, weightMap );
            runCtx->timer_pop();

            std::map< std::string, std::string > mapAttributes;
            if( err )
            {
                MB_SET_ERR( moab::MB_FAILURE, "Generating offline map with meshes failed." );
            }
            else
            {
                if( !runCtx->skip_io ) weightMap.Write( "outWeights.nc", mapAttributes );
            }
        }
    }
    else if( runCtx->meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m
#ifdef MOAB_HAVE_MPI
        MB_CHK_ERR( pcomm->check_all_shared_handles() );
#endif

        // print verbosely about the problem setting
        size_t velist[4] = {}, gvelist[4] = {};
        {
            moab::Range srcverts, srcelems;
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 0, srcverts ) );
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 2, srcelems ) );
            MB_CHK_ERR( moab::IntxUtils::fix_degenerate_quads( mbCore, runCtx->meshsets[0] ) );
            if( runCtx->enforceConvexity )
            {
                MB_CHK_ERR( moab::IntxUtils::enforce_convexity( mbCore, runCtx->meshsets[0], proc_id ) );
            }
            MB_CHK_ERR( areaAdaptor.positive_orientation( mbCore, runCtx->meshsets[0], radius_src ) );
            // if( !proc_id )
            //     outputFormatter.printf( 0, "The source set contains %lu vertices and %lu elements \n", srcverts.size(),
            //                             srcelems.size() );
            velist[0] = srcverts.size();
            velist[1] = srcelems.size();

            moab::Range tgtverts, tgtelems;
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 0, tgtverts ) );
            MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, tgtelems ) );
            MB_CHK_ERR( moab::IntxUtils::fix_degenerate_quads( mbCore, runCtx->meshsets[1] ) );
            if( runCtx->enforceConvexity )
            {
                MB_CHK_ERR( moab::IntxUtils::enforce_convexity( mbCore, runCtx->meshsets[1], proc_id ) );
            }
            MB_CHK_ERR( areaAdaptor.positive_orientation( mbCore, runCtx->meshsets[1], radius_dest ) );
            // if( !proc_id )
            //     outputFormatter.printf( 0, "The target set contains %lu vertices and %lu elements \n", tgtverts.size(),
            //                             tgtelems.size() );
            velist[2] = tgtverts.size();
            velist[3] = tgtelems.size();
        }
        //MB_CHK_ERR( mbCore->write_file( "source_mesh.h5m", nullptr, writeOptions, &runCtx->meshsets[0], 1 ) );
        //MB_CHK_ERR( mbCore->write_file( "target_mesh.h5m", nullptr, writeOptions, &runCtx->meshsets[1], 1 ) );

        // if( runCtx->nlayers && nprocs > 1 )
        // {
        //     remapper.ResetMeshSet( moab::Remapper::SourceMesh, runCtx->meshsets[3] );
        //     runCtx->meshes[0] = remapper.GetMesh( moab::Remapper::SourceMesh );  //  ?
        // }

        // First compute the covering set such that the target elements are fully covered by the
        // local source grid
        runCtx->timer_push( "construct covering set for intersection" );
        // if ghosting, no gnomonic
        if( runCtx->nlayers >= 1 ) runCtx->useGnomonicProjection = false;
        MB_CHK_ERR( remapper.ConstructCoveringSet( runCtx->epsrel, 1.0, 1.0, runCtx->boxeps, runCtx->rrmGrids,
                                                   runCtx->useGnomonicProjection, runCtx->nlayers ) );
        runCtx->timer_pop();

#ifdef MOAB_HAVE_MPI
        MPI_Reduce( velist, gvelist, 4, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD );
#else
        for( int i = 0; i < 4; i++ )
            gvelist[i] = velist[i];
#endif

        if( !proc_id && runCtx->print_diagnostics )
        {
            outputFormatter.printf( 0, "The source set contains %lu vertices and %lu elements \n", gvelist[0],
                                    gvelist[1] );
            outputFormatter.printf( 0, "The target set contains %lu vertices and %lu elements \n", gvelist[2],
                                    gvelist[3] );
        }

        // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
        runCtx->timer_push( "setup and compute mesh intersections" );
        MB_CHK_ERR( remapper.ComputeOverlapMesh( runCtx->kdtreeSearch, false ) );
        runCtx->timer_pop();

        // print some diagnostic checks to see if the overlap grid resolved the input meshes
        // correctly
        double dTotalOverlapArea = 0.0;
        if( runCtx->print_diagnostics )
        {
            moab::IntxAreaUtils areaAdaptorHuiller(
                moab::IntxAreaUtils::GaussQuadrature );  // lHuiller, GaussQuadrature
            double local_areas[3],
                global_areas[3];  // Array for Initial area, and through Method 1 and Method 2
            // local_areas[0] = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[1], radius_src );
            local_areas[0] = areaAdaptorHuiller.area_on_sphere( mbCore, runCtx->meshsets[0], radius_src );
            local_areas[1] = areaAdaptorHuiller.area_on_sphere( mbCore, runCtx->meshsets[1], radius_dest );
            local_areas[2] = areaAdaptorHuiller.area_on_sphere( mbCore, runCtx->meshsets[2], radius_src );

#ifdef MOAB_HAVE_MPI
            MPI_Allreduce( &local_areas[0], &global_areas[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
            global_areas[0] = local_areas[0];
            global_areas[1] = local_areas[1];
            global_areas[2] = local_areas[2];
#endif
            if( !proc_id )
            {
                outputFormatter.printf( 0,
                                        "initial area: source mesh = %12.14f, target mesh = "
                                        "%12.14f, overlap mesh = %12.14f\n",
                                        global_areas[0], global_areas[1], global_areas[2] );
                // outputFormatter.printf ( 0, " area with l'Huiller: %12.14f with Girard:
                // %12.14f\n", global_areas[2], global_areas[3] ); outputFormatter.printf ( 0, "
                // relative difference areas = %12.10e\n", fabs ( global_areas[2] - global_areas[3]
                // ) / global_areas[2] );
                outputFormatter.printf( 0, "relative error w.r.t source = %12.14e, and target = %12.14e\n",
                                        fabs( global_areas[0] - global_areas[2] ) / global_areas[0],
                                        fabs( global_areas[1] - global_areas[2] ) / global_areas[1] );
            }
            dTotalOverlapArea = global_areas[2];
        }

        if( runCtx->intxFilename.size() )
        {
            moab::EntityHandle writableOverlapSet;
            MB_CHK_SET_ERR( mbCore->create_meshset( moab::MESHSET_SET, writableOverlapSet ), "Can't create new set" );
            moab::EntityHandle meshOverlapSet = remapper.GetMeshSet( moab::Remapper::OverlapMesh );
            moab::Range ovEnts;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( meshOverlapSet, 2, ovEnts ), "Can't create new set" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( meshOverlapSet, 0, ovEnts ), "Can't create new set" );

#ifdef MOAB_HAVE_MPI
            // Do not remove ghosted entities if we still haven't computed weights
            // Remove ghosted entities from overlap set before writing the new mesh set to file
            if( nprocs > 1 )
            {
                moab::Range ghostedEnts;
                MB_CHK_ERR( remapper.GetOverlapAugmentedEntities( ghostedEnts ) );
                ovEnts = moab::subtract( ovEnts, ghostedEnts );
#ifdef MOAB_DBG
                if( !runCtx->skip_io )
                {
                    std::stringstream filename;
                    filename << "aug_overlap" << runCtx->pcomm->rank() << ".h5m";
                    MB_CHK_ERR( runCtx->mbcore->write_file( filename.str().c_str(), 0, 0, &meshOverlapSet, 1 ) );
                }
#endif
            }
#endif
            MB_CHK_SET_ERR( mbCore->add_entities( writableOverlapSet, ovEnts ), "adding local intx cells failed" );

#ifdef MOAB_HAVE_MPI
#ifdef MOAB_DBG
            if( nprocs > 1 && !runCtx->skip_io )
            {
                std::stringstream filename;
                filename << "writable_intx_" << runCtx->pcomm->rank() << ".h5m";
                MB_CHK_ERR( runCtx->mbcore->write_file( filename.str().c_str(), 0, 0, &writableOverlapSet, 1 ) );
            }
#endif
#endif

            size_t lastindex = runCtx->intxFilename.find_last_of( "." );
            sstr.str( "" );
            sstr << runCtx->intxFilename.substr( 0, lastindex ) << ".h5m";
            if( !runCtx->proc_id )
                std::cout << "Writing out the MOAB intersection mesh file to " << sstr.str() << std::endl;

            // Write out our computed intersection file
            if( !runCtx->skip_io )
            {
                MB_CHK_ERR( mbCore->write_file( sstr.str().c_str(), nullptr, writeOptions, &writableOverlapSet, 1 ) );
            }
        }

        if( runCtx->computeWeights )
        {
            runCtx->meshes[2] = remapper.GetMesh( moab::Remapper::OverlapMesh );
            if( !runCtx->proc_id ) std::cout << std::endl;

            runCtx->timer_push( "setup computation of weights" );
            // Call to generate the remapping weights with the tempest meshes
            moab::TempestOnlineMap* weightMap = new moab::TempestOnlineMap( &remapper );
            runCtx->timer_pop();

            runCtx->timer_push( "compute weights with TempestRemap" );
            MB_CHK_ERR( weightMap->GenerateRemappingWeights(
                runCtx->disc_methods[0],  // std::string strInputType
                runCtx->disc_methods[1],  // std::string strOutputType,
                runCtx->mapOptions,       // const GenerateOfflineMapAlgorithmOptions& options
                runCtx->doftag_names[0],  // const std::string& source_tag_name
                runCtx->doftag_names[1]   // const std::string& target_tag_name
                ) );
            runCtx->timer_pop();

            weightMap->PrintMapStatistics();

            // Invoke the CheckMap routine on the TempestRemap serial interface directly, if running
            // on a single process
            if( runCtx->fCheck )
            {
                const double dNormalTolerance = 1.0E-8;
                const double dStrictTolerance = 1.0E-12;
                weightMap->CheckMap( runCtx->fCheck, runCtx->fCheck, runCtx->fCheck && ( runCtx->ensureMonotonicity ),
                                     dNormalTolerance, dStrictTolerance, dTotalOverlapArea );
            }

            if( runCtx->outFilename.size() )
            {
                std::map< std::string, std::string > attrMap;
                attrMap["MOABversion"]   = std::string( MOAB_VERSION );
                attrMap["Title"]         = "MOAB-TempestRemap (mbtempest) Offline Regridding Weight Generator";
                attrMap["normalization"] = "ovarea";
                attrMap["remap_options"] = runCtx->mapOptions.strMethod;
                attrMap["domain_a"]      = runCtx->inFilenames[0];
                attrMap["domain_b"]      = runCtx->inFilenames[1];
                if( runCtx->intxFilename.size() ) attrMap["domain_aUb"] = runCtx->intxFilename;
                attrMap["map_aPb"]       = runCtx->outFilename;
                attrMap["methodorder_a"] = runCtx->disc_methods[0] + ":" + std::to_string( runCtx->disc_orders[0] ) +
                                           ":" + std::string( runCtx->doftag_names[0] );
                attrMap["concave_a"]     = runCtx->mapOptions.fSourceConcave ? "true" : "false";
                attrMap["methodorder_b"] = runCtx->disc_methods[1] + ":" + std::to_string( runCtx->disc_orders[1] ) +
                                           ":" + std::string( runCtx->doftag_names[1] );
                attrMap["concave_b"] = runCtx->mapOptions.fTargetConcave ? "true" : "false";
                attrMap["bubble"]    = runCtx->mapOptions.fNoBubble ? "false" : "true";
                attrMap["history"]   = historyStr;

                // Write the map file to disk in parallel using either HDF5 or SCRIP interface
                // in extra case; maybe need a better solution, just create it with the right meshset
                // from the beginning;
                if( !runCtx->skip_io )
                {
                    MB_CHK_ERR( weightMap->WriteParallelMap( runCtx->outFilename.c_str(), attrMap ) );
                }
            }

            if( runCtx->verifyWeights )
            {
                // Let us pick a sampling test function for solution evaluation
                // SH, SV, FH, USERVAR
                bool userVariable = false;
                moab::TempestOnlineMap::sample_function testFunction;
                if( !runCtx->variableToVerify.compare( "SH" ) )
                    testFunction = &sample_slow_harmonic;
                else if( !runCtx->variableToVerify.compare( "FH" ) )
                    testFunction = &sample_fast_harmonic;
                else if( !runCtx->variableToVerify.compare( "SV" ) )
                    testFunction = &sample_stationary_vortex;
                else
                {
                    userVariable = runCtx->variableToVerify.size() ? true : false;
                    testFunction = runCtx->variableToVerify.size() ? nullptr : sample_stationary_vortex;
                }

                moab::Tag srcAnalyticalFunction;
                moab::Tag tgtAnalyticalFunction;
                moab::Tag tgtProjectedFunction;
                if( testFunction )
                {
                    runCtx->timer_push( "describe a solution on source grid" );
                    MB_CHK_ERR( weightMap->DefineAnalyticalSolution( srcAnalyticalFunction, "AnalyticalSolnSrcExact",
                                                                     moab::Remapper::SourceMesh, testFunction ) );
                    runCtx->timer_pop();

                    runCtx->timer_push( "describe a solution on target grid" );

                    MB_CHK_ERR( weightMap->DefineAnalyticalSolution( tgtAnalyticalFunction, "AnalyticalSolnTgtExact",
                                                                     moab::Remapper::TargetMesh, testFunction,
                                                                     &tgtProjectedFunction, "ProjectedSolnTgt" ) );
                    // rval = mbCore->write_file ( "tgtWithSolnTag.h5m", nullptr, writeOptions,
                    // &runCtx->meshsets[1], 1 ); MB_CHK_ERR ( rval );
                    runCtx->timer_pop();
                }
                else
                {
                    MB_CHK_ERR( mbCore->tag_get_handle( runCtx->variableToVerify.c_str(), srcAnalyticalFunction ) );

                    MB_CHK_ERR( mbCore->tag_get_handle( "ProjectedSolnTgt", 1, moab::MB_TYPE_DOUBLE,
                                                        tgtProjectedFunction,
                                                        moab::MB_TAG_DENSE | moab::MB_TAG_CREAT ) );
                }

                if( !runCtx->skip_io )
                {
                    MB_CHK_ERR(
                        mbCore->write_file( "srcWithSolnTag.h5m", nullptr, writeOptions, &runCtx->meshsets[0], 1 ) );
                }

                runCtx->timer_push( "compute solution projection on target grid" );
                MB_CHK_ERR(
                    weightMap->ApplyWeights( srcAnalyticalFunction, tgtProjectedFunction, false, runCtx->cassType ) );
                runCtx->timer_pop();

                if( !runCtx->skip_io )
                {
                    MB_CHK_ERR(
                        mbCore->write_file( "tgtWithSolnTag2.h5m", nullptr, writeOptions, &runCtx->meshsets[1], 1 ) );
                }

                if( nprocs == 1 && runCtx->baselineFile.size() )
                {
                    // save the field from tgtWithSolnTag2 in a text file, and global ids for cells
                    moab::Range tgtCells;
                    MB_CHK_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, tgtCells ) );
                    std::vector< int > globIds;
                    globIds.resize( tgtCells.size() );
                    std::vector< double > vals;
                    vals.resize( tgtCells.size() );
                    moab::Tag projTag;
                    MB_CHK_ERR( mbCore->tag_get_handle( "ProjectedSolnTgt", projTag ) );
                    moab::Tag gid = mbCore->globalId_tag();
                    MB_CHK_ERR( mbCore->tag_get_data( gid, tgtCells, &globIds[0] ) );
                    MB_CHK_ERR( mbCore->tag_get_data( projTag, tgtCells, &vals[0] ) );
                    std::fstream fs;
                    fs.open( runCtx->baselineFile.c_str(), std::fstream::out );
                    fs << std::setprecision( 15 );  // maximum precision for doubles
                    for( size_t i = 0; i < tgtCells.size(); i++ )
                        fs << globIds[i] << " " << vals[i] << "\n";
                    fs.close();
                    // for good measure, save the source file too, with the tag AnalyticalSolnSrcExact
                    // it will be used later to test, along with a target file
                    if( !runCtx->skip_io )
                    {
                        MB_CHK_ERR( mbCore->write_file( "srcWithSolnTag.h5m", nullptr, writeOptions,
                                                        &runCtx->meshsets[0], 1 ) );
                    }
                }

                // compute error metrics if it is a known analytical functional
                if( !userVariable )
                {
                    runCtx->timer_push( "compute error metrics against analytical solution on target grid" );
                    std::map< std::string, double > errMetrics;
                    MB_CHK_ERR( weightMap->ComputeMetrics( moab::Remapper::TargetMesh, tgtAnalyticalFunction,
                                                           tgtProjectedFunction, errMetrics, true ) );
                    runCtx->timer_pop();
                }
            }

            delete weightMap;
        }
    }

    // Clean up
    remapper.clear();
    delete runCtx;
    delete mbCore;

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit( 0 );
}

static moab::ErrorCode CreateTempestMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    moab::ErrorCode rval = moab::MB_SUCCESS;
    int err;
    moab::DebugOutput& outputFormatter = ctx.outputFormatter;

    if( !ctx.proc_id )
    {
        outputFormatter.printf( 0, "Creating TempestRemap Mesh object ...\n" );
    }

    if( ctx.meshType == moab::TempestRemapper::OVERLAP_FILES )
    {
        ctx.timer_push( "create Tempest OverlapMesh" );
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapMesh( ctx.inFilenames[0], ctx.inFilenames[1], *tempest_mesh, ctx.outFilename, "NetCDF4",
                                   "exact", true );
        if( err )
            rval = moab::MB_FAILURE;
        else
            ctx.meshes.push_back( tempest_mesh );
        ctx.timer_pop();
    }
    else if( ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Load the meshes and validate
        ctx.meshsets.resize( 3 );
        ctx.meshes.resize( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet( moab::Remapper::OverlapMesh );

        ctx.timer_push( "load MOAB Source mesh" );
        // First the source
        MB_CHK_ERR(
            remapper.LoadMesh( moab::Remapper::SourceMesh, ctx.inFilenames[0], moab::TempestRemapper::DEFAULT ) );
        ctx.meshes[0] = remapper.GetMesh( moab::Remapper::SourceMesh );
        ctx.timer_pop();

        ctx.timer_push( "load MOAB Target mesh" );
        // Next the target
        MB_CHK_ERR(
            remapper.LoadMesh( moab::Remapper::TargetMesh, ctx.inFilenames[1], moab::TempestRemapper::DEFAULT ) );
        ctx.meshes[1] = remapper.GetMesh( moab::Remapper::TargetMesh );
        ctx.timer_pop();

        ctx.timer_push( "generate TempestRemap OverlapMesh" );
        // Now let us construct the overlap mesh, by calling TempestRemap interface directly
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapWithMeshes( *ctx.meshes[0], *ctx.meshes[1], *tempest_mesh, "" /*ctx.outFilename*/,
                                         "NetCDF4", "exact", false );
        ctx.timer_pop();
        if( err )
            rval = moab::MB_FAILURE;
        else
        {
            remapper.SetMesh( moab::Remapper::OverlapMesh, tempest_mesh );
            ctx.meshes[2] = remapper.GetMesh( moab::Remapper::OverlapMesh );
            // ctx.meshes.push_back(*tempest_mesh);
        }
    }
    else if( ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        ctx.meshsets.resize( 3 );
        ctx.meshes.resize( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet( moab::Remapper::OverlapMesh );

        const double radius_src  = 1.0 /*2.0*acos(-1.0)*/;
        const double radius_dest = 1.0 /*2.0*acos(-1.0)*/;

        std::vector< int > smetadata, tmetadata;
        //const char* additional_read_opts = ( ctx.n_procs > 1 ? "NO_SET_CONTAINING_PARENTS;" : "" );
        std::string additional_read_opts_src = get_file_read_options( ctx, ctx.inFilenames[0] );

        ctx.timer_push( "load MOAB Source mesh" );
        // Load the source mesh and validate
        MB_CHK_ERR( remapper.LoadNativeMesh( ctx.inFilenames[0], ctx.meshsets[0], smetadata,
                                             additional_read_opts_src.c_str() ) );
        if( smetadata.size() ) remapper.SetMeshType( moab::Remapper::SourceMesh, smetadata );
        ctx.timer_pop();

        ctx.timer_push( "preprocess MOAB Source mesh" );
        // Rescale the radius of both to compute the intersection
        MB_CHK_ERR( moab::IntxUtils::ScaleToRadius( ctx.mbcore, ctx.meshsets[0], radius_src ) );

        ctx.timer_pop();

        // Load the target mesh and validate
        std::string addititional_read_opts_tgt = get_file_read_options( ctx, ctx.inFilenames[1] );
        // addititional_read_opts_tgt += "PARALLEL_GHOSTS=2.0.3;PARALLEL_THIN_GHOST_LAYER;";

        ctx.timer_push( "load MOAB Target mesh" );
        MB_CHK_ERR( remapper.LoadNativeMesh( ctx.inFilenames[1], ctx.meshsets[1], tmetadata,
                                             addititional_read_opts_tgt.c_str() ) );
        if( tmetadata.size() ) remapper.SetMeshType( moab::Remapper::TargetMesh, tmetadata );
        ctx.timer_pop();

        ctx.timer_push( "preprocess MOAB Target mesh" );
        MB_CHK_ERR( moab::IntxUtils::ScaleToRadius( ctx.mbcore, ctx.meshsets[1], radius_dest ) );
        ctx.timer_pop();

        if( ctx.computeWeights )
        {
            ctx.timer_push( "convert MOAB meshes to TempestRemap meshes in memory" );
            // convert MOAB representation to TempestRemap's Mesh for source
            MB_CHK_ERR( remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh ) );
            ctx.meshes[0] = remapper.GetMesh( moab::Remapper::SourceMesh );

            // convert MOAB representation to TempestRemap's Mesh for target
            MB_CHK_ERR( remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh ) );
            ctx.meshes[1] = remapper.GetMesh( moab::Remapper::TargetMesh );
            ctx.timer_pop();
        }
    }
    else if( ctx.meshType == moab::TempestRemapper::ICO )
    {
        ctx.timer_push( "generate ICO mesh with TempestRemap" );
        err = GenerateICOMesh( *tempest_mesh, ctx.blockSize, ctx.computeDual, ctx.outFilename, "NetCDF4" );
        ctx.timer_pop();

        if( err )
            rval = moab::MB_FAILURE;
        else
            ctx.meshes.push_back( tempest_mesh );
    }
    else if( ctx.meshType == moab::TempestRemapper::RLL )
    {
        ctx.timer_push( "generate RLL mesh with TempestRemap" );
        err = GenerateRLLMesh( *tempest_mesh,                     // Mesh& meshOut,
                               ctx.blockSize * 2, ctx.blockSize,  // int nLongitudes, int nLatitudes,
                               0.0, 360.0,                        // double dLonBegin, double dLonEnd,
                               -90.0, 90.0,                       // double dLatBegin, double dLatEnd,
                               false, false, false,  // bool fGlobalCap, bool fFlipLatLon, bool fForceGlobal,
                               "" /*ctx.inFilename*/, "",
                               "",  // std::string strInputFile, std::string strInputFileLonName, std::string
                                    // strInputFileLatName,
                               ctx.outFilename, "NetCDF4",  // std::string strOutputFile, std::string strOutputFormat
                               true                         // bool fVerbose
        );
        ctx.timer_pop();

        if( err )
            rval = moab::MB_FAILURE;
        else
            ctx.meshes.push_back( tempest_mesh );
    }
    else  // default
    {
        ctx.timer_push( "generate CS mesh with TempestRemap" );
        err = GenerateCSMesh( *tempest_mesh, ctx.blockSize, ctx.outFilename, "NetCDF4" );
        ctx.timer_pop();
        if( err )
            rval = moab::MB_FAILURE;
        else
            ctx.meshes.push_back( tempest_mesh );
    }

    if( ctx.meshType != moab::TempestRemapper::OVERLAP_MOAB && !tempest_mesh )
    {
        std::cout << "Tempest Mesh is not a complete object; Quitting...";
        exit( -1 );
    }

    return rval;
}

#undef MOAB_DBG
///////////////////////////////////////////////
//         Test functions

double sample_slow_harmonic( double dLon, double dLat )
{
    return ( 2.0 + cos( dLat ) * cos( dLat ) * cos( 2.0 * dLon ) );
}

double sample_fast_harmonic( double dLon, double dLat )
{
    return ( 2.0 + pow( sin( 2.0 * dLat ), 16.0 ) * cos( 16.0 * dLon ) );
    // return (2.0 + pow(cos(2.0 * dLat), 16.0) * cos(16.0 * dLon));
}

double sample_constant( double /*dLon*/, double /*dLat*/ )
{
    return 1.0;
}

double sample_stationary_vortex( double dLon, double dLat )
{
    const double dLon0 = 0.0;
    const double dLat0 = 0.6;
    const double dR0   = 3.0;
    const double dD    = 5.0;
    const double dT    = 6.0;

    ///		Find the rotated longitude and latitude of a point on a sphere
    ///		with pole at (dLonC, dLatC).
    {
        double dSinC = sin( dLat0 );
        double dCosC = cos( dLat0 );
        double dCosT = cos( dLat );
        double dSinT = sin( dLat );

        double dTrm = dCosT * cos( dLon - dLon0 );
        double dX   = dSinC * dTrm - dCosC * dSinT;
        double dY   = dCosT * sin( dLon - dLon0 );
        double dZ   = dSinC * dSinT + dCosC * dTrm;

        dLon = atan2( dY, dX );
        if( dLon < 0.0 )
        {
            dLon += 2.0 * M_PI;
        }
        dLat = asin( dZ );
    }

    double dRho = dR0 * cos( dLat );
    double dVt  = 3.0 * sqrt( 3.0 ) / 2.0 / cosh( dRho ) / cosh( dRho ) * tanh( dRho );

    double dOmega;
    if( dRho == 0.0 )
    {
        dOmega = 0.0;
    }
    else
    {
        dOmega = dVt / dRho;
    }

    return ( 1.0 - tanh( dRho / dD * sin( dLon - dOmega * dT ) ) );
}

///////////////////////////////////////////////
