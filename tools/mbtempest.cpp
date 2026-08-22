/**
 * @file mbtempest.cpp
 * @brief MOAB-Tempest: A powerful mesh generation and remapping tool for climate and weather applications
 *
 * @section overview Overview
 * MOAB-Tempest is a command-line tool that provides mesh generation and conservative remapping capabilities
 * for climate and weather modeling. It combines the power of MOAB (Mesh-Oriented datABase) with the
 * TempestRemap library to enable high-performance, parallel mesh generation and remapping operations.
 *
 * @section features Key Features
 * - Generation of various spherical mesh types (Cubed-Sphere, RLL, Icosahedral)
 * - Support for high-order discretization methods (FV, CGLL, DGLL)
 * - Conservative remapping between different mesh types
 * - Parallel processing support via MPI
 * - Flexible I/O with support for multiple file formats
 * - Built-in analytical functions for testing and validation
 *
 * @section algorithms Supported Algorithms
 * - Mesh Generation:
 *   - Cubed-Sphere (CS) meshes
 *   - Regular Latitude-Longitude (RLL) meshes
 *   - Icosahedral (ICO) meshes
 *   - Overlap meshes for remapping
 * - Remapping Methods:
 *   - Finite Volume (FV)
 *   - Continuous Galerkin (CGLL)
 *   - Discontinuous Galerkin (DGLL)
 *   - Monotonic and high-order variants
 *
 * @section usage Basic Usage Examples
 * @code
 * # Generate a Cubed-Sphere mesh with resolution 25
 * ./mbtempest --type 0 --res 25 --file cubed_sphere_mesh.h5m
 *
 * # Generate a RLL mesh with resolution 90x180 (lon x lat)
 * ./mbtempest --type 1 --res 90 --file rll_mesh.h5m
 *
 * # Generate an Icosahedral mesh with resolution 25 (dual mesh)
 * ./mbtempest --type 2 --res 25 --dual --file icosahedral_dual_mesh.h5m
 *
 * # Compute overlap between two meshes
 * ./mbtempest --type 5 --load mesh1.h5m --load mesh2.h5m intx intersection_mesh.h5m
 *
 * # Generate a remapping weights file between two meshes: FV to FV (default)
 * ./mbtempest --type 5 --load source_mesh.h5m --load target_mesh.h5m --file weights.nc
 *
 * # Generate remapping weights file between two meshes: making it explicit (SE to FV)
 * ./mbtempest --type 5 --load source_mesh.h5m --load target_mesh.h5m \
 *             --order 4 --method cgll --global_id GLOBAL_DOFS \
 *             --order 1 --method fv --global_id GLOBAL_ID \
 *             --file weights_se_to_fv.nc
 * @endcode
 *
 * @section options Command Line Options
 * Run './mbtempest --help' for a complete list of available options.
 *
 * @section notes Notes
 * - For parallel execution, use MPI launcher (e.g., mpirun, mpiexec)
 * - Output formats: .h5m (MOAB), .nc (NetCDF), .exo (ExodusII)
 * - Requires MOAB and TempestRemap libraries
 *
 * @author MOAB Development Team
 * @date Created: 2023
 */

// standard C++ includes
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <memory>
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

#ifdef MOAB_HAVE_MPI
// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

/**
 * @brief Context class for MOAB-TempestRemap tool configuration and state management
 */
class ToolContext
{
  public:
    // Core components
    moab::Core* const mbcore;  ///< MOAB Core instance for mesh operations
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* const pcomm;  ///< Parallel communicator (nullptr in serial)
#endif
    const int proc_id;                  ///< MPI process rank (0 for serial)
    const int n_procs;                  ///< Total number of MPI processes (1 for serial)
    moab::DebugOutput outputFormatter;  ///< Formatter for debug output

    // Mesh and remapping configuration
    moab::TempestRemapper::TempestMeshType meshType{ moab::TempestRemapper::DEFAULT };  ///< Type of mesh to generate
    std::vector< std::string > inFilenames;      ///< Input filenames for source and target meshes
    std::vector< int > disc_orders;              ///< Discretization orders for source and target
    std::vector< std::string > disc_methods;     ///< Discretization methods (fv, cgll, dgll) for source and target
    std::vector< std::string > doftag_names;     ///< Degree of freedom tag names for source and target
    std::string outFilename{ "outputFile.nc" };  ///< Output filename for remapping results
    std::string intxFilename;                    ///< Intersection mesh filename (optional)
    std::string baselineFile;                    ///< Baseline file for verification (optional)
    std::string variableToVerify;                ///< Variable name for verification (optional)
    std::string fvMethod{ "none" };              ///< Finite volume method specification

    // Remapping options
    GenerateOfflineMapAlgorithmOptions mapOptions;  ///< Configuration for offline map generation
    moab::TempestOnlineMap::CAASType cassType{
        moab::TempestOnlineMap::CAAS_NONE };  ///< Conservative and accurate advection scheme type
    int ensureMonotonicity{ 0 };              ///< Monotonicity enforcement level (0=none, 1=basic, 2=full, 3=strict)
    bool rrmGrids{ false };                   ///< Flag to use RRM (Regional Refinement Meshes)
    bool kdtreeSearch{ true };                ///< Enable KD-tree for spatial searches
    bool fCheck{ false };                     ///< Enable additional checking during remapping
    bool fVolumetric{ false };                ///< Enable volumetric (3D) remapping
    bool useGnomonicProjection{ false };      ///< Use gnomonic projection for certain operations
    bool print_diagnostics{ false };          ///< Print detailed diagnostic information
    bool skip_intersection{ false };          ///< Skip intersection computation (for debugging)
    double boxeps{ 1e-7 };                    ///< Epsilon for bounding box checks
    double epsrel{ ReferenceTolerance };      ///< Relative tolerance for convergence

    // Mesh operations control
    bool skip_io{ false };             ///< Skip file I/O operations (for testing)
    bool computeDual{ false };         ///< Compute dual mesh
    bool computeWeights{ false };      ///< Compute interpolation weights
    bool verifyConservation{ false };  ///< Verify conservation properties
    bool verifyWeights{ false };       ///< Verify interpolation weights
    bool enforceConvexity{ false };    ///< Enforce convexity in mesh elements

    // Performance and debugging
    std::unique_ptr< moab::CpuTimer > timer;  ///< Timer for performance measurement
    double timer_ops{ 0.0 };                  ///< Operation timer value
    std::string opName;                       ///< Name of current operation being timed
    int nlayers{ 0 };                         ///< Number of ghost layers for parallel operations
    int blockSize{ 5 };                       ///< Block size for vectorized operations

    // Mesh data
    std::vector< Mesh* > meshes;                 ///< Collection of TempestRemap meshes
    std::vector< moab::EntityHandle > meshsets;  ///< MOAB entity sets for meshes

    /**
     * @brief Construct a new ToolContext object with MPI support
     * @param icore MOAB Core instance (must not be null)
     * @param p_pcomm Parallel communicator (must not be null in MPI mode)
     * @throw std::invalid_argument if icore is null or p_pcomm is null in MPI mode
     */
#ifdef MOAB_HAVE_MPI
    ToolContext( moab::Core* icore, moab::ParallelComm* p_pcomm )
        : mbcore( icore ), pcomm( p_pcomm ), proc_id( p_pcomm ? p_pcomm->rank() : 0 ),
          n_procs( p_pcomm ? p_pcomm->size() : 1 ), outputFormatter( std::cout, p_pcomm ? p_pcomm->rank() : 0, 0 )
    {
        if( !icore ) throw std::invalid_argument( "MOAB Core instance cannot be null" );
        if( !p_pcomm ) throw std::invalid_argument( "ParallelComm cannot be null in MPI mode" );
#else
    /**
     * @brief Construct a new ToolContext object (serial version)
     * @param icore MOAB Core instance (must not be null)
     * @throw std::invalid_argument if icore is null
     */
    explicit ToolContext( moab::Core* icore )
        : mbcore( icore ), proc_id( 0 ), n_procs( 1 ), outputFormatter( std::cout, 0, 0 )
    {
#endif
        // Initialize default values
        inFilenames.reserve( 2 );
        doftag_names = { "GLOBAL_ID", "GLOBAL_ID" };
        disc_orders  = { 1, 1 };
        disc_methods = { "fv", "fv" };

        // Initialize timer and output formatter
        timer = std::make_unique< moab::CpuTimer >();
        outputFormatter.set_prefix( "[MBTempest]: " );

        // Set default map options
        mapOptions.fNoConservation = false;
        mapOptions.fMonotone       = false;
        mapOptions.fNoCorrectAreas = false;
        mapOptions.fNoCheck        = false;
        mapOptions.nPin            = 1;
        mapOptions.nPout           = 1;
    }

    // Rule of Five - Delete copy/move operations as mbcore is const
    ~ToolContext()                               = default;
    ToolContext( const ToolContext& )            = delete;
    ToolContext& operator=( const ToolContext& ) = delete;
    ToolContext( ToolContext&& )                 = delete;
    ToolContext& operator=( ToolContext&& )      = delete;

    /**
     * @brief Start timing an operation
     * @param operation Name of the operation being timed
     */
    void timer_push( const std::string& operation )
    {
        timer_ops = timer->time_since_birth();
        opName    = operation;
    }

    /**
     * @brief Stop timing and log the operation duration
     */
    void timer_pop()
    {
        double locElapsed = timer->time_since_birth() - timer_ops;
        double avgElapsed = locElapsed;
        double maxElapsed = locElapsed;

#ifdef MOAB_HAVE_MPI
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, pcomm->comm() );
        MPI_Reduce( &locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, pcomm->comm() );
        avgElapsed /= n_procs;
#endif

        if( proc_id == 0 )
        {
            std::cout << "[LOG] Time taken to " << opName << ": max = " << maxElapsed << ", avg = " << avgElapsed
                      << "\n";
        }
        opName.clear();
    }

    /**
     * @brief Parse command line arguments
     * @param argc Argument count
     * @param argv Argument values
     * @return moab::ErrorCode indicating success or failure
     * @throw std::invalid_argument for invalid command line arguments
     */
    moab::ErrorCode ParseCLOptions( int argc, char** argv )
    {
        // Initialize variables for command line options
        int imeshType                  = 0;
        std::string expectedFName      = "output.exo";
        std::string expectedMethod     = "fv";
        std::string expectedFVMethod   = "none";
        std::string expectedDofTagName = "GLOBAL_ID";
        int expectedOrder              = 1;
        int useCAAS                    = 0;
        int nlayer_input               = -1;  // -1 means not set by user
        bool version_info              = false;

        // Print command line for debugging
        if( proc_id == 0 )
        {
            std::cout << "Command line options provided to mbtempest:\n  ";
            for( int i = 0; i < argc; ++i )
            {
                std::cout << argv[i] << " ";
            }
            std::cout << "\n" << std::endl;
        }

        // Create options object with description
        ProgOptions opts( "mbtempest - A mesh generation and remapping tool" );

        // Mesh generation options
        opts.addOpt< int >( "type,t",
                            "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP_FILES=3, "
                            "OVERLAP_MEMORY=4, OVERLAP_MOAB=5])",
                            &imeshType );

        opts.addOpt< int >( "res,r", "Resolution of the mesh (default=5)", &blockSize );

        opts.addOpt< void >( "dual,d", "Output the dual of the mesh (relevant only for ICO mesh type)", &computeDual );

        opts.addOpt< std::string >( "file,f", "Output computed mesh or remapping weights to specified filename",
                                    &outFilename );

        // Input/Output options
        opts.addOpt< std::string >(
            "load,l", "Input mesh filenames for source and target meshes. (relevant only when computing weights)",
            &expectedFName );

        opts.addOpt< void >( "advfront,a",
                             "Use the advancing front intersection instead of the Kd-tree based algorithm to compute "
                             "mesh intersections.",
                             &kdtreeSearch );

        opts.addOpt< std::string >( "intx,i", "Output TempestRemap intersection mesh filename", &intxFilename );

        opts.addOpt< void >(
            "weights,w",
            "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)",
            &computeWeights );

        // Discretization options
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

        // Advanced options
        opts.addOpt< std::string >( "fvmethod",
                                    "Sub-type method for FV-FV projections (invdist, delaunay, bilin, intbilin, "
                                    "intbilingb, none. Default: none)",
                                    &expectedFVMethod );

        opts.addOpt< void >(
            "noconserve", "Do not apply conservation to the resultant weights (relevant only when computing weights)",
            &mapOptions.fNoConservation );

        opts.addOpt< void >(
            "volumetric", "Apply a volumetric projection to compute the weights (relevant only when computing weights)",
            &fVolumetric );

        opts.addOpt< void >( "skip_intersection", "Skip mesh intersection computation.", &skip_intersection );

        opts.addOpt< void >( "skip_output", "For performance studies, skip all I/O operations.", &skip_io );

        opts.addOpt< void >( "gnomonic", "Use Gnomonic plane projections to compute coverage mesh.",
                             &useGnomonicProjection );

        opts.addOpt< void >( "enforce_convexity", "Check convexity of input meshes to compute mesh intersections",
                             &enforceConvexity );

        opts.addOpt< void >( "nobubble", "Do not use bubble on interior of spectral element nodes",
                             &mapOptions.fNoBubble );

        opts.addOpt< void >(
            "sparseconstraints",
            "Use sparse solver for constraints when we have high-valence (typical with high-res RLL mesh)",
            &mapOptions.fSparseConstraints );

        opts.addOpt< void >(
            "rrmgrids",
            "At least one of the meshes is a regionally refined grid (relevant to accelerate intersection computation)",
            &rrmGrids );

        opts.addOpt< void >( "checkmap", "Check the generated map for conservation and consistency", &fCheck );

        opts.addOpt< void >( "verify",
                             "Verify the accuracy of the maps by projecting analytical functions from source to target "
                             "grid by applying the maps",
                             &verifyWeights );

        opts.addOpt< std::string >( "var",
                                    "Tag name of the variable to use in the verification study (error metrics for user "
                                    "defined variables may not be available)",
                                    &variableToVerify );

        opts.addOpt< int >( "monotonicity", "Ensure monotonicity in the weight generation. Options=[0,1,2,3]",
                            &ensureMonotonicity );

        opts.addOpt< int >( "ghost",
                            "Number of ghost layers in coverage mesh (overrides automatic selection: 0 for FV order 1, "
                            "p+1 for FV order p>1)",
                            &nlayer_input );

        opts.addOpt< double >( "boxeps", "The tolerance for boxes (default=1e-7)", &boxeps );

        opts.addOpt< int >( "limiter", "Apply nonlinear filter after linear map application", &useCAAS );

        opts.addOpt< std::string >( "baseline", "Output baseline file", &baselineFile );

        opts.addOpt< void >( "manual", "Show documentation about usage with examples" );

        opts.addOpt< void >( "version", "Show version information", &version_info );

        // Parse command line
        opts.parseCommandLine( argc, argv );

        // Handle call for detailed information
        if( opts.numOptSet( "manual" ) > 0 )
        {
            if( this->proc_id == 0 )
            {
                this->printHelp( argv[0] );
            }
            exit( 0 );
        }

        if( version_info )
        {
            if( this->proc_id == 0 )
            {
                std::cout << "mbtempest is part of the MOAB library version " << std::string( MOAB_PACKAGE_VERSION )
                          << "\n";
            }
            exit( 0 );
        }

        // Process mesh type
        switch( imeshType )
        {
            case 0:
                this->meshType = moab::TempestRemapper::CS;
                break;
            case 1:
                this->meshType = moab::TempestRemapper::RLL;
                break;
            case 2:
                this->meshType = moab::TempestRemapper::ICO;
                break;
            case 3:
                this->meshType = moab::TempestRemapper::OVERLAP_FILES;
                break;
            case 4:
                this->meshType = moab::TempestRemapper::OVERLAP_MEMORY;
                break;
            case 5:
                this->meshType = moab::TempestRemapper::OVERLAP_MOAB;
                break;
            default:
                this->meshType = moab::TempestRemapper::DEFAULT;
                break;
        }

        // Process CAAS type
        switch( useCAAS )
        {
            case 1:
                this->cassType = moab::TempestOnlineMap::CAAS_GLOBAL;
                break;
            case 2:
                this->cassType = moab::TempestOnlineMap::CAAS_LOCAL;
                break;
            case 3:
                this->cassType = moab::TempestOnlineMap::CAAS_LOCAL_ADJACENT;
                break;
            case 4:
                this->cassType = moab::TempestOnlineMap::CAAS_QLT;
                break;
            default:
                this->cassType = moab::TempestOnlineMap::CAAS_NONE;
                break;
        }

        // Process input files if provided
        if( !expectedFName.empty() )
        {
            this->inFilenames = { expectedFName };
        }

        // Process discretization options through processMeshOptions to handle both single and multiple values
        // Set initial defaults that can be overridden by processMeshOptions
        this->fvMethod     = expectedFVMethod;
        this->disc_orders  = { expectedOrder, expectedOrder };
        this->disc_methods = { expectedMethod, expectedMethod };
        this->doftag_names = { expectedDofTagName, expectedDofTagName };

        // Let processMeshOptions handle all the discretization option processing
        this->processMeshOptions( opts );

        // Now use the processed values for map configuration
        this->mapOptions.nPin           = this->disc_orders[0];
        this->mapOptions.nPout          = this->disc_orders[1];
        this->mapOptions.fSourceConcave = false;
        this->mapOptions.fTargetConcave = false;
        this->mapOptions.strMethod      = "";

        // Configure map options with the processed values - this handles all remaining setup
        this->configureMapOptions( nlayer_input );

        // Print runtime parameters
        this->printRuntimeParameters();

        return moab::MB_SUCCESS;
    }

    /**
     * @brief Get the appropriate MOAB read options based on file extension and parallel configuration
     *
     * @param ctx Tool context containing parallel information
     * @param filename Input filename to determine read options
     * @return std::string MOAB read options string
     */
    std::string get_file_read_options( const std::string& filename )
    {
        // For serial execution, return default options
        if( n_procs <= 1 )
        {
            return "";
        }

        // Extract file extension
        const size_t last_dot = filename.find_last_of( "." );
        if( last_dot == std::string::npos )
        {
            return "";  // No extension found
        }

        const std::string extension = filename.substr( last_dot + 1 );

        // Handle H5M files
        if( extension == "h5m" )
        {
            return "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";
        }

        // Handle NetCDF files
        if( extension == "nc" )
        {
            // Default NetCDF options
#ifdef MOAB_HAVE_ZOLTAN
            std::string netcdf_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;";
#else
            std::string netcdf_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;";
#endif
            // Only rank 0 needs to determine the NetCDF file type
            if( proc_id == 0 )
            {
                NcFile ncFile( filename.c_str(), NcFile::ReadOnly );
                if( !ncFile.is_valid() )
                {
                    // Handle invalid file
                    return netcdf_options;
                }

                // Check for different NetCDF formats
                int format_flags = 0;
                for( int i = 0; i < ncFile.num_dims(); i++ )
                {
                    const std::string dim_name = ncFile.get_dim( i )->name();

                    if( dim_name == "grid_size" || dim_name == "grid_corners" || dim_name == "grid_rank" )
                    {
                        format_flags |= 1;  // SCRIP format
                    }
                    else if( dim_name == "nodeCount" || dim_name == "elementCount" || dim_name == "maxNodePElement" )
                    {
                        format_flags |= 2;  // ESMF format
                    }
                    else if( dim_name == "nCells" || dim_name == "nEdges" || dim_name == "nVertices" ||
                             dim_name == "vertexDegree" )
                    {
                        format_flags |= 4;  // MPAS format
                    }
                }

                // Apply format-specific options
                if( format_flags & 2 )
                {  // ESMF format
                    netcdf_options += "PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=;";
                }
                else if( format_flags & 1 )
                {                          // SCRIP format
                    netcdf_options += "";  // no extra options necessary for now
                }
                else if( format_flags & 4 )
                {  // MPAS format
                    netcdf_options += "PARALLEL_RESOLVE_SHARED_ENTS;NO_EDGES;NO_MIXED_ELEMENTS;VARIABLE=;";
                }
            }

            // Broadcast the options to all processes
#ifdef MOAB_HAVE_MPI
            int line_size = netcdf_options.size();
            MPI_Bcast( &line_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
            if( proc_id != 0 )
            {
                netcdf_options.resize( line_size );
            }
            MPI_Bcast( const_cast< char* >( netcdf_options.data() ), line_size, MPI_CHAR, 0, MPI_COMM_WORLD );
#endif

            return netcdf_options;
        }

        // Default options for other file types
        return "PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;";
    }

  private:
    /**
     * @brief Print detailed help message with usage examples
     * @param progName Program name
     */
    void printHelp( const char* progName ) const
    {
        if( this->proc_id != 0 ) return;

        std::cout << "MOAB-Tempest: A mesh generation and remapping tool\n"
                  << "==================================================\n\n"
                  << "Usage: " << progName << " [OPTIONS]\n\n"
                  << "Mesh Generation Options:\n"
                  << "  -t, --type TYPE       Type of mesh to generate (required for mesh generation):\n"
                  << "                           0 = Cubed-Sphere (CS)\n"
                  << "                           1 = Regular Latitude-Longitude (RLL)\n"
                  << "                           2 = Icosahedral (ICO)\n"
                  << "                           3 = TempestRemap overlap (thin interface))\n"
                  << "                           4 = MOAB with TempestRemap overlap in memory\n"
                  << "                           5 = Parallel handling of Overlap meshes with MOAB (recommended)\n\n"
                  << "  -r, --res N           Resolution (number of elements on edge, default: 10)\n"
                  << "  -f, --file FILE       Output filename (default: output.h5m)\n\n"
                  << "Discretization Options:\n"
                  << "  -m, --method METHOD   Discretization method (default: fv):\n"
                  << "                           fv   = Finite Volume\n"
                  << "                           cgll = Continuous Galerkin with Legendre-Gauss-Lobatto\n"
                  << "                           dgll = Discontinuous Galerkin with Legendre-Gauss-Lobatto\n\n"
                  << "  -o, --order N         Discretization order (default: 1, range: 1-4)\n\n"
                  << "Remapping Options:\n"
                  << "  --mono N              Monotonicity constraints (default: 0):\n"
                  << "                           0 = No monotonicity\n"
                  << "                           1 = Basic monotonicity\n"
                  << "                           2 = Full monotonicity with bounds\n"
                  << "                           3 = Strict monotonicity\n\n"
                  << "  --limiter TYPE        Nonlinear limiting (optional):\n"
                  << "                           none = No limiting (default)\n"
                  << "                           global = Global CAAS limiting\n"
                  << "                           local = Localized CAAS limiting\n"
                  << "                           qlt = Quasi-Local Tree-based limiting\n\n"
                  << "Input/Output Options:\n"
                  << "  -l, --load FILE       Load input mesh file (use twice for source and target)\n"
                  << "  -i, --global_id TAG   Global ID tag name (default: GLOBAL_ID)\n"
                  << "  --diagnostics         Print diagnostic information\n\n"
                  << "Miscellaneous Options:\n"
                  << "  --manual              Show this help message and exit\n"
                  << "  --version             Show version information\n\n"
                  << "Examples:\n"
                  << "  # Generate a cubed-sphere mesh with resolution 25\n"
                  << "  " << progName << " --type 0 --res 25 -f cs_mesh.h5m\n\n"
                  << "  # Generate a latitude-longitude mesh with resolution 180\n"
                  << "  " << progName << " --type 1 --res 180 -f rll_mesh.h5m\n\n"
                  << "  # Create a map between two meshes with order 4\n"
                  << "  " << progName << " --type 5 --load source_mesh.h5m --load target_mesh.h5m \\\n"
                  << "      --method cgll --order 4 --global_id GLOBAL_DOFS \\\n"
                  << "      --method fv --order 1 --limiter 1 --file map.nc\n";
    }

    /**
     * @brief Get mesh type as string
     * @return String representation of mesh type
     */
    std::string getMeshTypeName() const
    {
        switch( this->meshType )
        {
            case moab::TempestRemapper::CS:
                return "Cubed-Sphere";
            case moab::TempestRemapper::RLL:
                return "Latitude-Longitude";
            case moab::TempestRemapper::ICO:
                return "Icosahedral";
            case moab::TempestRemapper::OVERLAP_FILES:
                return "Overlap (files)";
            case moab::TempestRemapper::OVERLAP_MEMORY:
                return "Overlap (memory)";
            case moab::TempestRemapper::OVERLAP_MOAB:
                return "Overlap (MOAB)";
            default:
                return "Unknown";
        }
    }

    /**
     * @brief Process mesh options from command line
     * @param opts Program options
     * @param expectedFVMethod Expected finite volume method
     * @param nlayer_input Number of ghost layers
     */
    void processMeshOptions( ProgOptions& opts )
    {
        if( this->meshType <= moab::TempestRemapper::ICO ) return;

        // Process input files
        std::vector< std::string > inputFiles;
        opts.getOptAllArgs( "load,l", inputFiles );
        if( !inputFiles.empty() )
        {
            this->inFilenames = inputFiles;
            if( this->inFilenames.size() != 2 )
            {
                throw std::runtime_error( "Exactly two input filenames must be provided with -l/--load" );
            }
        }

        // Process discretization orders
        std::vector< int > orders;
        opts.getOptAllArgs( "order,o", orders );
        if( !orders.empty() )
        {
            this->disc_orders = orders;
            if( this->disc_orders.size() == 1 )
            {
                this->disc_orders.push_back( this->disc_orders[0] );
            }
            else if( this->disc_orders.size() != 2 )
            {
                throw std::runtime_error( "Must specify 1 or 2 values for order (source [target])" );
            }

            for( const auto& order : this->disc_orders )
            {
                if( order < 1 || order > 4 )
                {
                    throw std::runtime_error( "Discretization order must be between 1 and 4" );
                }
            }
        }

        // Process discretization methods
        std::vector< std::string > methods;
        opts.getOptAllArgs( "method,m", methods );
        if( !methods.empty() )
        {
            this->disc_methods = methods;
            if( this->disc_methods.size() == 1 )
            {
                // Use same method for both source and target
                this->disc_methods.push_back( this->disc_methods[0] );
            }
            else if( this->disc_methods.size() != 2 )
            {
                throw std::runtime_error( "Must specify 1 or 2 values for method (source [target])" );
            }

            // Validate method values
            for( const auto& method : this->disc_methods )
            {
                if( method != "fv" && method != "cgll" && method != "dgll" && method != "pcloud" )
                {
                    throw std::runtime_error( "Invalid method '" + method + "'. Must be one of: fv, cgll, dgll" );
                }
            }
        }

        // Process DOF tag names
        std::vector< std::string > tags;
        opts.getOptAllArgs( "global_id,i", tags );
        if( !tags.empty() )
        {
            this->doftag_names = tags;
            if( this->doftag_names.size() == 1 )
            {
                // Use same tag name for both source and target
                this->doftag_names.push_back( this->doftag_names[0] );
            }
            else if( this->doftag_names.size() != 2 )
            {
                throw std::runtime_error( "Must specify 1 or 2 values for DOF tag names (source [target])" );
            }
        }

        // Process output filename if specified
        std::string outFile;
        if( opts.getOpt( "file,f", &outFile ) )
        {
            this->outFilename = outFile;
        }
        // Note: configureMapOptions is now called from ParseCLOptions after processMeshOptions completes
    }

    /**
     * @brief Print all runtime parameters in a formatted way
     */
    void printRuntimeParameters() const
    {
        if( this->proc_id != 0 ) return;

        constexpr int width = 60;

        std::cout << std::string( width, '=' ) << "\n";
        std::cout << "  MOAB-TempestRemap Runtime Configuration " << "\n";
        std::cout << std::string( width, '=' );

        // Input files
        if( this->meshType == moab::TempestRemapper::OVERLAP_MOAB )
        {
            if( !this->inFilenames.empty() )
            {
                std::cout << "\n\nInput Files:";
                std::cout << "\n  Source mesh:          " << this->inFilenames[0];
                std::cout << "\n  Target mesh:          " << this->inFilenames[1];
            }

            std::cout << "\n\nOutput Files:";
            if( !skip_intersection )
                std::cout << "\n  Intersection mesh:    "
                          << ( this->computeWeights ? this->intxFilename : this->outFilename );
            if( computeWeights ) std::cout << "\n  Remap weights:        " << this->outFilename;
        }

        // Mesh configuration
        std::cout << "\n\nMesh Configuration:";
        std::cout << "\n  Mesh type:              " << this->getMeshTypeName();
        if( this->meshType <= moab::TempestRemapper::ICO )
            std::cout << "\n  Resolution:             " << this->blockSize;
        if( this->meshType == moab::TempestRemapper::ICO && computeDual )
            std::cout << "\n  Compute dual:           " << ( this->computeDual ? "Yes" : "No" );

        if( computeWeights )
        {
            std::cout << "\n  Gnomonic projection:    " << ( this->useGnomonicProjection ? "Yes" : "No" );
            std::cout << "\n  Intersection algorithm: " << ( this->kdtreeSearch ? "KdTree search" : "Advancing front" );

            // Discretization settings
            std::cout << "\n\nDiscretization:";
            std::cout << "\n  Source:             " << this->disc_methods[0] << " (order " << this->disc_orders[0]
                      << ")";
            std::cout << "\n  Target:             " << this->disc_methods[1] << " (order " << this->disc_orders[1]
                      << ")";

            // Remapping options
            std::cout << "\n\nRemapping Options:";
            std::cout << "\n  Method:             "
                      << ( this->mapOptions.strMethod.empty() ? "Default" : this->mapOptions.strMethod );
            std::cout << "\n  Monotonicity:       " << ( this->ensureMonotonicity ? "Yes" : "No" );
            std::cout << "\n  Volumetric:         " << ( this->fVolumetric ? "Yes" : "No" );
            std::cout << "\n  Check consistency:  " << ( this->fCheck ? "Yes" : "No" );
            std::cout << "\n  Skip intersection:  " << ( this->skip_intersection ? "Yes" : "No" );
        }

        // Parallel configuration
        std::cout << "\n\nParallel Configuration:";
        std::cout << "\n  MPI Processes:          " << this->n_procs;
        if( this->meshType > moab::TempestRemapper::ICO ) std::cout << "\n  Number of Ghost Layers: " << this->nlayers;

        std::cout << "\n\n" << std::string( width, '=' ) << "\n\n";
    }

    /**
     * @brief Configure map options based on command line parameters
     * @param nlayer_input Number of ghost layers
     */
    void configureMapOptions( int nlayer_input )
    {
        // Set polynomial orders with bounds checking
        this->mapOptions.nPin  = ( this->disc_orders.empty() ) ? 1 : this->disc_orders[0];
        this->mapOptions.nPout = ( this->disc_orders.size() > 1 ) ? this->disc_orders[1] : this->mapOptions.nPin;

        // Initialize flags
        this->mapOptions.fSourceConcave = false;
        this->mapOptions.fTargetConcave = false;
        this->mapOptions.strMethod.clear();

        // Configure finite volume method if specified
        if( this->fvMethod != "none" )
        {
            this->mapOptions.strMethod       = this->fvMethod + ";";
            this->mapOptions.fNoConservation = true;
        }

        // Configure monotonicity with validation
        this->ensureMonotonicity = std::max( 0, std::min( 3, this->ensureMonotonicity ) );  // Clamp to 0-3
        switch( this->ensureMonotonicity )
        {
            case 0:
                this->mapOptions.fMonotone = false;
                break;
            case 3:
                this->mapOptions.strMethod += "mono3;";
                this->mapOptions.fMonotone = true;
                break;
            case 2:
                this->mapOptions.strMethod += "mono2;";
                this->mapOptions.fMonotone = true;
                break;
            case 1:
            default:
                this->mapOptions.fMonotone = true;
                break;
        }

        // Set other options
        this->mapOptions.fNoCorrectAreas = false;
        this->mapOptions.fNoCheck        = !this->fCheck;

        // Add volumetric flag if needed
        if( this->fVolumetric )
        {
            this->mapOptions.strMethod += "volumetric;";
        }

        // Set number of ghost layers based on method and order.
        // FV order 1 needs 0 ghost layers; FV order p > 1 needs p+1 ghost layers.
        if( this->fvMethod == "delaunay" || this->fvMethod == "bilin" )
        {
            this->skip_intersection = true;
            this->nlayers           = 3;  // conservative
        }
        else
        {
            // order 1: no ghost layers
            // order p: p+1 layers (again, being conservative)
            this->nlayers = ( this->mapOptions.nPin > 1 ) ? this->mapOptions.nPin + 1 : 0;
        }

        // User-supplied value always overrides the internal default (even 0 is valid).
        if( nlayer_input >= 0 )
        {
            this->nlayers = nlayer_input;
        }

        // Configure output
        this->mapOptions.strOutputMapFile = this->outFilename;
        this->mapOptions.strOutputFormat  = "Netcdf4";
    }
};

// Forward declare some methods
static moab::ErrorCode CreateTempestMesh( ToolContext&, moab::TempestRemapper& remapper, Mesh* );
static inline constexpr double sample_constant( double dLon, double dLat ) noexcept;
static inline double sample_slow_harmonic( double dLon, double dLat ) noexcept;
static inline double sample_fast_harmonic( double dLon, double dLat ) noexcept;
static inline double sample_stationary_vortex( double dLon, double dLat ) noexcept;

/////////////////////////////////////////////////////////////

//#define MOAB_DBG
int main( int argc, char* argv[] )
{
    try
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
    moab::IntxAreaUtils areaAdaptor( moab::IntxAreaUtils::lHuiller );

    Mesh* tempest_mesh = new Mesh();
    MB_CHK_SET_ERR( CreateTempestMesh( *runCtx, remapper, tempest_mesh ), "Failed to create tempest mesh" );

    if( runCtx->meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Compute intersections with MOAB
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert( runCtx->meshes.size() == 3 );

#ifdef MOAB_HAVE_MPI
        MB_CHK_SET_ERR( pcomm->check_all_shared_handles(), "Failed to check all shared handles" );
#endif

        // Load the meshes and validate
        MB_CHK_SET_ERR( remapper.ConvertTempestMesh( moab::Remapper::SourceMesh ), "Failed to convert source mesh" );
        MB_CHK_SET_ERR( remapper.ConvertTempestMesh( moab::Remapper::TargetMesh ), "Failed to convert target mesh" );
        MB_CHK_SET_ERR( remapper.ConvertTempestMesh( moab::Remapper::OverlapMesh ), "Failed to convert overlap mesh" );
        if( !runCtx->skip_io )
        {
            MB_CHK_SET_ERR( mbCore->write_mesh( "tempest_intersection.h5m", &runCtx->meshsets[2], 1 ),
                            "Failed to write TempestRemap intersection mesh in MOAB format" );
        }

        // print verbosely about the problem setting
        size_t velist[6], gvelist[6];
        {
            moab::Range rintxverts, rintxelems;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 0, rintxverts ),
                            "Failed to get vertices" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 2, rintxelems ),
                            "Failed to get elements" );
            velist[0] = rintxverts.size();
            velist[1] = rintxelems.size();

            moab::Range bintxverts, bintxelems;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 0, bintxverts ),
                            "Failed to get vertices" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, bintxelems ),
                            "Failed to get elements" );
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
            MB_CHK_SET_ERR( mbintx->FindMaxEdges( runCtx->meshsets[0], runCtx->meshsets[1] ),
                            "Failed to find max edges" );

#ifdef MOAB_HAVE_MPI
            moab::Range local_verts;
            MB_CHK_SET_ERR( mbintx->build_processor_euler_boxes( runCtx->meshsets[1], local_verts ),
                            "Failed to build processor euler boxes" );

            runCtx->timer_pop();

            moab::EntityHandle covering_set;
            runCtx->timer_push( "communicate the mesh" );
            // we compute just intersection here, no need for extra ghost layers anyway
            // ghost layers are needed in coverage for bilinear map, which does not actually need intersection
            // this will be fixed in the future, bilinear map needs just coverage, not intersection
            // so I am not passing the ghost layer here, even though there is an option in runCtx for a ghost layer
            // NOTE: This is a communication-heavy kernel if mesh is distributed very differently
            MB_CHK_SET_ERR( mbintx->construct_covering_set( runCtx->meshsets[0], covering_set ),
                            "Failed to construct covering set" );
            runCtx->timer_pop();

            // print verbosely about the problem setting
            {
                moab::Range cintxverts, cintxelems;
                MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( covering_set, 0, cintxverts ),
                                "Failed to get vertices" );
                MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( covering_set, 2, cintxelems ),
                                "Failed to get elements" );
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
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( intxset, 2, intxelems ), "Failed to get elements" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( intxset, 0, intxverts, true ),
                            "Failed to get vertices" );
            outputFormatter.printf( 0, "The intersection set contains %lu elements and %lu vertices \n",
                                    intxelems.size(), intxverts.size() );

            double initial_sarea =
                areaAdaptor.area_on_sphere( mbCore, runCtx->meshsets[0],
                                            radius_src );  // use the target to compute the initial area
            double initial_tarea =
                areaAdaptor.area_on_sphere( mbCore, runCtx->meshsets[1],
                                            radius_dest );  // use the target to compute the initial area
            double intx_area = areaAdaptor.area_on_sphere( mbCore, intxset, radius_src );

            outputFormatter.printf( 0, "mesh areas: source = %12.10f, target = %12.10f, intersection = %12.10f \n",
                                    initial_sarea, initial_tarea, intx_area );
            outputFormatter.printf( 0, "relative error w.r.t source = %12.10e, target = %12.10e \n",
                                    fabs( intx_area - initial_sarea ) / initial_sarea,
                                    fabs( intx_area - initial_tarea ) / initial_tarea );
        }

        // Write out our computed intersection file
        if( !runCtx->skip_io )
        {
            MB_CHK_SET_ERR( mbCore->write_mesh( "moab_intersection.h5m", &intxset, 1 ),
                            "Failed to write the intersection" );
        }

        if( runCtx->computeWeights )
        {
            runCtx->timer_push( "compute weights with the Tempest meshes" );
            // Call to generate an offline map with the tempest meshes
            OfflineMap weightMap;
            if( GenerateOfflineMapWithMeshes( *runCtx->meshes[0], *runCtx->meshes[1], *runCtx->meshes[2],
                                              runCtx->disc_methods[0],  // std::string strInputType
                                              runCtx->disc_methods[1],  // std::string strOutputType,
                                              runCtx->mapOptions, weightMap ) != 0 )
                throw std::runtime_error( "Could not generate offline map with TempestRemap" );
            runCtx->timer_pop();

            std::map< std::string, std::string > mapAttributes;
            if( !runCtx->skip_io ) weightMap.Write( "outWeights.nc", mapAttributes );
        }
    }
    else if( runCtx->meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m
#ifdef MOAB_HAVE_MPI
        MB_CHK_SET_ERR( pcomm->check_all_shared_handles(), "Checking shared handles failed." );
#endif

        // print verbosely about the problem setting
        size_t velist[4] = { 0, 0, 0, 0 }, gvelist[4] = { 0, 0, 0, 0 };
        {
            moab::Range srcverts, srcelems;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 0, srcverts ),
                            "Failed to get vertices" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[0], 2, srcelems ),
                            "Failed to get elements" );
            MB_CHK_SET_ERR( moab::IntxUtils::fix_degenerate_quads( mbCore, runCtx->meshsets[0] ),
                            "Failed to fix degenerate quads" );
            if( runCtx->enforceConvexity )
            {
                MB_CHK_SET_ERR( moab::IntxUtils::enforce_convexity( mbCore, runCtx->meshsets[0], proc_id ),
                                "Failed to enforce convexity" );
            }
            MB_CHK_SET_ERR( areaAdaptor.positive_orientation( mbCore, runCtx->meshsets[0], radius_src ),
                            "Failed to enforce positive orientation" );
            velist[0] = srcverts.size();
            velist[1] = srcelems.size();

            moab::Range tgtverts, tgtelems;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 0, tgtverts ),
                            "Failed to get vertices" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, tgtelems ),
                            "Failed to get elements" );
            MB_CHK_SET_ERR( moab::IntxUtils::fix_degenerate_quads( mbCore, runCtx->meshsets[1] ),
                            "Failed to fix degenerate quads" );
            if( runCtx->enforceConvexity )
            {
                MB_CHK_SET_ERR( moab::IntxUtils::enforce_convexity( mbCore, runCtx->meshsets[1], proc_id ),
                                "Failed to enforce convexity" );
            }
            MB_CHK_SET_ERR( areaAdaptor.positive_orientation( mbCore, runCtx->meshsets[1], radius_dest ),
                            "Failed to enforce positive orientation" );
            velist[2] = tgtverts.size();
            velist[3] = tgtelems.size();
        }
        //MB_CHK_SET_ERR( mbCore->write_file( "source_mesh.h5m", nullptr, writeOptions, &runCtx->meshsets[0], 1 ), "Could not write source mesh" );
        //MB_CHK_SET_ERR( mbCore->write_file( "target_mesh.h5m", nullptr, writeOptions, &runCtx->meshsets[1], 1 ), "Could not write target mesh" );

        // if( runCtx->nlayers && nprocs > 1 )
        // {
        //     remapper.ResetMeshSet( moab::Remapper::SourceMesh, runCtx->meshsets[3] );
        //     runCtx->meshes[0] = remapper.GetMesh( moab::Remapper::SourceMesh );  //  ?
        // }

        // First compute the covering set such that the target elements are fully covered by the
        // local source grid
        runCtx->timer_push( "construct covering set for intersection" );
        // if ghosting, do not use gnomonic projection
        if( runCtx->nlayers > 0 ) runCtx->useGnomonicProjection = false;
        MB_CHK_SET_ERR( remapper.ConstructCoveringSet( runCtx->epsrel, 1.0, 1.0, runCtx->boxeps, runCtx->rrmGrids,
                                                       runCtx->useGnomonicProjection, runCtx->nlayers ),
                        "Failed to construct covering set" );
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

        if( runCtx->skip_intersection )
        {
            if( !proc_id ) outputFormatter.printf( 0, "Skipping mesh intersection computation.\n" );
        }
        else
        {
            // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
            runCtx->timer_push( "setup and compute mesh intersections" );
            MB_CHK_SET_ERR( remapper.ComputeOverlapMesh( runCtx->kdtreeSearch, false ),
                            "Failed to compute mesh intersections" );
            runCtx->timer_pop();
        }

        // print some diagnostic checks to see if the overlap grid resolved the input meshes
        // correctly
        // Compute ghost overlap elements once; reused for both area diagnostics and intx file write
        moab::Range ghostOverlapElems;
#ifdef MOAB_HAVE_MPI
        if( nprocs > 1 && !runCtx->skip_intersection )
            MB_CHK_SET_ERR( remapper.GetOverlapAugmentedEntities( ghostOverlapElems ),
                            "Failed to get ghost overlap entities" );
#endif

        double dTotalOverlapArea = 0.0;
        if( runCtx->print_diagnostics && !runCtx->skip_intersection )
        {
             // Areas for source, target, overlap meshes
            double local_areas[3]  = { 0, 0, 0 },
                   global_areas[3] = { 0, 0, 0 };

            // Helper: compute area of a meshset excluding cells with GRID_IMASK==0.
            // Both source and target SCRIP grids may have a land/sea mask; the intersection
            // only covers unmasked cells, so comparing full-mesh areas gives a misleading error.
            auto area_unmasked = [&]( moab::EntityHandle meshset, double radius ) -> double {
                moab::Tag imaskTag = 0;
                mbCore->tag_get_handle( "GRID_IMASK", imaskTag );
                if( !imaskTag ) return areaAdaptor.area_on_sphere( mbCore, meshset, radius );
                moab::Range cells;
                mbCore->get_entities_by_dimension( meshset, 2, cells );
                std::vector< int > masks( cells.size(), 1 );
                mbCore->tag_get_data( imaskTag, cells, masks.data() );
                moab::Range maskedCells;
                size_t idx = 0;
                for( auto it = cells.begin(); it != cells.end(); ++it, ++idx )
                    if( !masks[idx] ) maskedCells.insert( *it );
                moab::Range unmasked = moab::subtract( cells, maskedCells );
                moab::EntityHandle tmpSet;
                mbCore->create_meshset( moab::MESHSET_SET, tmpSet );
                mbCore->add_entities( tmpSet, unmasked );
                double area = areaAdaptor.area_on_sphere( mbCore, tmpSet, radius );
                mbCore->delete_entities( &tmpSet, 1 );
                return area;
            };

            local_areas[0] = area_unmasked( runCtx->meshsets[0], radius_src );
            local_areas[1] = area_unmasked( runCtx->meshsets[1], radius_dest );
            // Exclude ghost overlap elements from area sum to avoid double-counting after MPI_Allreduce
            {
                moab::Range ownedOverlapElems;
                MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[2], 2, ownedOverlapElems ),
                                "Failed to get overlap elements" );
                ownedOverlapElems = moab::subtract( ownedOverlapElems, ghostOverlapElems );
                moab::EntityHandle ownedOverlapSet;
                MB_CHK_SET_ERR( mbCore->create_meshset( moab::MESHSET_SET, ownedOverlapSet ),
                                "Can't create owned overlap meshset" );
                MB_CHK_SET_ERR( mbCore->add_entities( ownedOverlapSet, ownedOverlapElems ),
                                "Can't add owned overlap elements" );
                local_areas[2] = areaAdaptor.area_on_sphere( mbCore, ownedOverlapSet, radius_src );
                MB_CHK_SET_ERR( mbCore->delete_entities( &ownedOverlapSet, 1 ), "Can't delete temp meshset" );
            }

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
                outputFormatter.printf( 0, "relative error w.r.t source = %12.14e, and target = %12.14e\n",
                                        fabs( global_areas[0] - global_areas[2] ) / global_areas[0],
                                        fabs( global_areas[1] - global_areas[2] ) / global_areas[1] );
            }
            dTotalOverlapArea = global_areas[2];
        }

        if( runCtx->intxFilename.size() && !runCtx->skip_intersection )
        {
            moab::EntityHandle writableOverlapSet;
            MB_CHK_SET_ERR( mbCore->create_meshset( moab::MESHSET_SET, writableOverlapSet ), "Can't create new set" );
            moab::EntityHandle meshOverlapSet = remapper.GetMeshSet( moab::Remapper::OverlapMesh );
            moab::Range ovEnts;
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( meshOverlapSet, 2, ovEnts ), "Can't create new set" );
            MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( meshOverlapSet, 0, ovEnts ), "Can't create new set" );

#ifdef MOAB_HAVE_MPI
            // Exclude ghost overlap elements from the write: each ghost element is owned by another
            // rank and will be written from there. Including ghosts here causes duplicate entity
            // handles in the parallel HDF5 output and deadlocks the collective write.
            if( nprocs > 1 )
            {
                ovEnts = moab::subtract( ovEnts, ghostOverlapElems );
#ifdef MOAB_DBG
                if( !runCtx->skip_io )
                {
                    std::stringstream filename;
                    filename << "aug_overlap" << runCtx->pcomm->rank() << ".h5m";
                    MB_CHK_SET_ERR( mbCore->write_file( filename.str().c_str(), 0, 0, &meshOverlapSet, 1 ),
                                    "Failed to write the overlap set" );
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
                MB_CHK_SET_ERR( mbCore->write_file( filename.str().c_str(), 0, 0, &writableOverlapSet, 1 ),
                                "Failed to write the writable overlap set" );
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
                MB_CHK_SET_ERR( mbCore->write_file( sstr.str().c_str(), nullptr, writeOptions, &writableOverlapSet, 1 ),
                                "Failed to write the writable overlap set" );
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
            MB_CHK_SET_ERR( weightMap->GenerateRemappingWeights(
                                runCtx->disc_methods[0],  // std::string strInputType
                                runCtx->disc_methods[1],  // std::string strOutputType,
                                runCtx->mapOptions,       // const GenerateOfflineMapAlgorithmOptions& options
                                runCtx->doftag_names[0],  // const std::string& source_tag_name
                                runCtx->doftag_names[1]   // const std::string& target_tag_name
                                ),
                            "Failed to generate remapping weights" );
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

            if( runCtx->outFilename.size() && !runCtx->skip_io )
            {
                std::map< std::string, std::string > attrMap;
                attrMap["MOABversion"]   = std::string( MOAB_PACKAGE_VERSION );
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
                MB_CHK_SET_ERR( weightMap->WriteParallelMap( runCtx->outFilename.c_str(), attrMap ),
                                "Failed writing the parallel map to disk" );
            }

            if( runCtx->verifyWeights )
            {
                // Let us pick a sampling test function for solution evaluation
                // SH, SV, FH, C, USERVAR
                bool userVariable = false;
                moab::TempestOnlineMap::sample_function testFunction;
                if( !runCtx->variableToVerify.compare( "SH" ) )
                    testFunction = &sample_slow_harmonic;
                else if( !runCtx->variableToVerify.compare( "FH" ) )
                    testFunction = &sample_fast_harmonic;
                else if( !runCtx->variableToVerify.compare( "SV" ) )
                    testFunction = &sample_stationary_vortex;
                else if( !runCtx->variableToVerify.compare( "C" ) )
                    testFunction = &sample_constant;
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
                    // MB_CHK_SET_ERR( mbCore->tag_get_handle( runCtx->variableToVerify.c_str(), srcAnalyticalFunction ),
                    //                 "Failed to get analytical solution on source grid" );
                    MB_CHK_SET_ERR( weightMap->DefineAnalyticalSolution( srcAnalyticalFunction,
                                                                         "AnalyticalSolnSrcExact",
                                                                         moab::Remapper::SourceMesh, testFunction ),
                                    "Failed to define analytical solution on source grid" );
                    runCtx->timer_pop();

                    // runCtx->timer_push( "exchange solution on source grid" );
                    // moab::Range& srccovEnts = remapper.GetMeshEntities( moab::Remapper::CoveringMesh );
                    // MB_CHK_SET_ERR( pcomm->exchange_tags( srcAnalyticalFunction, srccovEnts ),
                    //                 "Failed to exchange analytical solution on source grid" );
                    // runCtx->timer_pop();

                    runCtx->timer_push( "describe a solution on target grid" );
                    MB_CHK_SET_ERR( weightMap->DefineAnalyticalSolution(
                                        tgtAnalyticalFunction, "AnalyticalSolnTgtExact", moab::Remapper::TargetMesh,
                                        testFunction, &tgtProjectedFunction, "ProjectedSolnTgt" ),
                                    "Failed to define analytical solution on target grid" );
                    runCtx->timer_pop();
                }
                else
                {
                    MB_CHK_SET_ERR( mbCore->tag_get_handle( runCtx->variableToVerify.c_str(), srcAnalyticalFunction ),
                                    "Failed to get analytical solution on source grid" );
                    MB_CHK_SET_ERR( mbCore->tag_get_handle( "ProjectedSolnTgt", 1, moab::MB_TYPE_DOUBLE,
                                                            tgtProjectedFunction,
                                                            moab::MB_TAG_DENSE | moab::MB_TAG_CREAT ),
                                    "Failed to get projected solution on target grid" );
                }

                // if( !runCtx->skip_io )
                {
                    MB_CHK_SET_ERR( mbCore->write_file( "srcWithSolnTag.h5m", nullptr, writeOptions,
                                                        &runCtx->meshsets[0], 1 ),
                                    "Failed to write the source mesh with solution tag" );
                }

                runCtx->timer_push( "compute solution projection on target grid" );
                MB_CHK_SET_ERR( weightMap->ApplyWeights( srcAnalyticalFunction, tgtProjectedFunction, false,
                                                         runCtx->cassType ),
                                "Failed to apply weights" );
                runCtx->timer_pop();

                // if( !runCtx->skip_io )
                {
                    MB_CHK_SET_ERR( mbCore->write_file( "tgtWithSolnTag2.h5m", nullptr, writeOptions,
                                                        &runCtx->meshsets[1], 1 ),
                                    "Failed to write the target mesh with projected solution tag" );
                }

                if( nprocs == 1 && runCtx->baselineFile.size() )
                {
                    // save the field from tgtWithSolnTag2 in a text file, and global ids for cells
                    moab::Range tgtEntities;
                    if( runCtx->disc_methods[1] == "pcloud" )
                    {
                        MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 0, tgtEntities ),
                                        "Failed to get entities by dimension" );
                    }
                    else
                    {
                        MB_CHK_SET_ERR( mbCore->get_entities_by_dimension( runCtx->meshsets[1], 2, tgtEntities ),
                                        "Failed to get entities by dimension" );
                    }
                    std::vector< int > globIds( tgtEntities.size() );
                    std::vector< double > vals( tgtEntities.size() );
                    moab::Tag projTag;
                    MB_CHK_SET_ERR( mbCore->tag_get_handle( "ProjectedSolnTgt", projTag ),
                                    "Failed to get projected solution tag" );
                    moab::Tag gid = mbCore->globalId_tag();
                    MB_CHK_SET_ERR( mbCore->tag_get_data( gid, tgtEntities, &globIds[0] ), "Failed to get global ids" );
                    MB_CHK_SET_ERR( mbCore->tag_get_data( projTag, tgtEntities, &vals[0] ),
                                    "Failed to get projected solution" );
                    std::fstream fs;
                    fs.open( runCtx->baselineFile.c_str(), std::fstream::out );
                    fs << std::setprecision( 15 );  // maximum precision for doubles
                    for( size_t i = 0; i < tgtEntities.size(); i++ )
                        fs << globIds[i] << " " << vals[i] << "\n";
                    fs.close();
                    // for good measure, save the source file too, with the tag AnalyticalSolnSrcExact
                    // it will be used later to test, along with a target file
                    if( !runCtx->skip_io )
                    {
                        MB_CHK_SET_ERR( mbCore->write_file( "srcWithSolnTag.h5m", nullptr, writeOptions,
                                                            &runCtx->meshsets[0], 1 ),
                                        "Failed to write the source mesh with solution tag" );
                    }
                }

                // compute error metrics if it is a known analytical functional
                if( !userVariable )
                {
                    runCtx->timer_push( "compute error metrics against analytical solution on target grid" );
                    std::map< std::string, double > errMetrics;
                    MB_CHK_SET_ERR( weightMap->ComputeMetrics( moab::Remapper::TargetMesh, tgtAnalyticalFunction,
                                                               tgtProjectedFunction, errMetrics, true ),
                                    "Failed to compute error metrics" );
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
    return 0;
    }
    catch( const std::exception& e )
    {
        std::cerr << "[mbtempest] Fatal error: " << e.what() << std::endl;
#ifdef MOAB_HAVE_MPI
        MPI_Abort( MPI_COMM_WORLD, 1 );
#endif
        return 1;
    }
    catch( ... )
    {
        std::cerr << "[mbtempest] Fatal: unknown exception caught" << std::endl;
#ifdef MOAB_HAVE_MPI
        MPI_Abort( MPI_COMM_WORLD, 1 );
#endif
        return 1;
    }
}

///////////////////////////////////////////////////////////////////////////////

// Helper functions for each mesh type
namespace
{

#define TR_CHK_SET_ERR( err, msg )                                                \
    if( err )                                                                     \
    {                                                                             \
        std::cout << "MOAB-TempestRemap Failure. ErrorCode (" << ( err ) << ") "; \
        MB_CHK_SET_ERR( moab::MB_FAILURE, msg );                                  \
    }

moab::ErrorCode handleOverlapMemory( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    using namespace moab;

    // resize the meshsets and meshes vectors
    ctx.meshsets.resize( 3 );
    ctx.meshes.resize( 3 );

    ctx.meshsets[0] = remapper.GetMeshSet( Remapper::SourceMesh );
    ctx.meshsets[1] = remapper.GetMeshSet( Remapper::TargetMesh );
    ctx.meshsets[2] = remapper.GetMeshSet( Remapper::OverlapMesh );

    // Load and process source mesh
    MB_CHK_SET_ERR( remapper.LoadMesh( Remapper::SourceMesh, ctx.inFilenames[0], TempestRemapper::DEFAULT ),
                    "Failed to load MOAB Source mesh" );

    // Load and process target mesh
    MB_CHK_SET_ERR( remapper.LoadMesh( Remapper::TargetMesh, ctx.inFilenames[1], TempestRemapper::DEFAULT ),
                    "Failed to load MOAB Target mesh" );

    // Generate overlap mesh
    TR_CHK_SET_ERR( GenerateOverlapWithMeshes( *ctx.meshes[0], *ctx.meshes[1], *tempest_mesh, "", "NetCDF4", "exact",
                                               false ),
                    "Failed to generate TempestRemap OverlapMesh" );

    remapper.SetMesh( Remapper::OverlapMesh, tempest_mesh );
    ctx.meshes[2] = remapper.GetMesh( Remapper::OverlapMesh );

    return moab::MB_SUCCESS;
}

moab::ErrorCode handleOverlapMOAB( ToolContext& ctx, moab::TempestRemapper& remapper )
{
    using namespace moab;

    // resize the meshsets and meshes vectors
    ctx.meshsets.resize( 3 );
    ctx.meshes.resize( 3 );

    ctx.meshsets[0] = remapper.GetMeshSet( Remapper::SourceMesh );
    ctx.meshsets[1] = remapper.GetMeshSet( Remapper::TargetMesh );
    ctx.meshsets[2] = remapper.GetMeshSet( Remapper::OverlapMesh );

    constexpr double radius_src  = 1.0;
    constexpr double radius_dest = 1.0;

    // Load and process target mesh
    {
        std::vector< int > metadata;
        std::string additional_read_opts_tgt = ctx.get_file_read_options( ctx.inFilenames[1] );
        if( ctx.n_procs > 1 && ctx.disc_methods[1].compare( "fv" ) != 0 )  // target discretization is cgll or dgll
        {
            // auto pcomm = new ParallelComm( ctx.mbcore, MPI_COMM_WORLD );
            // add one ghost layer to the target mesh
            // additional_read_opts_tgt = additional_read_opts_tgt +  "PARALLEL_GHOSTS=3.0.2;PARALLEL_THIN_GHOST_LAYER;SKIP_AUGMENT_WITH_GHOSTS;PRINT_PARALLEL;";
            // additional_read_opts_tgt = additional_read_opts_tgt +  "PARALLEL_COMM=1;";
            // additional_read_opts_tgt = additional_read_opts_tgt +  "PARALLEL_GHOSTS=3.0.1;";
            // additional_read_opts_tgt = additional_read_opts_tgt +  "PARALLEL_COMM=" + std::to_string(ctx.pcomm->get_id()) + ";";
        }

        MB_CHK_SET_ERR( remapper.LoadNativeMesh( ctx.inFilenames[1], ctx.meshsets[1], metadata,
                                                 additional_read_opts_tgt.c_str() ),
                        "Failed to load MOAB Target mesh" );

#ifdef MOAB_HAVE_MPI
        if( ctx.n_procs > 1 && ctx.disc_methods[1].compare( "fv" ) != 0 &&
            false )  // target discretization is cgll or dgll
        {
            Range beforeGhost, afterGhost;
            ctx.mbcore->get_entities_by_dimension( ctx.meshsets[1], 2, beforeGhost );

            ctx.pcomm->set_debug_verbosity( 5 );
            MB_CHK_SET_ERR( ctx.pcomm->exchange_ghost_cells( 2, 0, 1, 0, true, true, &ctx.meshsets[1] ),
                            "Failed to exchange ghost cells for MOAB Target mesh" );
            ctx.pcomm->set_debug_verbosity( 0 );

            ctx.mbcore->get_entities_by_dimension( ctx.meshsets[1], 2, afterGhost );
            std::cout << ctx.proc_id << ": N(before) = " << beforeGhost.size() << ", N(after) = " << afterGhost.size()
                      << std::endl;

            std::vector< Tag > taglist;
            taglist.push_back( ctx.mbcore->globalId_tag() );
            Tag gdofTag;
            MB_CHK_SET_ERR( ctx.mbcore->tag_get_handle( "GLOBAL_DOFS", gdofTag ),
                            "Failed to get global dofs tag for MOAB Target mesh" );
            taglist.push_back( gdofTag );
            MB_CHK_SET_ERR( ctx.pcomm->exchange_tags( taglist, taglist, afterGhost ),
                            "Failed to exchange global dofs for MOAB Target mesh" );
            // std::set< unsigned int > commprocs;
            // MB_CHK_SET_ERR( ctx.pcomm->get_comm_procs( commprocs ),
            //                 "Failed to get commprocs for MOAB Target mesh" );
            // if (ctx.proc_id == 0)
            // {
            //     std::cout << ctx.proc_id << ": commprocs = [";
            //     for( auto p : commprocs ) std::cout << p << ", ";
            //     std::cout << "]\n";

            //     std::cout << ctx.proc_id << ": N(after) = " << afterGhost.size() << std::endl;
            //     for (auto eh: afterGhost)
            //     {
            //         std::cout << ctx.mbcore->type_from_handle(eh) << ": " << eh << std::endl;
            //     }
            // }
        }
#endif

        if( !metadata.empty() )
        {
            remapper.SetMeshType( Remapper::TargetMesh, metadata );
        }

        MB_CHK_SET_ERR( IntxUtils::ScaleToRadius( ctx.mbcore, ctx.meshsets[1], radius_dest ),
                        "Failed to preprocess MOAB Target mesh" );
    }

    // Load and process source mesh
    {
        std::vector< int > metadata;
        auto additional_read_opts_src = ctx.get_file_read_options( ctx.inFilenames[0] );
#ifdef MOAB_HAVE_MPI
        if( ctx.n_procs > 1 )
        {
            // auto pcomm = new ParallelComm( ctx.mbcore, MPI_COMM_WORLD );
            additional_read_opts_src =
                additional_read_opts_src + "PARALLEL_COMM=" + std::to_string( ctx.pcomm->get_id() ) + ";";
        }
#endif
        MB_CHK_SET_ERR( remapper.LoadNativeMesh( ctx.inFilenames[0], ctx.meshsets[0], metadata,
                                                 additional_read_opts_src.c_str() ),
                        "Failed to load MOAB Source mesh" );

        if( !metadata.empty() )
        {
            remapper.SetMeshType( Remapper::SourceMesh, metadata );
        }

        MB_CHK_SET_ERR( IntxUtils::ScaleToRadius( ctx.mbcore, ctx.meshsets[0], radius_src ),
                        "Failed to preprocess MOAB Source mesh" );
    }

    if( ctx.computeWeights )
    {
        // Convert MOAB to TempestRemap meshes
        MB_CHK_SET_ERR( remapper.ConvertMeshToTempest( Remapper::SourceMesh ),
                        "Failed to convert MOAB Source mesh to TempestRemap mesh" );
        ctx.meshes[0] = remapper.GetMesh( Remapper::SourceMesh );

        MB_CHK_SET_ERR( remapper.ConvertMeshToTempest( Remapper::TargetMesh ),
                        "Failed to convert MOAB Target mesh to TempestRemap mesh" );
        ctx.meshes[1] = remapper.GetMesh( Remapper::TargetMesh );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode handleOverlapFiles( ToolContext& ctx, Mesh* tempest_mesh )
{
    ctx.timer_push( "create Tempest OverlapMesh" );
    TR_CHK_SET_ERR( GenerateOverlapMesh( ctx.inFilenames[0], ctx.inFilenames[1], *tempest_mesh, ctx.outFilename,
                                         "NetCDF4", "exact", true ),
                    "Failed to create Tempest OverlapMesh" );
    ctx.timer_pop();

    // Add the overlap mesh to the list of meshes
    ctx.meshes.push_back( tempest_mesh );
    return moab::MB_SUCCESS;
}

/**
 * @brief Convert a generated TempestRemap mesh to MOAB format and write as h5m file.
 *
 * When the output filename has a .h5m extension, the mesh is converted from TempestRemap
 * format to MOAB format in memory and written as a native MOAB HDF5 file. This allows
 * the generated meshes to be loaded by mbtempest type 5 (OVERLAP_MOAB) workflows.
 * The TempestRemap format file is still written (with .g extension) for compatibility.
 */
moab::ErrorCode convertAndWriteMOABMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    // Check if output filename has .h5m extension
    const std::string& outFile = ctx.outFilename;
    const size_t dot           = outFile.find_last_of( "." );
    if( dot == std::string::npos ) return moab::MB_SUCCESS;

    const std::string ext = outFile.substr( dot + 1 );
    if( ext != "h5m" ) return moab::MB_SUCCESS;

    // Register the TempestRemap mesh with the remapper as SourceMesh
    remapper.SetMesh( moab::Remapper::SourceMesh, tempest_mesh, false );

    // Convert TempestRemap mesh to MOAB format
    ctx.timer_push( "convert TempestRemap mesh to MOAB format" );
    MB_CHK_SET_ERR( remapper.ConvertTempestMesh( moab::Remapper::SourceMesh ),
                    "Failed to convert TempestRemap mesh to MOAB format" );
    ctx.timer_pop();

    // Fix degenerate quads: RLL meshes from TempestRemap have polar cells stored as
    // 4-node quads with duplicate vertices. Convert these to proper triangles so the
    // intersection algorithm can handle them correctly.
    moab::EntityHandle meshSet = remapper.GetMeshSet( moab::Remapper::SourceMesh );
    MB_CHK_SET_ERR( moab::IntxUtils::fix_degenerate_quads( ctx.mbcore, meshSet ),
                    "Failed to fix degenerate quads in converted mesh" );
    ctx.timer_push( "write MOAB mesh to h5m file" );
    MB_CHK_SET_ERR( ctx.mbcore->write_file( outFile.c_str(), nullptr, nullptr, &meshSet, 1 ),
                    "Failed to write MOAB mesh to h5m file" );
    ctx.timer_pop();

    if( !ctx.proc_id )
        ctx.outputFormatter.printf( 0, "Wrote MOAB mesh to %s\n", outFile.c_str() );

    return moab::MB_SUCCESS;
}

moab::ErrorCode handleICOMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    std::string trFilename = ctx.outFilename;
    const size_t dot       = trFilename.find_last_of( "." );
    if( dot != std::string::npos && trFilename.substr( dot + 1 ) == "h5m" )
        trFilename = trFilename.substr( 0, dot ) + ".g";

    ctx.timer_push( "generate ICO mesh with TempestRemap" );
    TR_CHK_SET_ERR( GenerateICOMesh( *tempest_mesh, ctx.blockSize, ctx.computeDual, trFilename, "NetCDF4" ),
                    "Failed to generate ICO mesh with TempestRemap" );
    ctx.timer_pop();

    // Add the ICO mesh to the list of meshes
    ctx.meshes.push_back( tempest_mesh );

    MB_CHK_SET_ERR( convertAndWriteMOABMesh( ctx, remapper, tempest_mesh ),
                    "Failed to convert and write MOAB mesh" );

    return moab::MB_SUCCESS;
}

moab::ErrorCode handleRLLMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    std::string trFilename = ctx.outFilename;
    const size_t dot       = trFilename.find_last_of( "." );
    if( dot != std::string::npos && trFilename.substr( dot + 1 ) == "h5m" )
        trFilename = trFilename.substr( 0, dot ) + ".g";

    ctx.timer_push( "generate RLL mesh with TempestRemap" );
    TR_CHK_SET_ERR( GenerateRLLMesh( *tempest_mesh,                     // Mesh& meshOut,
                                     ctx.blockSize * 2, ctx.blockSize,  // int nLongitudes, int nLatitudes,
                                     0.0, 360.0,                        // double dLonBegin, double dLonEnd,
                                     -90.0, 90.0,                       // double dLatBegin, double dLatEnd,
                                     false, false, false,  // bool fGlobalCap, bool fFlipLatLon, bool fForceGlobal,
                                     "" /*ctx.inFilename*/,
                                     "",             // std::string strInputFile, std::string strInputFileLonName
                                     "",             // std::string strInputFileLatName
                                     trFilename,     // std::string strOutputFile
                                     "NetCDF4",      // std::string strOutputFormat
                                     true            // bool fVerbose
                                     ),
                    "Failed to generate RLL mesh with TempestRemap" );
    ctx.timer_pop();

    // Add the RLL mesh to the list of meshes
    ctx.meshes.push_back( tempest_mesh );

    MB_CHK_SET_ERR( convertAndWriteMOABMesh( ctx, remapper, tempest_mesh ),
                    "Failed to convert and write MOAB mesh" );

    return moab::MB_SUCCESS;
}

moab::ErrorCode handleCSMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    // Generate the TempestRemap Exodus mesh (always written as .g for TempestRemap compatibility)
    std::string trFilename = ctx.outFilename;
    const size_t dot       = trFilename.find_last_of( "." );
    if( dot != std::string::npos && trFilename.substr( dot + 1 ) == "h5m" )
        trFilename = trFilename.substr( 0, dot ) + ".g";

    ctx.timer_push( "generate CS mesh with TempestRemap" );
    TR_CHK_SET_ERR( GenerateCSMesh( *tempest_mesh, ctx.blockSize, trFilename, "NetCDF4" ),
                    "Failed to generate CS mesh with TempestRemap" );
    ctx.timer_pop();

    // Add the CS mesh to the list of meshes
    ctx.meshes.push_back( tempest_mesh );

    // Convert and write as MOAB h5m if requested
    MB_CHK_SET_ERR( convertAndWriteMOABMesh( ctx, remapper, tempest_mesh ),
                    "Failed to convert and write MOAB mesh" );

    return moab::MB_SUCCESS;
}

}  // namespace

/**
 * @brief Creates a TempestRemap mesh based on the provided context and mesh type
 *
 * @param ctx Tool context containing configuration and state
 * @param remapper TempestRemap instance for mesh operations
 * @param tempest_mesh Output parameter for the created mesh
 * @return moab::ErrorCode Status of the operation
 */
static moab::ErrorCode CreateTempestMesh( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    using namespace moab;
    using RemapperType = moab::TempestRemapper;

    auto& outputFormatter = ctx.outputFormatter;

    try
    {
        switch( ctx.meshType )
        {
            case RemapperType::OVERLAP_FILES:
                if( !ctx.proc_id ) outputFormatter.printf( 0, "Creating TempestRemap overlap mesh ...\n" );
                return handleOverlapFiles( ctx, tempest_mesh );

            case RemapperType::OVERLAP_MEMORY:
                if( !ctx.proc_id )
                    outputFormatter.printf( 0, "Convert MOAB overlap files to TempestRemap format in-memory ...\n" );
                return handleOverlapMemory( ctx, remapper, tempest_mesh );

            case RemapperType::OVERLAP_MOAB:
                if( !ctx.proc_id )
                    outputFormatter.printf( 0, "Convert MOAB meshes to TempestRemap format in-memory ...\n" );
                return handleOverlapMOAB( ctx, remapper );

            case RemapperType::ICO:
                if( !ctx.proc_id ) outputFormatter.printf( 0, "Creating TempestRemap ICO mesh ...\n" );
                return handleICOMesh( ctx, remapper, tempest_mesh );

            case RemapperType::RLL:
                if( !ctx.proc_id ) outputFormatter.printf( 0, "Creating TempestRemap RLL mesh ...\n" );
                return handleRLLMesh( ctx, remapper, tempest_mesh );

            default:  // Default to CS mesh
                if( !ctx.proc_id ) outputFormatter.printf( 0, "Creating TempestRemap CS mesh ...\n" );
                return handleCSMesh( ctx, remapper, tempest_mesh );
        }
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error in CreateTempestMesh: " << e.what() << "\n";
        return MB_FAILURE;
    }
}

#undef MOAB_DBG

///////////////////////////////////////////////
/**
 * @brief Sample functions for testing remapping operations
 *
 * These functions provide analytical test cases with different spatial patterns
 * for verifying remapping accuracy and performance.
 */

// Constants for sample functions
namespace
{
// Constants for sample_stationary_vortex
constexpr double VORTEX_LON0 = 0.0;
constexpr double VORTEX_LAT0 = 0.6;
constexpr double VORTEX_R0   = 3.0;
constexpr double VORTEX_D    = 5.0;
constexpr double VORTEX_T    = 6.0;

}  // namespace

/**
 * @brief Constant sample function
 *
 * @return double Always returns 1.0
 */
static inline constexpr double sample_constant( double /*dLon*/, double /*dLat*/ ) noexcept
{
    return 1.0;
}

/**
 * @brief Sample function with slow harmonic variation
 *
 * @param dLon Longitude in radians
 * @param dLat Latitude in radians
 * @return double Function value at (dLon, dLat)
 */
static inline double sample_slow_harmonic( double dLon, double dLat ) noexcept
{
    const double cosLat = std::cos( dLat );
    return 2.0 + cosLat * cosLat * std::cos( 2.0 * dLon );
}

/**
 * @brief Sample function with fast harmonic variation
 *
 * @param dLon Longitude in radians
 * @param dLat Latitude in radians
 * @return double Function value at (dLon, dLat)
 */
static inline double sample_fast_harmonic( double dLon, double dLat ) noexcept
{
    const double sin2Lat = std::sin( 2.0 * dLat );
    return 2.0 +
           sin2Lat * sin2Lat * sin2Lat * sin2Lat * sin2Lat * sin2Lat * sin2Lat * sin2Lat * std::cos( 16.0 * dLon );
}

/**
 * @brief Sample function representing a stationary vortex
 *
 * @param dLon Longitude in radians
 * @param dLat Latitude in radians
 * @return double Function value at (dLon, dLat)
 */
static inline double sample_stationary_vortex( double dLon, double dLat ) noexcept
{
    // Find the rotated longitude and latitude of a point on a sphere
    // with pole at (dLonC, dLatC)
    const double dSinC = std::sin( VORTEX_LAT0 );
    const double dCosC = std::cos( VORTEX_LAT0 );
    const double dSinT = std::sin( dLat );
    const double dCosT = std::cos( dLat );

    const double dTrm = dCosT * std::cos( dLon - VORTEX_LON0 );
    const double dX   = dSinC * dTrm - dCosC * dSinT;
    const double dY   = dCosT * std::sin( dLon - VORTEX_LON0 );
    const double dZ   = dSinC * dSinT + dCosC * dTrm;

    // Calculate new longitude and latitude in rotated coordinate system
    double dNewLon = std::atan2( dY, dX );
    if( dNewLon < 0.0 )
    {
        dNewLon += 2.0 * M_PI;
    }
    const double dNewLat = std::asin( dZ );

    // Calculate vortex profile
    const double dRho = VORTEX_R0 * std::cos( dNewLat );
    const double dVt  = 3.0 * std::sqrt( 3.0 ) / 2.0 / std::cosh( dRho ) / std::cosh( dRho ) * std::tanh( dRho );

    // Calculate angular velocity (avoid division by zero)
    const double dOmega = ( dRho == 0.0 ) ? 0.0 : ( dVt / dRho );

    // Return the final vortex profile
    return ( 1.0 - std::tanh( dRho / VORTEX_D * std::sin( dNewLon - dOmega * VORTEX_T ) ) );
}

///////////////////////////////////////////////
