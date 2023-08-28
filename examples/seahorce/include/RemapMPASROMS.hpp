#ifndef __remap_mpas_roms_hpp__
#define __remap_mpas_roms_hpp__

#include "moab/MOABConfig.h"

#include "moab/Core.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"

// 3D settings
constexpr int mpas_zreflevels = 60;
constexpr int mpas_zlevels    = 60;
constexpr int roms_zlevels    = 100;
constexpr int nvars           = 2;

// tag name data
const char* mpas_twod_tagnames[nvars]       = { "salinity", "temperature" };
const char* mpas_threed_cum_tagnames[nvars] = { "salinity_3d", "temperature_3d" };
const char* mpas_threed_tagnames[nvars]     = { "Salinity3d", "Temperature3d" };
const char* roms_twod_tagnames[nvars]       = { "Salinity2DROMS", "Temperature2DROMS" };
const char* roms_threed_tagnames[nvars]     = { "Salinity3dROMS", "Temperature3dROMS" };

// write the map file to disk; comment out to just compute in-memory
#define VERTICAL_INTERPOLATION
// #define VERTICAL_INTERPOLANT_LINEAR
#define WRITE_MAP_FILE

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "GridElements.h"
#include "OfflineMap.h"
#endif

// Error check routines and utility macros
#define dbgprint( MSG )                                           \
    do                                                            \
    {                                                             \
        if( context.proc_id == 0 ) std::cout << MSG << std::endl; \
    } while( false )

#define dbgprintall( MSG )                                                \
    do                                                                    \
    {                                                                     \
        std::cout << "[" << context.proc_id << "]: " << MSG << std::endl; \
    } while( false )

#define println( MSG )                                    \
    do                                                    \
    {                                                     \
        if( proc_id == 0 ) std::cout << MSG << std::endl; \
    } while( false )

// get number of arguments with __NARG__
#define __NARG__( ... )  __NARG_I_( __VA_ARGS__, __RSEQ_N() )
#define __NARG_I_( ... ) __ARG_N( __VA_ARGS__ )
#define __ARG_N( _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, \
                 _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, _41, _42,  \
                 _43, _44, _45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62,  \
                 _63, N, ... )                                                                                        \
    N
#define __RSEQ_N()                                                                                                    \
    63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36,   \
        35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, \
        7, 6, 5, 4, 3, 2, 1, 0

// general definition for any function name
#define _VFUNC_( name, n ) name##n
#define _VFUNC( name, n )  _VFUNC_( name, n )
#define VFUNC( func, ... ) _VFUNC( func, __NARG__( __VA_ARGS__ ) )( __VA_ARGS__ )

// #define FOO( ... )         VFUNC( FOO, __VA_ARGS__ )

#define runchk( ... ) VFUNC( runchk, __VA_ARGS__ )

#define runchk1( CODE )     \
    do                      \
    {                       \
        MB_CHK_ERR( CODE ); \
    } while( false )

#define runchk2( CODE, MSG )        \
    do                              \
    {                               \
        moab::ErrorCode err = CODE; \
        MB_CHK_SET_ERR( err, MSG ); \
    } while( false )

#define runchk_cont( CODE, MSG )                               \
    do                                                         \
    {                                                          \
        moab::ErrorCode err = CODE;                            \
        if( err ) std::cout << "Error:: " << MSG << std::endl; \
        MB_CHK_ERR_CONT( err );                                \
    } while( false )

/// @brief class RuntimeContext
/// The RunttimeContext stores and manages the MPAS-ROMS coupler specific
/// to multiscale modeling. The class also provides other utility functions
/// to profile operations etc
struct RuntimeContext
{
#ifdef MOAB_HAVE_TEMPESTREMAP
    Mesh meshInput;
    Mesh meshOutput;
    Mesh meshOverlap;
    OfflineMap weightMap;
#endif

    std::string mpas_filename{ "mpas_grid.h5m" };  /// input MPAS grid file
    std::string roms_filename{ "roms_grid.h5m" };  /// input ROMS grid file
    int dimension{ 2 };                            /// dimensionality of the problem: 2 or 3
    bool ensureMonotonicity{ false };  /// flag indicating use of monotone approximations for projections (TempestRemap)
    bool computeTR{ false };           /// flag indicating use of TempestRemap conservative schemes for projections
    bool computeShepard{ false };      /// flag indicating use of Shepard approximations for projections
    bool computeMBA{ false };          /// flag indicating use of multilevel B-spline approximations for projections
    bool nearestNeighbor{ false };     /// flag indicating use of nearest neighbor method for projections
    bool normalize{ false };           /// normalize the dataset to original source data integral
                                       /// (more relevant for 2D; not implemented for 3D at the moment)
    bool use_3dprojection{ false };    /// perform 2D projections (if false), else compute 3D projections
    bool generateExtrusions{ false };  /// flag to indicate extruded 3D mesh computation
                                       ///(as opposed to loading from disk)
    bool threetwooneD{ false };        /// Perform 3D projections as 2Dx1D or 1Dx2D
                                       /// (tensor product computations) as opposed to full 3D
    bool oneDfirst{ false };           /// flag to specify that 1Dx2D should be performed as opposed
                                       /// to 2Dx1D (only relevant when threetwooneD is true)
    std::string strMethod{ "" };       /// remapping method used to other scalar fields (temperature, salinity)
    std::string bathymetryMethod{ "bilin" };  /// remapping method used for Bathymetry field
    int bathymetryOrder{ 1 };                 /// order of the reconstruction for Bathymetry field
    int fieldOrder{ 3 };         /// order of the reconstruction for other scalar fields (temperature, salinity)
    int proc_id{ 1 };            /// process identifier
    int num_procs{ 1 };          /// total number of processes
    double last_counter{ 0.0 };  /// last time counter between push/pop timer

    // MOAB objects
    moab::Interface* moab_interface{ nullptr };
    moab::ParallelComm* parallel_communicator{ nullptr };
    moab::EntityHandle partnset{ 0 };
    moab::EntityHandle mpasset, mpas_covering_set, romsset;

    double mpas_zref_heights[mpas_zreflevels];

    /// @brief Constructor: allocate MOAB interface and communicator, and initialize
    /// other data members with some default values
    RuntimeContext( MPI_Comm comm = MPI_COMM_WORLD )
    {
        moab::ErrorCode err;
        // Create the moab instance
        moab_interface = new( std::nothrow ) moab::Core;
        if( NULL == moab_interface ) exit( 1 );

        // Create sets for the mesh and partition.  Then pass these to the load_file functions to populate the mesh.
        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, partnset ), "Creating partition set failed" );

        // Create the parallel communicator object with the partition handle associated with MOAB
        parallel_communicator = moab::ParallelComm::get_pcomm( moab_interface, partnset, &comm );

        proc_id   = parallel_communicator->rank();
        num_procs = parallel_communicator->size();

        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, mpasset ), "Can't create new set" );
        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, mpas_covering_set ), "Can't create new set" );
        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, romsset ), "Can't create new set" );
    }

    /// @brief Destructor: deallocate MOAB interface and communicator
    ~RuntimeContext()
    {
        delete parallel_communicator;
        delete moab_interface;
    }

    /// @brief Parse the runtime command line options
    /// @param argc - number of command line arguments
    /// @param argv - command line arguments as string list
    void ParseCLOptions( int argc, char* argv[] )
    {
        ProgOptions opts;

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

        // only MBA is right now tested with 3D projections?
        if( use_3dprojection ) computeMBA = true;
    }

    /// @brief Method that prints out the runtime parameters in use for provenance
    void describe()
    {
        // Print out the input parameters
        if( proc_id == 0 )
        {
            println( "********** Remap MPAS-to-ROMS **********\n" );

            println( " -- Runtime Parameters -- " );
            println( "   MPAS mesh file: " << mpas_filename );
            println( "   ROMS mesh file: " << roms_filename );

            println( "        Dimension: " << dimension );
            println( "        Algorithm: " << ( nearestNeighbor ? "Nearest Neighbor mapping"
                                                : computeTR     ? "TempestRemap Conservative mapping"
                                                : computeMBA    ? "Multilevel B-spline Approximation"
                                                                : "Shepard interpolant" ) );
            if( computeTR ) println( "           Method: " << strMethod );
            println( " Bathymetry Order: " << bathymetryOrder );
            println( "      Field Order: " << fieldOrder );
            println( std::endl );
        }
    }

    /// @brief Measure and start the timer to profile a task
    /// @param operation String name of the task being measured
    inline void timer_push( std::string operation )
    {
        mTimerOps = mTimer.time_since_birth();
        mOpName   = operation;
    }

    /// @brief Stop the timer and store the elapsed duration
    /// @param nruns Optional argument used to average the measured time
    void timer_pop( const int nruns = 1 )
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

    /// @brief Return the last elapsed time
    /// @return last_counter from timer_pop was called
    inline double last_elapsed() const
    {
        return last_counter;
    }

  private:
    moab::CpuTimer mTimer;
    double mTimerOps{ 0.0 };
    std::string mOpName;
};

#endif  // __remap_mpas_roms_hpp__