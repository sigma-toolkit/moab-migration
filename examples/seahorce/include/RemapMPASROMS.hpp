#ifndef __remap_mpas_roms_hpp__
#define __remap_mpas_roms_hpp__

// MOAB config include
#include "moab/MOABConfig.h"
#include "ExampleConfig.hpp"
#include "ExampleErrorHandler.hpp"

// MOAB includes
#include "moab/Core.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"

#include <Eigen/Dense>

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "GridElements.h"
#include "OfflineMap.h"
#endif

enum RemappingMethod
{
    DefaultRemappingMethod         = -1,
    NearestNeighborInterpolant     = 0,
    TempestRemapFV                 = 1,
    TempestRemapBilinear           = 2,
    TempestRemapInvDist            = 3,
    TempestRemapDelaunay           = 4,
    TempestRemapIntegratedBilinear = 5,
    DelaunayInterpolant            = 6,
    ShepardInterpolant             = 7,
    TrilinearTensor2D1D            = 8,
    TrilinearTensor1D2D            = 9,
    MultilevelBsplineApproximation = 10
};

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
    bool computeTRMaps{ false };       /// flag indicating use of TempestRemap conservative schemes for projections
    // bool computeShepardInterpolant{ false };      /// flag indicating use of Shepard approximations for projections
    bool computeMBAInterpolant{ false };  /// flag indicating use of multilevel B-spline approximations for projections
    // bool computeNNInterpolant{ false };  /// flag indicating use of nearest neighbor method for projections
    // bool computeDelaunayInterpolant{ false };  ///flag indicating use of delaunay natural neighbor interpolant
    bool normalize{ false };           /// normalize the dataset to original source data integral
                                       /// (more relevant for 2D; not implemented for 3D at the moment)
    bool use_3dprojection{ false };    /// perform 2D projections (if false), else compute 3D projections
    bool generateExtrusions{ false };  /// flag to indicate extruded 3D mesh computation
                                       ///(as opposed to loading from disk)
    bool threetwooneD{ false };        /// Perform 3D projections as 2Dx1D or 1Dx2D
                                       /// (tensor product computations) as opposed to full 3D
    bool oneDfirst{ false };           /// flag to specify that 1Dx2D should be performed as opposed
                                       /// to 2Dx1D (only relevant when threetwooneD is true)
    // std::string strMethod{ "" };       /// remapping method used to other scalar fields (temperature, salinity)
    // std::string bathymetryMethod{ "" };  /// remapping method used for Bathymetry field
    // int bathymetryOrder{ 3 };                 /// order of the reconstruction for Bathymetry field
    // int fieldOrder{ 3 };         /// order of the reconstruction for other scalar fields (temperature, salinity)
    int proc_id{ 1 };            /// process identifier
    int num_procs{ 1 };          /// total number of processes
    double last_counter{ 0.0 };  /// last time counter between push/pop timer
    bool useCAAS{ false };

    std::map< std::string, std::pair< RemappingMethod, int > > field_methods;

    // MOAB objects
    moab::Interface* moab_interface{ nullptr };
    moab::ParallelComm* parallel_communicator{ nullptr };
    moab::EntityHandle partnset{ 0 };
    moab::EntityHandle mpasset, mpas_covering_set, romsset;

    // Eigen::Matrix< int, Eigen::Dynamic, 4 > dual_mpas_tetrahedron;
    std::vector< int > nvertcache;
    std::vector< int > tetrahedraconn;
    Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > vertex_to_element;
    std::vector< double > dual_mpas_tetrahedron_centroids;
    // moab::Range mpas_elems, mpas_verts;
    // moab::Range roms_elems, roms_verts;
    // moab::Range mpas3d_elems, mpas3d_verts, mpas3d_dual_elems;
    // moab::Range roms3d_elems, roms3d_verts;

    double mpas_zref_heights[mpas_zreflevels];

    /// @brief Constructor: allocate MOAB interface and communicator, and initialize
    /// other data members with some default values
    RuntimeContext( MPI_Comm comm = MPI_COMM_WORLD )
    {
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

    static RemappingMethod GetMethod( std::string methodName )
    {
        /*
            DefaultRemappingMethod         = -1,
            NearestNeighborInterpolant     = 0,
            TempestRemapFV                 = 1,
            TempestRemapBilinear           = 2,
            TempestRemapInvDist            = 3,
            DelaunayInterpolant            = 4,
            ShepardInterpolant             = 5,
            TrilinearTensor2D1D            = 6,
            TrilinearTensor1D2D            = 7,
            MultilevelBsplineApproximation = 8
        */
        if( !methodName.compare( "NearestNeighborInterpolant" ) )
            return NearestNeighborInterpolant;
        else if( !methodName.compare( "TempestRemapFV" ) )
            return TempestRemapFV;
        else if( !methodName.compare( "TempestRemapBilinear" ) )
            return TempestRemapBilinear;
        else if( !methodName.compare( "TempestRemapInvDist" ) )
            return TempestRemapInvDist;
        else if( !methodName.compare( "TempestRemapDelaunay" ) )
            return TempestRemapDelaunay;
        else if( !methodName.compare( "TempestRemapIntegratedBilinear" ) )
            return TempestRemapIntegratedBilinear;
        else if( !methodName.compare( "DelaunayInterpolant" ) )
            return DelaunayInterpolant;
        else if( !methodName.compare( "ShepardInterpolant" ) )
            return ShepardInterpolant;
        else if( !methodName.compare( "TrilinearTensor2D1D" ) )
            return TrilinearTensor2D1D;
        else if( !methodName.compare( "TrilinearTensor2D1D" ) )
            return TrilinearTensor2D1D;
        else if( !methodName.compare( "MultilevelBsplineApproximation" ) )
            return MultilevelBsplineApproximation;
        else
            return DefaultRemappingMethod;
    }

    static std::string GetMethod( RemappingMethod methodEnum )
    {
        /*
            DefaultRemappingMethod         = -1,
            NearestNeighborInterpolant     = 0,
            TempestRemapFV                 = 1,
            TempestRemapBilinear           = 2,
            TempestRemapInvDist            = 3,
            DelaunayInterpolant            = 4,
            ShepardInterpolant             = 5,
            TrilinearTensor2D1D            = 6,
            TrilinearTensor1D2D            = 7,
            MultilevelBsplineApproximation = 8
        */
        if( methodEnum == NearestNeighborInterpolant )
            return "NearestNeighborInterpolant";
        else if( methodEnum == TempestRemapFV )
            return "TempestRemapFV";
        else if( methodEnum == TempestRemapBilinear )
            return "TempestRemapBilinear";
        else if( methodEnum == TempestRemapInvDist )
            return "TempestRemapInvDist";
        else if( methodEnum == TempestRemapDelaunay )
            return "TempestRemapDelaunay";
        else if( methodEnum == TempestRemapIntegratedBilinear )
            return "TempestRemapIntegratedBilinear";
        else if( methodEnum == DelaunayInterpolant )
            return "DelaunayInterpolant";
        else if( methodEnum == ShepardInterpolant )
            return "ShepardInterpolant";
        else if( methodEnum == TrilinearTensor2D1D )
            return "TrilinearTensor2D1D";
        else if( methodEnum == TrilinearTensor2D1D )
            return "TrilinearTensor2D1D";
        else if( methodEnum == MultilevelBsplineApproximation )
            return "MultilevelBsplineApproximation";
        else
            return "DefaultRemappingMethod";
    }

    /// @brief Function to split a given option as: field:method:order into appropriate pieces
    /// @param str Original string containing all options
    /// @param strings Output list of split strings
    /// @param separator Separator character (default: ':')
    void split_option( std::string str, std::vector< std::string >& strings, char separator = ':' )
    {
        strings.clear();
        int startIndex = 0, endIndex = 0;
        for( size_t i = 0; i <= str.size(); i++ )
        {
            // If we reached the end of the word or the end of the input.
            if( str[i] == separator || i == str.size() )
            {
                endIndex = i;
                std::string temp;
                temp.append( str, startIndex, endIndex - startIndex );
                strings.push_back( temp );
                startIndex = endIndex + 1;
            }
        }
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
        opts.addOpt< void >( "mono", "Ensure monotonicity in the weight generation (only for TR-FV)",
                             &ensureMonotonicity );
        opts.addOpt< void >( "caas", "Use CAAS limiter", &useCAAS );
        opts.addOpt< void >( "321D", "Compute three-dimensional projections using a 2Dx1D approach", &threetwooneD );
        opts.addOpt< void >( "1D2D", "Use 1Dx2D as opposed to 2Dx1D for 321D projection", &oneDfirst );
        opts.addOpt< void >(
            "normalize",
            "Re-normalize interpolant to preserve global field integral (only 2D and requires mesh intersection)",
            &normalize );

        std::string fmethodorder = "";
        opts.addOpt< std::string >( "methodorder",
                                    "Format: field:method, field: { Bathymetry, Temperature, Salinity }, "
                                    "Method: { NearestNeighborInterpolant, "
                                    "TempestRemapFV, "
                                    "TempestRemapBilinear, "
                                    "TempestRemapInvDist, "
                                    "TempestRemapDelaunay, "
                                    "TempestRemapIntegratedBilinear, "
                                    "DelaunayInterpolant, "
                                    "ShepardInterpolant, "
                                    "TrilinearTensor2D1D, "
                                    "TrilinearTensor1D2D, "
                                    "MultilevelBsplineApproximation }",
                                    &fmethodorder );

        opts.parseCommandLine( argc, argv );

        {
            field_methods["Bathymetry"]  = std::make_pair< RemappingMethod, int >( TempestRemapBilinear, 1 );
            field_methods["Salinity"]    = std::make_pair< RemappingMethod, int >( TempestRemapBilinear, 1 );
            field_methods["Temperature"] = std::make_pair< RemappingMethod, int >( TempestRemapBilinear, 1 );

            std::vector< std::string > fieldmethods;
            opts.getOptAllArgs( "methodorder", fieldmethods );

            for( auto fmethod : fieldmethods )
            {
                std::vector< std::string > optionStorage;
                split_option( fmethod, optionStorage );
                assert( optionStorage.size() > 1 );
                std::string fieldname = optionStorage[0];
                if( ( !fieldname.compare( "Bathymetry" ) || !fieldname.compare( "Salinity" ) ||
                      !fieldname.compare( "Temperature" ) ) &&
                    optionStorage.size() > 1 )
                {
                    std::string tmpmethod   = optionStorage[1];
                    RemappingMethod rmethod = GetMethod( tmpmethod );
                    int methodorder         = 1;
                    if( optionStorage.size() > 2 )  // order of the method
                        methodorder = atoi( optionStorage[2].c_str() );

                    field_methods[fieldname] = std::make_pair( rmethod, methodorder );
                    if( rmethod == TempestRemapFV || rmethod == TempestRemapBilinear ||
                        rmethod == TempestRemapInvDist || rmethod == TempestRemapDelaunay ||
                        rmethod == TempestRemapIntegratedBilinear )
                        computeTRMaps = true;
                    if( rmethod == MultilevelBsplineApproximation ) computeMBAInterpolant = true;
                }
                else
                {
                    std::cout << "Error: Ignoring specification for non-standard field: " << fieldname << std::endl;
                }
            }
        }

        if( dimension == 3 ) use_3dprojection = true;
    }

    /// @brief Method that prints out the runtime parameters in use for provenance
    void describe()
    {
        // Print out the input parameters
        if( proc_id == 0 )
        {
            std::stringstream sstr;
            sstr << "********** Remap MPAS-to-ROMS **********" << std::endl << std::endl;

            sstr << " -- Runtime Parameters -- " << std::endl;
            sstr << "   MPAS mesh file: " << mpas_filename << std::endl;
            sstr << "   ROMS mesh file: " << roms_filename << std::endl;

            sstr << "        Dimension: " << dimension << std::endl;
            // sstr << "        Algorithm: "
            //  << ( computeNNInterpolant    ? "Nearest Neighbor mapping"
            //           : computeTRProjection   ? "TempestRemap Conservative mapping"
            //   : computeMBAInterpolant ? "Multilevel B-spline Approximation"
            //                                   : "Shepard interpolant" )
            //      << std::endl;
            // if( computeTRProjection ) sstr << "           Method: " << strMethod << std::endl;
            sstr << "   Field: Bathymetry, Method: " << GetMethod( field_methods["Bathymetry"].first )
                 << ", Order: " << field_methods["Bathymetry"].second << std::endl;
            sstr << "   Field: Salinity, Method: " << GetMethod( field_methods["Salinity"].first )
                 << ", Order: " << field_methods["Salinity"].second << std::endl;
            sstr << "   Field: Temperature, Method: " << GetMethod( field_methods["Temperature"].first )
                 << ", Order: " << field_methods["Temperature"].second << std::endl;
            sstr << std::endl;

            std::cout << sstr.str() << std::endl;
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
    void timer_pop()
    {
        double locElapsed = mTimer.time_since_birth() - mTimerOps;
        if( proc_id == 0 )
        {
            std::cout << "[LOG] Time taken to " << mOpName.c_str() << " := " << locElapsed << std::endl;
            last_counter = locElapsed;
        }
        mOpName.clear();
    }

    /// @brief Return the last elapsed time
    /// @return last_counter from timer_pop was called
    inline double last_elapsed() const
    {
        return last_counter;
    }

    /// @brief Computes the field projection for the variable according to method parameters requested by user
    /// @param varProjectSrc Name of the source variable tag
    /// @param varProjectDst Name of the target variable tag after projection
    /// @param srcelems Source element list
    /// @param dstelems Target element list
    /// @param constantoffset If internal normalization is to be performed
    /// @param src3delems If 3D, source element list of extruded cells
    /// @param dst3delems If 3D, target element list of extruded cells
    /// @return Error code
    moab::ErrorCode ComputeFieldProjections( std::string varProjectSrc,
                                             std::string varProjectDst,
                                             std::vector< moab::EntityHandle >& srcelems,
                                             std::vector< moab::EntityHandle >& dstelems,
                                             const double constantoffset                   = 0.0,
                                             std::vector< moab::EntityHandle >* src3delems = nullptr,
                                             std::vector< moab::EntityHandle >* dst3delems = nullptr );

  private:
    moab::CpuTimer mTimer;
    double mTimerOps{ 0.0 };
    std::string mOpName;
};

// Forward declarations

#endif  // __remap_mpas_roms_hpp__