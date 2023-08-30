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

    Eigen::Matrix< int, Eigen::Dynamic, 4 > dual_mpas_tetrahedron;
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
          std::stringstream sstr;
          sstr << "********** Remap MPAS-to-ROMS **********" << std::endl << std::endl;

          sstr << " -- Runtime Parameters -- " << std::endl;
          sstr << "   MPAS mesh file: " << mpas_filename << std::endl;
          sstr << "   ROMS mesh file: " << roms_filename << std::endl;

          sstr << "        Dimension: " << dimension << std::endl;
          sstr << "        Algorithm: "
               << ( nearestNeighbor ? "Nearest Neighbor mapping"
                    : computeTR     ? "TempestRemap Conservative mapping"
                    : computeMBA    ? "Multilevel B-spline Approximation"
                                    : "Shepard interpolant" )
               << std::endl;
          if( computeTR ) sstr << "           Method: " << strMethod << std::endl;
          sstr << " Bathymetry Order: " << bathymetryOrder << std::endl;
          sstr << "      Field Order: " << fieldOrder << std::endl;
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

  private:
    moab::CpuTimer mTimer;
    double mTimerOps{ 0.0 };
    std::string mOpName;
};

// Forward declarations

moab::ErrorCode ComputeFieldProjections( moab::Interface* mbi,
                                         RuntimeContext& context,
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

#endif  // __remap_mpas_roms_hpp__