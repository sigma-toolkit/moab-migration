/**
 * @file iMOAB.cpp
 * @brief Implementation of iMOAB C API for MOAB mesh database operations
 *
 * This file implements the iMOAB interface - a C API for MOAB (Mesh Oriented datABase)
 * that provides simplified access to MOAB functionality for coupled applications,
 * particularly in climate and multiphysics simulations.
 */

/**
 * @defgroup iMOAB iMOAB C Interface
 * @brief C API for MOAB mesh database operations in coupled applications
 *
 * The iMOAB interface provides a simplified C API for accessing MOAB functionality,
 * designed for use in coupled climate and multiphysics applications. It handles
 * mesh I/O, parallel communication, remapping, and data transfer between components.
 *
 * @{
 */

/**
 * @defgroup iMOABInit Initialization and Finalization
 * @brief Functions for initializing and finalizing the iMOAB library
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABApp Application Management
 * @brief Functions for registering and managing application instances
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABIO Mesh I/O Operations
 * @brief Functions for reading and writing mesh files
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABQuery Mesh Query Operations
 * @brief Functions for querying mesh information (vertices, elements, blocks, BCs)
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABTag Tag Operations
 * @brief Functions for defining and accessing tag data on mesh entities
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABParallel Parallel Communication
 * @brief Functions for parallel mesh operations and data transfer
 * @ingroup iMOAB
 */

/**
 * @defgroup iMOABRemap Remapping and Intersection
 * @brief Functions for mesh intersection and field remapping (TempestRemap integration)
 * @ingroup iMOAB
 */

/** @} */  // end of iMOAB group

#include "moab/MOABConfig.h"
#include "moab/Core.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "moab/ParCommGraph.hpp"
#include "moab/ParallelMergeMesh.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#endif
#include "DebugOutput.hpp"
#include "moab/iMOAB.h"

/* this is needed because of direct access to hdf5/mhdf */
#ifdef MOAB_HAVE_HDF5
#include "mhdf.h"
#include <H5Tpublic.h>
#endif

#include "moab/CartVect.hpp"
#include "MBTagConventions.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/MergeMesh.hpp"

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "STLStringHelper.h"
#include "moab/IntxMesh/IntxUtils.hpp"

#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"
#endif

// C++ includes
#include <cassert>
#include <sstream>
#include <iostream>

using namespace moab;

// #define VERBOSE

// global variables ; should they be organized in a structure, for easier references?
// or how do we keep them global?

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
struct TempestMapAppData
{
    moab::TempestRemapper* remapper;
    std::map< std::string, moab::TempestOnlineMap* > weightMaps;
    iMOAB_AppID pid_src;
    iMOAB_AppID pid_dest;
    int num_src_ghost_layers;  // number of ghost layers
    int num_tgt_ghost_layers;  // number of ghost layers
};
#endif

struct appData
{
    EntityHandle file_set;     // primary file set containing all entities of interest
    int global_id;             // external component id, unique for application
    std::string name;          // name of the application
    Range all_verts;           // local vertices would be all_verts if no ghosting was required
    Range local_verts;         // it could include shared, but not owned at the interface
    Range owned_verts;         // owned_verts <= local_verts <= all_verts
    Range ghost_vertices;      // locally ghosted from other processors
    Range primary_elems;       // all primary entities (owned + ghosted)
    Range owned_elems;         // only owned entities
    Range ghost_elems;         // only ghosted entities (filtered)
    int dimension;             // 2 or 3, dimension of primary elements (redundant?)
    long num_global_elements;  // reunion of all elements in primary_elements; either from hdf5
                               // reading or from reduce
    long num_global_vertices;  // reunion of all nodes, after sharing is resolved; it could be
                               // determined from hdf5 reading
    int num_ghost_layers;      // number of ghost layers
    Range mat_sets;
    std::map< int, int > matIndex;  // map from global block id to index in mat_sets
    Range neu_sets;
    Range diri_sets;
    std::map< std::string, Tag > tagMap;
    std::vector< Tag > tagList;
    bool point_cloud;
    bool is_fortran;

#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm;
    // constructor for this ParCommGraph takes the joint comm and the MPI groups for each
    // application
    std::map< int, ParCommGraph* > pgraph;  // map from context () to the parcommgraph*
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
    EntityHandle secondary_file_set;  // secondary file set (typically a child set like covering mesh)
                                      // so we assume only one covering set for all maps on this intx app
                                      // we can have multiple weightMaps, but only one coverage set for all maps
    TempestMapAppData tempestData;
    std::map< std::string, std::string > metadataMap;
#endif
};

struct GlobalContext
{
    // are there reasons to have multiple moab inits? Is ref count needed?
    Interface* MBI;
    // we should also have the default tags stored, initialized
    Tag material_tag, neumann_tag, dirichlet_tag,
        globalID_tag;  // material, neumann, dirichlet,  globalID
    int refCountMB;
    int iArgc;
    iMOAB_String* iArgv;

    std::map< std::string, int > appIdMap;  // from app string (uppercase) to app id
    std::map< int, appData > appDatas;      // the same order as pcomms
    int globalrank, worldprocs;
    bool MPI_initialized;

    GlobalContext()
    {
        MBI        = 0;
        refCountMB = 0;
    }
};

static struct GlobalContext context;

/**
 * @brief Initialize iMOAB library and create MOAB instance.
 * @ingroup iMOABInit
 *
 * @details Initializes the iMOAB interface by creating a MOAB Core instance and setting up
 * standard tags for material sets, boundary conditions, and global IDs. Uses reference counting
 * to support multiple initialization calls. Detects MPI initialization status if built with MPI support.
 *
 * @par Initialization Steps:
 * -# Store command-line arguments (argc/argv) for later use
 * -# Create MOAB Core instance on first initialization (refCountMB == 0)
 * -# Retrieve standard tag handles: MATERIAL_SET, NEUMANN_SET, DIRICHLET_SET, GLOBAL_ID
 * -# Store tag handles in global context for efficient access
 * -# Detect MPI initialization and store world size/rank if MPI is active
 * -# Increment reference count to support nested initialization
 *
 * @param[in] argc  Number of command-line arguments (can be 0)
 * @param[in] argv  Array of command-line argument strings (can be NULL if argc is 0)
 *
 * @pre MPI must be initialized before calling if using parallel features
 * @post MOAB Core instance created and ready for use
 * @post Standard tags available via context.material_tag, context.neumann_tag, etc.
 *
 * @note Uses reference counting - safe to call multiple times
 * @note Fortran interface should use iMOAB_InitializeFortran() instead
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_Finalize()
 */
ErrCode iMOAB_Initialize( int argc, iMOAB_String* argv )
{
    if( argc ) IMOAB_CHECKPOINTER( argv, 1 );

    // Store command-line arguments for potential later use
    context.iArgc = argc;
    context.iArgv = argv;  // Shallow copy - caller must maintain argv lifetime

    // Create MOAB instance only on first initialization
    if( 0 == context.refCountMB )
    {
        context.MBI = new( std::nothrow ) moab::Core;

        // Retrieve standard MOAB tags that are commonly used across applications
        const char* const shared_set_tag_names[] = { MATERIAL_SET_TAG_NAME, NEUMANN_SET_TAG_NAME,
                                                     DIRICHLET_SET_TAG_NAME, GLOBAL_ID_TAG_NAME };
        // Tag purposes: materials (blocks), surface-BC (Neumann), vertex-BC (Dirichlet), global-id
        Tag gtags[4];
        for( int i = 0; i < 4; i++ )
        {
            MB_CHK_ERR(
                context.MBI->tag_get_handle( shared_set_tag_names[i], 1, MB_TYPE_INTEGER, gtags[i], MB_TAG_ANY ) );
        }

        // Cache standard tag handles in global context for efficient repeated access
        context.material_tag  = gtags[0];  // Material/block identification
        context.neumann_tag   = gtags[1];  // Surface boundary conditions
        context.dirichlet_tag = gtags[2];  // Vertex boundary conditions
        context.globalID_tag  = gtags[3];  // Global entity identification

        // Check if MPI is initialized and cache world communicator info
        context.MPI_initialized = false;
#ifdef MOAB_HAVE_MPI
        int flagInit;
        MPI_Initialized( &flagInit );

        if( flagInit && !context.MPI_initialized )
        {
            MPI_Comm_size( MPI_COMM_WORLD, &context.worldprocs );
            MPI_Comm_rank( MPI_COMM_WORLD, &context.globalrank );
            context.MPI_initialized = true;
        }
#endif
    }

    // Increment reference count to track active users
    context.refCountMB++;
    return moab::MB_SUCCESS;
}

/**
 * @brief Fortran-compatible wrapper for iMOAB_Initialize.
 * @ingroup iMOABInit
 *
 * @details Calls iMOAB_Initialize with zero arguments, suitable for Fortran applications
 * that don't pass command-line arguments through the interface.
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_Initialize()
 */
ErrCode iMOAB_InitializeFortran()
{
    return iMOAB_Initialize( 0, 0 );
}

/**
 * @brief Finalize iMOAB library and cleanup MOAB instance.
 * @ingroup iMOABInit
 *
 * @details Decrements reference count and deletes MOAB Core instance when count reaches zero.
 * Safe to call multiple times - must be called once for each Initialize call.
 *
 * @post Reference count decremented
 * @post MOAB instance deleted when refCountMB reaches 0
 * @post All application data should be cleaned up before final Finalize
 *
 * @warning Must call iMOAB_DeregisterApplication for all applications before final Finalize
 * @note Uses reference counting to support nested Initialize/Finalize pairs
 *
 * @return moab::MB_SUCCESS on success
 *
 * @see iMOAB_Initialize()
 */
ErrCode iMOAB_Finalize()
{
    // Decrement reference count
    context.refCountMB--;

    // Delete MOAB instance only when last user finalizes
    if( 0 == context.refCountMB )
    {
        delete context.MBI;
    }

    return MB_SUCCESS;
}

/**
 * @brief Fast string-to-integer hash function for generating application IDs.
 *
 * @details Uses Fowler-Noll-Vo (FNV) hash algorithm to generate unique integer IDs from
 * application names combined with component IDs. Non-cryptographic but fast and robust.
 *
 * @param[in] str        Application name string to hash
 * @param[in] identifier Optional component ID to incorporate into hash (default 0)
 *
 * @return Positive integer hash value (masked to ensure non-negative)
 *
 * @see https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
 */
static int apphash( const iMOAB_String str, int identifier = 0 )
{
    std::string appstr( str );
    // FNV-1a: Fast non-cryptographic hash by Fowler, Noll, and Vo
    unsigned int h = 2166136261u;  // FNV offset basis

    // Hash the application name string
    for( char c : appstr )
    {
        h ^= static_cast< unsigned char >( c );  // XOR with byte
        h *= 16777619u;                          // Multiply by FNV prime
    }

    // Incorporate component identifier into the hash for uniqueness
    if( identifier )
    {
        h ^= static_cast< unsigned int >( identifier & 0xFFFFFFFF );
        h *= 16777619u;  // FNV prime again
    }
    return static_cast< int >( h & 0x7FFFFFFF );  // Mask to ensure positive value
}

/**
 * @brief Register a new application instance with iMOAB.
 * @ingroup iMOABApp
 *
 * @details Creates a new application context with unique ID, allocates a mesh set for data storage,
 * and sets up parallel communicator if MPI is enabled. Each application represents a distinct
 * component in a coupled simulation (e.g., atmosphere, ocean, land).
 *
 * @par Registration Process:
 * -# Validate application name is unique (not already registered)
 * -# Generate unique application ID using FNV hash of name + component ID
 * -# Create MOAB mesh set for storing application's mesh data
 * -# Initialize application data structure (appData) with default values
 * -# Create ParallelComm instance for parallel operations if MPI enabled
 * -# Store application context in global map indexed by application ID
 *
 * @param[in]  app_name  Unique name for this application instance
 * @param[in]  comm      MPI communicator for this application (MPI builds only)
 * @param[in]  compid    External component ID (must be positive)
 * @param[out] pid       Generated application ID (output parameter)
 *
 * @pre iMOAB_Initialize must have been called
 * @pre app_name must be unique (not already registered)
 * @pre compid must be positive
 *
 * @post Application registered and ready for mesh loading
 * @post Application ID stored in pid
 * @post Mesh set created and stored in context.appDatas[pid]
 * @post ParallelComm created if MPI enabled
 *
 * @note Application ID is computed as hash(app_name, compid) for uniqueness
 * @note Fortran interface should use iMOAB_RegisterApplicationFortran()
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if already registered or compid invalid
 *
 * @see iMOAB_DeregisterApplication()
 */
ErrCode iMOAB_RegisterApplication( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
                                   MPI_Comm* comm,
#endif
                                   int* compid,
                                   iMOAB_AppID pid )
{
    IMOAB_CHECKPOINTER( app_name, 1 );
#ifdef MOAB_HAVE_MPI
    IMOAB_CHECKPOINTER( comm, 2 );
    IMOAB_CHECKPOINTER( compid, 3 );
#else
    IMOAB_CHECKPOINTER( compid, 2 );
#endif

    // will create a parallel comm for this application too, so there will be a
    // mapping from *pid to file set and to parallel comm instances
    std::string name( app_name );

    if( context.appIdMap.find( name ) != context.appIdMap.end() )
    {
        std::cout << " application " << name << " already registered \n";
        return moab::MB_FAILURE;
    }

    // trivial hash: just use the id that the user provided and assume it is unique
    // *pid = *compid;
    // hash: compute a unique hash as a combination of the appname and the component id
    *pid = apphash( app_name, *compid );
    // store map of application name and the ID we just generated
    context.appIdMap[name] = *pid;
    int rankHere           = 0;
#ifdef MOAB_HAVE_MPI
    MPI_Comm_rank( *comm, &rankHere );
#endif
    if( !rankHere )
        std::cout << " application " << name << " with ID = " << *pid << " and external id: " << *compid
                  << "  is registered now \n";
    if( *compid <= 0 )
    {
        std::cout << " convention for external application is to have its id positive \n";
        return moab::MB_FAILURE;
    }

    // create now the file set that will be used for loading the model in
    EntityHandle file_set;
    MB_CHK_SET_ERR( context.MBI->create_meshset( MESHSET_SET, file_set ), "can't create file set" );

    appData app_data;
    app_data.file_set  = file_set;
    app_data.global_id = *compid;  // will be used mostly for par comm graph
    app_data.name      = name;     // save the name of application

#ifdef MOAB_HAVE_TEMPESTREMAP
    app_data.tempestData.remapper             = nullptr;  // Only allocate as needed
    app_data.tempestData.num_src_ghost_layers = 0;
    app_data.tempestData.num_tgt_ghost_layers = 0;
#endif

    // set some default values
    app_data.num_ghost_layers = 0;
    app_data.point_cloud      = false;
    app_data.is_fortran       = false;
#ifdef MOAB_HAVE_TEMPESTREMAP
    app_data.secondary_file_set = app_data.file_set;
#endif

#ifdef MOAB_HAVE_MPI
    if( *comm ) app_data.pcomm = new ParallelComm( context.MBI, *comm );
#endif
    context.appDatas[*pid] = app_data;  // Store application data indexed by generated ID
    return moab::MB_SUCCESS;
}

/**
 * @brief Fortran-compatible wrapper for iMOAB_RegisterApplication.
 * @ingroup iMOABApp
 *
 * @details Registers application from Fortran code by converting Fortran MPI communicator to C
 * and calling C-style registration. Sets is_fortran flag to enable proper MPI handle conversions.
 *
 * @param[in]  app_name  Unique name for this application instance
 * @param[in]  comm      Fortran MPI communicator (integer handle, MPI builds only)
 * @param[in]  compid    External component ID (must be positive)
 * @param[out] pid       Generated application ID (output parameter)
 *
 * @note Automatically converts Fortran MPI_Comm to C MPI_Comm using MPI_Comm_f2c
 * @note Sets is_fortran flag in application data for future MPI handle conversions
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE otherwise
 *
 * @see iMOAB_RegisterApplication(), MPI_Comm_f2c()
 */
ErrCode iMOAB_RegisterApplicationFortran( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
                                          int* comm,
#endif
                                          int* compid,
                                          iMOAB_AppID pid )
{
    IMOAB_CHECKPOINTER( app_name, 1 );
#ifdef MOAB_HAVE_MPI
    IMOAB_CHECKPOINTER( comm, 2 );
    IMOAB_CHECKPOINTER( compid, 3 );
#else
    IMOAB_CHECKPOINTER( compid, 2 );
#endif

    ErrCode err;
    assert( app_name != nullptr );
    std::string name( app_name );

#ifdef MOAB_HAVE_MPI
    MPI_Comm ccomm;
    if( comm )
    {
        // Convert from Fortran communicator (integer) to C communicator
        // See MPI standard: http://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node361.htm
        ccomm = MPI_Comm_f2c( (MPI_Fint)*comm );
    }
#endif

    // Call C-style registration function with converted communicator
    err = iMOAB_RegisterApplication( app_name,
#ifdef MOAB_HAVE_MPI
                                     &ccomm,
#endif
                                     compid, pid );

    // Mark application as Fortran-based for proper MPI handle conversions in future calls
    context.appDatas[*pid].is_fortran = true;

    return err;
}

/**
 * @brief Deregister an application instance and cleanup associated resources.
 * @ingroup iMOABApp
 *
 * @details Removes application from iMOAB, deletes mesh entities, frees parallel communicator,
 * and cleans up all associated data structures. Must be called before iMOAB_Finalize().
 *
 * @par Cleanup Process:
 * -# Validate application ID exists
 * -# Delete all mesh entities in application's file set
 * -# Delete communication graphs (ParCommGraph instances)
 * -# Delete parallel communicator if MPI enabled
 * -# Delete TempestRemap objects if present
 * -# Delete application's mesh set
 * -# Remove application from global context maps
 *
 * @param[in] pid  Application ID to deregister
 *
 * @pre Application must have been registered via iMOAB_RegisterApplication
 * @post All mesh data deleted
 * @post Parallel communicator deleted
 * @post Application ID invalid and cannot be reused
 *
 * @warning Do not access application ID after deregistration
 * @note Safe to deregister in any order relative to other applications
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if application not found
 *
 * @see iMOAB_RegisterApplication()
 */
ErrCode iMOAB_DeregisterApplication( iMOAB_AppID pid )
{
    // Look up application in global context
    auto appIterator = context.appDatas.find( *pid );

    // Validate application exists before attempting deletion
    if( appIterator == context.appDatas.end() ) return MB_FAILURE;

    // we found the application
    appData& data = appIterator->second;
    int rankHere  = 0;
#ifdef MOAB_HAVE_MPI
    rankHere = data.pcomm->rank();
#endif
    if( !rankHere )
        std::cout << " application with ID: " << *pid << " global id: " << data.global_id << " name: " << data.name
                  << " is de-registered now \n";

    EntityHandle fileSet = data.file_set;
    // get all entities part of the file set
    Range fileents;
    MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( fileSet, fileents, /*recursive */ true ),
                    "can't get file entities" );
    fileents.insert( fileSet );
    MB_CHK_SET_ERR( context.MBI->get_entities_by_type( fileSet, MBENTITYSET, fileents ),
                    "can't get file entities" );  // append all mesh sets

#ifdef MOAB_HAVE_TEMPESTREMAP
    if( data.tempestData.remapper ) delete data.tempestData.remapper;
    if( data.tempestData.weightMaps.size() ) data.tempestData.weightMaps.clear();
#endif

#ifdef MOAB_HAVE_MPI
    // NOTE: we could get the pco also with the following workflow.
    // ParallelComm * pcomm = ParallelComm::get_pcomm(context.MBI, *pid);

    auto& pargs = data.pgraph;
    // free the parallel comm graphs associated with this app
    for( auto mt = pargs.begin(); mt != pargs.end(); ++mt )
    {
        ParCommGraph* pgr = mt->second;
        if( pgr != nullptr )
        {
            delete pgr;
            pgr = nullptr;
        }
    }
    // now free the ParallelComm resources
    if( data.pcomm )
    {
        delete data.pcomm;
        data.pcomm = nullptr;
    }
#endif

    // delete first all except vertices
    Range vertices = fileents.subset_by_type( MBVERTEX );
    Range noverts  = subtract( fileents, vertices );

    MB_CHK_SET_ERR( context.MBI->delete_entities( noverts ), "can't delete entities" );
    // now retrieve connected elements that still exist (maybe in other sets, pids?)
    Range adj_ents_left;
    MB_CHK_SET_ERR( context.MBI->get_adjacencies( vertices, 1, false, adj_ents_left, Interface::UNION ),
                    "can't get 1D adjacencies" );
    MB_CHK_SET_ERR( context.MBI->get_adjacencies( vertices, 2, false, adj_ents_left, Interface::UNION ),
                    "can't get 2D adjacencies" );
    MB_CHK_SET_ERR( context.MBI->get_adjacencies( vertices, 3, false, adj_ents_left, Interface::UNION ),
                    "can't get 3D adjacencies" );

    if( !adj_ents_left.empty() )
    {
        Range conn_verts;
        MB_CHK_SET_ERR( context.MBI->get_connectivity( adj_ents_left, conn_verts ), "can't get connectivity" );
        vertices = subtract( vertices, conn_verts );
    }

    MB_CHK_SET_ERR( context.MBI->delete_entities( vertices ), "can't delete vertices" );

    for( auto mit = context.appIdMap.begin(); mit != context.appIdMap.end(); ++mit )
    {
        if( *pid == mit->second )
        {
#ifdef MOAB_HAVE_MPI
            if( data.pcomm )
            {
                delete data.pcomm;
                data.pcomm = nullptr;
            }
#endif
            context.appIdMap.erase( mit );
            break;
        }
    }

    // now we can finally delete the application data itself
    context.appDatas.erase( appIterator );

    return moab::MB_SUCCESS;
}

/**
 * @brief Fortran-compatible wrapper for iMOAB_DeregisterApplication.
 * @ingroup iMOABApp
 *
 * @details Deregisters application from Fortran code by clearing Fortran flag and calling
 * C-style deregistration.
 *
 * @param[in] pid  Application ID to deregister
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if application not found
 *
 * @see iMOAB_DeregisterApplication()
 */
ErrCode iMOAB_DeregisterApplicationFortran( iMOAB_AppID pid )
{
    // Clear Fortran-specific flag before passing to C-style deregistration
    context.appDatas[*pid].is_fortran = false;

    // Release all data structure allocations via C-style function
    return iMOAB_DeregisterApplication( pid );
}

/**
 * @brief Read mesh file header information without loading full mesh.
 * @ingroup iMOABIO
 *
 * @details Quickly extracts metadata from HDF5 mesh file including vertex count, element count,
 * spatial dimension, and partition count. Useful for memory planning before full load.
 *
 * @par Implementation:
 * - Opens HDF5 file in read-only mode using mhdf library
 * - Reads file summary without loading actual mesh data
 * - Counts elements by type (edges, faces, regions)
 * - Extracts partition information from parallel decomposition tags
 *
 * @param[in]  filename             Path to HDF5 mesh file
 * @param[out] num_global_vertices  Total number of vertices in mesh
 * @param[out] num_global_elements  Total number of elements in mesh (all types)
 * @param[out] num_dimension        Spatial dimension (2D or 3D)
 * @param[out] num_parts            Number of parallel partitions
 *
 * @pre HDF5 support must be enabled (MOAB_HAVE_HDF5)
 * @pre File must be valid HDF5 mesh file
 *
 * @post Output parameters populated with file metadata
 * @post File remains closed (no mesh data loaded)
 *
 * @note Lightweight operation - does not load mesh into memory
 * @note All output parameters are optional (can be NULL)
 * @warning Returns failure if HDF5 support not enabled
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE on error
 *
 * @see iMOAB_LoadMesh()
 */
ErrCode iMOAB_ReadHeaderInfo( const iMOAB_String filename,
                              int* num_global_vertices,
                              int* num_global_elements,
                              int* num_dimension,
                              int* num_parts )
{
    IMOAB_CHECKPOINTER( filename, 1 );
    IMOAB_ASSERT( strlen( filename ), "Invalid filename length." );

#ifdef MOAB_HAVE_HDF5
    std::string filen( filename );

    int edges   = 0;
    int faces   = 0;
    int regions = 0;
    if( num_global_vertices ) *num_global_vertices = 0;
    if( num_global_elements ) *num_global_elements = 0;
    if( num_dimension ) *num_dimension = 0;
    if( num_parts ) *num_parts = 0;

    mhdf_FileHandle file;
    mhdf_Status status;
    unsigned long max_id;
    struct mhdf_FileDesc* data;

    file = mhdf_openFile( filen.c_str(), 0, &max_id, -1, &status );

    if( mhdf_isError( &status ) )
    {
        fprintf( stderr, "%s: %s\n", filename, mhdf_message( &status ) );
        return moab::MB_FAILURE;
    }

    data = mhdf_getFileSummary( file, H5T_NATIVE_ULONG, &status,
                                1 );  // will use extra set info; will get parallel partition tag info too!

    if( mhdf_isError( &status ) )
    {
        fprintf( stderr, "%s: %s\n", filename, mhdf_message( &status ) );
        return moab::MB_FAILURE;
    }

    if( num_dimension ) *num_dimension = data->nodes.vals_per_ent;
    if( num_global_vertices ) *num_global_vertices = (int)data->nodes.count;

    for( int i = 0; i < data->num_elem_desc; i++ )
    {
        struct mhdf_ElemDesc* el_desc = &( data->elems[i] );
        struct mhdf_EntDesc* ent_d    = &( el_desc->desc );

        if( 0 == strcmp( el_desc->type, mhdf_EDGE_TYPE_NAME ) )
        {
            edges += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_TRI_TYPE_NAME ) )
        {
            faces += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_QUAD_TYPE_NAME ) )
        {
            faces += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_POLYGON_TYPE_NAME ) )
        {
            faces += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_TET_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_PYRAMID_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_PRISM_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mdhf_KNIFE_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mdhf_HEX_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_POLYHEDRON_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }

        if( 0 == strcmp( el_desc->type, mhdf_SEPTAHEDRON_TYPE_NAME ) )
        {
            regions += ent_d->count;
        }
    }

    if( num_parts ) *num_parts = data->numEntSets[0];

    // is this required?
    if( edges > 0 )
    {
        if( num_dimension ) *num_dimension = 1;  // I don't think it will ever return 1
        if( num_global_elements ) *num_global_elements = edges;
    }

    if( faces > 0 )
    {
        if( num_dimension ) *num_dimension = 2;
        if( num_global_elements ) *num_global_elements = faces;
    }

    if( regions > 0 )
    {
        if( num_dimension ) *num_dimension = 3;
        if( num_global_elements ) *num_global_elements = regions;
    }

    mhdf_closeFile( file, &status );

    free( data );

#else
    std::cout << filename
              << ": Please reconfigure with HDF5. Cannot retrieve header information for file "
                 "formats other than a h5m file.\n";
    if( num_global_vertices ) *num_global_vertices = 0;
    if( num_global_elements ) *num_global_elements = 0;
    if( num_dimension ) *num_dimension = 0;
    if( num_parts ) *num_parts = 0;
#endif

    return moab::MB_SUCCESS;
}

/**
 * @brief Load mesh file into application's mesh set.
 * @ingroup iMOABIO
 *
 * @details Loads mesh from file (HDF5, Exodus, VTK, etc.) into application's mesh set with optional
 * ghost layer exchange for parallel simulations. Automatically configures parallel reading options.
 *
 * @par Loading Process:
 * -# Construct read options string from user options + automatic parallel settings
 * -# Add PARALLEL_COMM option for parallel HDF5/NetCDF files
 * -# Add ghost layer options (PARALLEL_GHOSTS, PARTITION_DISTRIBUTE) if requested
 * -# Load mesh into application's file set using MOAB::load_file()
 * -# Update mesh information (vertex/element counts, ranges, etc.)
 * -# Exchange ghost layers if num_ghost_layers > 0
 *
 * @param[in] pid              Application ID
 * @param[in] filename         Path to mesh file to load
 * @param[in] read_options     Optional reader options (semicolon-separated, can be NULL)
 * @param[in] num_ghost_layers Number of ghost element layers to exchange (0 for none)
 *
 * @pre Application must be registered via iMOAB_RegisterApplication
 * @pre File must exist and be readable
 * @pre For parallel: MPI must be initialized
 *
 * @post Mesh loaded into application's file set
 * @post Mesh info updated (vertices, elements, ranges)
 * @post Ghost layers exchanged if requested
 *
 * @note Supported formats: HDF5 (.h5m), Exodus (.exo), VTK (.vtk), NetCDF (.nc), etc.
 * @note PARALLEL_COMM option automatically added for parallel HDF5/NetCDF
 * @note Ghost exchange requires PARALLEL_GHOSTS and PARTITION_DISTRIBUTE options
 * @warning Do not manually specify PARALLEL_COMM in read_options (will cause error)
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_WriteMesh(), iMOAB_UpdateMeshInfo()
 */
ErrCode iMOAB_LoadMesh( iMOAB_AppID pid,
                        const iMOAB_String filename,
                        const iMOAB_String read_options,
                        int* num_ghost_layers )
{
    IMOAB_CHECKPOINTER( filename, 2 );
    IMOAB_ASSERT( strlen( filename ), "Invalid filename length." );
    IMOAB_CHECKPOINTER( num_ghost_layers, 4 );

    // make sure we use the file set and pcomm associated with the *pid
    std::ostringstream newopts;
    if( read_options ) newopts << read_options;

#ifdef MOAB_HAVE_MPI

    if( context.MPI_initialized )
    {
        if( context.worldprocs > 1 )
        {
            std::string opts( ( read_options ? read_options : "" ) );
            std::string pcid( "PARALLEL_COMM=" );
            std::size_t found = opts.find( pcid );

            if( found != std::string::npos )
            {
                std::cerr << " cannot specify PARALLEL_COMM option, it is implicit \n";
                return moab::MB_FAILURE;
            }

            // in serial, apply PARALLEL_COMM option only for h5m files; it does not work for .g
            // files (used in test_remapping)
            std::string filen( filename );
            std::string::size_type idx = filen.rfind( '.' );

            if( idx != std::string::npos )
            {
                ParallelComm* pco     = context.appDatas[*pid].pcomm;
                std::string extension = filen.substr( idx + 1 );
                if( ( extension == std::string( "h5m" ) ) || ( extension == std::string( "nc" ) ) )
                    newopts << ";;PARALLEL_COMM=" << pco->get_id();
            }

            if( *num_ghost_layers >= 1 )
            {
                // if we want ghosts, we will want additional entities, the last .1
                // because the addl ents can be edges, faces that are part of the neumann sets
                std::string pcid2( "PARALLEL_GHOSTS=" );
                std::size_t found2 = opts.find( pcid2 );

                if( found2 != std::string::npos )
                {
                    std::cout << " PARALLEL_GHOSTS option is already specified, ignore passed "
                                 "number of layers \n";
                }
                else
                {
                    // dimension of primary entities is 3 here, but it could be 2 for climate
                    // meshes; we would need to pass PARALLEL_GHOSTS explicitly for 2d meshes, for
                    // example:  ";PARALLEL_GHOSTS=2.0.1"
                    newopts << ";PARALLEL_GHOSTS=3.0." << *num_ghost_layers << ".3";
                }
            }
        }
    }
#else
    IMOAB_ASSERT( *num_ghost_layers == 0, "Cannot provide ghost layers in serial." );
#endif

    // Now let us actually load the MOAB file with the appropriate read options
    MB_CHK_SET_ERR( context.MBI->load_file( filename, &context.appDatas[*pid].file_set, newopts.str().c_str() ),
                    "can't load file" );

#ifdef VERBOSE
    // some debugging stuff
    std::ostringstream outfile;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.appDatas[*pid].pcomm;
    int rank          = pco->rank();
    int nprocs        = pco->size();
    outfile << "TaskMesh_n" << nprocs << "." << rank << ".h5m";
#else
    outfile << "TaskMesh_n1.0.h5m";
#endif
    // the mesh contains ghosts too, but they are not part of mat/neumann set
    // write in serial the file, to see what tags are missing
    MB_CHK_SET_ERR( context.MBI->write_file( outfile.str().c_str() ),
                    "can't write file" );  // everything on current task, written in serial
#endif

    // Update ghost layer information
    context.appDatas[*pid].num_ghost_layers = *num_ghost_layers;

    // Update mesh information
    return iMOAB_UpdateMeshInfo( pid );
}

static ErrCode internal_WriteMesh( iMOAB_AppID pid,
                                   const iMOAB_String filename,
                                   const iMOAB_String write_options,
                                   bool primary_set = true )
{
    IMOAB_CHECKPOINTER( filename, 2 );
    IMOAB_ASSERT( strlen( filename ), "Invalid filename length." );

    appData& data        = context.appDatas[*pid];
    EntityHandle fileSet = ( primary_set ? data.file_set : 0 );

    std::ostringstream newopts;
#ifdef MOAB_HAVE_MPI
    std::string write_opts( ( write_options ? write_options : "" ) );
    std::string pcid( "PARALLEL_COMM=" );

    if( write_opts.find( pcid ) != std::string::npos )
    {
        std::cerr << " cannot specify PARALLEL_COMM option, it is implicit \n";
        return moab::MB_FAILURE;
    }

    // if write in parallel, add pc option, to be sure about which ParallelComm instance is used
    std::string pw( "PARALLEL=WRITE_PART" );
    if( write_opts.find( pw ) != std::string::npos )
    {
        ParallelComm* pco = data.pcomm;
        newopts << "PARALLEL_COMM=" << pco->get_id() << ";";
    }

#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
    if( !primary_set )
    {
        if( data.tempestData.remapper != nullptr )
            fileSet = data.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
        else if( data.file_set != data.secondary_file_set )
            fileSet = data.secondary_file_set;
        else
            MB_CHK_SET_ERR( moab::MB_FAILURE, "Invalid secondary file set handle" );
    }
#endif

    // append user write options to the one we have built
    if( write_options ) newopts << write_options;

    std::vector< Tag > copyTagList = data.tagList;
    // append Global ID and Parallel Partition
    std::string gid_name_tag( "GLOBAL_ID" );

    // export global id tag, we need it always
    if( data.tagMap.find( gid_name_tag ) == data.tagMap.end() )
    {
        Tag gid = context.MBI->globalId_tag();
        copyTagList.push_back( gid );
    }
    // also Parallel_Partition PARALLEL_PARTITION
    std::string pp_name_tag( "PARALLEL_PARTITION" );

    // write parallel part tag too, if it exists
    if( data.tagMap.find( pp_name_tag ) == data.tagMap.end() )
    {
        Tag ptag;
        context.MBI->tag_get_handle( pp_name_tag.c_str(), ptag );
        if( ptag ) copyTagList.push_back( ptag );
    }

    // Now let us actually write the file to disk with appropriate options
    if( primary_set )
    {
        MB_CHK_ERR( context.MBI->write_file( filename, 0, newopts.str().c_str(), &fileSet, 1, copyTagList.data(),
                                             copyTagList.size() ) );
    }
    else
    {
        MB_CHK_ERR( context.MBI->write_file( filename, 0, newopts.str().c_str(), &fileSet, 1 ) );
    }

    return moab::MB_SUCCESS;
}

/**
 * @brief Write application's mesh to file.
 * @ingroup iMOABIO
 *
 * @details Writes mesh data from application's mesh set to file in various formats (HDF5, Exodus, VTK, etc.).
 * Automatically includes GLOBAL_ID and PARALLEL_PARTITION tags for parallel simulations.
 *
 * @par Writing Process:
 * -# Construct write options from user options + automatic parallel settings
 * -# Add PARALLEL=WRITE_PART for parallel HDF5 output if needed
 * -# Append GLOBAL_ID tag if not already in tag list
 * -# Append PARALLEL_PARTITION tag if available and not in tag list
 * -# Write mesh using MOAB::write_file() with appropriate options
 *
 * @param[in] pid           Application ID
 * @param[in] filename      Output file path
 * @param[in] write_options Optional writer options (semicolon-separated, can be NULL)
 *
 * @pre Application must be registered and mesh loaded
 * @pre For parallel: output directory must be writable by all processes
 *
 * @post Mesh written to file with requested format
 * @post GLOBAL_ID and PARALLEL_PARTITION tags included if available
 *
 * @note Supported formats: HDF5 (.h5m), Exodus (.exo), VTK (.vtk), etc.
 * @note For parallel writes, use PARALLEL=WRITE_PART option
 * @note Tags in application's tagList are automatically written
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_LoadMesh(), iMOAB_WriteLocalMesh()
 */
ErrCode iMOAB_WriteMesh( iMOAB_AppID pid, const iMOAB_String filename, const iMOAB_String write_options )
{
    return internal_WriteMesh( pid, filename, write_options );
}

/**
 * @brief Write each process's local mesh to separate file for debugging.
 * @ingroup iMOABIO
 *
 * @details Creates per-process HDF5 file with naming convention: prefix_nprocs_rank.h5m.
 * Useful for debugging parallel mesh distribution and ghost layer issues.
 *
 * @param[in] pid     Application ID
 * @param[in] prefix  Output file prefix (rank and process count appended automatically)
 *
 * @pre Application must be registered and mesh loaded
 *
 * @post Each process writes its local mesh to: prefix_<nprocs>_<rank>.h5m
 *
 * @note Output files show per-process mesh partition (useful for parallel debugging)
 * @note Files include both owned and ghost entities
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_WriteMesh()
 */
ErrCode iMOAB_WriteLocalMesh( iMOAB_AppID pid, iMOAB_String prefix )
{
    IMOAB_CHECKPOINTER( prefix, 2 );
    IMOAB_ASSERT( strlen( prefix ), "Invalid prefix string length." );

    // Construct filename with parallel decomposition info: prefix_nprocs_rank.h5m
    std::ostringstream file_name;
    int rank = 0, size = 1;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = context.appDatas[*pid].pcomm;
    rank                = pcomm->rank();
    size                = pcomm->size();
#endif
    file_name << prefix << "_" << size << "_" << rank << ".h5m";

    // Write this process's local mesh partition to its own file
    MB_CHK_ERR( context.MBI->write_file( file_name.str().c_str(), 0, 0, &context.appDatas[*pid].file_set, 1 ) );

    return moab::MB_SUCCESS;
}

/**
 * @brief Update application's mesh information after loading or modification.
 * @ingroup iMOABQuery
 *
 * @details Refreshes all cached mesh metadata including vertex/element counts, ranges, owned/ghost
 * entities, and material/boundary condition sets. Must be called after mesh loading or ghosting.
 *
 * @par Update Process:
 * -# Clear all existing cached ranges (vertices, elements, sets)
 * -# Query all vertices from file set
 * -# Determine mesh dimension (3D → 2D → 1D → 0D) by checking for elements
 * -# Extract primary elements of detected dimension
 * -# Separate owned vs ghost entities based on parallel ownership
 * -# Query material sets (MATERIAL_SET tag)
 * -# Query boundary condition sets (NEUMANN_SET, DIRICHLET_SET tags)
 * -# Separate local vs owned vertices based on parallel partitioning
 *
 * @param[in] pid  Application ID
 *
 * @pre Application must be registered and mesh loaded
 *
 * @post Mesh dimension determined and cached
 * @post All vertex/element ranges updated
 * @post Owned/ghost entity ranges computed
 * @post Material and BC sets identified
 *
 * @note Call this after iMOAB_LoadMesh or ghost layer exchange
 * @note Automatically called by iMOAB_LoadMesh
 * @warning Clears all existing cached mesh information
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_LoadMesh(), iMOAB_GetMeshInfo()
 */
ErrCode iMOAB_UpdateMeshInfo( iMOAB_AppID pid )
{
    // this will include ghost elements info
    appData& data        = context.appDatas[*pid];
    EntityHandle fileSet = data.file_set;

    // first clear all data ranges; this can be called after ghosting
    data.all_verts.clear();
    data.primary_elems.clear();
    data.local_verts.clear();
    data.owned_verts.clear();
    data.ghost_vertices.clear();
    data.owned_elems.clear();
    data.ghost_elems.clear();
    data.mat_sets.clear();
    data.neu_sets.clear();
    data.diri_sets.clear();

    // Let us get all the vertex entities
    MB_CHK_SET_ERR( context.MBI->get_entities_by_type( fileSet, MBVERTEX, data.all_verts, true ),
                    "can't get vertices" );

    // Let us check first entities of dimension = 3
    data.dimension = 3;
    MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( fileSet, data.dimension, data.primary_elems, true ),
                    "can't get primary elements" );

    if( data.primary_elems.empty() )
    {
        // Now 3-D elements. Let us check entities of dimension = 2
        data.dimension = 2;
        MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( fileSet, data.dimension, data.primary_elems, true ),
                        "can't get primary elements" );

        if( data.primary_elems.empty() )
        {
            // Now 3-D/2-D elements. Let us check entities of dimension = 1
            data.dimension = 1;
            MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( fileSet, data.dimension, data.primary_elems, true ),
                            "can't get primary elements" );

            if( data.primary_elems.empty() )
            {
                // no elements of dimension 1 or 2 or 3; it could happen for point clouds
                data.dimension = 0;
            }
        }
    }

    // check if the current mesh is just a point cloud
    data.point_cloud = ( ( data.primary_elems.size() == 0 && data.all_verts.size() > 0 ) || data.dimension == 0 );

#ifdef MOAB_HAVE_MPI

    if( context.MPI_initialized )
    {
        ParallelComm* pco = context.appDatas[*pid].pcomm;

        // filter ghost vertices, from local
        MB_CHK_SET_ERR( pco->filter_pstatus( data.all_verts, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.local_verts ),
                        "can't filter ghost vertices" );

        // Store handles for all ghosted entities
        data.ghost_vertices = subtract( data.all_verts, data.local_verts );

        // filter ghost elements, from local
        MB_CHK_SET_ERR( pco->filter_pstatus( data.primary_elems, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.owned_elems ),
                        "can't filter ghost elements" );

        data.ghost_elems = subtract( data.primary_elems, data.owned_elems );
        // now update global number of primary cells and global number of vertices
        // determine first number of owned vertices
        // Get local owned vertices
        MB_CHK_SET_ERR( pco->filter_pstatus( data.all_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &data.owned_verts ),
                        "can't filter ghost vertices" );
        int local[2], global[2];
        local[0] = data.owned_verts.size();
        local[1] = data.owned_elems.size();
        MPI_Allreduce( local, global, 2, MPI_INT, MPI_SUM, pco->comm() );
        MB_CHK_SET_ERR( iMOAB_SetGlobalInfo( pid, &( global[0] ), &( global[1] ) ), "can't set global info" );
    }
    else
    {
        data.local_verts = data.all_verts;
        data.owned_elems = data.primary_elems;
    }

#else

    data.local_verts = data.all_verts;
    data.owned_elems = data.primary_elems;

#endif

    // Get the references for some standard internal tags such as material blocks, BCs, etc
    MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.material_tag ), 0, 1,
                                                               data.mat_sets, Interface::UNION ),
                    "can't get material sets" );

    MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.neumann_tag ), 0, 1,
                                                               data.neu_sets, Interface::UNION ),
                    "can't get neumann sets" );

    MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.dirichlet_tag ), 0, 1,
                                                               data.diri_sets, Interface::UNION ),
                    "can't get dirichlet sets" );

    return moab::MB_SUCCESS;
}

/**
 * @brief Query mesh statistics including vertices, elements, and boundary conditions.
 * @ingroup iMOABQuery
 *
 * @details Retrieves counts of vertices, elements, material blocks, and boundary condition sets.
 * Each output array has 3 components: [owned, ghost, total]. Useful for memory allocation
 * and understanding parallel mesh distribution.
 *
 * @par Array Format:
 * All output arrays use 3-element format: [owned, ghost, total]
 * - owned: Entities owned by this process
 * - ghost: Ghost/halo entities from neighboring processes
 * - total: Sum of owned and ghost entities
 *
 * @param[in]  pid                   Application ID
 * @param[out] num_visible_vertices  Vertex counts [local(non-ghost), ghost, total]
 * @param[out] num_visible_elements  Element counts [owned, ghost, total]
 * @param[out] num_visible_blocks    Material block counts [owned, ghost, total]
 * @param[out] num_visible_surfaceBC Surface BC counts (total face count across all BC sets)
 * @param[out] num_visible_vertexBC  Vertex BC counts (total vertex count across all BC sets)
 *
 * @pre Application must be registered and mesh loaded
 * @pre iMOAB_UpdateMeshInfo must have been called
 *
 * @post Output arrays populated with mesh statistics
 *
 * @note All parameters are optional (can be NULL if not needed)
 * @note Surface BC count includes all faces in Neumann sets
 * @note Vertex BC count includes all vertices in Dirichlet sets
 * @note Material blocks queried via MATERIAL_SET tag
 * @note BCs queried via NEUMANN_SET and DIRICHLET_SET tags
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_UpdateMeshInfo(), iMOAB_LoadMesh()
 */
ErrCode iMOAB_GetMeshInfo( iMOAB_AppID pid,
                           int* num_visible_vertices,
                           int* num_visible_elements,
                           int* num_visible_blocks,
                           int* num_visible_surfaceBC,
                           int* num_visible_vertexBC )
{
    appData& data        = context.appDatas[*pid];
    EntityHandle fileSet = data.file_set;

    // this will include ghost elements
    // first clear all data ranges; this can be called after ghosting
    if( num_visible_elements )
    {
        num_visible_elements[2] = static_cast< int >( data.primary_elems.size() );
        // separate ghost and local/owned primary elements
        num_visible_elements[0] = static_cast< int >( data.owned_elems.size() );
        num_visible_elements[1] = static_cast< int >( data.ghost_elems.size() );
    }
    if( num_visible_vertices )
    {
        num_visible_vertices[2] = static_cast< int >( data.all_verts.size() );
        num_visible_vertices[1] = static_cast< int >( data.ghost_vertices.size() );
        // local are those that are not ghosts; they may include shared, not owned vertices
        num_visible_vertices[0] = num_visible_vertices[2] - num_visible_vertices[1];
    }

    if( num_visible_blocks )
    {
        MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.material_tag ), 0,
                                                                   1, data.mat_sets, Interface::UNION ),
                        "can't get material sets" );

        num_visible_blocks[2] = data.mat_sets.size();
        num_visible_blocks[0] = num_visible_blocks[2];
        num_visible_blocks[1] = 0;
    }

    if( num_visible_surfaceBC )
    {
        MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.neumann_tag ), 0, 1,
                                                                   data.neu_sets, Interface::UNION ),
                        "can't get neumann sets" );

        num_visible_surfaceBC[2] = 0;
        // count how many faces are in each neu set, and how many regions are
        // adjacent to them;
        int numNeuSets = (int)data.neu_sets.size();

        for( int i = 0; i < numNeuSets; i++ )
        {
            Range subents;
            EntityHandle nset = data.neu_sets[i];
            MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( nset, data.dimension - 1, subents ),
                            "can't get neumann sets" );

            for( Range::iterator it = subents.begin(); it != subents.end(); ++it )
            {
                EntityHandle subent = *it;
                Range adjPrimaryEnts;
                MB_CHK_SET_ERR( context.MBI->get_adjacencies( &subent, 1, data.dimension, false, adjPrimaryEnts ),
                                "can't get adjacencies" );

                num_visible_surfaceBC[2] += (int)adjPrimaryEnts.size();
            }
        }

        num_visible_surfaceBC[0] = num_visible_surfaceBC[2];
        num_visible_surfaceBC[1] = 0;
    }

    if( num_visible_vertexBC )
    {
        MB_CHK_SET_ERR( context.MBI->get_entities_by_type_and_tag( fileSet, MBENTITYSET, &( context.dirichlet_tag ), 0,
                                                                   1, data.diri_sets, Interface::UNION ),
                        "can't get dirichlet sets" );

        num_visible_vertexBC[2] = 0;
        int numDiriSets         = (int)data.diri_sets.size();

        for( int i = 0; i < numDiriSets; i++ )
        {
            Range verts;
            EntityHandle diset = data.diri_sets[i];
            MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( diset, 0, verts ), "can't get dirichlet sets" );

            num_visible_vertexBC[2] += (int)verts.size();
        }

        num_visible_vertexBC[0] = num_visible_vertexBC[2];
        num_visible_vertexBC[1] = 0;
    }

    return moab::MB_SUCCESS;
}

/**
 * @brief Retrieve global IDs for all vertices in application's mesh.
 * @ingroup iMOABQuery
 *
 * @details Extracts GLOBAL_ID tag values for all vertices. Global IDs are unique across all
 * processes and are used for parallel mesh operations and entity identification.
 *
 * @param[in]  pid              Application ID
 * @param[in]  vertices_length  Number of vertices (must match actual vertex count)
 * @param[out] global_vertex_ID Array of global IDs (size = vertices_length)
 *
 * @pre Application must be registered and mesh loaded
 * @pre vertices_length must equal total vertex count from iMOAB_GetMeshInfo
 * @pre global_vertex_ID array must be pre-allocated
 *
 * @post global_vertex_ID populated with GLOBAL_ID tag values
 *
 * @note Global IDs are 1-indexed and unique across all processes
 * @note Array must be sized correctly or assertion will fail
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 *
 * @see iMOAB_GetMeshInfo(), iMOAB_GetVertexOwnership()
 */
ErrCode iMOAB_GetVertexID( iMOAB_AppID pid, int* vertices_length, iMOAB_GlobalID* global_vertex_ID )
{
    IMOAB_CHECKPOINTER( vertices_length, 2 );
    IMOAB_CHECKPOINTER( global_vertex_ID, 3 );

    const Range& verts = context.appDatas[*pid].all_verts;

    // Validate array size matches actual vertex count to prevent buffer overruns
    IMOAB_ASSERT( *vertices_length == static_cast< int >( verts.size() ), "Invalid vertices length provided" );

    // Extract GLOBAL_ID tag values for all vertices in order
    return context.MBI->tag_get_data( context.globalID_tag, verts, global_vertex_ID );
}

/**
 * @brief Retrieve parallel ownership information for all vertices.
 * @ingroup iMOABQuery
 *
 * @details Queries which MPI rank owns each vertex. In parallel simulations, each vertex is owned
 * by exactly one process, though it may be ghosted on others.
 *
 * @param[in]  pid                      Application ID
 * @param[in]  vertices_length          Number of vertices (must match actual count)
 * @param[out] visible_global_rank_ID   Array of owning ranks (size = vertices_length)
 *
 * @pre Application must be registered and mesh loaded
 * @pre For parallel: mesh must have parallel ownership information
 * @pre visible_global_rank_ID array must be pre-allocated
 *
 * @post visible_global_rank_ID[i] = rank that owns vertex i
 *
 * @note In serial runs, all vertices owned by rank 0
 * @note Ghost vertices report the rank of their owner, not local rank
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if count mismatch
 *
 * @see iMOAB_GetVertexID(), iMOAB_GetElementOwnership()
 */
ErrCode iMOAB_GetVertexOwnership( iMOAB_AppID pid, int* vertices_length, int* visible_global_rank_ID )
{
    assert( vertices_length && *vertices_length );
    assert( visible_global_rank_ID );

    Range& verts = context.appDatas[*pid].all_verts;

#ifdef MOAB_HAVE_MPI
    // Parallel: query actual ownership from ParallelComm
    ParallelComm* pco = context.appDatas[*pid].pcomm;

    int i = 0;
    for( Range::iterator vit = verts.begin(); vit != verts.end(); ++vit, i++ )
    {
        // Get owning rank for this vertex (may be local or remote)
        MB_CHK_SET_ERR( pco->get_owner( *vit, visible_global_rank_ID[i] ), "can't get owner" );
    }

    // Validate we processed exactly the expected number of vertices
    if( i != *vertices_length )
    {
        return moab::MB_INVALID_SIZE;  // Array size mismatch
    }

#else
    // Serial: all vertices owned by rank 0
    if( (int)verts.size() != *vertices_length )
    {
        return moab::MB_INVALID_SIZE;  // Array size mismatch
    }  // warning array allocation problem

    int i = 0;
    for( Range::iterator vit = verts.begin(); vit != verts.end(); ++vit, i++ )
    {
        visible_global_rank_ID[i] = 0;
    }  // all vertices are owned by processor 0, as this is serial run

#endif

    return moab::MB_SUCCESS;
}

/**
 * @brief Retrieve coordinates for all vertices in interleaved format.
 * @ingroup iMOABQuery
 *
 * @details Extracts 3D coordinates (x,y,z) for all vertices in interleaved format:
 * [x0,y0,z0, x1,y1,z1, ...]. Coordinates are returned in same order as vertex IDs.
 *
 * @param[in]  pid           Application ID
 * @param[in]  coords_length Length of coordinates array (must equal 3 * num_vertices)
 * @param[out] coordinates   Interleaved coordinate array (size = 3 * num_vertices)
 *
 * @pre Application must be registered and mesh loaded
 * @pre coords_length must equal 3 * vertex count
 * @pre coordinates array must be pre-allocated
 *
 * @post coordinates filled with interleaved x,y,z values
 *
 * @note Format: [x0,y0,z0, x1,y1,z1, x2,y2,z2, ...]
 * @note Always returns 3D coordinates even for 2D meshes (z=0)
 * @warning Array size must be exact or function returns failure
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if size mismatch
 *
 * @see iMOAB_GetVertexID(), iMOAB_GetMeshInfo()
 */
ErrCode iMOAB_GetVisibleVerticesCoordinates( iMOAB_AppID pid, int* coords_length, double* coordinates )
{
    Range& verts = context.appDatas[*pid].all_verts;

    // Validate array size: need 3 coordinates (x,y,z) per vertex
    if( *coords_length != 3 * (int)verts.size() )
    {
        return moab::MB_INVALID_SIZE;  // Size mismatch
    }

    // Extract coordinates in interleaved format (deep copy)
    MB_CHK_SET_ERR( context.MBI->get_coords( verts, coordinates ), "can't get coordinates" );

    return moab::MB_SUCCESS;
}

/**
 * @brief Retrieve global IDs for all material blocks in mesh.
 * @ingroup iMOABQuery
 *
 * @details Extracts MATERIAL_SET tag values for all material blocks. Blocks are identified by
 * unique integer IDs and contain elements with similar material properties.
 *
 * @param[in]  pid              Application ID
 * @param[in]  block_length     Number of blocks (must match actual count)
 * @param[out] global_block_IDs Array of material set IDs (size = block_length)
 *
 * @pre Application must be registered and mesh loaded
 * @pre block_length must match block count from iMOAB_GetMeshInfo
 * @pre global_block_IDs array must be pre-allocated
 *
 * @post global_block_IDs populated with MATERIAL_SET tag values
 * @post Internal matIndex map populated for fast block lookup
 *
 * @note Material blocks are mesh sets containing elements of same material
 * @note Creates internal index map for subsequent block queries
 * @note Block IDs are user-defined and may not be consecutive
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if count mismatch
 *
 * @see iMOAB_GetBlockInfo(), iMOAB_GetMeshInfo()
 */
ErrCode iMOAB_GetBlockID( iMOAB_AppID pid, int* block_length, iMOAB_GlobalID* global_block_IDs )
{
    Range& matSets = context.appDatas[*pid].mat_sets;

    // Validate array size matches actual block count
    if( *block_length != (int)matSets.size() )
    {
        return moab::MB_INVALID_SIZE;  // Size mismatch
    }

    // Extract MATERIAL_SET tag values for all material blocks
    MB_CHK_SET_ERR( context.MBI->tag_get_data( context.material_tag, matSets, global_block_IDs ),
                    "can't get material IDs" );

    // Build internal index map: block_id -> array_index for fast lookup
    std::map< int, int >& matIdx = context.appDatas[*pid].matIndex;
    for( unsigned i = 0; i < matSets.size(); i++ )
    {
        matIdx[global_block_IDs[i]] = i;  // Cache block index for future queries
    }

    return moab::MB_SUCCESS;
}

/**
 * @brief Retrieve information about a specific material block.
 * @ingroup iMOABQuery
 *
 * @details Queries element count and connectivity size for specified material block.
 * Assumes all elements in block have same topology type.
 *
 * @param[in]  pid                   Application ID
 * @param[in]  global_block_ID       Material block ID to query
 * @param[out] vertices_per_element  Number of vertices per element in this block
 * @param[out] num_elements_in_block Total number of elements in this block
 *
 * @pre Application must be registered and mesh loaded
 * @pre iMOAB_GetBlockID must have been called first (to populate matIndex)
 * @pre global_block_ID must be valid block ID from iMOAB_GetBlockID
 *
 * @post vertices_per_element set to connectivity size
 * @post num_elements_in_block set to element count
 *
 * @note All elements in block assumed to have same topology
 * @note vertices_per_element depends on element type (4=tet, 8=hex, etc.)
 * @warning Returns failure if block ID not found in matIndex map
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if block not found
 *
 * @see iMOAB_GetBlockID(), iMOAB_GetBlockElementConnectivities()
 */
ErrCode iMOAB_GetBlockInfo( iMOAB_AppID pid,
                            iMOAB_GlobalID* global_block_ID,
                            int* vertices_per_element,
                            int* num_elements_in_block )
{
    assert( global_block_ID );

    std::map< int, int >& matMap      = context.appDatas[*pid].matIndex;
    std::map< int, int >::iterator it = matMap.find( *global_block_ID );

    if( it == matMap.end() )
    {
        return moab::MB_FAILURE;
    }  // error in finding block with id

    int blockIndex          = matMap[*global_block_ID];
    EntityHandle matMeshSet = context.appDatas[*pid].mat_sets[blockIndex];
    Range blo_elems;
    MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( matMeshSet, blo_elems ), "can't get block elements" );

    if( blo_elems.empty() ) return moab::MB_FAILURE;

    EntityType type = context.MBI->type_from_handle( blo_elems[0] );

    if( !blo_elems.all_of_type( type ) ) return moab::MB_FAILURE;

    const EntityHandle* conn = nullptr;
    int num_verts            = 0;
    MB_CHK_SET_ERR( context.MBI->get_connectivity( blo_elems[0], conn, num_verts ), "can't get connectivity" );

    *vertices_per_element  = num_verts;
    *num_elements_in_block = (int)blo_elems.size();

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetVisibleElementsInfo( iMOAB_AppID pid,
                                      int* num_visible_elements,
                                      iMOAB_GlobalID* element_global_IDs,
                                      int* ranks,
                                      iMOAB_GlobalID* block_IDs )
{
    appData& data = context.appDatas[*pid];
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.appDatas[*pid].pcomm;
#endif

    MB_CHK_SET_ERR( context.MBI->tag_get_data( context.globalID_tag, data.primary_elems, element_global_IDs ),
                    "can't get global IDs" );

    int i = 0;

    for( Range::iterator eit = data.primary_elems.begin(); eit != data.primary_elems.end(); ++eit, ++i )
    {
#ifdef MOAB_HAVE_MPI
        MB_CHK_SET_ERR( pco->get_owner( *eit, ranks[i] ), "can't get owner" );

#else
        /* everything owned by task 0 */
        ranks[i] = 0;
#endif
    }

    for( Range::iterator mit = data.mat_sets.begin(); mit != data.mat_sets.end(); ++mit )
    {
        EntityHandle matMeshSet = *mit;
        Range elems;
        MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( matMeshSet, elems ), "can't get block elements" );

        int valMatTag;
        MB_CHK_SET_ERR( context.MBI->tag_get_data( context.material_tag, &matMeshSet, 1, &valMatTag ),
                        "can't get material tag" );

        for( Range::iterator eit = elems.begin(); eit != elems.end(); ++eit )
        {
            EntityHandle eh = *eit;
            int index       = data.primary_elems.index( eh );

            if( -1 == index )
            {
                return moab::MB_FAILURE;
            }

            if( -1 >= *num_visible_elements )
            {
                return moab::MB_FAILURE;
            }

            block_IDs[index] = valMatTag;
        }
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetBlockElementConnectivities( iMOAB_AppID pid,
                                             iMOAB_GlobalID* global_block_ID,
                                             int* connectivity_length,
                                             int* element_connectivity )
{
    assert( global_block_ID );      // ensure global block ID argument is not null
    assert( connectivity_length );  // ensure connectivity length argument is not null

    appData& data                     = context.appDatas[*pid];
    std::map< int, int >& matMap      = data.matIndex;
    std::map< int, int >::iterator it = matMap.find( *global_block_ID );

    if( it == matMap.end() )
    {
        return moab::MB_FAILURE;
    }  // error in finding block with id

    int blockIndex          = matMap[*global_block_ID];
    EntityHandle matMeshSet = data.mat_sets[blockIndex];
    std::vector< EntityHandle > elems;

    MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( matMeshSet, elems ), "can't get block elements" );

    if( elems.empty() ) return moab::MB_FAILURE;

    std::vector< EntityHandle > vconnect;
    MB_CHK_SET_ERR( context.MBI->get_connectivity( &elems[0], elems.size(), vconnect ), "can't get connectivity" );

    if( *connectivity_length != (int)vconnect.size() )
    {
        return moab::MB_FAILURE;
    }  // mismatched sizes

    for( int i = 0; i < *connectivity_length; i++ )
    {
        int inx = data.all_verts.index( vconnect[i] );

        if( -1 == inx )
        {
            return moab::MB_FAILURE;
        }  // error, vertex not in local range

        element_connectivity[i] = inx;
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetElementConnectivity( iMOAB_AppID pid,
                                      iMOAB_LocalID* elem_index,
                                      int* connectivity_length,
                                      int* element_connectivity )
{
    assert( elem_index );           // ensure element index argument is not null
    assert( connectivity_length );  // ensure connectivity length argument is not null

    appData& data = context.appDatas[*pid];
    assert( ( *elem_index >= 0 ) && ( *elem_index < (int)data.primary_elems.size() ) );

    int num_nodes;
    const EntityHandle* conn;

    EntityHandle eh = data.primary_elems[*elem_index];

    MB_CHK_SET_ERR( context.MBI->get_connectivity( eh, conn, num_nodes ), "can't get connectivity" );

    if( *connectivity_length < num_nodes )
    {
        return moab::MB_FAILURE;
    }  // wrong number of vertices

    for( int i = 0; i < num_nodes; i++ )
    {
        int index = data.all_verts.index( conn[i] );

        if( -1 == index )
        {
            return moab::MB_FAILURE;
        }

        element_connectivity[i] = index;
    }

    *connectivity_length = num_nodes;

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetElementOwnership( iMOAB_AppID pid,
                                   iMOAB_GlobalID* global_block_ID,
                                   int* num_elements_in_block,
                                   int* element_ownership )
{
    assert( global_block_ID );        // ensure global block ID argument is not null
    assert( num_elements_in_block );  // ensure number of elements in block argument is not null

    std::map< int, int >& matMap = context.appDatas[*pid].matIndex;

    std::map< int, int >::iterator it = matMap.find( *global_block_ID );

    if( it == matMap.end() )
    {
        return moab::MB_FAILURE;
    }  // error in finding block with id

    int blockIndex          = matMap[*global_block_ID];
    EntityHandle matMeshSet = context.appDatas[*pid].mat_sets[blockIndex];
    Range elems;

    MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( matMeshSet, elems ), "can't get block elements" );

    if( elems.empty() ) return moab::MB_FAILURE;

    if( *num_elements_in_block != (int)elems.size() )
    {
        return moab::MB_FAILURE;
    }  // bad memory allocation

    int i = 0;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.appDatas[*pid].pcomm;
#endif

    for( Range::iterator vit = elems.begin(); vit != elems.end(); vit++, i++ )
    {
#ifdef MOAB_HAVE_MPI
        MB_CHK_SET_ERR( pco->get_owner( *vit, element_ownership[i] ), "can't get owner" );
#else
        element_ownership[i] = 0; /* owned by 0 */
#endif
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetElementID( iMOAB_AppID pid,
                            iMOAB_GlobalID* global_block_ID,
                            int* num_elements_in_block,
                            iMOAB_GlobalID* global_element_ID,
                            iMOAB_LocalID* local_element_ID )
{
    assert( global_block_ID );        // ensure global block ID argument is not null
    assert( num_elements_in_block );  // ensure number of elements in block argument is not null

    appData& data                = context.appDatas[*pid];
    std::map< int, int >& matMap = data.matIndex;

    std::map< int, int >::iterator it = matMap.find( *global_block_ID );

    if( it == matMap.end() )
    {
        return moab::MB_FAILURE;
    }  // error in finding block with id

    int blockIndex          = matMap[*global_block_ID];
    EntityHandle matMeshSet = data.mat_sets[blockIndex];
    Range elems;
    MB_CHK_SET_ERR( context.MBI->get_entities_by_handle( matMeshSet, elems ), "can't get block elements" );

    if( elems.empty() ) return moab::MB_FAILURE;

    if( *num_elements_in_block != (int)elems.size() )
    {
        return moab::MB_FAILURE;
    }  // bad memory allocation

    MB_CHK_SET_ERR( context.MBI->tag_get_data( context.globalID_tag, elems, global_element_ID ),
                    "can't get global IDs" );

    // check that elems are among primary_elems in data
    for( int i = 0; i < *num_elements_in_block; i++ )
    {
        local_element_ID[i] = data.primary_elems.index( elems[i] );

        if( -1 == local_element_ID[i] )
        {
            return moab::MB_FAILURE;
        }  // error, not in local primary elements
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetPointerToSurfaceBC( iMOAB_AppID pid,
                                     int* surface_BC_length,
                                     iMOAB_LocalID* local_element_ID,
                                     int* reference_surface_ID,
                                     int* boundary_condition_value )
{
    // we have to fill bc data for neumann sets;/

    // it was counted above, in GetMeshInfo
    appData& data  = context.appDatas[*pid];
    int numNeuSets = (int)data.neu_sets.size();

    int index = 0;  // index [0, surface_BC_length) for the arrays returned

    for( int i = 0; i < numNeuSets; i++ )
    {
        Range subents;
        EntityHandle nset = data.neu_sets[i];
        MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( nset, data.dimension - 1, subents ),
                        "can't get subentities" );

        int neuVal;
        MB_CHK_SET_ERR( context.MBI->tag_get_data( context.neumann_tag, &nset, 1, &neuVal ), "can't get neumann tag" );

        for( Range::iterator it = subents.begin(); it != subents.end(); ++it )
        {
            EntityHandle subent = *it;
            Range adjPrimaryEnts;
            MB_CHK_SET_ERR( context.MBI->get_adjacencies( &subent, 1, data.dimension, false, adjPrimaryEnts ),
                            "can't get adjacencies" );

            // get global id of the primary ents, and side number of the quad/subentity
            // this is moab ordering
            for( Range::iterator pit = adjPrimaryEnts.begin(); pit != adjPrimaryEnts.end(); pit++ )
            {
                EntityHandle primaryEnt = *pit;
                // get global id
                /*int globalID;
                rval = context.MBI->tag_get_data(gtags[3], &primaryEnt, 1, &globalID);
                if (MB_SUCCESS!=rval)
                  return moab::MB_FAILURE;
                global_element_ID[index] = globalID;*/

                // get local element id
                local_element_ID[index] = data.primary_elems.index( primaryEnt );

                if( -1 == local_element_ID[index] )
                {
                    return moab::MB_FAILURE;
                }  // did not find the element locally

                int side_number, sense, offset;
                MB_CHK_SET_ERR( context.MBI->side_number( primaryEnt, subent, side_number, sense, offset ),
                                "can't get side number" );

                reference_surface_ID[index]     = side_number + 1;  // moab is from 0 to 5, it needs 1 to 6
                boundary_condition_value[index] = neuVal;
                index++;
            }
        }
    }

    if( index != *surface_BC_length )
    {
        return moab::MB_FAILURE;
    }  // error in array allocations

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetPointerToVertexBC( iMOAB_AppID pid,
                                    int* vertex_BC_length,
                                    iMOAB_LocalID* local_vertex_ID,
                                    int* boundary_condition_value )
{
    // it was counted above, in GetMeshInfo
    appData& data   = context.appDatas[*pid];
    int numDiriSets = (int)data.diri_sets.size();
    int index       = 0;  // index [0, *vertex_BC_length) for the arrays returned

    for( int i = 0; i < numDiriSets; i++ )
    {
        Range verts;
        EntityHandle diset = data.diri_sets[i];
        MB_CHK_SET_ERR( context.MBI->get_entities_by_dimension( diset, 0, verts ), "can't get vertices" );

        int diriVal;
        MB_CHK_SET_ERR( context.MBI->tag_get_data( context.dirichlet_tag, &diset, 1, &diriVal ),
                        "can't get dirichlet tag" );

        for( Range::iterator vit = verts.begin(); vit != verts.end(); ++vit )
        {
            EntityHandle vt = *vit;
            /*int vgid;
            rval = context.MBI->tag_get_data(gtags[3], &vt, 1, &vgid);
            if (MB_SUCCESS!=rval)
              return moab::MB_FAILURE;
            global_vertext_ID[index] = vgid;*/
            local_vertex_ID[index] = data.all_verts.index( vt );

            if( -1 == local_vertex_ID[index] )
            {
                return moab::MB_FAILURE;
            }  // vertex was not found

            boundary_condition_value[index] = diriVal;
            index++;
        }
    }

    if( *vertex_BC_length != index )
    {
        return moab::MB_FAILURE;
    }  // array allocation issue

    return moab::MB_SUCCESS;
}

// Utility function
static void split_tag_names( std::string input_names,
                             std::string& separator,
                             std::vector< std::string >& list_tag_names )
{
    size_t pos = 0;
    std::string token;
    while( ( pos = input_names.find( separator ) ) != std::string::npos )
    {
        token = input_names.substr( 0, pos );
        if( !token.empty() ) list_tag_names.push_back( token );
        // std::cout << token << std::endl;
        input_names.erase( 0, pos + separator.length() );
    }
    if( !input_names.empty() )
    {
        // if leftover something, or if not ended with delimiter
        list_tag_names.push_back( input_names );
    }
    return;
}

ErrCode iMOAB_DefineTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* tag_type,
                                int* components_per_entity,
                                int* tag_index )
{
    // we have 6 types of tags supported so far
    // check if tag type is valid
    if( *tag_type < 0 || *tag_type > 5 ) return moab::MB_FAILURE;

    DataType tagDataType;
    TagType tagType;
    void* defaultVal        = nullptr;
    int* defInt             = new int[*components_per_entity];
    double* defDouble       = new double[*components_per_entity];
    EntityHandle* defHandle = new EntityHandle[*components_per_entity];

    for( int i = 0; i < *components_per_entity; i++ )
    {
        defInt[i]    = 0;
        defDouble[i] = -1e+10;
        defHandle[i] = static_cast< EntityHandle >( 0 );
    }

    switch( *tag_type )
    {
        case 0:
            tagDataType = MB_TYPE_INTEGER;
            tagType     = MB_TAG_DENSE;
            defaultVal  = defInt;
            break;

        case 1:
            tagDataType = MB_TYPE_DOUBLE;
            tagType     = MB_TAG_DENSE;
            defaultVal  = defDouble;
            break;

        case 2:
            tagDataType = MB_TYPE_HANDLE;
            tagType     = MB_TAG_DENSE;
            defaultVal  = defHandle;
            break;

        case 3:
            tagDataType = MB_TYPE_INTEGER;
            tagType     = MB_TAG_SPARSE;
            defaultVal  = defInt;
            break;

        case 4:
            tagDataType = MB_TYPE_DOUBLE;
            tagType     = MB_TAG_SPARSE;
            defaultVal  = defDouble;
            break;

        case 5:
            tagDataType = MB_TYPE_HANDLE;
            tagType     = MB_TAG_SPARSE;
            defaultVal  = defHandle;
            break;

        default: {
            delete[] defInt;
            delete[] defDouble;
            delete[] defHandle;
            return moab::MB_FAILURE;
        }  // error
    }

    // split storage names if separated list
    std::string tag_name( tag_storage_name );
    // first separate the names of the tags
    // we assume that there are separators ":" between the tag names
    std::vector< std::string > tagNames;
    std::string separator( ":" );
    split_tag_names( tag_name, separator, tagNames );

    ErrorCode rval           = moab::MB_SUCCESS;  // assume success already :)
    appData& data            = context.appDatas[*pid];
    int already_defined_tags = 0;

    Tag tagHandle;
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        rval = context.MBI->tag_get_handle( tagNames[i].c_str(), *components_per_entity, tagDataType, tagHandle,
                                            tagType | MB_TAG_EXCL | MB_TAG_CREAT, defaultVal );

        if( MB_ALREADY_ALLOCATED == rval )
        {
            std::map< std::string, Tag >& mTags        = data.tagMap;
            std::map< std::string, Tag >::iterator mit = mTags.find( tagNames[i].c_str() );

            if( mit == mTags.end() )
            {
                // add it to the map
                mTags[tagNames[i]] = tagHandle;
                // push it to the list of tags, too
                *tag_index = (int)data.tagList.size();
                data.tagList.push_back( tagHandle );
            }
            rval = MB_SUCCESS;
            already_defined_tags++;
        }
        else if( MB_SUCCESS == rval )
        {
            data.tagMap[tagNames[i]] = tagHandle;
            *tag_index               = (int)data.tagList.size();
            data.tagList.push_back( tagHandle );
        }
        else
        {
            rval = moab::MB_FAILURE;  // some tags were not created
        }
    }
    // we don't need default values anymore, avoid leaks
    int rankHere = 0;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco = context.appDatas[*pid].pcomm;
    rankHere          = pco->rank();
#endif
    if( !rankHere && already_defined_tags )
        std::cout << " application with ID: " << *pid << " global id: " << data.global_id << " name: " << data.name
                  << " has " << already_defined_tags << " already defined tags out of " << tagNames.size()
                  << " tags \n";
    delete[] defInt;
    delete[] defDouble;
    delete[] defHandle;
    return rval;
}

ErrCode iMOAB_SetIntTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* num_tag_storage_length,
                                int* ent_type,
                                int* tag_storage_data )
{
    std::string tag_name( tag_storage_name );

    // Get the application data
    appData& data = context.appDatas[*pid];
    // check if tag is defined
    if( data.tagMap.find( tag_name ) == data.tagMap.end() ) return moab::MB_FAILURE;

    Tag tag = data.tagMap[tag_name];

    int tagLength = 0;
    MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLength ) );

    DataType dtype;
    MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );

    if( dtype != MB_TYPE_INTEGER )
    {
        MB_CHK_SET_ERR( moab::MB_FAILURE, "The tag is not of integer type." );
    }

    // set it on a subset of entities, based on type and length
    // if *entity_type = 0, then use vertices; else elements
    Range* ents_to_set  = ( *ent_type == 0 ? &data.all_verts : &data.primary_elems );
    int nents_to_be_set = *num_tag_storage_length / tagLength;

    if( nents_to_be_set > (int)ents_to_set->size() )
    {
        return moab::MB_FAILURE;
    }  // to many entities to be set or too few

    // now set the tag data
    MB_CHK_ERR( context.MBI->tag_set_data( tag, *ents_to_set, tag_storage_data ) );

    return moab::MB_SUCCESS;  // no error
}

ErrCode iMOAB_GetIntTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* num_tag_storage_length,
                                int* ent_type,
                                int* tag_storage_data )
{
    std::string tag_name( tag_storage_name );

    appData& data = context.appDatas[*pid];

    if( data.tagMap.find( tag_name ) == data.tagMap.end() )
    {
        return moab::MB_FAILURE;
    }  // tag not defined

    Tag tag = data.tagMap[tag_name];

    int tagLength = 0;
    MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLength ) );

    DataType dtype;
    MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );

    if( dtype != MB_TYPE_INTEGER )
    {
        MB_CHK_SET_ERR( moab::MB_FAILURE, "The tag is not of integer type." );
    }

    // set it on a subset of entities, based on type and length
    // if *entity_type = 0, then use vertices; else elements
    Range* ents_to_get = ( *ent_type == 0 ? &data.all_verts : &data.primary_elems );
    int nents_to_get   = *num_tag_storage_length / tagLength;

    if( nents_to_get > (int)ents_to_get->size() )
    {
        return moab::MB_FAILURE;
    }  // to many entities to get, or too little

    // now set the tag data
    MB_CHK_ERR( context.MBI->tag_get_data( tag, *ents_to_get, tag_storage_data ) );

    return moab::MB_SUCCESS;  // no error
}

ErrCode iMOAB_SetDoubleTagStorage( iMOAB_AppID pid,
                                   const iMOAB_String tag_storage_names,
                                   int* num_tag_storage_length,
                                   int* ent_type,
                                   double* tag_storage_data )
{
    std::string tag_names( tag_storage_names );
    // exactly the same code as for int tag :) maybe should check the type of tag too
    std::vector< std::string > tagNames;
    std::vector< Tag > tagHandles;
    std::string separator( ":" );
    split_tag_names( tag_names, separator, tagNames );

    appData& data = context.appDatas[*pid];
    // set it on a subset of entities, based on type and length
    // if *entity_type = 0, then use vertices; else elements
    Range* ents_to_set = ( *ent_type == 0 ? &data.all_verts : &data.primary_elems );

    int nents_to_be_set = (int)( *ents_to_set ).size();
    int position        = 0;
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        if( data.tagMap.find( tagNames[i] ) == data.tagMap.end() )
        {
            return moab::MB_FAILURE;
        }  // some tag not defined yet in the app

        Tag tag = data.tagMap[tagNames[i]];

        int tagLength = 0;
        MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLength ) );

        DataType dtype;
        MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );

        if( dtype != MB_TYPE_DOUBLE )
        {
            return moab::MB_FAILURE;
        }

        // set it on the subset of entities, based on type and length
        if( position + tagLength * nents_to_be_set > *num_tag_storage_length )
            return moab::MB_FAILURE;  // too many entity values to be set

        MB_CHK_ERR( context.MBI->tag_set_data( tag, *ents_to_set, &tag_storage_data[position] ) );
        // increment position to next tag
        position = position + tagLength * nents_to_be_set;
    }
    return moab::MB_SUCCESS;  // no error
}

ErrCode iMOAB_SetDoubleTagStorageWithGid( iMOAB_AppID pid,
                                          const iMOAB_String tag_storage_names,
                                          int* num_tag_storage_length,
                                          int* ent_type,
                                          double* tag_storage_data,
                                          int* globalIds )
{
    std::string tag_names( tag_storage_names );
    // exactly the same code as for int tag :) maybe should check the type of tag too
    std::vector< std::string > tagNames;
    std::vector< Tag > tagHandles;
    std::string separator( ":" );
    split_tag_names( tag_names, separator, tagNames );

    appData& data = context.appDatas[*pid];
    // set it on a subset of entities, based on type and length
    // if *entity_type = 0, then use vertices; else elements
    Range* ents_to_set  = ( *ent_type == 0 ? &data.all_verts : &data.primary_elems );
    int nents_to_be_set = (int)( *ents_to_set ).size();

    Tag gidTag = context.MBI->globalId_tag();
    std::vector< int > gids;
    gids.resize( nents_to_be_set );
    MB_CHK_ERR( context.MBI->tag_get_data( gidTag, *ents_to_set, &gids[0] ) );

    // so we will need to set the tags according to the global id passed;
    // so the order in tag_storage_data is the same as the order in globalIds, but the order
    // in local range is gids
    std::map< int, EntityHandle > eh_by_gid;
    int i = 0;
    for( Range::iterator it = ents_to_set->begin(); it != ents_to_set->end(); ++it, ++i )
    {
        eh_by_gid[gids[i]] = *it;
    }

    /// @todo Allow for tags of different length
    size_t nbLocalVals = *num_tag_storage_length / tagNames.size();  // assumes all tags have the same length?

    // check global ids to have different values
    std::set< int > globalIdsSet( globalIds, globalIds + nbLocalVals );
    if( globalIdsSet.size() < nbLocalVals )
    {
        std::cout << "iMOAB_SetDoubleTagStorageWithGid: for pid:" << *pid << " tags[0]:" << tagNames[0]
                  << " global ids passed are not unique, major error\n";
        std::cout << " nbLocalVals:" << nbLocalVals << " globalIdsSet.size():" << globalIdsSet.size()
                  << " first global id:" << globalIds[0] << "\n";
        return moab::MB_FAILURE;
    }

    std::vector< int > tagLengths( tagNames.size() );
    std::vector< Tag > tagList;
#ifdef MOAB_HAVE_MPI
    size_t total_tag_len = 0;
#endif
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        if( data.tagMap.find( tagNames[i] ) == data.tagMap.end() )
        {
            MB_SET_ERR( moab::MB_FAILURE, "tag missing" );
        }  // some tag not defined yet in the app

        Tag tag = data.tagMap[tagNames[i]];
        tagList.push_back( tag );

        int tagLength = 0;
        MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLength ) );

#ifdef MOAB_HAVE_MPI
        total_tag_len += tagLength;
#endif
        tagLengths[i] = tagLength;
        DataType dtype;
        MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );

        if( dtype != MB_TYPE_DOUBLE )
        {
            MB_SET_ERR( moab::MB_FAILURE, "tag not double type" );
        }
    }
    bool serial = true;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco  = context.appDatas[*pid].pcomm;
    unsigned num_procs = pco->size();
    if( num_procs > 1 ) serial = false;
#endif

    if( serial )
    {
        // we do not assume anymore that the number of entities has to match
        // we will set only what matches, and skip entities that do not have corresponding global ids
        //assert( total_tag_len * nents_to_be_set - *num_tag_storage_length == 0 );
        // tags are unrolled, we loop over global ids first, then careful about tags
        for( int i = 0; i < nents_to_be_set; i++ )
        {
            int gid                                       = globalIds[i];
            std::map< int, EntityHandle >::iterator mapIt = eh_by_gid.find( gid );
            if( mapIt == eh_by_gid.end() ) continue;
            EntityHandle eh = mapIt->second;
            // now loop over tags
            int indexInTagValues = 0;  //
            for( size_t j = 0; j < tagList.size(); j++ )
            {
                indexInTagValues += i * tagLengths[j];
                MB_CHK_ERR( context.MBI->tag_set_data( tagList[j], &eh, 1, &tag_storage_data[indexInTagValues] ) );
                // advance the pointer/index
                indexInTagValues += ( nents_to_be_set - i ) * tagLengths[j];  // at the end of tag data
            }
        }
    }
#ifdef MOAB_HAVE_MPI
    else  // it can be not serial only if pco->size() > 1, parallel
    {
        // in this case, we have to use 2 crystal routers, to send data to the processor that needs it
        // we will create first a tuple to rendevous points, then from there send to the processor that requested it
        // it is a 2-hop global gather scatter
        // we do not expect the sizes to match
        //assert( nbLocalVals * tagNames.size() - *num_tag_storage_length == 0 );
        TupleList TLsend;
        TLsend.initialize( 2, 0, 0, total_tag_len, nbLocalVals );  //  to proc, marker(gid), total_tag_len doubles
        TLsend.enableWriteAccess();
        // the processor id that processes global_id is global_id / num_ents_per_proc

        int indexInRealLocal = 0;
        for( size_t i = 0; i < nbLocalVals; i++ )
        {
            // to proc, marker, element local index, index in el
            int marker              = globalIds[i];
            int to_proc             = marker % num_procs;
            int n                   = TLsend.get_n();
            TLsend.vi_wr[2 * n]     = to_proc;  // send to processor
            TLsend.vi_wr[2 * n + 1] = marker;
            int indexInTagValues    = 0;
            // tag data collect by number of tags
            for( size_t j = 0; j < tagList.size(); j++ )
            {
                indexInTagValues += i * tagLengths[j];
                for( int k = 0; k < tagLengths[j]; k++ )
                {
                    TLsend.vr_wr[indexInRealLocal++] = tag_storage_data[indexInTagValues + k];
                }
                indexInTagValues += ( nbLocalVals - i ) * tagLengths[j];
            }
            TLsend.inc_n();
        }

        //assert( nbLocalVals * total_tag_len - indexInRealLocal == 0 );
        // send now requests, basically inform the rendez-vous point who needs a particular global id
        // send the data to the other processors:
        ( pco->proc_config().crystal_router() )->gs_transfer( 1, TLsend, 0 );
        TupleList TLreq;
        TLreq.initialize( 2, 0, 0, 0, nents_to_be_set );
        TLreq.enableWriteAccess();
        for( int i = 0; i < nents_to_be_set; i++ )
        {
            // to proc, marker
            int marker             = gids[i];
            int to_proc            = marker % num_procs;
            int n                  = TLreq.get_n();
            TLreq.vi_wr[2 * n]     = to_proc;  // send to processor
            TLreq.vi_wr[2 * n + 1] = marker;
            // tag data collect by number of tags
            TLreq.inc_n();
        }

        // perform communication with crystal router
        pco->proc_config().crystal_router()->gs_transfer( 1, TLreq, 0 );

        // we know now that process TLreq.vi_wr[2 * n] needs tags for gid TLreq.vi_wr[2 * n + 1]
        // we should first order by global id, and then build the new TL with send to proc, global id and
        // tags for it
        // sort by global ids the tuple lists
        moab::TupleList::buffer sort_buffer;
        sort_buffer.buffer_init( TLreq.get_n() );
        TLreq.sort( 1, &sort_buffer );
        sort_buffer.reset();
        sort_buffer.buffer_init( TLsend.get_n() );
        TLsend.sort( 1, &sort_buffer );
        sort_buffer.reset();
        // now send the tag values to the proc that requested it
        // in theory, for a full  partition, TLreq  and TLsend should have the same size, and
        // each dof should have exactly one target proc. Is that true or not in general ?
        // how do we plan to use this? Is it better to store the comm graph for future

        // start copy from comm graph settle
        TupleList TLBack;
        TLBack.initialize( 3, 0, 0, total_tag_len, 0 );  // to proc, marker, tag from proc , tag values
        TLBack.enableWriteAccess();

        int n1 = TLreq.get_n();
        int n2 = TLsend.get_n();

        int indexInTLreq  = 0;
        int indexInTLsend = 0;  // advance both, according to the marker
        if( n1 > 0 && n2 > 0 )
        {

            while( indexInTLreq < n1 && indexInTLsend < n2 )  // if any is over, we are done
            {
                int currentValue1 = TLreq.vi_rd[2 * indexInTLreq + 1];
                int currentValue2 = TLsend.vi_rd[2 * indexInTLsend + 1];
                if( currentValue1 < currentValue2 )
                {
                    // we have a big problem; basically, we are saying that
                    // dof currentValue is on one model and not on the other
                    indexInTLreq++;
                    continue;
                }

                if( currentValue1 > currentValue2 )
                {
                    indexInTLsend++;
                    continue;
                }

                int size1 = 1;
                while( indexInTLreq + size1 < n1 && currentValue1 == TLreq.vi_rd[2 * ( indexInTLreq + size1 ) + 1] )
                    size1++;
                int size2 = 1;
                while( indexInTLsend + size2 < n2 && currentValue2 == TLsend.vi_rd[2 * ( indexInTLsend + size2 ) + 1] )
                    size2++;

                // must be found in both lists, find the start and end indices
                for( int i1 = 0; i1 < size1; i1++ )
                {
                    for( int i2 = 0; i2 < size2; i2++ )
                    {
                        // send the info back to components
                        int n = TLBack.get_n();
                        TLBack.reserve();
                        TLBack.vi_wr[3 * n] = TLreq.vi_rd[2 * ( indexInTLreq + i1 )];  // send back to the proc marker
                                                                                       // came from, info from comp2
                        TLBack.vi_wr[3 * n + 1] = currentValue1;  // initial value (resend, just for verif ?)
                        TLBack.vi_wr[3 * n + 2] = TLsend.vi_rd[2 * ( indexInTLsend + i2 )];  // from proc on comp2
                        // also fill tag values
                        for( size_t k = 0; k < total_tag_len; k++ )
                        {
                            TLBack.vr_rd[total_tag_len * n + k] =
                                TLsend.vr_rd[total_tag_len * indexInTLsend + k];  // deep copy of tag values
                        }
                    }
                }
                indexInTLreq += size1;
                indexInTLsend += size2;
            }
        }

        // invoke crystal router
        pco->proc_config().crystal_router()->gs_transfer( 1, TLBack, 0 );

        // end copy from comm graph
        // after we are done sending, we need to set those tag values, in a reverse process compared to send
        n1             = TLBack.get_n();
        double* ptrVal = &TLBack.vr_rd[0];  //
        for( int i = 0; i < n1; i++ )
        {
            const int gid    = TLBack.vi_rd[3 * i + 1];  // marker
            const auto mapIt = eh_by_gid.find( gid );
            if( mapIt == eh_by_gid.end() ) continue;

            EntityHandle eh = mapIt->second;
            // now loop over tags
            for( size_t j = 0; j < tagList.size(); j++ )
            {
                MB_CHK_ERR( context.MBI->tag_set_data( tagList[j], &eh, 1, (void*)ptrVal ) );
                // advance the pointer/index
                ptrVal += tagLengths[j];  // at the end of tag data per call
            }
        }
    }
#endif
    return MB_SUCCESS;
}

#ifdef MOAB_HAVE_TEMPESTREMAP
// Helper: get the 2D entities of the TempestRemap CoveringMesh attached to this app.
// Returns an empty range (and success) if no remapper is attached.
static ErrCode get_coverage_entities( iMOAB_AppID pid, moab::Range& covEnts )
{
    appData& data = context.appDatas[*pid];
    covEnts.clear();
    if( data.tempestData.remapper == nullptr ) return moab::MB_SUCCESS;
    covEnts = data.tempestData.remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
    return moab::MB_SUCCESS;
}

/**
 * \brief Return the number of 2D elements in the TempestRemap CoveringMesh of this app
 *        and (optionally) their global IDs and centroid coordinates.
 *
 * Used to populate analytic source fields directly on the coverage entities of a
 * dual-map intersection app, bypassing the cmpAtm->cplAtm->cplDualMap migration
 * chain — useful for BFB testing of dual-map CAAS where every rank must see an
 * identical input set independent of how the source mesh was partitioned.
 *
 * Call once with gids=NULL, centroids=NULL to query num_cov_elems, then allocate
 * and call again to fill the buffers.
 *
 * \param[in]    pid             Application ID with an attached TempestRemap remapper
 *                               (e.g. the dual-map intersection app).
 * \param[inout] num_cov_elems   On input: ignored (or capacity); on output: number of
 *                               2D entities in the covering mesh on this rank.
 * \param[out]   gids            Optional. If non-null, filled with the global IDs of
 *                               the covering entities (size num_cov_elems).
 * \param[out]   centroids       Optional. If non-null, filled with the arithmetic mean
 *                               of vertex coords per element (size 3*num_cov_elems,
 *                               x,y,z interleaved).
 * \return ErrCode               MB_SUCCESS on success.
 */
ErrCode iMOAB_GetCoverageMeshInfo( iMOAB_AppID pid, int* num_cov_elems, int* gids, double* centroids )
{
    moab::Range covEnts;
    MB_CHK_ERR( get_coverage_entities( pid, covEnts ) );
    *num_cov_elems = static_cast< int >( covEnts.size() );
    if( covEnts.empty() ) return moab::MB_SUCCESS;

    if( gids != nullptr )
    {
        std::vector< int > tmpGids( covEnts.size() );
        MB_CHK_ERR( context.MBI->tag_get_data( context.globalID_tag, covEnts, tmpGids.data() ) );
        std::copy( tmpGids.begin(), tmpGids.end(), gids );
    }

    if( centroids != nullptr )
    {
        const moab::EntityHandle* conn;
        int numNodes;
        std::vector< double > coords;
        int i = 0;
        for( moab::Range::iterator it = covEnts.begin(); it != covEnts.end(); ++it, ++i )
        {
            MB_CHK_ERR( context.MBI->get_connectivity( *it, conn, numNodes ) );
            coords.resize( 3 * numNodes );
            MB_CHK_ERR( context.MBI->get_coords( conn, numNodes, coords.data() ) );
            double cx = 0.0, cy = 0.0, cz = 0.0;
            for( int v = 0; v < numNodes; v++ )
            {
                cx += coords[3 * v + 0];
                cy += coords[3 * v + 1];
                cz += coords[3 * v + 2];
            }
            centroids[3 * i + 0] = cx / numNodes;
            centroids[3 * i + 1] = cy / numNodes;
            centroids[3 * i + 2] = cz / numNodes;
        }
    }
    return moab::MB_SUCCESS;
}

/**
 * \brief Set a double tag on every entity of the TempestRemap CoveringMesh of this app.
 *
 * Values are written in the same order as iMOAB_GetCoverageMeshInfo returns gids
 * (i.e. moab::Range iteration order over the covering set entities). The tag must
 * already be defined on the app via iMOAB_DefineTagStorage.
 *
 * This bypasses the partition-dependent migration path and is useful for BFB tests
 * that must populate the source field on the actual coverage cells consumed by
 * iMOAB_ApplyScalarProjectionWeights.
 *
 * \param[in] pid                       Application ID with an attached TempestRemap remapper.
 * \param[in] tag_storage_name          Name of the (already-defined) double tag.
 * \param[in] num_tag_storage_length    Total number of values supplied (= num_cov_elems *
 *                                      components_per_entity).
 * \param[in] tag_storage_data          The values to write, in covering-Range order.
 * \return ErrCode                      MB_SUCCESS on success.
 */
ErrCode iMOAB_SetDoubleTagStorageOnCoverage( iMOAB_AppID pid,
                                             const iMOAB_String tag_storage_name,
                                             int* num_tag_storage_length,
                                             double* tag_storage_data )
{
    moab::Range covEnts;
    MB_CHK_ERR( get_coverage_entities( pid, covEnts ) );
    if( covEnts.empty() ) return moab::MB_SUCCESS;

    std::string tagName( tag_storage_name );
    moab::Tag tag = nullptr;
    if( moab::MB_SUCCESS != context.MBI->tag_get_handle( tagName.c_str(), tag ) || !tag ) return moab::MB_FAILURE;

    int tagLen = 0;
    MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLen ) );
    moab::DataType dtype;
    MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );
    if( dtype != moab::MB_TYPE_DOUBLE ) return moab::MB_FAILURE;

    const int needed = tagLen * static_cast< int >( covEnts.size() );
    if( *num_tag_storage_length < needed ) return moab::MB_FAILURE;

    MB_CHK_ERR( context.MBI->tag_set_data( tag, covEnts, tag_storage_data ) );
    return moab::MB_SUCCESS;
}
#endif  // MOAB_HAVE_TEMPESTREMAP

ErrCode iMOAB_GetDoubleTagStorage( iMOAB_AppID pid,
                                   const iMOAB_String tag_storage_names,
                                   int* num_tag_storage_length,
                                   int* ent_type,
                                   double* tag_storage_data )
{
    // exactly the same code, except tag type check
    std::string tag_names( tag_storage_names );
    // exactly the same code as for int tag :) maybe should check the type of tag too
    std::vector< std::string > tagNames;
    std::vector< Tag > tagHandles;
    std::string separator( ":" );
    split_tag_names( tag_names, separator, tagNames );

    appData& data = context.appDatas[*pid];

    // set it on a subset of entities, based on type and length
    Range* ents_to_get = nullptr;

    if( *ent_type == 0 )  // vertices
    {
        ents_to_get = &data.all_verts;
    }
    else if( *ent_type == 1 )
    {
        ents_to_get = &data.primary_elems;
    }
    int nents_to_get = (int)ents_to_get->size();
    int position     = 0;
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        if( data.tagMap.find( tagNames[i] ) == data.tagMap.end() )
        {
            return moab::MB_FAILURE;
        }  // tag not defined

        Tag tag = data.tagMap[tagNames[i]];

        int tagLength = 0;
        MB_CHK_ERR( context.MBI->tag_get_length( tag, tagLength ) );

        DataType dtype;
        MB_CHK_ERR( context.MBI->tag_get_data_type( tag, dtype ) );

        if( dtype != MB_TYPE_DOUBLE )
        {
            return moab::MB_FAILURE;
        }

        if( position + nents_to_get * tagLength > *num_tag_storage_length )
            return moab::MB_FAILURE;  // too many entity values to get

        MB_CHK_ERR( context.MBI->tag_get_data( tag, *ents_to_get, &tag_storage_data[position] ) );
        position = position + nents_to_get * tagLength;
    }

    return moab::MB_SUCCESS;  // no error
}

ErrCode iMOAB_SynchronizeTags( iMOAB_AppID pid, int* num_tag, int* tag_indices, int* ent_type )
{
#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    Range ent_exchange;
    std::vector< Tag > tags;

    for( int i = 0; i < *num_tag; i++ )
    {
        if( tag_indices[i] < 0 || tag_indices[i] >= (int)data.tagList.size() )
        {
            return moab::MB_FAILURE;
        }  // error in tag index

        tags.push_back( data.tagList[tag_indices[i]] );
    }

    if( *ent_type == 0 )
    {
        ent_exchange = data.all_verts;
    }
    else if( *ent_type == 1 )
    {
        ent_exchange = data.primary_elems;
    }
    else
    {
        return moab::MB_FAILURE;
    }  // unexpected type

    ParallelComm* pco = context.appDatas[*pid].pcomm;

    MB_CHK_ERR( pco->exchange_tags( tags, tags, ent_exchange ) );

#else
    /* do nothing if serial */
    UNUSED( pid );
    UNUSED( num_tag );
    UNUSED( tag_indices );
    UNUSED( ent_type );
#endif

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_ReduceTagsMax( iMOAB_AppID pid, int* tag_index, int* ent_type )
{
#ifdef MOAB_HAVE_MPI
    appData& data = context.appDatas[*pid];
    Range ent_exchange;

    if( *tag_index < 0 || *tag_index >= (int)data.tagList.size() )
    {
        return moab::MB_FAILURE;
    }  // error in tag index

    Tag tagh = data.tagList[*tag_index];

    if( *ent_type == 0 )
    {
        ent_exchange = data.all_verts;
    }
    else if( *ent_type == 1 )
    {
        ent_exchange = data.primary_elems;
    }
    else
    {
        return moab::MB_FAILURE;
    }  // unexpected type

    ParallelComm* pco = context.appDatas[*pid].pcomm;
    // we could do different MPI_Op; do not bother now, we will call from fortran
    MB_CHK_ERR( pco->reduce_tags( tagh, MPI_MAX, ent_exchange ) );

#else
    /* do nothing if serial */
    UNUSED( pid );
    UNUSED( tag_index );
    UNUSED( ent_type );
#endif
    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetNeighborElements( iMOAB_AppID pid,
                                   iMOAB_LocalID* local_index,
                                   int* num_adjacent_elements,
                                   iMOAB_LocalID* adjacent_element_IDs )
{
    assert( local_index && *local_index >= 0 );

    // one neighbor for each subentity of dimension-1
    MeshTopoUtil mtu( context.MBI );
    appData& data   = context.appDatas[*pid];
    EntityHandle eh = data.primary_elems[*local_index];
    Range adjs;
    MB_CHK_SET_ERR( mtu.get_bridge_adjacencies( eh, data.dimension - 1, data.dimension, adjs ),
                    "Getting bridge adjacencies failed" );

    if( *num_adjacent_elements < (int)adjs.size() )
    {
        return moab::MB_FAILURE;
    }  // not dimensioned correctly

    *num_adjacent_elements = (int)adjs.size();

    for( int i = 0; i < *num_adjacent_elements; i++ )
    {
        adjacent_element_IDs[i] = data.primary_elems.index( adjs[i] );
    }

    return moab::MB_SUCCESS;
}

#if 0

ErrCode iMOAB_GetNeighborVertices ( iMOAB_AppID pid, iMOAB_LocalID* local_vertex_ID, int* num_adjacent_vertices, iMOAB_LocalID* adjacent_vertex_IDs )
{
    return moab::MB_SUCCESS;
}

#endif

ErrCode iMOAB_CreateVertices( iMOAB_AppID pid, int* coords_len, int* dim, double* coordinates )
{
    appData& data = context.appDatas[*pid];

    if( !data.local_verts.empty() )  // we should have no vertices in the app
    {
        return moab::MB_FAILURE;
    }

    int nverts = *coords_len / *dim;

    MB_CHK_ERR( context.MBI->create_vertices( coordinates, nverts, data.local_verts ) );

    MB_CHK_ERR( context.MBI->add_entities( data.file_set, data.local_verts ) );

    // also add the vertices to the all_verts range
    data.all_verts.merge( data.local_verts );
    return moab::MB_SUCCESS;
}

ErrCode iMOAB_CreateElements( iMOAB_AppID pid,
                              int* num_elem,
                              int* type,
                              int* num_nodes_per_element,
                              int* connectivity,
                              int* block_ID )
{
    // Create elements
    appData& data = context.appDatas[*pid];

    ReadUtilIface* read_iface;
    MB_CHK_ERR( context.MBI->query_interface( read_iface ) );

    EntityType mbtype = (EntityType)( *type );
    EntityHandle actual_start_handle;
    EntityHandle* array = nullptr;
    MB_CHK_ERR(
        read_iface->get_element_connect( *num_elem, *num_nodes_per_element, mbtype, 1, actual_start_handle, array ) );

    // fill up with actual connectivity from input; assume the vertices are in order, and start
    // vertex is the first in the current data vertex range
    EntityHandle firstVertex = data.local_verts[0];

    for( int j = 0; j < *num_elem * ( *num_nodes_per_element ); j++ )
    {
        array[j] = connectivity[j] + firstVertex - 1;
    }  // assumes connectivity uses 1 based array (from fortran, mostly)

    Range new_elems( actual_start_handle, actual_start_handle + *num_elem - 1 );

    MB_CHK_ERR( context.MBI->add_entities( data.file_set, new_elems ) );

    data.primary_elems.merge( new_elems );

    // add to adjacency
    MB_CHK_ERR( read_iface->update_adjacencies( actual_start_handle, *num_elem, *num_nodes_per_element, array ) );

    // organize all new elements in block, with the given block ID; if the block set is not
    // existing, create  a new mesh set;
    Range sets;
    int set_no            = *block_ID;
    const void* setno_ptr = &set_no;

    EntityHandle block_set;
    ErrorCode rval = context.MBI->get_entities_by_type_and_tag( data.file_set, MBENTITYSET, &context.material_tag,
                                                                &setno_ptr, 1, sets );

    if( MB_FAILURE == rval || sets.empty() )
    {
        // create a new set, with this block ID
        MB_CHK_ERR( context.MBI->create_meshset( MESHSET_SET, block_set ) );

        MB_CHK_ERR( context.MBI->tag_set_data( context.material_tag, &block_set, 1, &set_no ) );

        // add the material set to file set
        MB_CHK_ERR( context.MBI->add_entities( data.file_set, &block_set, 1 ) );
    }
    else
    {
        block_set = sets[0];
    }  // first set is the one we want

    /// add the new ents to the clock set
    MB_CHK_ERR( context.MBI->add_entities( block_set, new_elems ) );

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_SetGlobalInfo( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems )
{
    appData& data            = context.appDatas[*pid];
    data.num_global_vertices = *num_global_verts;
    data.num_global_elements = *num_global_elems;
    return moab::MB_SUCCESS;
}

ErrCode iMOAB_GetGlobalInfo( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems )
{
    appData& data = context.appDatas[*pid];
    if( nullptr != num_global_verts )
    {
        *num_global_verts = data.num_global_vertices;
    }
    if( nullptr != num_global_elems )
    {
        *num_global_elems = data.num_global_elements;
    }

    return moab::MB_SUCCESS;
}

#ifdef MOAB_HAVE_MPI

// this makes sense only for parallel runs
ErrCode iMOAB_ResolveSharedEntities( iMOAB_AppID pid, int* num_verts, int* marker )
{
    appData& data     = context.appDatas[*pid];
    ParallelComm* pco = context.appDatas[*pid].pcomm;
    EntityHandle cset = data.file_set;
    int dum_id        = 0;
    ErrorCode rval;
    if( data.primary_elems.empty() )
    {
        // skip actual resolve, assume vertices are distributed already ,
        // no need to share them
    }
    else
    {
        // create an integer tag for resolving ; maybe it can be a long tag in the future
        // (more than 2 B vertices;)

        Tag stag;
        MB_CHK_ERR( context.MBI->tag_get_handle( "__sharedmarker", 1, MB_TYPE_INTEGER, stag,
                                                 MB_TAG_CREAT | MB_TAG_DENSE, &dum_id ) );

        if( *num_verts > (int)data.local_verts.size() )
        {
            return moab::MB_FAILURE;
        }  // we are not setting the size

        MB_CHK_ERR( context.MBI->tag_set_data( stag, data.local_verts, (void*)marker ) );  // assumes integer tag

        MB_CHK_ERR( pco->resolve_shared_ents( cset, -1, -1, &stag ) );

        MB_CHK_ERR( context.MBI->tag_delete( stag ) );
    }
    // provide partition tag equal to rank
    Tag part_tag;
    dum_id = -1;
    rval   = context.MBI->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_tag,
                                          MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );

    if( part_tag == nullptr || ( ( rval != MB_SUCCESS ) && ( rval != MB_ALREADY_ALLOCATED ) ) )
    {
        std::cout << " can't get par part tag.\n";
        return moab::MB_FAILURE;
    }

    int rank = pco->rank();
    MB_CHK_ERR( context.MBI->tag_set_data( part_tag, &cset, 1, &rank ) );

    return moab::MB_SUCCESS;
}

// this assumes that this was not called before
ErrCode iMOAB_DetermineGhostEntities( iMOAB_AppID pid, int* ghost_dim, int* num_ghost_layers, int* bridge_dim )
{
    // verify we have valid ghost layers input specified. If invalid, exit quick.
    if( num_ghost_layers && *num_ghost_layers <= 0 )
    {
        return moab::MB_SUCCESS;
    }  // nothing to do

    appData& data     = context.appDatas[*pid];
    ParallelComm* pco = context.appDatas[*pid].pcomm;

    // maybe we should be passing this too; most of the time we do not need additional ents collective call
    constexpr int addl_ents = 0;
    MB_CHK_ERR( pco->exchange_ghost_cells( *ghost_dim, *bridge_dim, 1,  // get only one layer of ghost entities
                                           addl_ents, true, true, &data.file_set ) );
    for( int i = 2; i <= *num_ghost_layers; i++ )
    {
        MB_CHK_ERR( pco->correct_thin_ghost_layers() );                     // correct for thin layers
        MB_CHK_ERR( pco->exchange_ghost_cells( *ghost_dim, *bridge_dim, i,  // iteratively get one extra layer
                                               addl_ents, true, true, &data.file_set ) );
    }

    // Update ghost layer information
    data.num_ghost_layers = *num_ghost_layers;

    // now re-establish all mesh info; will reconstruct mesh info, based solely on what is in the file set
    return iMOAB_UpdateMeshInfo( pid );
}

/**
 * @brief Transfer mesh data from sender component to receiver component.
 * @ingroup iMOABParallel
 *
 * @details Implements distributed mesh transfer protocol that moves mesh entities (elements, vertices,
 * and associated data) from one set of MPI processes to another. Transfer is orchestrated through a
 * communication graph that maps sender processes to receiver processes based on partitioning strategy.
 *
 * @par Algorithm Workflow:
 * -# <b>COMMUNICATOR SETUP:</b> Resolve MPI communicator and group handles (handles Fortran F2C conversions).
 *    Extract sender group from MOAB ParallelComm. Create ParCommGraph for sender→receiver mapping.
 * -# <b>ENTITY SELECTION:</b> Determine owned entities to transfer: owned_elems or local_verts (point cloud mode).
 *    Store ParCommGraph in context.appDatas[pid].pgraph[rcompid].
 * -# <b>PARTITIONING STRATEGY:</b>
 *    - method == 0 (TRIVIAL): Gather global entity counts, compute block distribution, broadcast via send_graph
 *    - method != 0 (COMPUTED): Use graph/geometric partitioning considering connectivity and spatial distribution
 * -# <b>MESH TRANSFER:</b> Pack mesh entities and data, use ParCommGraph to route to receivers via send_mesh_parts
 *
 * @par Data Transferred:
 * - <b>Mesh entities:</b> elements (triangles, quads, etc.) and vertices
 * - <b>Geometric data:</b> vertex coordinates, element connectivity
 * - <b>Metadata:</b> material sets, boundary conditions, tags
 * - <b>Global IDs:</b> for entity identification across processes
 *
 * @param[in] pid                Application ID for sender component
 * @param[in] joint_communicator MPI communicator spanning both sender and receiver groups
 * @param[in] receivingGroup     MPI group of receiving processes
 * @param[in] rcompid            Receiver component ID (key for communication graph storage)
 * @param[in] method             Partitioning method: 0=trivial, !=0=computed
 *
 * @note <b>Communication:</b> Collective (MPI_Allgather) and point-to-point (non-blocking sends)
 * @warning Cleans up MPI groups on early return to avoid resource leaks
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE on error
 */
ErrCode iMOAB_SendMesh( iMOAB_AppID pid,
                        MPI_Comm* joint_communicator,
                        MPI_Group* receivingGroup,
                        int* rcompid,
                        int* method )
{
    // Validate input parameters to ensure they are not null pointers
    assert( joint_communicator != nullptr );
    assert( receivingGroup != nullptr );
    assert( rcompid != nullptr );

    int ierr;
    // Get application data and parallel communicator for this component
    appData& data     = context.appDatas[*pid];
    ParallelComm* pco = context.appDatas[*pid].pcomm;

    // Handle Fortran-to-C MPI handle conversion if needed
    // Fortran passes MPI handles as integers, C++ expects MPI_Comm/MPI_Group objects
    MPI_Comm global = ( data.is_fortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) )
                                        : *joint_communicator );
    MPI_Group recvGroup =
        ( data.is_fortran ? MPI_Group_f2c( *reinterpret_cast< MPI_Fint* >( receivingGroup ) ) : *receivingGroup );

    // Extract sender communicator from MOAB's parallel communicator
    // This represents the group of processes that will send mesh data
    MPI_Comm sender = pco->comm();
    // Extract sender group from the sender communicator to understand sender process topology
    MPI_Group senderGroup;
    ierr = MPI_Comm_group( sender, &senderGroup );
    if( ierr != 0 ) return moab::MB_FAILURE;

    // Create communication graph to manage sender->receiver process mapping
    // This graph will determine which entities go to which receiver processes
    ParCommGraph* cgraph =
        new ParCommGraph( global, senderGroup, recvGroup, context.appDatas[*pid].global_id, *rcompid );

    // Store the communication graph in the application's graph map for later use
    // The key is the receiver component ID (*rcompid)
    context.appDatas[*pid].pgraph[*rcompid] = cgraph;

    // Get current sender rank for potential debugging/logging
    int sender_rank = -1;
    MPI_Comm_rank( sender, &sender_rank );

    // Determine which entities to send: primary elements or vertices (for point clouds)
    std::vector< int > number_elems_per_part;
    Range owned = context.appDatas[*pid].owned_elems;
    if( owned.size() == 0 )
    {
        // If no elements, this must be a point cloud - send vertices instead
        owned = context.appDatas[*pid].local_verts;
    }

    // Choose partitioning strategy based on method parameter
    if( *method == 0 )  // Trivial partitioning: simple block distribution
    {
        // Count local entities on this sender process
        int local_owned_elem = (int)owned.size();
        int size             = pco->size();
        int rank             = pco->rank();

        // Prepare array to hold entity counts from all sender processes
        number_elems_per_part.resize( size );
        number_elems_per_part[rank] = local_owned_elem;
        // Gather entity counts from all sender processes using MPI_Allgather
        // This allows each sender to know the total entity count for load balancing
#if ( MPI_VERSION >= 2 )
        // Use "in place" option for efficiency (MPI 2.0+)
        ierr = MPI_Allgather( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &number_elems_per_part[0], 1, MPI_INT, sender );
#else
        // Fallback for older MPI versions
        {
            std::vector< int > all_tmp( size );
            ierr = MPI_Allgather( &number_elems_per_part[rank], 1, MPI_INT, &all_tmp[0], 1, MPI_INT, sender );
            number_elems_per_part = all_tmp;
        }
#endif

        if( ierr != 0 )
        {
            return moab::MB_FAILURE;
        }

        // Compute trivial partition: each receiver gets roughly equal share of entities
        // This creates a simple block distribution based on entity counts
        MB_CHK_ERR( cgraph->compute_trivial_partition( number_elems_per_part ) );

        // Send the partition mapping to all receiver processes
        // This tells receivers which senders will send them data
        MB_CHK_ERR( cgraph->send_graph( global ) );
    }
    else  // Advanced partitioning: graph-based or geometric partitioning
    {
        // Use sophisticated partitioning that considers entity connectivity or spatial distribution
        // This can lead to better load balancing and reduced communication
        MB_CHK_ERR( cgraph->compute_partition( pco, owned, *method ) );

        // Send the computed partition mapping to receiver processes
        // This includes more complex sender->receiver mappings than trivial partition
        MB_CHK_ERR( cgraph->send_graph_partition( pco, global ) );
    }
    // Execute the actual mesh transfer based on the communication graph
    // This packs mesh entities and sends them to appropriate receiver processes
    // pco is needed for MOAB operations (packing), not for MPI communication
    MB_CHK_ERR( cgraph->send_mesh_parts( global, pco, owned ) );

    // Clean up temporary MPI group to prevent memory leaks
    MPI_Group_free( &senderGroup );
    return moab::MB_SUCCESS;
}

/**
 * @brief Receive mesh data on target component from sending component.
 * @ingroup iMOABParallel
 *
 * @details Implements receiver side of distributed mesh transfer protocol. Receives mesh entities
 * and data from multiple sender processes, consolidates them into local mesh, and handles entity
 * deduplication based on global IDs.
 *
 * @par Algorithm Workflow:
 * -# <b>COMMUNICATOR SETUP:</b> Resolve MPI communicator and group handles (Fortran F2C conversions).
 *    Extract receiver group from MOAB ParallelComm. Create ParCommGraph with reverse mapping.
 * -# <b>COMMUNICATION GRAPH RECEPTION:</b> Receive communication graph via receive_comm_graph.
 *    Parse packed graph array to extract sender ranks for this receiver process.
 * -# <b>SENDER IDENTIFICATION:</b> Determine current receiver's rank index. Walk through packed
 *    graph array to find all sender ranks that contribute. Build senders_local list.
 * -# <b>MESH RECEPTION:</b> Invoke receive_mesh to pull mesh parts from identified sender processes.
 *    Deposit mesh entities into data.file_set (local mesh container).
 * -# <b>POST-PROCESSING:</b> Optional vertex merging for entities with identical GLOBAL_ID from
 *    different senders. Re-establish mesh information and update local mesh metadata.
 *
 * @par Data Received:
 * - <b>Mesh entities:</b> elements and vertices with geometric data
 * - <b>Connectivity:</b> element-to-vertex relationships
 * - <b>Tags:</b> material properties, boundary conditions, custom data
 * - <b>Global IDs:</b> for entity identification and deduplication
 *
 * @param[in] pid                Application ID for receiver component
 * @param[in] joint_communicator MPI communicator spanning both sender and receiver groups
 * @param[in] sendingGroup       MPI group of sending processes
 * @param[in] scompid            Sender component ID (key for communication graph storage)
 *
 * @note <b>Communication:</b> Collective (receive graph) and point-to-point (receive mesh data)
 * @warning Vertex merging only performed if receiving from 2+ senders
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE on error
 */
ErrCode iMOAB_ReceiveMesh( iMOAB_AppID pid, MPI_Comm* joint_communicator, MPI_Group* sendingGroup, int* scompid )
{
    // Validate input parameters to ensure they are not null pointers
    assert( joint_communicator != nullptr );
    assert( sendingGroup != nullptr );
    assert( scompid != nullptr );

    ErrorCode rval;
    // Get application data and parallel communicator for this component
    appData& data          = context.appDatas[*pid];
    ParallelComm* pco      = context.appDatas[*pid].pcomm;
    MPI_Comm receive       = pco->comm();
    EntityHandle local_set = data.file_set;

    // Handle Fortran-to-C MPI handle conversion if needed
    MPI_Comm global = ( data.is_fortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) )
                                        : *joint_communicator );
    MPI_Group sendGroup =
        ( data.is_fortran ? MPI_Group_f2c( *reinterpret_cast< MPI_Fint* >( sendingGroup ) ) : *sendingGroup );

    // Extract receiver group from the receiver communicator to understand receiver process topology
    MPI_Group receiverGroup;
    int ierr = MPI_Comm_group( receive, &receiverGroup );CHK_MPI_ERR( ierr );

    // Create communication graph with reverse mapping (sender->receiver)
    // This graph will manage the reception of mesh data from sender processes
    ParCommGraph* cgraph =
        new ParCommGraph( global, sendGroup, receiverGroup, *scompid, context.appDatas[*pid].global_id );

    // Store the communication graph in the application's graph map for later use
    // The key is the sender component ID (*scompid)
    context.appDatas[*pid].pgraph[*scompid] = cgraph;

    // Get current receiver rank for process identification
    int receiver_rank = -1;
    MPI_Comm_rank( receive, &receiver_rank );

    // Receive the communication graph from sender processes
    // This graph tells this receiver which senders will send it data
    std::vector< int > pack_array;
    MB_CHK_ERR( cgraph->receive_comm_graph( global, pco, pack_array ) );

    // Determine this receiver's index within the receiver group
    int current_receiver = cgraph->receiver( receiver_rank );

    // Parse the packed communication graph to find which senders contribute to this receiver
    std::vector< int > senders_local;
    size_t n = 0;

    // Walk through the packed graph array to find sender ranks for this receiver
    while( n < pack_array.size() )
    {
        if( current_receiver == pack_array[n] )
        {
            // Found this receiver's entry - extract all contributing sender ranks
            for( int j = 0; j < pack_array[n + 1]; j++ )
            {
                senders_local.push_back( pack_array[n + 2 + j] );
            }
            break;
        }
        // Move to next receiver's entry in the packed array
        n = n + 2 + pack_array[n + 1];
    }

#ifdef VERBOSE
    std::cout << " receiver " << current_receiver << " at rank " << receiver_rank << " will receive from "
              << senders_local.size() << " tasks: ";

    for( int k = 0; k < (int)senders_local.size(); k++ )
    {
        std::cout << " " << senders_local[k];
    }

    std::cout << "\n";
#endif

    // Validate that we have senders contributing to this receiver
    if( senders_local.empty() )
    {
        std::cout << " we do not have any senders for receiver rank " << receiver_rank << "\n";
    }

    // Receive mesh data from all identified sender processes
    // This unpacks mesh entities and deposits them into the local mesh set
    MB_CHK_ERR( cgraph->receive_mesh( global, pco, local_set, senders_local ) );

    // Post-process received mesh data to handle potential vertex duplication
    // Vertices with identical GLOBAL_ID from different senders need to be merged
    Tag idtag;
    MB_CHK_ERR( context.MBI->tag_get_handle( "GLOBAL_ID", idtag ) );

    // Get all entities in the local mesh set
    Range local_ents;
    MB_CHK_ERR( context.MBI->get_entities_by_handle( local_set, local_ents ) );

    // Only perform vertex merging for regular meshes (not point clouds)
    if( !local_ents.all_of_type( MBVERTEX ) )
    {
        // Merge vertices only if we received data from multiple senders
        if( (int)senders_local.size() >= 2 )  // Need to remove duplicate vertices
        {

            // Separate vertices from elements for processing
            Range local_verts = local_ents.subset_by_type( MBVERTEX );
            Range local_elems = subtract( local_ents, local_verts );

            // Temporarily remove vertices from local set for merging
            MB_CHK_ERR( context.MBI->remove_entities( local_set, local_verts ) );

#ifdef VERBOSE
            std::cout << "current_receiver " << current_receiver << " local verts: " << local_verts.size() << "\n";
#endif
            // Merge vertices with identical GLOBAL_ID values
            MergeMesh mm( context.MBI );
            MB_CHK_ERR( mm.merge_using_integer_tag( local_verts, idtag ) );

            // Get the merged vertices back from element connectivity
            Range new_verts;
            MB_CHK_ERR( context.MBI->get_connectivity( local_elems, new_verts ) );

#ifdef VERBOSE
            std::cout << "after merging: new verts: " << new_verts.size() << "\n";
#endif
            // Add the merged vertices back to the local set
            MB_CHK_ERR( context.MBI->add_entities( local_set, new_verts ) );
        }
    }
    else
        // Mark this as a point cloud if we only have vertices
        data.point_cloud = true;

    if( !data.point_cloud )
    {
        // For regular meshes, resolve shared entities (vertices) across process boundaries
        MB_CHK_ERR( pco->resolve_shared_ents( local_set, -1, -1, &idtag ) );
    }
    else
    {
        // For point clouds, set partition tag to current rank for visualization
        Tag densePartTag;
        rval = context.MBI->tag_get_handle( "partition", densePartTag );
        if( nullptr != densePartTag && MB_SUCCESS == rval )
        {
            Range local_verts;
            MB_CHK_ERR( context.MBI->get_entities_by_dimension( local_set, 0, local_verts ) );
            std::vector< int > vals;
            int rank = pco->rank();
            vals.resize( local_verts.size(), rank );
            MB_CHK_ERR( context.MBI->tag_set_data( densePartTag, local_verts, &vals[0] ) );
        }
    }
    // Set the parallel partition tag to identify which process owns this mesh set
    Tag part_tag;
    int dum_id = -1;
    rval       = context.MBI->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_tag,
                                              MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );

    if( part_tag == nullptr || ( ( rval != MB_SUCCESS ) && ( rval != MB_ALREADY_ALLOCATED ) ) )
    {
        std::cout << " can't get par part tag.\n";
        return moab::MB_FAILURE;
    }

    // Tag the local set with the current process rank
    int rank = pco->rank();
    MB_CHK_ERR( context.MBI->tag_set_data( part_tag, &local_set, 1, &rank ) );

    // Ensure that the GLOBAL_ID tag is properly defined for entity identification
    int tagtype = 0;  // dense, integer
    int numco   = 1;  // size
    int tagindex;     // not used
    MB_CHK_ERR( iMOAB_DefineTagStorage( pid, "GLOBAL_ID", &tagtype, &numco, &tagindex ) );

    // Update mesh information with current data statistics
    MB_CHK_ERR( iMOAB_UpdateMeshInfo( pid ) );

    // Clean up temporary MPI group to prevent memory leaks
    MPI_Group_free( &receiverGroup );

    return moab::MB_SUCCESS;
}

/**
 * @brief Send element tag data from calling component to another component.
 * @ingroup iMOABParallel
 *
 * @details Transfers tag data (scalar or vector values) associated with mesh elements from one
 * component to another using pre-established communication graph. Transfer based on communication
 * patterns determined by previous mesh transfer or communication graph computation operations.
 *
 * @par Algorithm Workflow:
 * -# <b>COMMUNICATION GRAPH VALIDATION:</b> Look up ParCommGraph for specified context_id,
 *    validate existence and initialization, extract graph and parallel communicator.
 * -# <b>ENTITY SELECTION:</b> Determine target entities based on component type (point cloud: local_verts,
 *    regular mesh: owned_elems). Handle special cases (TempestRemap coverage mesh, intersection coverage sets).
 * -# <b>TAG PROCESSING:</b> Parse tag names using colon separator for multiple tags, resolve tag handles,
 *    validate existence. Support scalar and vector tag data.
 * -# <b>DATA TRANSFER:</b> Use ParCommGraph to determine routing, pack tag values associated with selected
 *    entities, execute non-blocking point-to-point transfers via send_tag_values.
 *
 * @par Data Transferred:
 * - <b>Tag values:</b> Scalar or vector data associated with mesh entities
 * - <b>Entity mappings:</b> Which entities the tag values belong to
 * - <b>Tag metadata:</b> Tag type, size, and component information
 * - <b>Global IDs:</b> For entity identification across processes
 *
 * @param[in] pid                 Application ID for sender component
 * @param[in] tag_storage_name    Tag name(s) to send (colon-separated for multiple tags)
 * @param[in] joint_communicator  MPI communicator spanning sender and receiver groups
 * @param[in] context_id          Communication graph context ID
 *
 * @note <b>Use cases:</b> Transfer solution data, send material properties/BCs, exchange computed quantities
 * @warning Communication graph must exist before calling (establish via mesh transfer or compute graph)
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if graph missing or tags not found
 */
ErrCode iMOAB_SendElementTag( iMOAB_AppID pid,
                              const iMOAB_String tag_storage_name,
                              MPI_Comm* joint_communicator,
                              int* context_id )
{
    // Get application data and look up the communication graph for this context
    appData& data                               = context.appDatas[*pid];
    std::map< int, ParCommGraph* >::iterator mt = data.pgraph.find( *context_id );
    if( mt == data.pgraph.end() )
    {
        // Communication graph not found - this means mesh transfer hasn't been performed yet
        std::cout << " no par com graph for context_id:" << *context_id << " available contexts:";
        for( auto mit = data.pgraph.begin(); mit != data.pgraph.end(); mit++ )
            std::cout << "  " << mit->first;
        std::cout << "\n";
        return moab::MB_FAILURE;
    }

    // Extract communication graph and parallel communicator
    ParCommGraph* cgraph = mt->second;
    ParallelComm* pco    = context.appDatas[*pid].pcomm;

    // Handle Fortran-to-C MPI handle conversion if needed
    MPI_Comm global = ( data.is_fortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) )
                                        : *joint_communicator );

    // Determine target entities: vertices for point clouds, elements for regular meshes
    Range owned = ( data.point_cloud ? data.local_verts : data.owned_elems );

#ifdef MOAB_HAVE_TEMPESTREMAP
    // Handle TempestRemap coverage mesh entities for intersection applications
    if( data.tempestData.remapper != nullptr )  // This is the case for intersection operations
    {
        EntityHandle cover_set = data.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( cover_set, 2, owned ) );
    }
#else
    // Handle coverage set entities for intersection applications (non-TempestRemap case)
    // This covers cases where elements have been instantiated in coverage sets during intersection
    EntityHandle cover_set = cgraph->get_cover_set();  // Non-null only for intersection applications
    if( 0 != cover_set ) MB_CHK_ERR( context.MBI->get_entities_by_dimension( cover_set, 2, owned ) );
#endif

    // Parse tag names from the input string (multiple tags separated by colons)
    std::string tag_name( tag_storage_name );

    // Parse multiple tag names separated by colons
    std::vector< std::string > tagNames;
    std::vector< Tag > tagHandles;
    std::string separator( ":" );
    split_tag_names( tag_name, separator, tagNames );

    // Resolve tag handles for each specified tag name
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        Tag tagHandle;
        ErrorCode rval = context.MBI->tag_get_handle( tagNames[i].c_str(), tagHandle );
        if( MB_SUCCESS != rval || nullptr == tagHandle )
        {
            MB_CHK_SET_ERR( moab::MB_FAILURE,
                            "can't get tag handle with name: " << tagNames[i].c_str() << " at index " << i );
        }
        tagHandles.push_back( tagHandle );
    }

    // Execute the tag data transfer using the communication graph
    // This packs tag values and sends them to appropriate receiver processes
    // pco is needed for MOAB operations (packing), not for MPI communication
    MB_CHK_ERR( cgraph->send_tag_values( global, pco, owned, tagHandles ) );

    return moab::MB_SUCCESS;
}

/**
 * @brief Receive element tag data from another component.
 * @ingroup iMOABParallel
 *
 * @details Receives tag data (scalar or vector values) associated with mesh elements from another
 * component using pre-established communication graph. Received tag values are applied to
 * corresponding entities in local mesh, enabling data exchange between coupled components.
 *
 * @par Algorithm Workflow:
 * -# <b>COMMUNICATION GRAPH VALIDATION:</b> Look up ParCommGraph for specified context_id,
 *    validate existence and initialization, extract graph and parallel communicator.
 * -# <b>ENTITY SELECTION:</b> Determine target entities based on component type (point cloud: local_verts,
 *    regular mesh: owned_elems). Handle special cases (coverage set, TempestRemap coverage mesh).
 * -# <b>TAG PROCESSING:</b> Parse tag names using colon separator for multiple tags, resolve tag handles,
 *    validate existence. Support scalar and vector tag data.
 * -# <b>DATA RECEPTION:</b> Use ParCommGraph to receive tag data via receive_tag_values, unpack received
 *    tag values and associate with local entities, execute non-blocking point-to-point receives.
 *
 * @par Data Received:
 * - <b>Tag values:</b> Scalar or vector data from sender component
 * - <b>Entity mappings:</b> Which entities the tag values belong to
 * - <b>Tag metadata:</b> Tag type, size, and component information
 * - <b>Global IDs:</b> For entity identification and matching
 *
 * @param[in] pid                 Application ID for receiver component
 * @param[in] tag_storage_name    Tag name(s) to receive (colon-separated for multiple tags)
 * @param[in] joint_communicator  MPI communicator spanning sender and receiver groups
 * @param[in] context_id          Communication graph context ID
 *
 * @note <b>Use cases:</b> Receive solution data, get material properties/BCs, accept computed quantities
 * @warning Communication graph must exist before calling (establish via mesh transfer or compute graph)
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if graph missing or tags not found
 */
ErrCode iMOAB_ReceiveElementTag( iMOAB_AppID pid,
                                 const iMOAB_String tag_storage_name,
                                 MPI_Comm* joint_communicator,
                                 int* context_id )
{
    // Get application data and look up the communication graph for this context
    appData& data                               = context.appDatas[*pid];
    std::map< int, ParCommGraph* >::iterator mt = data.pgraph.find( *context_id );
    if( mt == data.pgraph.end() )
    {
        // Communication graph not found - this means mesh transfer hasn't been performed yet
        std::cout << " no par com graph for context_id:" << *context_id << " available contexts:";
        for( auto mit = data.pgraph.begin(); mit != data.pgraph.end(); mit++ )
            std::cout << "  " << mit->first;
        std::cout << "\n";
        return moab::MB_FAILURE;
    }

    // Extract communication graph and parallel communicator
    ParCommGraph* cgraph = mt->second;
    ParallelComm* pco    = context.appDatas[*pid].pcomm;

    // Handle Fortran-to-C MPI handle conversion if needed
    MPI_Comm global = ( data.is_fortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) )
                                        : *joint_communicator );

    // Determine target entities: vertices for point clouds, elements for regular meshes
    Range owned = ( data.point_cloud ? data.local_verts : data.owned_elems );

    // Handle coverage set entities for intersection applications
    // This covers cases where elements have been instantiated in coverage sets during intersection
    EntityHandle cover_set = cgraph->get_cover_set();
    if( 0 != cover_set ) MB_CHK_ERR( context.MBI->get_entities_by_dimension( cover_set, 2, owned ) );

#ifdef MOAB_HAVE_TEMPESTREMAP
    // Handle TempestRemap coverage mesh entities for intersection applications
    if( data.tempestData.remapper != nullptr )  // This is the case for intersection operations
    {
        cover_set = data.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( cover_set, 2, owned ) );
    }
#endif

    // Parse tag names from the input string (multiple tags separated by colons)
    std::string tag_name( tag_storage_name );

    // Parse multiple tag names separated by colons
    std::vector< std::string > tagNames;
    std::vector< Tag > tagHandles;
    std::string separator( ":" );
    split_tag_names( tag_name, separator, tagNames );

    // Resolve tag handles for each specified tag name
    for( size_t i = 0; i < tagNames.size(); i++ )
    {
        Tag tagHandle;
        ErrorCode rval = context.MBI->tag_get_handle( tagNames[i].c_str(), tagHandle );
        if( MB_SUCCESS != rval || nullptr == tagHandle )
        {
            MB_CHK_SET_ERR( moab::MB_FAILURE,
                            " can't get tag handle for tag named:" << tagNames[i].c_str() << " at index " << i );
        }
        tagHandles.push_back( tagHandle );
    }

#ifdef VERBOSE
    std::cout << pco->rank() << ". Looking to receive data for tags: " << tag_name
              << " and file set = " << ( data.file_set ) << "\n";
#endif
    // Execute the tag data reception using the communication graph
    // This receives tag values from sender processes and applies them to local entities
    // pco is needed for MOAB operations (unpacking), not for MPI communication
    MB_CHK_SET_ERR( cgraph->receive_tag_values( global, pco, owned, tagHandles ), "failed to receive tag values" );

#ifdef VERBOSE
    std::cout << pco->rank() << ". Looking to receive data for tags: " << tag_name << "\n";
#endif

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_FreeSenderBuffers( iMOAB_AppID pid, int* context_id )
{
    // Find the communication graph that holds the send buffer information
    // This function is called on the sender side only to clean up after data transfer
    appData& data                               = context.appDatas[*pid];
    std::map< int, ParCommGraph* >::iterator mt = data.pgraph.find( *context_id );
    if( mt == data.pgraph.end() ) return moab::MB_FAILURE;  // Communication graph not found

    // Release all non-blocking send buffers associated with this communication context
    // This frees memory and completes any pending asynchronous send operations
    mt->second->release_send_buffers();
    return moab::MB_SUCCESS;
}

/**
 * @brief Compute communication graph between two components using rendezvous algorithm.
 * @ingroup iMOABParallel
 *
 * @details Implements sophisticated rendezvous-based algorithm to determine which processes from
 * component 1 need to communicate with which processes from component 2. Uses global entity IDs
 * (or DOFs) as "rendezvous points" to discover communication patterns between components.
 *
 * @par Algorithm Workflow:
 * -# <b>COMMUNICATOR SETUP:</b> Resolve MPI communicator and group handles (Fortran F2C conversions),
 *    create bidirectional ParCommGraph instances for both components, store graphs in pgraph maps.
 * -# <b>ENTITY ID COLLECTION:</b> Extract entity IDs based on type parameter:
 *    - type=1: GLOBAL_DOFS from MBQUAD elements (spectral elements)
 *    - type=2: GLOBAL_ID from MBVERTEX entities (vertex-based coupling)
 *    - type=3: GLOBAL_ID from 2D elements (finite volume meshes)
 *    - Handle TempestRemap coverage mesh entities
 * -# <b>RENDEZVOUS ALGORITHM:</b> Create TupleList for each component with (target_proc, entity_id) pairs.
 *    Use hash-based distribution (entity_id % numProcs). Execute crystal router (gs_transfer) to send
 *    entity IDs to rendezvous processes. Sort received data by entity ID for efficient matching.
 * -# <b>COMMUNICATION DISCOVERY:</b> Synchronously iterate through both sorted TupleLists, find matching
 *    entity IDs between components (intersection), record which processes share each entity, build
 *    reverse communication maps.
 * -# <b>GRAPH CONSTRUCTION:</b> Send communication information back to original processes via crystal
 *    router. Each process learns communication partners. Store patterns in ParCommGraph.
 *
 * @par Data Processed:
 * - <b>Entity IDs:</b> Global identifiers serving as rendezvous points
 * - <b>Process mappings:</b> Which processes own which entities
 * - <b>Communication patterns:</b> Sender→receiver relationships
 *
 * @param[in] pid1               Application ID for component 1 (or negative if not on these ranks)
 * @param[in] pid2               Application ID for component 2 (or negative if not on these ranks)
 * @param[in] joint_communicator MPI communicator spanning both component groups
 * @param[in] group1             MPI group for component 1 processes
 * @param[in] group2             MPI group for component 2 processes
 * @param[in] type1              Entity type for component 1 (1=DOFs, 2=vertices, 3=elements)
 * @param[in] type2              Entity type for component 2 (1=DOFs, 2=vertices, 3=elements)
 * @param[in] comp1              Component 1 external ID
 * @param[in] comp2              Component 2 external ID
 *
 * @note <b>Complexity:</b> O(N log N) sorting, O(N/P) communication per process
 * @note <b>Use cases:</b> Coupling different mesh types, establishing comm for data transfer, mesh intersection prep
 * @warning Both components must have consistent global entity IDs for rendezvous to work
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 */
//#define VERBOSE
ErrCode iMOAB_ComputeCommGraph( iMOAB_AppID pid1,
                                iMOAB_AppID pid2,
                                MPI_Comm* joint_communicator,
                                MPI_Group* group1,
                                MPI_Group* group2,
                                int* type1,
                                int* type2,
                                int* comp1,
                                int* comp2 )
{
    // Validate input parameters
    assert( joint_communicator );
    assert( group1 );
    assert( group2 );
    int localRank = 0, numProcs = 1;

    // Determine if either component is using Fortran (affects MPI handle conversion)
    bool isFortran = false;
    if( *pid1 >= 0 ) isFortran = isFortran || context.appDatas[*pid1].is_fortran;
    if( *pid2 >= 0 ) isFortran = isFortran || context.appDatas[*pid2].is_fortran;

    // Handle Fortran-to-C MPI handle conversion if needed
    MPI_Comm global =
        ( isFortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) ) : *joint_communicator );
    MPI_Group srcGroup = ( isFortran ? MPI_Group_f2c( *reinterpret_cast< MPI_Fint* >( group1 ) ) : *group1 );
    MPI_Group tgtGroup = ( isFortran ? MPI_Group_f2c( *reinterpret_cast< MPI_Fint* >( group2 ) ) : *group2 );

    // Get process information for the joint communicator
    MPI_Comm_rank( global, &localRank );
    MPI_Comm_size( global, &numProcs );

    // Create bidirectional communication graphs for both components

    // Create communication graphs for both components (bidirectional)
    ParCommGraph* cgraph     = nullptr;
    ParCommGraph* cgraph_rev = nullptr;

    if( *pid1 >= 0 )
    {
        // Create communication graph for component 1 -> component 2
        appData& data = context.appDatas[*pid1];
        auto mt       = data.pgraph.find( *comp2 );
        if( mt != data.pgraph.end() ) data.pgraph.erase( mt );  // Remove existing graph if present
        cgraph                                 = new ParCommGraph( global, srcGroup, tgtGroup, *comp1, *comp2 );
        context.appDatas[*pid1].pgraph[*comp2] = cgraph;  // Store with target component ID as key
    }

    if( *pid2 >= 0 )
    {
        // Create reverse communication graph for component 2 -> component 1
        appData& data = context.appDatas[*pid2];
        auto mt       = data.pgraph.find( *comp1 );
        if( mt != data.pgraph.end() ) data.pgraph.erase( mt );  // Remove existing graph if present
        cgraph_rev                             = new ParCommGraph( global, tgtGroup, srcGroup, *comp2, *comp1 );
        context.appDatas[*pid2].pgraph[*comp1] = cgraph_rev;  // Store with source component ID as key
    }

    // Initialize TupleLists for rendezvous algorithm
    // Each tuple contains: (target_process, entity_id) pairs for hash-based distribution
    TupleList TLcomp1;
    TLcomp1.initialize( 2, 0, 0, 0, 0 );  // 2 integers: target_proc, entity_id
    TupleList TLcomp2;
    TLcomp2.initialize( 2, 0, 0, 0, 0 );  // 2 integers: target_proc, entity_id

    // Enable write access for building the tuple lists
    TLcomp1.enableWriteAccess();

    // Get tag handles for entity identification
    // GLOBAL_DOFS: for spectral elements (type 1)
    // GLOBAL_ID: for vertices (type 2) and 2D elements (type 3)
    Tag gdsTag;
    int lenTagType1 = 1;
    if( 1 == *type1 || 1 == *type2 )
    {
        // Get GLOBAL_DOFS tag for spectral elements (usually 16 DOFs per element)
        MB_CHK_ERR( context.MBI->tag_get_handle( "GLOBAL_DOFS", gdsTag ) );
        MB_CHK_ERR( context.MBI->tag_get_length( gdsTag, lenTagType1 ) );  // Usually 16 DOFs per element
    }
    Tag gidTag = context.MBI->globalId_tag();  // Standard GLOBAL_ID tag

    // Collect entity IDs from component 1 for rendezvous algorithm
    std::vector< int > valuesComp1;
    if( *pid1 >= 0 )
    {
        appData& data1     = context.appDatas[*pid1];
        EntityHandle fset1 = data1.file_set;

        // Handle TempestRemap coverage mesh for intersection applications
#ifdef MOAB_HAVE_TEMPESTREMAP
        if( data1.tempestData.remapper != nullptr )  // This is the case for intersection operations
            fset1 = data1.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
#endif

        Range ents_of_interest;
        if( *type1 == 1 )  // Spectral elements: get GLOBAL_DOFS from MBQUAD elements
        {
            assert( gdsTag );
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset1, MBQUAD, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() * lenTagType1 );
            MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, ents_of_interest, &valuesComp1[0] ) );
        }
        else if( *type1 == 2 )  // Vertex-based coupling: get GLOBAL_ID from MBVERTEX entities
        {
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset1, MBVERTEX, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp1[0] ) );
        }
        else if( *type1 == 3 )  // Finite volume meshes: get GLOBAL_ID from 2D elements
        {
            MB_CHK_ERR( context.MBI->get_entities_by_dimension( fset1, 2, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp1[0] ) );
        }
        else
        {
            MB_CHK_ERR( MB_FAILURE );  // Only types 1, 2, or 3 are supported
        }

        // Build tuple list for component 1: (target_process, entity_id) pairs
        // Use hash-based distribution: entity_id % numProcs determines target process
        std::set< int > uniq( valuesComp1.begin(), valuesComp1.end() );  // Remove duplicates
        TLcomp1.resize( uniq.size() );
        for( std::set< int >::iterator sit = uniq.begin(); sit != uniq.end(); sit++ )
        {
            int marker               = *sit;               // Entity ID
            int to_proc              = marker % numProcs;  // Hash-based target process
            int n                    = TLcomp1.get_n();
            TLcomp1.vi_wr[2 * n]     = to_proc;  // Target process
            TLcomp1.vi_wr[2 * n + 1] = marker;   // Entity ID
            TLcomp1.inc_n();
        }
    }

    // Execute crystal router transfer for component 1 entity IDs
    // This sends entity IDs to rendezvous processes based on hash distribution
    ProcConfig pc( global );                            // Process configuration for crystal router
    pc.crystal_router()->gs_transfer( 1, TLcomp1, 0 );  // Transfer to joint tasks with markers

    // Sort component 1 tuple list by entity ID (key 1) for efficient matching
#ifdef VERBOSE
    std::stringstream ff1;
    ff1 << "TLcomp1_" << localRank << ".txt";
    TLcomp1.print_to_file( ff1.str().c_str() );  // Debug output before sorting
#endif
    moab::TupleList::buffer sort_buffer;
    sort_buffer.buffer_init( TLcomp1.get_n() );
    TLcomp1.sort( 1, &sort_buffer );  // Sort by entity ID (column 1)
    sort_buffer.reset();
#ifdef VERBOSE
    TLcomp1.print_to_file( ff1.str().c_str() );  // Debug output after sorting
#endif
    // Now process component 2 in the same way as component 1
    TLcomp2.enableWriteAccess();

    // Collect entity IDs from component 2 for rendezvous algorithm
    std::vector< int > valuesComp2;
    if( *pid2 >= 0 )
    {
        appData& data2     = context.appDatas[*pid2];
        EntityHandle fset2 = data2.file_set;

        // Handle TempestRemap coverage mesh for intersection applications
#ifdef MOAB_HAVE_TEMPESTREMAP
        if( data2.tempestData.remapper != nullptr )  // This is the case for intersection operations
            fset2 = data2.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
#endif

        Range ents_of_interest;
        if( *type2 == 1 )  // Spectral elements: get GLOBAL_DOFS from MBQUAD elements
        {
            assert( gdsTag );
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset2, MBQUAD, ents_of_interest ) );
            valuesComp2.resize( ents_of_interest.size() * lenTagType1 );
            MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, ents_of_interest, &valuesComp2[0] ) );
        }
        else if( *type2 == 2 )  // Vertex-based coupling: get GLOBAL_ID from MBVERTEX entities
        {
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset2, MBVERTEX, ents_of_interest ) );
            valuesComp2.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp2[0] ) );
        }
        else if( *type2 == 3 )  // Finite volume meshes: get GLOBAL_ID from 2D elements
        {
            MB_CHK_ERR( context.MBI->get_entities_by_dimension( fset2, 2, ents_of_interest ) );
            valuesComp2.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp2[0] ) );
        }
        else
        {
            MB_CHK_ERR( MB_FAILURE );  // Only types 1, 2, or 3 are supported
        }
        // Build tuple list for component 2: (target_process, entity_id) pairs
        // Use hash-based distribution: entity_id % numProcs determines target process
        std::set< int > uniq( valuesComp2.begin(), valuesComp2.end() );  // Remove duplicates
        TLcomp2.resize( uniq.size() );
        for( std::set< int >::iterator sit = uniq.begin(); sit != uniq.end(); sit++ )
        {
            int marker               = *sit;               // Entity ID
            int to_proc              = marker % numProcs;  // Hash-based target process
            int n                    = TLcomp2.get_n();
            TLcomp2.vi_wr[2 * n]     = to_proc;  // Target process
            TLcomp2.vi_wr[2 * n + 1] = marker;   // Entity ID
            TLcomp2.inc_n();
        }
    }

    // Execute crystal router transfer for component 2 entity IDs
    pc.crystal_router()->gs_transfer( 1, TLcomp2, 0 );  // Transfer to joint tasks with markers

    // Sort component 2 tuple list by entity ID (key 1) for efficient matching
#ifdef VERBOSE
    std::stringstream ff2;
    ff2 << "TLcomp2_" << localRank << ".txt";
    TLcomp2.print_to_file( ff2.str().c_str() );  // Debug output before sorting
#endif
    sort_buffer.buffer_reserve( TLcomp2.get_n() );
    TLcomp2.sort( 1, &sort_buffer );  // Sort by entity ID (column 1)
    sort_buffer.reset();
#ifdef VERBOSE
    TLcomp2.print_to_file( ff2.str().c_str() );  // Debug output after sorting
#endif
    // RENDEZVOUS ALGORITHM: Find matching entity IDs between components
    // Now we need to send back communication information from the rendezvous point
    // Loop synchronously over both sorted tuple lists to find matching entity IDs
    // Build new tuple lists to send communication information back to original processes

    // Tuple lists for sending communication info back to components
    // Format: (target_process, entity_id, source_process_from_other_component)
    TupleList TLBackToComp1;
    TLBackToComp1.initialize( 3, 0, 0, 0, 0 );  // 3 integers: target_proc, entity_id, source_proc_from_comp2
    TLBackToComp1.enableWriteAccess();

    TupleList TLBackToComp2;
    TLBackToComp2.initialize( 3, 0, 0, 0, 0 );  // 3 integers: target_proc, entity_id, source_proc_from_comp1
    TLBackToComp2.enableWriteAccess();

    // Get sizes of both tuple lists for synchronized iteration
    int n1 = TLcomp1.get_n();
    int n2 = TLcomp2.get_n();

    // Synchronized iteration through both sorted tuple lists to find matching entity IDs
    int indexInTLComp1 = 0;
    int indexInTLComp2 = 0;
    if( n1 > 0 && n2 > 0 )
    {
        while( indexInTLComp1 < n1 && indexInTLComp2 < n2 )  // Continue until either list is exhausted
        {
            // Get current entity IDs from both components
            int currentValue1 = TLcomp1.vi_rd[2 * indexInTLComp1 + 1];  // Entity ID from comp1
            int currentValue2 = TLcomp2.vi_rd[2 * indexInTLComp2 + 1];  // Entity ID from comp2

            if( currentValue1 < currentValue2 )
            {
                // Entity ID exists in comp1 but not in comp2 - skip comp1 entry
                // This is normal for non-overlapping meshes
                indexInTLComp1++;
                continue;
            }
            if( currentValue1 > currentValue2 )
            {
                // Entity ID exists in comp2 but not in comp1 - skip comp2 entry
                // This is normal for non-overlapping meshes
                indexInTLComp2++;
                continue;
            }
            // Found matching entity ID! Count consecutive entries with same ID in both lists
            int size1 = 1;
            int size2 = 1;
            while( indexInTLComp1 + size1 < n1 && currentValue1 == TLcomp1.vi_rd[2 * ( indexInTLComp1 + size1 ) + 1] )
                size1++;  // Count consecutive entries with same entity ID in comp1
            while( indexInTLComp2 + size2 < n2 && currentValue2 == TLcomp2.vi_rd[2 * ( indexInTLComp2 + size2 ) + 1] )
                size2++;  // Count consecutive entries with same entity ID in comp2

            // Create communication pairs for all combinations of processes that share this entity ID
            for( int i1 = 0; i1 < size1; i1++ )
            {
                for( int i2 = 0; i2 < size2; i2++ )
                {
                    // Send communication info back to component 1 processes
                    int n = TLBackToComp1.get_n();
                    TLBackToComp1.reserve();
                    TLBackToComp1.vi_wr[3 * n] =
                        TLcomp1.vi_rd[2 * ( indexInTLComp1 + i1 )];  // Target process from comp1
                    TLBackToComp1.vi_wr[3 * n + 1] = currentValue1;  // Shared entity ID
                    TLBackToComp1.vi_wr[3 * n + 2] =
                        TLcomp2.vi_rd[2 * ( indexInTLComp2 + i2 )];  // Source process from comp2

                    // Send communication info back to component 2 processes
                    n = TLBackToComp2.get_n();
                    TLBackToComp2.reserve();
                    TLBackToComp2.vi_wr[3 * n] =
                        TLcomp2.vi_rd[2 * ( indexInTLComp2 + i2 )];  // Target process from comp2
                    TLBackToComp2.vi_wr[3 * n + 1] = currentValue1;  // Shared entity ID
                    TLBackToComp2.vi_wr[3 * n + 2] =
                        TLcomp1.vi_rd[2 * ( indexInTLComp1 + i1 )];  // Source process from comp1
                }
            }

            // Advance indices past all entries with the current entity ID
            indexInTLComp1 += size1;
            indexInTLComp2 += size2;
        }
    }
    // Send communication information back to original component processes
    pc.crystal_router()->gs_transfer( 1, TLBackToComp1, 0 );  // Send to component 1 processes
    pc.crystal_router()->gs_transfer( 1, TLBackToComp2, 0 );  // Send to component 2 processes

    // Process communication information on component 1 processes
    if( *pid1 >= 0 )
    {
        // Sort received communication information by entity ID for efficient processing
#ifdef VERBOSE
        std::stringstream f1;
        f1 << "TLBack1_" << localRank << ".txt";
        TLBackToComp1.print_to_file( f1.str().c_str() );  // Debug output before sorting
#endif
        sort_buffer.buffer_reserve( TLBackToComp1.get_n() );
        TLBackToComp1.sort( 1, &sort_buffer );  // Sort by entity ID (column 1)
        sort_buffer.reset();
#ifdef VERBOSE
        TLBackToComp1.print_to_file( f1.str().c_str() );  // Debug output after sorting
#endif

        // Set up communication graph for component 1 based on discovered communication patterns
        // This establishes which processes component 1 needs to communicate with
        cgraph->settle_comm_by_ids( *comp1, TLBackToComp1, valuesComp1 );
    }
    // Process communication information on component 2 processes
    if( *pid2 >= 0 )
    {
        // Sort received communication information by source process for efficient processing
#ifdef VERBOSE
        std::stringstream f2;
        f2 << "TLBack2_" << localRank << ".txt";
        TLBackToComp2.print_to_file( f2.str().c_str() );  // Debug output before sorting
#endif
        sort_buffer.buffer_reserve( TLBackToComp2.get_n() );
        TLBackToComp2.sort( 2, &sort_buffer );  // Sort by source process (column 2)
        sort_buffer.reset();
#ifdef VERBOSE
        TLBackToComp2.print_to_file( f2.str().c_str() );  // Debug output after sorting
#endif

        // Set up reverse communication graph for component 2 based on discovered communication patterns
        // This establishes which processes component 2 needs to communicate with
        cgraph_rev->settle_comm_by_ids( *comp2, TLBackToComp2, valuesComp2 );
    }

    // Communication graph computation complete - both components now know their communication partners
    return moab::MB_SUCCESS;
}

/**
 * @brief Merge duplicate vertices in a mesh and update connectivity information.
 * @ingroup iMOABParallel
 *
 * @details This function consolidates duplicate vertices that may exist after mesh operations
 * such as intersection or migration. Vertices with identical GLOBAL_ID values are
 * merged, and connectivity of elements is updated accordingly. The function performs
 * both local (within-process) and parallel (across-process) vertex merging.
 *
 * @par Algorithm Workflow:
 * -# <b>TAG COLLECTION:</b> Collect critical tags that need to be preserved during merge:
 *    - GLOBAL_ID: Required for vertex identification (fatal error if missing)
 *    - area: Optional - element area values to preserve
 *    - frac: Optional - fraction values for conservation
 * -# <b>LOCAL VERTEX DEDUPLICATION:</b> Remove duplicate vertices within each process using
 *    spatial tolerance (1.0e-9). IntxUtils::remove_duplicate_vertices merges vertices and
 *    updates connectivity, preserving tag values from first occurrence.
 * -# <b>MATERIAL SET CLEANUP:</b> Clear all elements from material sets and re-populate with
 *    valid elements from the current file set. Ensures material set integrity after topology changes.
 * -# <b>MESH INFO UPDATE:</b> Refresh internal mesh information (vertex counts, element counts,
 *    owned/ghosted entity information) and recalculate mesh statistics.
 * -# <b>PARALLEL VERTEX MERGING:</b> ParallelMergeMesh handles vertices shared across process
 *    boundaries. Identifies vertices with identical GLOBAL_ID on different processes and merges
 *    them while updating parallel communication graph.
 * -# <b>GLOBAL ID ASSIGNMENT:</b> Reassign consistent global IDs to vertices after merge.
 *    Only updates vertex global IDs (cells already have correct IDs).
 * -# <b>PARTITION TAG SETUP:</b> Create/retrieve PARALLEL_PARTITION tag for visualization
 *    indicating which process owns each mesh set.
 *
 * @par Data Processed:
 * - <b>Vertices:</b> Spatial coordinates and associated tags
 * - <b>Elements:</b> Connectivity information referencing vertices
 * - <b>Tags:</b> GLOBAL_ID (required), area, frac (optional)
 * - <b>Material sets:</b> Collections of elements with shared properties
 *
 * @par Common Use Cases:
 * - After mesh intersection: Overlapping meshes create duplicate vertices
 * - After mesh migration: Vertices copied to multiple processes need merging
 * - After remapping: Source/target vertices may coincide
 * - Quality improvement: Remove numerical duplicates from mesh operations
 *
 * @param[in] pid  Application ID for the mesh
 *
 * @pre GLOBAL_ID tag must exist on vertices
 * @post Duplicate vertices merged both locally and across processes
 * @post Element connectivity updated to reference merged vertices
 * @post PARALLEL_PARTITION tag set for visualization
 *
 * @note <b>Parallelism:</b> Performs both local (per-process) and parallel (across-process) merging
 * @warning Modifies mesh topology and connectivity - may invalidate cached entity handles
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE if GLOBAL_ID tag missing
 */
ErrCode iMOAB_MergeVertices( iMOAB_AppID pid )
{
    // =========================================================================
    // Validate inputs and initialize local variables
    // =========================================================================

    // Get application data and parallel communicator for this component
    // data contains mesh info, ranges, and pcomm holds parallel communication state
    appData& data     = context.appDatas[*pid];
    ParallelComm* pco = data.pcomm;

    // Initialize return value to success state
    ErrorCode rval = MB_SUCCESS;

    // =========================================================================
    // STEP 1: Collect tags that need to be preserved during vertex merging
    // =========================================================================

    // Initialize vector to hold tags that should be preserved during merge
    // These tags contain important data that must be maintained when vertices are merged
    // The merge operation will preserve tag values from the first vertex in each duplicate set
    std::vector< Tag > tagsList;
    Tag tag;

    // GLOBAL_ID is mandatory - used to identify duplicate vertices across processes
    // Without this tag, we cannot determine which vertices should be merged
    tag = context.MBI->globalId_tag();
    if( !tag )
    {
        // Fatal error: GLOBAL_ID tag is required for merge operations
        // Return failure as merge cannot proceed without vertex identification capability
        return moab::MB_FAILURE;
    }
    // Add GLOBAL_ID to the list of tags to preserve during merge
    tagsList.push_back( tag );

    // Area tag (optional) - contains element area values for conservation checks
    // Used in remapping applications to ensure conservation properties
    rval = context.MBI->tag_get_handle( "area", tag );
    if( tag && MB_SUCCESS == rval )
    {
        // Area tag found, include it in preserved tags
        tagsList.push_back( tag );
    }

    // frac tag (optional) - contains fractional coverage for remapping conservation
    // Used to track fractional ownership in intersection-based remapping
    rval = context.MBI->tag_get_handle( "frac", tag );
    if( tag && MB_SUCCESS == rval )
    {
        // Fraction tag found, include it in preserved tags
        tagsList.push_back( tag );
    }

    // ========================================================================
    // STEP 2: Remove duplicate vertices within this process (local merge)
    // ========================================================================

    // Set spatial tolerance for considering vertices identical (meters on unit sphere)
    // This tolerance accounts for numerical precision issues in floating-point coordinates
    // 1.0e-9 is a typical tolerance for spherical meshes near unit radius
    const double tol = 1.0e-9;

    // Perform local vertex deduplication using MOAB's IntxUtils
    // This function:
    // 1. Identifies spatially coincident vertices within tolerance
    // 2. Updates element connectivity to reference merged vertices
    // 3. Preserves tag values from the first vertex in each duplicate set
    // 4. Removes duplicate vertices from the mesh
    MB_CHK_ERR( IntxUtils::remove_duplicate_vertices( context.MBI, data.file_set, tol,
                                                      tagsList ) );  // Exit on failure to merge local vertices

    // ========================================================================
    // STEP 3: Clean up material sets that may reference deleted elements
    // ========================================================================

    // Material sets organize elements by material properties (e.g., land, ocean, ice)
    // After vertex merge, some elements may have been deleted or connectivity modified
    // Update material sets to maintain mesh data integrity

    // Retrieve all material sets currently defined on the file set
    // Material sets are identified by having a MATERIAL_SET_TAG_NAME tag
    MB_CHK_ERR( context.MBI->get_entities_by_type_and_tag(
        data.file_set, MBENTITYSET, &( context.material_tag ), nullptr, 1, data.mat_sets,
        Interface::UNION ) );  // Exit if we can't retrieve material sets

    // Update material sets if any were found
    if( !data.mat_sets.empty() )
    {
        // For simplicity, we only clean the first material set
        // In production code, all material sets should be systematically updated
        // This is a limitation that could be improved in future versions
        EntityHandle matSet = data.mat_sets[0];
        Range elems;

        // Get current 2D elements in this material set (may include stale references)
        MB_CHK_ERR(
            context.MBI->get_entities_by_dimension( matSet, 2, elems ) );  // Exit if we can't get material set contents

        // Remove all current elements from the material set
        // This clears potentially stale element references
        MB_CHK_ERR( context.MBI->remove_entities( matSet, elems ) );  // Exit on failure to remove stale elements

        // Get only valid 2D elements from the current file set
        // This excludes any elements that may have been deleted during merge
        elems.clear();  // Reset range for reuse
        MB_CHK_ERR(
            context.MBI->get_entities_by_dimension( data.file_set, 2, elems ) );  // Exit if we can't get valid elements

        // Add back only the valid elements to the material set
        // This ensures material set consistency after topology changes
        MB_CHK_ERR( context.MBI->add_entities( matSet, elems ) );  // Exit on failure to update material set
    }

    // ========================================================================
    // STEP 4: Update internal mesh information after topology changes
    // ========================================================================

    // Recalculate mesh statistics and update internal data structures
    // This includes vertex/element counts, owned/ghosted entity information
    // iMOAB_UpdateMeshInfo() refreshes all cached mesh information
    MB_CHK_ERR( iMOAB_UpdateMeshInfo( pid ) );  // Exit on failure to update mesh info

    // ========================================================================
    // STEP 5: Perform parallel merge of vertices shared across processes
    // ========================================================================

    // Create parallel merge mesh object with the same spatial tolerance
    // ParallelMergeMesh handles vertices shared across process boundaries
    ParallelMergeMesh pmm( pco, tol );

    // Perform parallel vertex merging across all processes
    // This identifies vertices with identical GLOBAL_ID on different processes
    // and merges them while updating the parallel communication graph
    // Parameters:
    //   - data.file_set: The mesh set to operate on
    //   - false: Skip local merge (already done in step 2)
    //   - 2: Update 2D elements after merge (dimension of elements to fix connectivity for)
    MB_CHK_ERR( pmm.merge( data.file_set,
                           /* do not do local merge*/ false,
                           /*  2d cells*/ 2 ) );  // Exit on failure of parallel merge

    // ========================================================================
    // STEP 6: Assign consistent global IDs to vertices after merge
    // ========================================================================

    // After merging, vertices need new globally unique IDs to maintain consistency
    // Only update vertices (dimension 0) - elements already have correct IDs
    // assign_global_ids() coordinates across all processes to ensure uniqueness
    MB_CHK_ERR( pco->assign_global_ids( data.file_set, /*dim*/ 0 ) );  // Exit on failure to assign global IDs

    // ========================================================================
    // STEP 7: Set partition tag for visualization and debugging
    // ========================================================================

    // Set up PARALLEL_PARTITION tag to indicate which process owns each mesh set
    // This is essential for visualization tools like ParaView/VisIt to show
    // domain decomposition and parallel partitioning

    Tag part_tag;     // Handle for the partition tag
    int dum_id = -1;  // Default value if tag doesn't exist yet

    // Create or retrieve the partition tag (sparse storage, one int per set)
    // Sparse storage is efficient since we only need one value per mesh set
    rval = context.MBI->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_tag,
                                        MB_TAG_CREAT | MB_TAG_SPARSE, &dum_id );

    // Validate that partition tag was created or retrieved successfully
    if( part_tag == nullptr || ( ( MB_SUCCESS != rval ) && ( MB_ALREADY_ALLOCATED != rval ) ) )
    {
        // Error creating/retrieving partition tag
        std::cout << "Failed to create/retrieve PARALLEL_PARTITION tag.\n";
        return moab::MB_FAILURE;
    }

    // Tag the file set with this process's rank to show ownership
    // This allows visualization tools to distinguish between different process domains
    int rank = pco->rank();
    MB_CHK_ERR(
        context.MBI->tag_set_data( part_tag, &data.file_set, 1, &rank ) );  // Exit on failure to set partition data

    // ========================================================================
    // FUNCTION EPILOGUE: Return success status
    // ========================================================================

    // All merge operations completed successfully
    return moab::MB_SUCCESS;
}

#ifdef MOAB_HAVE_TEMPESTREMAP

/**
 * @brief Establishes parallel communication graph for coverage mesh data exchange in coupled simulations.
 *
 * @details In multi-physics coupling (e.g., atmosphere-ocean), components run on different MPI task sets
 * with non-identical mesh decompositions. This function builds optimized communication patterns
 * based on intersection results, enabling targeted data exchange instead of all-to-all communication.
 *
 * @par Algorithm - Two-Hop Crystal Router:
 * -# <b>DISCOVERY PHASE:</b> Intersection/coupler tasks analyze coverage mesh elements, identify which
 *    source tasks originally owned them (via orig_sending_processor tag), and send ID requirements
 *    back to original source tasks via crystal router.
 * -# <b>GRAPH CONSTRUCTION:</b> Source tasks receive requirements and create new ParCommGraph with
 *    coverage-specific send patterns. Receiver tasks create corresponding receive-side graph.
 *    Both graphs stored with unique context_id for later data operations.
 *
 * @par Workflow:
 * - Validate application IDs and retrieve existing ParCommGraphs
 * - On receiver tasks: analyze intersection/coverage elements, build idsFromProcs map, create
 *   receiver-side graph, pack requirements into TupleList
 * - Crystal router exchange: receivers send requirements to senders
 * - On sender tasks: create sender-side graph from received requirements
 *
 * @par Data Structures:
 * - <b>ParCommGraph:</b> Manages point-to-point MPI communication with sender/receiver ID mappings
 * - <b>TupleList:</b> Bulk MPI communication structure (vi_wr/vi_rd for integers, vr_wr/vr_rd for reals)
 * - <b>Tags:</b> GLOBAL_ID (element IDs), orig_sending_processor (source task), SourceParent (parent ID)
 *
 * @param[in] joint_communicator  MPI communicator spanning all participating tasks (component + coupler)
 * @param[in] pid_src            Source component application ID
 * @param[in] pid_migr           Migrated/coverage mesh application ID (on coupler)
 * @param[in] pid_intx           Intersection mesh application ID
 * @param[in] source_id          External ID of source component
 * @param[in] migration_id       External ID of migration/coupler component
 * @param[in] context_id         Unique identifier for this coverage graph
 *
 * @pre Existing ParCommGraph between pid_src and pid_migr must exist
 * @pre Coverage mesh computed with orig_sending_processor tags
 * @pre Global IDs assigned consistently across all tasks
 *
 * @post New ParCommGraphs created on both sender and receiver sides with context_id
 * @post Graphs contain precise send/receive ID mappings for coverage data exchange
 *
 * @note <b>Performance:</b> O(N_cov * log(P)) time complexity, O(N_cov) communication volume
 * @warning <b>Collective operation:</b> Must be called on all tasks in joint_communicator
 *
 * @return moab::MB_SUCCESS on success, moab::MB_FAILURE on error
 */
ErrCode iMOAB_CoverageGraph( MPI_Comm* joint_communicator,
                             iMOAB_AppID pid_src,
                             iMOAB_AppID pid_migr,
                             iMOAB_AppID pid_intx,
                             int* source_id,
                             int* migration_id,
                             int* context_id )
{
    // =========================================================================
    // STEP 1: Input validation and retrieve existing communication graphs
    // =========================================================================

    assert( joint_communicator != nullptr );

    std::vector< int > srcSenders, receivers;  // Process lists from existing graphs
    ParCommGraph* sendGraph = nullptr;         // Original send graph (source → migration)
    bool is_fortran_context = false;           // Handle Fortran/C communicator conversion

    // Validate application IDs exist before accessing application data
    if( pid_src && *pid_src > 0 && context.appDatas.find( *pid_src ) == context.appDatas.end() )
        MB_CHK_SET_ERR( moab::MB_FAILURE, "Invalid source application ID specified: " << *pid_src );
    if( pid_migr && *pid_migr > 0 && context.appDatas.find( *pid_migr ) == context.appDatas.end() )
        MB_CHK_SET_ERR( moab::MB_FAILURE, "Invalid migration/coverage application ID specified: " << *pid_migr );
    if( pid_intx && *pid_intx > 0 && context.appDatas.find( *pid_intx ) == context.appDatas.end() )
        MB_CHK_SET_ERR( moab::MB_FAILURE, "Invalid intersection application ID specified: " << *pid_intx );

    // Retrieve existing send graph from source application
    // This was created during initial migration and defines original communication pattern
    if( *pid_src >= 0 )
    {
        appData& dataSrc       = context.appDatas[*pid_src];
        int default_context_id = *migration_id;  // Original context uses migration_id as key
        assert( dataSrc.global_id == *source_id );
        is_fortran_context = dataSrc.is_fortran || is_fortran_context;
        if( dataSrc.pgraph.find( default_context_id ) != dataSrc.pgraph.end() )
            sendGraph = dataSrc.pgraph[default_context_id];
        else
            MB_CHK_SET_ERR( moab::MB_FAILURE, "Could not find source ParCommGraph with default migration context" );

        // Extract sender/receiver process lists from existing graph
        srcSenders = sendGraph->senders();
        receivers  = sendGraph->receivers();
#ifdef VERBOSE
        std::cout << "senders: " << srcSenders.size() << " first sender: " << srcSenders[0] << std::endl;
#endif
    }

    // Retrieve existing receive graph from migration/coupler application
    // This graph will be copied and modified for coverage-specific communication
    ParCommGraph* recvGraph = nullptr;
    if( *pid_migr >= 0 )
    {
        appData& dataMigr      = context.appDatas[*pid_migr];
        is_fortran_context     = dataMigr.is_fortran || is_fortran_context;
        int default_context_id = *source_id;  // Original context uses source_id as key
        assert( dataMigr.global_id == *migration_id );
        if( dataMigr.pgraph.find( default_context_id ) != dataMigr.pgraph.end() )
            recvGraph = dataMigr.pgraph[default_context_id];
        else
            MB_CHK_SET_ERR( moab::MB_FAILURE,
                            "Could not find coverage receiver ParCommGraph with default migration context" );
        // Extract sender/receiver lists from migration perspective
        srcSenders = recvGraph->senders();
        receivers  = recvGraph->receivers();
#ifdef VERBOSE
        std::cout << "receivers: " << receivers.size() << " first receiver: " << receivers[0] << std::endl;
#endif
    }

    // =========================================================================
    // STEP 2: Initialize TupleList for coverage ID exchange and get rank
    // =========================================================================

    if( *pid_intx >= 0 ) is_fortran_context = context.appDatas[*pid_intx].is_fortran || is_fortran_context;

    // TupleList for packing coverage element IDs to send back to source tasks
    // Format: (destination_proc, global_id)
    TupleList TLcovIDs;
    TLcovIDs.initialize( 2, 0, 0, 0, 0 );  // 2 integers per tuple, dynamic growth
    TLcovIDs.enableWriteAccess();

    // Handle Fortran vs C communicator conversion
    MPI_Comm global = ( is_fortran_context ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( joint_communicator ) )
                                           : *joint_communicator );
    int currentRankInJointComm = -1;CHK_MPI_ERR( MPI_Comm_rank( global, &currentRankInJointComm ) );

    // Get global ID tag for element identification
    Tag gidTag = context.MBI->globalId_tag();
    if( !gidTag ) return moab::MB_TAG_NOT_FOUND;

    // =========================================================================
    // STEP 3: On receiver tasks, analyze coverage mesh and build ID requirements
    // =========================================================================

    // Check if this rank is a receiver (intersection/coupler task)
    if( find( receivers.begin(), receivers.end(), currentRankInJointComm ) != receivers.end() )
    {
        appData& dataIntx      = context.appDatas[*pid_intx];
        EntityHandle intx_set  = dataIntx.file_set;
        EntityHandle cover_set = dataIntx.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );

        // Get tag that identifies which source task originally sent each element
        Tag orgSendProcTag;
        MB_CHK_ERR( context.MBI->tag_get_handle( "orig_sending_processor", orgSendProcTag ) );
        if( !orgSendProcTag ) return moab::MB_TAG_NOT_FOUND;

        // Get intersection elements (if they exist)
        Range intx_cells;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( intx_set, 2, intx_cells ) );

        // Build map: source_processor → set_of_required_global_IDs
        // This tells each source task which elements we need from it
        std::map< int, std::set< int > > idsFromProcs;
        // Process intersection elements to find required source element IDs
        if( intx_cells.size() )
        {
            // SourceParent tag contains global ID of source element that created this intx element
            Tag parentTag;
            MB_CHK_ERR( context.MBI->tag_get_handle( "SourceParent", parentTag ) );
            if( !parentTag ) return moab::MB_TAG_NOT_FOUND;

            for( auto it = intx_cells.begin(); it != intx_cells.end(); ++it )
            {
                EntityHandle intx_cell = *it;

                // Get original processor that sent this element during migration
                int origProc;
                MB_CHK_ERR( context.MBI->tag_get_data( orgSendProcTag, &intx_cell, 1, &origProc ) );

                // Skip ghost elements (origProc == -1)
                if( origProc < 0 ) continue;

                // Get global ID of source parent element
                int gidCell;
                MB_CHK_ERR( context.MBI->tag_get_data( parentTag, &intx_cell, 1, &gidCell ) );

                // Record that we need gidCell from origProc
                idsFromProcs[origProc].insert( gidCell );
            }
        }

        // Process coverage mesh elements (superset of intersection elements)
        // Coverage includes all source elements needed, not just those in active intersection
        {
            Range cover_cells;
            MB_CHK_ERR( context.MBI->get_entities_by_dimension( cover_set, 2, cover_cells ) );

            Tag gidTag = context.MBI->globalId_tag();

            for( auto it = cover_cells.begin(); it != cover_cells.end(); ++it )
            {
                const EntityHandle cov_cell = *it;

                // Get original source processor for this coverage element
                int origProc;
                MB_CHK_ERR( context.MBI->tag_get_data( orgSendProcTag, &cov_cell, 1, &origProc ) );

                // Skip unassigned/ghost elements
                if( origProc < 0 ) continue;

                // Get global ID and add to requirements for this source processor
                int gidCell;
                MB_CHK_ERR( context.MBI->tag_get_data( gidTag, &cov_cell, 1, &gidCell ) );
                assert( gidCell > 0 );
                idsFromProcs[origProc].insert( gidCell );
            }
        }

#ifdef VERBOSE
        std::ofstream dbfile;
        std::stringstream outf;
        outf << "idsFromProc_0" << currentRankInJointComm << ".txt";
        dbfile.open( outf.str().c_str() );
        dbfile << "Writing this to a file.\n";

        dbfile << " map size:" << idsFromProcs.size()
               << std::endl;  // on the receiver side, these show how much data to receive
        // from the sender (how many ids, and how much tag data later; we need to size up the
        // receiver buffer) arrange in tuples , use map iterators to send the ids
        for( std::map< int, std::set< int > >::iterator mt = idsFromProcs.begin(); mt != idsFromProcs.end(); mt++ )
        {
            std::set< int >& setIds = mt->second;
            dbfile << "from id: " << mt->first << " receive " << setIds.size() << " cells \n";
            int counter = 0;
            for( std::set< int >::iterator st = setIds.begin(); st != setIds.end(); st++ )
            {
                int valueID = *st;
                dbfile << " " << valueID;
                counter++;
                if( counter % 10 == 0 ) dbfile << "\n";
            }
            dbfile << "\n";
        }
        dbfile.close();
#endif

        // Create new receiver-side ParCommGraph for coverage-specific communication
        if( nullptr != recvGraph )
        {
            ParCommGraph* newRecvGraph = new ParCommGraph( *recvGraph );  // Copy original structure
            newRecvGraph->set_context_id( *context_id );                  // Differentiate from original graph
            newRecvGraph->SetReceivingAfterCoverage( idsFromProcs );      // Store coverage requirements
            newRecvGraph->set_cover_set( cover_set );                     // Associate with coverage mesh set

            // Store in application's graph map; replace any existing graph for this context.
            auto& pgraphMap = context.appDatas[*pid_migr].pgraph;
            auto existing   = pgraphMap.find( *context_id );
            if( existing != pgraphMap.end() )
            {
                delete existing->second;
                pgraphMap.erase( existing );
            }
            pgraphMap[*context_id] = newRecvGraph;
        }

        // Pack requirements into TupleList for crystal router communication
        // Format: (destination_processor, required_global_id)
        for( std::map< int, std::set< int > >::iterator mit = idsFromProcs.begin(); mit != idsFromProcs.end(); ++mit )
        {
            int procToSendTo       = mit->first;  // Source task that owns these IDs
            std::set< int >& idSet = mit->second;
            for( std::set< int >::iterator sit = idSet.begin(); sit != idSet.end(); ++sit )
            {
                int n = TLcovIDs.get_n();
                TLcovIDs.reserve();
                TLcovIDs.vi_wr[2 * n]     = procToSendTo;  // Destination processor
                TLcovIDs.vi_wr[2 * n + 1] = *sit;          // Required global ID
            }
        }
    }

    // =========================================================================
    // STEP 4: Crystal router exchange - send ID requirements to source tasks
    // =========================================================================

    ProcConfig pc( global );
    // Collective operation: receivers send requirements → senders receive them
    pc.crystal_router()->gs_transfer( 1, TLcovIDs, 0 );

    // =========================================================================
    // STEP 5: On sender tasks, create sender-side graph from received requirements
    // =========================================================================

    if( nullptr != sendGraph )
    {
        appData& dataSrc = context.appDatas[*pid_src];

        // Create new sender-side ParCommGraph for coverage communication
        ParCommGraph* newSendGraph = new ParCommGraph( *sendGraph );  // Copy original structure
        newSendGraph->set_context_id( *context_id );                  // Unique context for this graph
        dataSrc.pgraph[*context_id] = newSendGraph;

        // Process received TLcovIDs to build send patterns: receiver → vector<ids_to_send>
        MB_CHK_ERR( newSendGraph->settle_send_graph( TLcovIDs ) );
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_DumpCommGraph( iMOAB_AppID pid, int* context_id, int* is_sender, const iMOAB_String prefix )
{
    assert( prefix && strlen( prefix ) );

    ParCommGraph* cgraph = context.appDatas[*pid].pgraph[*context_id];
    std::string prefix_str( prefix );

    if( nullptr != cgraph )
        cgraph->dump_comm_information( prefix_str, *is_sender );
    else
    {
        std::cout << " cannot find ParCommGraph on app with pid " << *pid << " name: " << context.appDatas[*pid].name
                  << " context: " << *context_id << "\n";
    }
    return moab::MB_SUCCESS;
}

#endif  // #ifdef MOAB_HAVE_TEMPESTREMAP

#endif  // #ifdef MOAB_HAVE_MPI

#ifdef MOAB_HAVE_TEMPESTREMAP

#ifdef MOAB_HAVE_NETCDF

/**
 * @brief Internal helper: redistributes area values from trivial to arbitrary mesh distribution.
 *
 * @details NetCDF remapping files store areas in "trivial" distribution (contiguous blocks: rank 0 owns
 * IDs 1 to N/P, rank 1 owns N/P+1 to 2N/P, etc.), but MOAB meshes use spatial/graph partitioning
 * with non-contiguous global IDs. This function bridges the gap using two-hop crystal router.
 *
 * @par Algorithm - Two-Hop Rendezvous:
 * <b>TRIVIAL DISTRIBUTION:</b> Rank r owns global IDs [r*nL+1, (r+1)*nL] where nL=N/size @n
 * <b>MESH DISTRIBUTION:</b> Ranks own arbitrary non-contiguous IDs from spatial decomposition
 *
 * -# <b>HOP 1 (REQUESTS):</b> Each rank sends requests to trivial owners for its local element IDs
 *    - For each local element with gid, compute trivial owner: (gid-1)/nL
 *    - Pack tuple: (owner_rank, gid, local_index)
 *    - Crystal router sends requests to trivial owners
 * -# <b>HOP 2 (RESPONSES):</b> Trivial owners look up areas and send back to requesters
 *    - For each request, compute local index: gid - 1 - rank*nL
 *    - Look up area value: trvArea[local_index]
 *    - Pack response: (requester_rank, gid, original_index, area_value)
 *    - Crystal router sends responses back
 *    - Requesters apply area values to mesh elements using original_index
 *
 * @par Serial Optimization:
 * Direct mapping when size==1 (no communication needed)
 *
 * @par Data Structures:
 * - <b>TupleList:</b> Bulk MPI communication (vi_wr/vi_rd for ints, vr_wr/vr_rd for doubles)
 * - <b>Crystal router:</b> Collective all-to-all with O(log P) communication stages
 *
 * @param[in] pid       Application ID for mesh
 * @param[in] N         Total number of elements (global size)
 * @param[in] trvArea   Area values in trivial distribution order
 *                      (size = N/size for most ranks, N/size + N%size for last rank)
 *
 * @pre Global IDs numbered 1 to N (1-indexed)
 * @pre "aream" tag exists on mesh elements
 * @pre trvArea on each rank contains its trivial partition slice
 *
 * @todo Auto-create "aream" tag if missing instead of failing
 * @todo Support variable-length tags, not just single double values
 *
 * @note <b>Performance:</b> O(n_local + log P) time, O(N) communication, 2 crystal router passes
 * @warning <b>Collective operation:</b> All ranks must call (crystal router is collective)
 *
 * @return moab::MB_SUCCESS on success, error code otherwise
 */
static ErrCode set_aream_from_trivial_distribution( iMOAB_AppID pid, int N, std::vector< double >& trvArea )
{
    // =========================================================================
    // STEP 1: Setup - Get parallel context and retrieve tags
    // =========================================================================

    appData& data = context.appDatas[*pid];
    int size = 1, rank = 0;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = context.appDatas[*pid].pcomm;
    size                = pcomm->size();
    rank                = pcomm->rank();
#endif

    // Get "aream" tag for storing area values on mesh elements
    /// @todo Auto-create if missing instead of failing
    Tag areaTag;
    MB_CHK_ERR( context.MBI->tag_get_handle( "aream", areaTag ) );

    // Get local mesh elements and their global IDs
    const Range& ents_to_set = data.primary_elems;
    size_t nents_to_be_set   = ents_to_set.size();

    Tag gidTag = context.MBI->globalId_tag();
    std::vector< int > globalIds( nents_to_be_set );
    MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_to_set, &globalIds[0] ) );

    // =========================================================================
    // STEP 2: Serial path - direct mapping without communication
    // =========================================================================

    const bool serial = ( size == 1 );
    if( serial )
    {
        // Single process: directly map global IDs to trivial array indices
        // Since there's only one rank, trivial and mesh distributions are identical
        for( size_t i = 0; i < nents_to_be_set; i++ )
        {
            int gid        = globalIds[i];  // Global ID (1-indexed)
            int indexInVal = gid - 1;       // Convert to 0-indexed array position
            assert( indexInVal < N );       // Safety: prevent out-of-bounds
            EntityHandle eh = ents_to_set[i];
            // Direct lookup and assignment
            MB_CHK_ERR( context.MBI->tag_set_data( areaTag, &eh, 1, &trvArea[indexInVal] ) );
        }
    }
#ifdef MOAB_HAVE_MPI
    // =========================================================================
    // STEP 3: Parallel path - HOP 1: Build and send requests
    // =========================================================================
    else
    {
        // Compute trivial distribution parameters
        int nL = N / size;  // Base number of elements per rank in trivial distribution

        // Initialize request TupleList: (destination_rank, global_id, local_index)
        /// @todo Extend to support variable-length tags, not just single double values
        TupleList TLreq;
        TLreq.initialize( 3, 0, 0, 0, nents_to_be_set );  // 3 integers per tuple
        TLreq.enableWriteAccess();

        // Build requests: for each local element, determine trivial owner and request area
        for( size_t i = 0; i < nents_to_be_set; i++ )
        {
            int marker  = globalIds[i];         // Global ID we need area for
            int to_proc = ( marker - 1 ) / nL;  // Rank that owns this ID in trivial layout

            // Handle last rank which may own extra elements (N % size)
            if( to_proc == size ) to_proc = size - 1;

            // Pack request tuple
            int n                  = TLreq.get_n();
            TLreq.vi_wr[3 * n]     = to_proc;  // Destination: trivial owner rank
            TLreq.vi_wr[3 * n + 1] = marker;   // Global ID being requested
            TLreq.vi_wr[3 * n + 2] = i;        // Local index (for response routing)
            TLreq.inc_n();
        }

        // Crystal router: send requests to trivial owners (collective operation)
        ( pcomm->proc_config().crystal_router() )->gs_transfer( 1, TLreq, 0 );

        // =====================================================================
        // STEP 4: Process received requests and build responses
        // =====================================================================

        int sizeBack = TLreq.get_n();  // Number of requests received by this rank

        // Initialize response TupleList: (requester_rank, global_id, orig_index, area_value)
        TupleList TLBack;
        TLBack.initialize( 3, 0, 0, 1, sizeBack );  // 3 ints + 1 double per tuple
        TLBack.enableWriteAccess();

        // Process each request: lookup area value and prepare response
        for( int i = 0; i < sizeBack; i++ )
        {
            int from_proc  = TLreq.vi_wr[3 * i];      // Who requested this
            int marker     = TLreq.vi_wr[3 * i + 1];  // Which global ID
            int orig_index = TLreq.vi_wr[3 * i + 2];  // Their local index

            // Compute local index in trvArea for this rank's trivial partition
            // Formula: global_id - 1 - rank_offset where rank_offset = rank * nL
            int index       = marker - 1 - rank * nL;
            double area_val = trvArea[index];  // Lookup area value

            // Pack response tuple
            TLBack.vi_wr[3 * i]     = from_proc;   // Destination: original requester
            TLBack.vi_wr[3 * i + 1] = marker;      // Echo back global ID
            TLBack.vi_wr[3 * i + 2] = orig_index;  // Echo back their local index
            TLBack.vr_wr[i]         = area_val;    // Area value payload
            TLBack.inc_n();
        }

        // =====================================================================
        // STEP 5: Crystal router - send responses back to requesters
        // =====================================================================

        ( pcomm->proc_config().crystal_router() )->gs_transfer( 1, TLBack, 0 );

        // =====================================================================
        // STEP 6: Apply received area values to local mesh elements
        // =====================================================================

        int n1 = TLBack.get_n();  // Number of responses received (should equal nents_to_be_set)

        for( int i = 0; i < n1; i++ )
        {
            // Extract response data
            // int gid         = TLBack.vi_rd[3 * i + 1];  // Global ID (unused, we have orig_index)
            int origIndex   = TLBack.vi_rd[3 * i + 2];  // Local element index from step 3
            double area_val = TLBack.vr_rd[i];          // Area value from trivial owner

            // Get element handle and set area tag
            EntityHandle eh = ents_to_set[origIndex];
            MB_CHK_ERR( context.MBI->tag_set_data( areaTag, &eh, 1, &area_val ) );
        }
    }
#endif

    return MB_SUCCESS;
}

ErrCode iMOAB_LoadMapFile( iMOAB_AppID pid_source,
                           iMOAB_AppID pid_target,
                           iMOAB_AppID pid_intersection,
                           int* srctype,
                           int* tgttype,
                           int* arearead,
                           const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
                           const iMOAB_String remap_weights_filename )
{
    assert( srctype && tgttype );

    // get the local degrees of freedom, from the pid_cpl and type of mesh
    // Get the source and target data and pcomm objects
    appData& data_source = context.appDatas[*pid_source];
    appData& data_target = context.appDatas[*pid_target];
    appData& data_intx   = context.appDatas[*pid_intersection];
    // do we need to check for tempest remap?
    // if we do not have it, do not do anything anyway
    //#ifdef MOAB_HAVE_TEMPESTREMAP
    TempestMapAppData& tdata = data_intx.tempestData;

    // check if the remapped context is null; we need to fix that, if so
    // what if we read a map and compute a map, on a particular iMOAB app a2o for example?
    // we compute an intx map and read a bilinear map
    // this is rather wrong, we need to fix it   FIXME
    if( tdata.remapper == nullptr )
    {
        // do not compute coverage anymore in advance;
        // need to initialize the coverage set creation?
        // or should we really have just one enclosing coverage set for all maps in here ?
        // Now allocate and initialize the remapper object
#ifdef MOAB_HAVE_MPI
        ParallelComm* pco_intx = data_intx.pcomm;
        tdata.remapper         = new moab::TempestRemapper( context.MBI, pco_intx );
#else
        tdata.remapper = new moab::TempestRemapper( context.MBI );
#endif
        tdata.remapper->meshValidate     = true;
        tdata.remapper->constructEdgeMap = true;

        // Do not create new filesets; Use the sets from our respective applications
        tdata.remapper->initialize( false );
        tdata.remapper->GetMeshSet( moab::Remapper::SourceMesh )  = data_source.file_set;
        tdata.remapper->GetMeshSet( moab::Remapper::TargetMesh )  = data_target.file_set;
        tdata.remapper->GetMeshSet( moab::Remapper::OverlapMesh ) = data_intx.file_set;
        // create a unique coverage set; it will be
        moab::EntityHandle covering_set_new;
        MB_CHK_SET_ERR( context.MBI->create_meshset( moab::MESHSET_SET, covering_set_new ), "Can't create new set" );
        tdata.remapper->GetMeshSet( moab::Remapper::CoveringMesh ) = covering_set_new;
    }

    // Setup loading of weights onto TempestOnlineMap
    // Set the context for the remapping weights computation
    tdata.weightMaps[std::string( solution_weights_identifier )] = new moab::TempestOnlineMap( tdata.remapper );

    // Now allocate and initialize the remapper object
    moab::TempestOnlineMap* weightMap = tdata.weightMaps[std::string( solution_weights_identifier )];
    assert( weightMap != nullptr );

    // EntityHandle source_set   = data_source.file_set;
    // EntityHandle covering_set = tdata.remapper->GetMeshSet( Remapper::CoveringMesh );
    EntityHandle target_set = data_target.file_set;  // default: row based partition

    int src_elem_dof_length = 1, tgt_elem_dof_length = 1;  // default=1: FV - element average DoF value
    // tags of interest are either GLOBAL_DOFS (SE) or GLOBAL_ID (FV)
    Tag gdsTag = nullptr;
    if( *srctype == 1 || *tgttype == 1 )
    {
        MB_CHK_ERR( context.MBI->tag_get_handle( "GLOBAL_DOFS", gdsTag ) );
        assert( gdsTag );
    }
    //#endif
    // Find the DoF tag length
    // if it fails, usually it is 16
    // first check if we need to query the source
    if( *srctype == 1 )  // spectral element
        MB_CHK_ERR( context.MBI->tag_get_length( gdsTag, src_elem_dof_length ) );
    // find the DoF tag length
    // next check if we need to query the target
    if( *tgttype == 1 )  // spectral element
        MB_CHK_ERR( context.MBI->tag_get_length( gdsTag, tgt_elem_dof_length ) );

    Tag gidTag = context.MBI->globalId_tag();
    std::vector< int > tgtDofValues;  // srcDofValues,

    // populate first tuple
    // will be filled with entities on coupler, from which we will get the DOFs, based on type
    Range tgt_ents_of_interest;

    if( *tgttype == 1 )  // spectral element
    {
        MB_CHK_ERR( context.MBI->tag_get_length( gdsTag, tgt_elem_dof_length ) );
        MB_CHK_ERR( context.MBI->get_entities_by_type( target_set, MBQUAD, tgt_ents_of_interest ) );
        tgtDofValues.resize( tgt_ents_of_interest.size() * tgt_elem_dof_length );
        MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, tgt_ents_of_interest, &tgtDofValues[0] ) );
    }
    else if( *tgttype == 2 )
    {
        // vertex global ids
        MB_CHK_ERR( context.MBI->get_entities_by_type( target_set, MBVERTEX, tgt_ents_of_interest ) );
        tgtDofValues.resize( tgt_ents_of_interest.size() * tgt_elem_dof_length );
        MB_CHK_ERR( context.MBI->tag_get_data( gidTag, tgt_ents_of_interest, &tgtDofValues[0] ) );
    }
    else if( *tgttype == 3 )  // for FV meshes, just get the global id of cell
    {
        // element global ids
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( target_set, 2, tgt_ents_of_interest ) );
        tgtDofValues.resize( tgt_ents_of_interest.size() * tgt_elem_dof_length );
        MB_CHK_ERR( context.MBI->tag_get_data( gidTag, tgt_ents_of_interest, &tgtDofValues[0] ) );
    }
    else
    {
        MB_CHK_ERR( MB_FAILURE );  // we know only type 1 or 2 or 3
    }

    // pass tgt ordered dofs, and unique
    // we need to read area_b and set aream tag on target cells, too
    std::vector< int > sortTgtDofs( tgtDofValues.begin(), tgtDofValues.end() );
    std::sort( sortTgtDofs.begin(), sortTgtDofs.end() );
    sortTgtDofs.erase( std::unique( sortTgtDofs.begin(), sortTgtDofs.end() ), sortTgtDofs.end() );  // remove duplicates

    ///   the tag should be created already in the e3sm workflow; if not, create it here
    Tag areaTag;
    ErrorCode rval =
        context.MBI->tag_get_handle( "aream", 1, MB_TYPE_DOUBLE, areaTag, MB_TAG_DENSE | MB_TAG_EXCL | MB_TAG_CREAT );
    if( MB_ALREADY_ALLOCATED == rval )
    {
#ifdef MOAB_HAVE_MPI
        if( 0 == data_intx.pcomm->rank() )
#endif
            std::cout << " aream tag already defined \n ";
    }

    std::vector< double > trvAreaA, trvAreaB;  // passed by reference
    int nA, nB;                                // passed by reference, so returned
    MB_CHK_SET_ERR( weightMap->ReadParallelMap( remap_weights_filename, sortTgtDofs, *arearead, trvAreaA, nA, trvAreaB,
                                                nB ),
                    "reading map from disk failed" );
    // trivially distributed areaAs and areaBs will need to be set on their correct source and target cells, as an aream tag

    // if we are on target mesh (row based partition)
    tdata.pid_src  = pid_source;
    tdata.pid_dest = pid_target;

    /// @todo Perform filter based on available source/target nnz data in the map when coverage set exists.
    /// The idea is to look at the source coverage and target meshes and figure out only relevant elements
    /// that need to participate in the meshes.

    //tdata.remapper->SetMeshSet( Remapper::SourceMesh, source_set, &srcc_ents_of_interest );
    // we have read the area A from map file, and we will set it as a aream double tag on the source set, knowing that we
    // read it trivially, with a trivial distribution by the global DOFs
    // local , private method:
    if( 1 == *arearead || 3 == *arearead )
        MB_CHK_SET_ERR( set_aream_from_trivial_distribution( pid_source, nA, trvAreaA ),
                        " fail to set aream on source " );
    if( 2 == *arearead || 3 == *arearead )
        MB_CHK_SET_ERR( set_aream_from_trivial_distribution( pid_target, nB, trvAreaB ),
                        " fail to set aream on target " );

    //tdata.remapper->SetMeshSet( Remapper::CoveringMesh, covering_set, &src_ents_of_interest );
    weightMap->SetSourceNDofsPerElement( src_elem_dof_length );
    //weightMap->set_col_dc_dofs( srcDofValues );  // will set col_dtoc_dofmap

    tdata.remapper->SetMeshSet( Remapper::TargetMesh, target_set, &tgt_ents_of_interest );
    weightMap->SetDestinationNDofsPerElement( tgt_elem_dof_length );
    weightMap->set_row_dc_dofs( tgtDofValues );  // will set row_dtoc_dofmap

    /// @todo Ideally, we should get this metadata from remap_weights_filename and propagate it
    std::string metadataStr = std::string( remap_weights_filename ) + ";FV:1:GLOBAL_ID;FV:1:GLOBAL_ID";

    data_intx.metadataMap[std::string( solution_weights_identifier )] = metadataStr;

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_WriteMapFile( iMOAB_AppID pid_intersection,
                            const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
                            const iMOAB_String remap_weights_filename )
{
    assert( solution_weights_identifier && strlen( solution_weights_identifier ) );
    assert( remap_weights_filename && strlen( remap_weights_filename ) );

    // Get the source and target data and pcomm objects
    appData& data_intx       = context.appDatas[*pid_intersection];
    TempestMapAppData& tdata = data_intx.tempestData;

    // Get the handle to the remapper object
    assert( tdata.remapper != nullptr );

    // Now get online weights object and ensure it is valid
    moab::TempestOnlineMap* weightMap = tdata.weightMaps[std::string( solution_weights_identifier )];
    assert( weightMap != nullptr );

    std::string filename = std::string( remap_weights_filename );

    std::string metadataStr = data_intx.metadataMap[std::string( solution_weights_identifier )];
    std::map< std::string, std::string > attrMap;
    attrMap["title"]         = "MOAB-TempestRemap Online Regridding Weight Generator";
    attrMap["normalization"] = "ovarea";
    attrMap["map_aPb"]       = filename;

    // const std::string delim = ";";
    // size_t pos = 0, index = 0;
    // std::vector< std::string > stringAttr( 3 );
    // // use find() function to get the position of the delimiters
    // while( ( pos = metadataStr.find( delim ) ) != std::string::npos )
    // {
    //     std::string token1 = metadataStr.substr( 0, pos );  // store the substring
    //     if( token1.size() > 0 || index == 0 ) stringAttr[index++] = token1;
    //     std::cout << "\t Found token: "  << token1 << std::endl;
    //     metadataStr.erase(
    //         0, pos + delim.length() ); /* erase() function store the current positon and move to next token. */
    // }
    // std::cout << "\t Found token: " << metadataStr << " -- and index = " << index << std::endl;
    // stringAttr[index] = metadataStr;  // store the last token of the string.
    // assert( index == 2 );
    // attrMap["remap_options"] = stringAttr[0];
    // attrMap["methodorder_b"] = stringAttr[1];
    // attrMap["methodorder_a"] = stringAttr[2];
    attrMap["concave_a"]   = "false";  // defaults
    attrMap["concave_b"]   = "false";  // defaults
    attrMap["bubble"]      = "true";   // defaults
    attrMap["MOABversion"] = std::string( MOAB_PACKAGE_VERSION_STRING );

    // Write the map file to disk in parallel using either HDF5 or SCRIP interface
    MB_CHK_ERR( weightMap->WriteParallelMap( filename, attrMap ) );

    return moab::MB_SUCCESS;
}
//#define VERBOSE
#ifdef MOAB_HAVE_MPI
ErrCode iMOAB_MigrateMapMesh( iMOAB_AppID pid1,
                              iMOAB_AppID pid2,
                              MPI_Comm* jointcomm,
                              MPI_Group* groupA,
                              MPI_Group* groupB,
                              int* type,
                              int* comp1,
                              int* comp2 )
{
    assert( jointcomm );
    assert( groupA );
    assert( groupB );
    bool is_fortran = false;
    if( *pid1 >= 0 ) is_fortran = context.appDatas[*pid1].is_fortran || is_fortran;
    if( *pid2 >= 0 ) is_fortran = context.appDatas[*pid2].is_fortran || is_fortran;

    MPI_Comm joint_communicator =
        ( is_fortran ? MPI_Comm_f2c( *reinterpret_cast< MPI_Fint* >( jointcomm ) ) : *jointcomm );

    int localRank = 0, numProcs = 1;

    // Get the local rank and number of processes in the joint communicator
    MPI_Comm_rank( joint_communicator, &localRank );
    MPI_Comm_size( joint_communicator, &numProcs );

    // each model has a list of global ids that will need to be sent by gs to rendezvous the other
    // model on the joint_communicator
    TupleList TLcomp1;
    TLcomp1.initialize( 2, 0, 0, 0, 0 );  // to proc, marker
    TupleList TLcomp2;
    TLcomp2.initialize( 2, 0, 0, 0, 0 );  // to proc, marker

    // will push_back a new tuple, if needed
    TLcomp1.enableWriteAccess();

    // tags of interest are either GLOBAL_DOFS or GLOBAL_ID
    Tag gdsTag;

    // find the values on first cell
    int lenTagType1 = 1;
    if( *type == 1 )
    {
        MB_CHK_ERR( context.MBI->tag_get_handle( "GLOBAL_DOFS", gdsTag ) );
        MB_CHK_ERR( context.MBI->tag_get_length( gdsTag, lenTagType1 ) );  // usually it is 16
    }
    Tag gidTag = context.MBI->globalId_tag();

    std::vector< int > valuesComp1;

    // populate first tuple
    Range ents_of_interest;  // will be filled with entities on pid1, that need to be distributed,
                             // rearranged in split_ranges map
    if( *pid1 >= 0 )
    {
        appData& data1     = context.appDatas[*pid1];
        EntityHandle fset1 = data1.file_set;

        if( *type == 1 )
        {
            assert( gdsTag );
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset1, MBQUAD, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() * lenTagType1 );
            MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, ents_of_interest, &valuesComp1[0] ) );
        }
        else if( *type == 2 )
        {
            MB_CHK_ERR( context.MBI->get_entities_by_type( fset1, MBVERTEX, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp1[0] ) );  // just global ids
        }
        else if( *type == 3 )  // for FV meshes, just get the global id of cell
        {
            MB_CHK_ERR( context.MBI->get_entities_by_dimension( fset1, 2, ents_of_interest ) );
            valuesComp1.resize( ents_of_interest.size() );
            MB_CHK_ERR( context.MBI->tag_get_data( gidTag, ents_of_interest, &valuesComp1[0] ) );  // just global ids
        }
        else
        {
            MB_CHK_ERR( MB_FAILURE );  // we know only type 1 or 2 or 3
        }
        // now fill the tuple list with info and markers
        // because we will send only the ids, order and compress the list
        std::set< int > uniq( valuesComp1.begin(), valuesComp1.end() );
        TLcomp1.resize( uniq.size() );
        for( std::set< int >::iterator sit = uniq.begin(); sit != uniq.end(); sit++ )
        {
            // to proc, marker, element local index, index in el
            int marker               = *sit;
            int to_proc              = marker % numProcs;
            int n                    = TLcomp1.get_n();
            TLcomp1.vi_wr[2 * n]     = to_proc;  // send to processor
            TLcomp1.vi_wr[2 * n + 1] = marker;
            TLcomp1.inc_n();
        }
    }

    ProcConfig pc( joint_communicator );  // proc config does the crystal router
    pc.crystal_router()->gs_transfer( 1, TLcomp1,
                                      0 );  // communication towards joint tasks, with markers
    // sort by value (key 1)
#ifdef VERBOSE
    std::stringstream ff1;
    ff1 << "TLcomp1_" << localRank << ".txt";
    TLcomp1.print_to_file( ff1.str().c_str() );  // it will append!
#endif
    moab::TupleList::buffer sort_buffer;
    sort_buffer.buffer_init( TLcomp1.get_n() );
    TLcomp1.sort( 1, &sort_buffer );
    sort_buffer.reset();
#ifdef VERBOSE
    // after sorting
    TLcomp1.print_to_file( ff1.str().c_str() );  // it will append!
#endif
    // do the same, for the other component, number2, with type2
    // start copy
    TLcomp2.enableWriteAccess();

    //moab::TempestOnlineMap* weightMap = nullptr;  // declare it outside, but it will make sense only for *pid2 >= 0
    // we know that :) (or *pid2 >= 0, it means we are on the coupler PEs, read map exists, and coupler procs exist)
    // populate second tuple with ids  from read map: we need row_gdofmap and col_gdofmap
    std::vector< int > valuesComp2;
    if( *pid2 >= 0 )  // we are now on coupler, map side
    {
        appData& data2           = context.appDatas[*pid2];
        TempestMapAppData& tdata = data2.tempestData;
        // could be more than one map, read from file
        // get all column dofs, and create a union
        // assert( tdata.weightMaps.size() == 1 );
        // maybe we need to check it is the map we expect
        for( auto mapIt = tdata.weightMaps.begin(); mapIt != tdata.weightMaps.end(); ++mapIt )
        {
            moab::TempestOnlineMap* weightMap = mapIt->second;
            std::vector< int > valueDofs;
            MB_CHK_ERR( weightMap->fill_col_ids( valueDofs ) );
            valuesComp2.insert( valuesComp2.end(), valueDofs.begin(), valueDofs.end() );
        }
        // std::vector<int> ids_of_interest;
        // do a deep copy of the ids of interest: row ids
        // we are interested in col ids, source
        // new method from moab::TempestOnlineMap

        // now fill the tuple list with info and markers
        std::set< int > uniq( valuesComp2.begin(), valuesComp2.end() );
        TLcomp2.resize( uniq.size() );
        for( std::set< int >::iterator sit = uniq.begin(); sit != uniq.end(); sit++ )
        {
            // to proc, marker, element local index, index in el
            int marker               = *sit;
            int to_proc              = marker % numProcs;
            int n                    = TLcomp2.get_n();
            TLcomp2.vi_wr[2 * n]     = to_proc;  // send to processor
            TLcomp2.vi_wr[2 * n + 1] = marker;
            TLcomp2.inc_n();
        }
    }
    pc.crystal_router()->gs_transfer( 1, TLcomp2,
                                      0 );  // communication towards joint tasks, with markers
    // sort by value (key 1)
    // in the rendez-vous approach, markers meet at meet point (rendez-vous) processor marker % nprocs
#ifdef VERBOSE
    std::stringstream ff2;
    ff2 << "TLcomp2_" << localRank << ".txt";
    TLcomp2.print_to_file( ff2.str().c_str() );
#endif
    sort_buffer.buffer_reserve( TLcomp2.get_n() );
    TLcomp2.sort( 1, &sort_buffer );
    sort_buffer.reset();
    // end copy
#ifdef VERBOSE
    TLcomp2.print_to_file( ff2.str().c_str() );
#endif
    // need to send back the info, from the rendezvous point, for each of the values
    /* so go over each value, on local process in joint communicator,

    now have to send back the info needed for communication;
     loop in in sync over both TLComp1 and TLComp2, in local process;
      So, build new tuple lists, to send synchronous communication
      populate them at the same time, based on marker, that is indexed
    */

    TupleList TLBackToComp1;
    TLBackToComp1.initialize( 3, 0, 0, 0, 0 );  // to proc, marker, from proc on comp2,
    TLBackToComp1.enableWriteAccess();

    int n1 = TLcomp1.get_n();
    int n2 = TLcomp2.get_n();

    int indexInTLComp1 = 0;
    int indexInTLComp2 = 0;  // advance both, according to the marker
    if( n1 > 0 && n2 > 0 )
    {

        while( indexInTLComp1 < n1 && indexInTLComp2 < n2 )  // if any is over, we are done
        {
            int currentValue1 = TLcomp1.vi_rd[2 * indexInTLComp1 + 1];
            int currentValue2 = TLcomp2.vi_rd[2 * indexInTLComp2 + 1];
            if( currentValue1 < currentValue2 )
            {
                // we have a big problem; basically, we are saying that
                // dof currentValue is on one model and not on the other
                // std::cout << " currentValue1:" << currentValue1 << " missing in comp2" << "\n";
                indexInTLComp1++;
                continue;
            }
            if( currentValue1 > currentValue2 )
            {
                // std::cout << " currentValue2:" << currentValue2 << " missing in comp1" << "\n";
                indexInTLComp2++;
                continue;
            }
            int size1 = 1;
            int size2 = 1;
            while( indexInTLComp1 + size1 < n1 && currentValue1 == TLcomp1.vi_rd[2 * ( indexInTLComp1 + size1 ) + 1] )
                size1++;
            while( indexInTLComp2 + size2 < n2 && currentValue2 == TLcomp2.vi_rd[2 * ( indexInTLComp2 + size2 ) + 1] )
                size2++;
            // must be found in both lists, find the start and end indices
            for( int i1 = 0; i1 < size1; i1++ )
            {
                for( int i2 = 0; i2 < size2; i2++ )
                {
                    // send the info back to components
                    int n = TLBackToComp1.get_n();
                    TLBackToComp1.reserve();
                    TLBackToComp1.vi_wr[3 * n] =
                        TLcomp1.vi_rd[2 * ( indexInTLComp1 + i1 )];  // send back to the proc marker
                                                                     // came from, info from comp2
                    TLBackToComp1.vi_wr[3 * n + 1] = currentValue1;  // initial value (resend?)
                    TLBackToComp1.vi_wr[3 * n + 2] = TLcomp2.vi_rd[2 * ( indexInTLComp2 + i2 )];  // from proc on comp2
                }
            }
            indexInTLComp1 += size1;
            indexInTLComp2 += size2;
        }
    }
    pc.crystal_router()->gs_transfer( 1, TLBackToComp1, 0 );  // communication towards original tasks, with info about

    TupleList TLv;  // vertices
    TupleList TLc;  // cells if needed (not type 2)

    if( *pid1 >= 0 )
    {
        // we are on original comp 1 tasks
        // before ordering
        // now for each value in TLBackToComp1.vi_rd[3*i+1], on current proc, we know the
        // processors it communicates with
#ifdef VERBOSE
        std::stringstream f1;
        f1 << "TLBack1_" << localRank << ".txt";
        TLBackToComp1.print_to_file( f1.str().c_str() );
#endif
        // so we are now on pid1, we know now each cell where it has to go
        int n = TLBackToComp1.get_n();
        std::map< int, std::set< int > > uniqueIDs;
        for( int i = 0; i < n; i++ )
        {
            int to_proc  = TLBackToComp1.vi_wr[3 * i + 2];
            int globalId = TLBackToComp1.vi_wr[3 * i + 1];
            uniqueIDs[to_proc].insert( globalId );
        }
        // gidTag is gid tag
        // gdsTag is GLOBAL_DOFS , used only for spectral (type 1)
        // (*type 1 is spectral, right now not used in E3SM)

        std::map< int, Range > splits;
        for( size_t i = 0; i < ents_of_interest.size(); i++ )
        {
            EntityHandle ent = ents_of_interest[i];
            for( int j = 0; j < lenTagType1; j++ )
            {
                int marker = valuesComp1[i * lenTagType1 + j];
                for( auto mit = uniqueIDs.begin(); mit != uniqueIDs.end(); mit++ )
                {
                    int proc                = mit->first;
                    std::set< int >& setIds = mit->second;
                    if( setIds.find( marker ) != setIds.end() )
                    {
                        splits[proc].insert( ent );
                    }
                }
            }
        }

        std::map< int, Range > verts_to_proc;
        int numv = 0, numc = 0;
        for( auto it = splits.begin(); it != splits.end(); it++ )
        {
            int to_proc = it->first;
            Range verts;
            if( *type != 2 )
            {
                MB_CHK_ERR( context.MBI->get_connectivity( it->second, verts ) );
                numc += (int)it->second.size();
            }
            else
                verts = it->second;
            verts_to_proc[to_proc] = verts;
            numv += (int)verts.size();
        }
        // first vertices:
        TLv.initialize( 2, 0, 0, 3, numv );  // to proc, GLOBAL ID, 3 real coordinates
        TLv.enableWriteAccess();
        // use the global id of vertices for connectivity
        for( auto it = verts_to_proc.begin(); it != verts_to_proc.end(); it++ )
        {
            int to_proc  = it->first;
            Range& verts = it->second;
            for( Range::iterator vit = verts.begin(); vit != verts.end(); ++vit )
            {
                EntityHandle v   = *vit;
                int n            = TLv.get_n();  // current size of tuple list
                TLv.vi_wr[2 * n] = to_proc;      // send to processor

                MB_CHK_ERR( context.MBI->tag_get_data( gidTag, &v, 1, &( TLv.vi_wr[2 * n + 1] ) ) );
                MB_CHK_ERR( context.MBI->get_coords( &v, 1, &( TLv.vr_wr[3 * n] ) ) );
                TLv.inc_n();  // increment tuple list size
            }
        }
        if( *type != 2 )
        {
            // to proc, ID cell, gdsTag, nbv, id conn,
            int size_tuple =
                2 + ( ( *type != 1 ) ? 0 : lenTagType1 ) + 1 + 10;  // 10 is the max number of vertices in cell

            std::vector< int > gdvals;

            TLc.initialize( size_tuple, 0, 0, 0, numc );  // to proc, GLOBAL ID, 3 real coordinates
            TLc.enableWriteAccess();
            for( auto it = splits.begin(); it != splits.end(); it++ )
            {
                int to_proc  = it->first;
                Range& cells = it->second;
                for( Range::iterator cit = cells.begin(); cit != cells.end(); ++cit )
                {
                    EntityHandle cell         = *cit;
                    int n                     = TLc.get_n();  // current size of tuple list
                    TLc.vi_wr[size_tuple * n] = to_proc;
                    int current_index         = 2;
                    MB_CHK_ERR( context.MBI->tag_get_data( gidTag, &cell, 1, &( TLc.vi_wr[size_tuple * n + 1] ) ) );
                    if( 1 == *type )
                    {
                        MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, &cell, 1,
                                                               &( TLc.vi_wr[size_tuple * n + current_index] ) ) );
                        current_index += lenTagType1;
                    }
                    // now get connectivity
                    const EntityHandle* conn = NULL;
                    int nnodes               = 0;
                    MB_CHK_ERR( context.MBI->get_connectivity( cell, conn, nnodes ) );
                    // fill nnodes:
                    TLc.vi_wr[size_tuple * n + current_index] = nnodes;
                    MB_CHK_ERR( context.MBI->tag_get_data( gidTag, conn, nnodes,
                                                           &( TLc.vi_wr[size_tuple * n + current_index + 1] ) ) );
                    TLc.inc_n();  // increment tuple list size
                }
            }
        }
    }
    else if( *pid2 >= 0 )  // TLv and TLc should be able to receive if *pid2 >= 0
                           // this case will not happen if pid1 and pid2 are both on the coupler side
                           // we need to cover the case if map migrate was used directly
                           // in one hop projection; right now, we prefer 2 hop projection
    {
        TLv.initialize( 2, 0, 0, 3, 0 );  // no vertices here yet, for sure
        TLv.enableWriteAccess();          // to be able to receive stuff, even if nothing is here yet, on this task
        if( *type != 2 )                  // for point cloud, we do not need to initialize TLc (for cells)
        {
            // we still need to initialize the tuples with the right size, as in form_tuples_to_migrate_mesh
            int size_tuple = 2 + ( ( *type != 1 ) ? 0 : lenTagType1 ) + 1 +
                             10;  // 10 is the max number of vertices in cell; kind of arbitrary
            TLc.initialize( size_tuple, 0, 0, 0, 0 );
            TLc.enableWriteAccess();
        }
    }
    pc.crystal_router()->gs_transfer( 1, TLv, 0 );  // communication towards coupler tasks, with mesh vertices
    if( *type != 2 ) pc.crystal_router()->gs_transfer( 1, TLc, 0 );  // those are cells

    if( *pid2 >= 0 )  // will receive the mesh, on coupler pes!, the coverage mesh !
    {
        appData& dataIntx        = context.appDatas[*pid2];
        TempestMapAppData& tdata = dataIntx.tempestData;
        Range primary_ents;                  // vertices for type 2, cells of dim 2 for type 1 or 3
        std::vector< int > values_entities;  // will be the size of primary_ents3 * lenTagType1
        EntityHandle fset3 = tdata.remapper->GetMeshSet( Remapper::CoveringMesh );

        // When an intersection-based covering mesh already exists, ReceiveElementTag and
        // ApplyWeights must share the same entity handles.  Skip populating fset3 with MAP
        // cells so we don't corrupt the set used by both callers.
        Range intx_cov_cells;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( fset3, 2, intx_cov_cells ) );
        const bool has_intersection_coverage = !intx_cov_cells.empty();

        if( !has_intersection_coverage )
        {
            // start copy
            std::map< int, EntityHandle > vertexMap;  //
            Range verts;
            // always form vertices and add them to the fset3;
            int n = TLv.get_n();
            EntityHandle vertex;
            for( int i = 0; i < n; i++ )
            {
                int gid = TLv.vi_rd[2 * i + 1];
                if( vertexMap.find( gid ) == vertexMap.end() )
                {
                    // need to form this vertex
                    MB_CHK_ERR( context.MBI->create_vertex( &( TLv.vr_rd[3 * i] ), vertex ) );
                    vertexMap[gid] = vertex;
                    verts.insert( vertex );
                    MB_CHK_ERR( context.MBI->tag_set_data( gidTag, &vertex, 1, &gid ) );
                }
            }
            MB_CHK_ERR( context.MBI->add_entities( fset3, verts ) );
            if( 2 == *type )
            {
                values_entities.resize( verts.size() );  // just get the ids of vertices
                MB_CHK_ERR( context.MBI->tag_get_data( gidTag, verts, &values_entities[0] ) );
                primary_ents = verts;
                //return MB_SUCCESS;
            }
            else
            {
                n = TLc.get_n();
                int size_tuple =
                    2 + ( ( *type != 1 ) ? 0 : lenTagType1 ) + 1 + 10;  // 10 is the max number of vertices in cell

                EntityHandle new_element;

                std::map< int, EntityHandle >
                    cellMap;  // do not create one if it already exists, maybe from other processes
                for( int i = 0; i < n; i++ )
                {
                    // int from_proc  = TLc.vi_rd[size_tuple * i];
                    int globalIdEl = TLc.vi_rd[size_tuple * i + 1];
                    if( cellMap.find( globalIdEl ) == cellMap.end() )  // need to create the cell
                    {
                        int current_index = 2;
                        if( 1 == *type ) current_index += lenTagType1;
                        int nnodes = TLc.vi_rd[size_tuple * i + current_index];
                        std::vector< EntityHandle > conn;
                        conn.resize( nnodes );
                        for( int j = 0; j < nnodes; j++ )
                        {
                            conn[j] = vertexMap[TLc.vi_rd[size_tuple * i + current_index + j + 1]];
                        }
                        //
                        EntityType entType = MBQUAD;
                        if( nnodes > 4 ) entType = MBPOLYGON;
                        if( nnodes < 4 ) entType = MBTRI;
                        MB_CHK_SET_ERR( context.MBI->create_element( entType, &conn[0], nnodes, new_element ),
                                        "can't create new element " );
                        primary_ents.insert( new_element );
                        cellMap[globalIdEl] = new_element;
                        MB_CHK_SET_ERR( context.MBI->tag_set_data( gidTag, &new_element, 1, &globalIdEl ),
                                        "can't set global id tag on cell " );
                        if( 1 == *type )
                        {
                            // set the gds tag
                            MB_CHK_SET_ERR( context.MBI->tag_set_data( gdsTag, &new_element, 1,
                                                                       &( TLc.vi_rd[size_tuple * i + 2] ) ),
                                            "can't set gds tag on cell " );
                        }
                    }
                }
                MB_CHK_ERR( context.MBI->add_entities( fset3, primary_ents ) );
                if( 1 == *type )
                {
                    values_entities.resize( lenTagType1 * primary_ents.size() );
                    MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, primary_ents, &values_entities[0] ) );
                }
                else  // *type == 3
                {
                    values_entities.resize( primary_ents.size() );  // just get the global ids !
                    MB_CHK_ERR( context.MBI->tag_get_data( gidTag, primary_ents, &values_entities[0] ) );
                }
            }
        }  // end if( !has_intersection_coverage )
        else
        {
            // Reuse the intersection covering mesh as the entity set for weight application;
            // rebuild values_entities from it so col_dtoc_dofmap maps those entities to MAP columns.
            primary_ents = intx_cov_cells;
            if( 1 == *type )
            {
                values_entities.resize( lenTagType1 * primary_ents.size() );
                MB_CHK_ERR( context.MBI->tag_get_data( gdsTag, primary_ents, &values_entities[0] ) );
            }
            else
            {
                values_entities.resize( primary_ents.size() );
                MB_CHK_ERR( context.MBI->tag_get_data( gidTag, primary_ents, &values_entities[0] ) );
            }
        }

        int ndofPerEl = 1;
        if( 1 == *type ) ndofPerEl = (int)( sqrt( lenTagType1 ) );

        tdata.remapper->SetMeshSet( Remapper::CoveringMesh, fset3, &primary_ents );
        // dump covering mesh in a file, to look at it
        // should be one covering mesh per task, should cover the target mesh set
#ifdef VERBOSE
        std::stringstream fcov;
        fcov << "MapCover_" << numProcs << "_" << localRank << ".h5m";
        context.MBI->write_file( fcov.str().c_str(), 0, 0, &fset3, 1 );
        EntityHandle fset2 = tdata.remapper->GetMeshSet( Remapper::TargetMesh );
        std::stringstream ftarg;
        ftarg << "TargMap_" << numProcs << "_" << localRank << ".h5m";
        context.MBI->write_file( ftarg.str().c_str(), 0, 0, &fset2, 1 );
#endif
        for( auto mapIt = tdata.weightMaps.begin(); mapIt != tdata.weightMaps.end(); ++mapIt )
        {
            moab::TempestOnlineMap* weightMap = mapIt->second;
            weightMap->SetSourceNDofsPerElement( ndofPerEl );
            weightMap->set_col_dc_dofs( values_entities );  // will set col_dtoc_dofmap
        }
    }

    // ---------------------------------------------------------------------
    // BfB guardrail: reconcile the loaded maps against the source mesh.
    //
    // A weight matrix column the migrated coverage could not supply means the
    // map references a source DoF that is not in the source mesh. There are
    // two causes, which we must NOT conflate:
    //   (a) the source mesh is a masked subset of the map's source grid
    //       (e.g. a land-only ELM mesh vs a full-grid r05->ne30pg2 map). The
    //       missing columns are GLOBALLY absent (no rank owns them), so the
    //       rendezvous correctly delivered every cell that exists; dropping
    //       them is decomposition-independent and BfB-safe (it matches MCT,
    //       whose land AV has no ocean points either).
    //   (b) the source mesh is complete but a referenced column still went
    //       undelivered -> the map was generated against a different mesh, or
    //       the migration lost a cell. Silently dropping would break BfB, so
    //       we fail loudly here, at setup, instead of deep inside CAAS.
    //
    // We distinguish them with a single integer reduction (total source cells
    // vs the map's declared n_a) -- no global set, no Allgatherv, O(1) memory.
    {
        long localSrcCells = ( *pid1 >= 0 ) ? (long)ents_of_interest.size() : 0;
        long globalSrcCells = 0;
        MPI_Allreduce( &localSrcCells, &globalSrcCells, 1, MPI_LONG, MPI_SUM, joint_communicator );

        int localBadMap = 0, badGid = -1, badNa = -1;
        long droppedCols = 0;
        if( *pid2 >= 0 )
        {
            TempestMapAppData& tdata2 = context.appDatas[*pid2].tempestData;
            for( auto mapIt = tdata2.weightMaps.begin(); mapIt != tdata2.weightMaps.end(); ++mapIt )
            {
                moab::TempestOnlineMap* weightMap = mapIt->second;
                int exGid    = -1;
                int nAbsent  = weightMap->CountAbsentColumns( exGid );
                if( nAbsent <= 0 ) continue;
                const int na = weightMap->GlobalSourceDofCount();
                if( na > 0 && globalSrcCells >= (long)na )
                {
                    // case (b): complete source but a column is missing -> error
                    if( !localBadMap ) { localBadMap = 1; badGid = exGid; badNa = na; }
                }
                else
                {
                    // case (a): masked source -> BfB-safe drop of absent columns
                    droppedCols += weightMap->DropAbsentColumns();
                }
            }
        }

        int globalBadMap = 0;
        MPI_Allreduce( &localBadMap, &globalBadMap, 1, MPI_INT, MPI_MAX, joint_communicator );
        if( globalBadMap )
        {
            if( localBadMap )
                fprintf( stderr,
                         "FATAL: iMOAB_MigrateMapMesh on rank %d: a weight map references source "
                         "DoF GID=%d that the migrated coverage could not supply, yet the source "
                         "mesh is complete (global source cells=%ld >= map n_a=%d). The map was "
                         "generated against a different source mesh, or the migration dropped a "
                         "cell. Refusing to proceed (silently dropping would break BfB).\n",
                         localRank, badGid, globalSrcCells, badNa );
            return moab::MB_FAILURE;
        }

        long globalDropped = 0;
        MPI_Reduce( &droppedCols, &globalDropped, 1, MPI_LONG, MPI_SUM, 0, joint_communicator );
        if( localRank == 0 && globalDropped > 0 )
            std::cout << " iMOAB_MigrateMapMesh: source mesh is a masked subset of the map grid; "
                      << "dropped " << globalDropped << " globally-absent source-column references "
                      << "across all maps and ranks (BfB-safe).\n";
    }

    // compute par comm graph
    int ierr = iMOAB_ComputeCommGraph( pid1, pid2, jointcomm, groupA, groupB, type, type, comp1, comp2 );
    if( ierr != 0 ) return moab::MB_FAILURE;

    return moab::MB_SUCCESS;
}
#undef VERBOSE
#endif  // #ifdef MOAB_HAVE_MPI

#endif  // #ifdef MOAB_HAVE_NETCDF

#define USE_API
static ErrCode ComputeSphereRadius( iMOAB_AppID pid, double* radius )
{
    moab::CartVect pos;

    Range& verts = context.appDatas[*pid].all_verts;
    *radius      = 1.0;
    if( verts.empty() ) return moab::MB_SUCCESS;
    moab::EntityHandle firstVertex = ( verts[0] );

    // coordinate data
    MB_CHK_ERR( context.MBI->get_coords( &( firstVertex ), 1, (double*)&( pos[0] ) ) );

    // compute the distance from origin
    /// @todo We could do this in a loop to verify if the pid represents a spherical mesh
    *radius = pos.length();
    return moab::MB_SUCCESS;
}

ErrCode iMOAB_SetMapGhostLayers( iMOAB_AppID pid, int* n_src_ghost_layers, int* n_tgt_ghost_layers )
{
    appData& data = context.appDatas[*pid];

    // Set the number of ghost layers
    data.tempestData.num_src_ghost_layers = *n_src_ghost_layers;  // number of ghost layers for source
    data.tempestData.num_tgt_ghost_layers = *n_tgt_ghost_layers;  // number of ghost layers for source

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_ComputeCoverageMesh( iMOAB_AppID pid_src, iMOAB_AppID pid_tgt, iMOAB_AppID pid_intx )
{
    // Default constant parameters
    constexpr bool validate        = true;
    constexpr bool meshCleanup     = true;
    bool gnomonic                  = true;
    constexpr double defaultradius = 1.0;
    constexpr double boxeps        = 1.e-10;

    // Other constant parameters
    const double epsrel  = ReferenceTolerance;  // ReferenceTolerance is defined in Defines.h in tempestremap source ;
    double radius_source = 1.0;
    double radius_target = 1.0;

    // Error code definitions
    ErrCode ierr;

    // Get the source and target data and pcomm objects
    appData& data_src  = context.appDatas[*pid_src];
    appData& data_tgt  = context.appDatas[*pid_tgt];
    appData& data_intx = context.appDatas[*pid_intx];
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco_intx = context.appDatas[*pid_intx].pcomm;
#endif

    // Mesh intersection has already been computed; Return early.
    TempestMapAppData& tdata = data_intx.tempestData;
    if( tdata.remapper != nullptr ) return moab::MB_SUCCESS;  // nothing to do

    int rank = 0;
#ifdef MOAB_HAVE_MPI
    if( pco_intx )
    {
        rank = pco_intx->rank();
        MB_CHK_ERR( pco_intx->check_all_shared_handles() );
    }
#endif

    moab::DebugOutput outputFormatter( std::cout, rank, 0 );
    outputFormatter.set_prefix( "[iMOAB_ComputeCoverageMesh]: " );

    ierr = iMOAB_UpdateMeshInfo( pid_src );MB_CHK_ERR( ierr );
    ierr = iMOAB_UpdateMeshInfo( pid_tgt );MB_CHK_ERR( ierr );

    // Rescale the radius of both to compute the intersection
    MB_CHK_ERR( ComputeSphereRadius( pid_src, &radius_source ) );
    MB_CHK_ERR( ComputeSphereRadius( pid_tgt, &radius_target ) );
#ifdef VERBOSE
    if( !rank )
        outputFormatter.printf( 0, "Radius of spheres: source = %12.14f, and target = %12.14f\n", radius_source,
                                radius_target );
#endif

    /* Let make sure that the radius match for source and target meshes. If not, rescale now and
     * unscale later. */
    if( fabs( radius_source - radius_target ) > 1e-10 )
    { /* the radii are different */
        MB_CHK_ERR( IntxUtils::ScaleToRadius( context.MBI, data_src.file_set, defaultradius ) );
        MB_CHK_ERR( IntxUtils::ScaleToRadius( context.MBI, data_tgt.file_set, defaultradius ) );
    }

    // Online workflows use the adaptive default (Van Oosterom-Strackee), which is
    // robust for sliver and polar cells and does not depend on TempestRemap.
    IntxAreaUtils areaAdaptor( IntxAreaUtils::DEFAULT_AREA_METHOD );

    if( meshCleanup )
    {
        // Address issues for source mesh first
        // fixes to enforce positive orientation of the vertices (outward normal)
        MB_CHK_ERR(
            areaAdaptor.positive_orientation( context.MBI, data_src.file_set, defaultradius /*radius_source*/ ) );

        // fixes to clean up any degenerate quadrangular elements present in the mesh (RLL specifically?)
        MB_CHK_ERR( IntxUtils::fix_degenerate_quads( context.MBI, data_src.file_set ) );

        // fixes to enforce convexity in case concave elements are present
        MB_CHK_ERR( moab::IntxUtils::enforce_convexity( context.MBI, data_src.file_set, rank ) );

        // Address issues for target mesh first
        // fixes to enforce positive orientation of the vertices (outward normal)
        MB_CHK_ERR(
            areaAdaptor.positive_orientation( context.MBI, data_tgt.file_set, defaultradius /*radius_target*/ ) );

        // fixes to clean up any degenerate quadrangular elements present in the mesh (RLL specifically?)
        MB_CHK_ERR( IntxUtils::fix_degenerate_quads( context.MBI, data_tgt.file_set ) );

        // fixes to enforce convexity in case concave elements are present
        MB_CHK_ERR( moab::IntxUtils::enforce_convexity( context.MBI, data_tgt.file_set, rank ) );
    }

    // print verbosely about the problem setting
#ifdef VERBOSE
    {
        moab::Range rintxverts, rintxelems;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_src.file_set, 0, rintxverts ) );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_src.file_set, data_src.dimension, rintxelems ) );

        moab::Range bintxverts, bintxelems;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_tgt.file_set, 0, bintxverts ) );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_tgt.file_set, data_tgt.dimension, bintxelems ) );

        if( !rank )
        {
            outputFormatter.printf( 0, "The source set contains %zu vertices and %zu elements\n", rintxverts.size(),
                                    rintxelems.size() );
            outputFormatter.printf( 0, "The target set contains %zu vertices and %zu elements\n", bintxverts.size(),
                                    bintxelems.size() );
        }
    }
#endif
    // use_kdtree_search = ( srctgt_areas_glb[0] < srctgt_areas_glb[1] );
    // advancing front method will fail unless source mesh has no holes (area = 4 * pi) on unit sphere
    // use_kdtree_search = true;

    data_intx.dimension = data_tgt.dimension;
    // set the context for the source and destination applications
    tdata.pid_src  = pid_src;
    tdata.pid_dest = pid_tgt;

    // Now allocate and initialize the remapper object
#ifdef MOAB_HAVE_MPI
    tdata.remapper = new moab::TempestRemapper( context.MBI, pco_intx );
#else
    tdata.remapper = new moab::TempestRemapper( context.MBI );
#endif
    tdata.remapper->meshValidate     = validate;
    tdata.remapper->constructEdgeMap = true;

    // Do not create new filesets; Use the sets from our respective applications
    tdata.remapper->initialize( false );
    tdata.remapper->GetMeshSet( moab::Remapper::SourceMesh )  = data_src.file_set;
    tdata.remapper->GetMeshSet( moab::Remapper::TargetMesh )  = data_tgt.file_set;
    tdata.remapper->GetMeshSet( moab::Remapper::OverlapMesh ) = data_intx.file_set;

    // First, compute the covering source set.
    if( tdata.num_src_ghost_layers >= 1 ) gnomonic = false;  // do not use gnomonic when we need ghost layers;
    MB_CHK_SET_ERR( tdata.remapper->ConstructCoveringSet( epsrel, 1.0, 1.0, boxeps, false, gnomonic,
                                                          tdata.num_src_ghost_layers ),
                    "failed to compute covering set" );

#ifdef MOAB_HAVE_TEMPESTREMAP
    // set the reference to the covering set in the source PID
    data_intx.secondary_file_set = tdata.remapper->GetMeshSet( moab::Remapper::CoveringMesh );
#endif

    return moab::MB_SUCCESS;
}

#ifdef MOAB_HAVE_TEMPESTREMAP
ErrCode iMOAB_WriteCoverageMesh( iMOAB_AppID pid, const iMOAB_String prefix )
{
    IMOAB_CHECKPOINTER( prefix, 2 );
    IMOAB_ASSERT( strlen( prefix ), "Invalid prefix." );

    appData& data        = context.appDatas[*pid];
    EntityHandle fileSet = data.file_set;

    if( data.tempestData.remapper != nullptr )
        fileSet = data.tempestData.remapper->GetMeshSet( Remapper::CoveringMesh );
    else if( data.file_set != data.secondary_file_set )
        fileSet = data.secondary_file_set;
    else
        MB_CHK_SET_ERR( moab::MB_FAILURE, "Invalid secondary file set handle" );

    std::ostringstream file_name;
    int rank = 0, size = 1;

#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = data.pcomm;
    rank                = pcomm->rank();
    size                = pcomm->size();
#endif

    file_name << prefix << "_" << size << "_" << rank << ".h5m";
    // Now let us actually write the file to disk with appropriate options
    MB_CHK_ERR( context.MBI->write_file( file_name.str().c_str(), 0, 0, &fileSet, 1 ) );
    return moab::MB_SUCCESS;
}
#endif

ErrCode iMOAB_ComputeMeshIntersectionOnSphere( iMOAB_AppID pid_src, iMOAB_AppID pid_tgt, iMOAB_AppID pid_intx )
{
    // Default constant parameters
    constexpr bool validate          = true;
    constexpr bool use_kdtree_search = true;
    constexpr double defaultradius   = 1.0;

    // Get the source and target data and pcomm objects
    appData& data_src  = context.appDatas[*pid_src];
    appData& data_tgt  = context.appDatas[*pid_tgt];
    appData& data_intx = context.appDatas[*pid_intx];
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco_intx = context.appDatas[*pid_intx].pcomm;
#endif

    // Mesh intersection has already been computed; Return early.
    TempestMapAppData& tdata = data_intx.tempestData;

    if( tdata.remapper == nullptr )
    {
        // user has not called the coverage mesh computation routine -- so explicitly call it now
        // this check supports the traditional workflow of directly computing mesh intersection
        // and letting this routine compute coverage mesh as needed
        MB_CHK_ERR( iMOAB_ComputeCoverageMesh( pid_src, pid_tgt, pid_intx ) );
    }

    int rank = 0;
#ifdef MOAB_HAVE_MPI
    rank = pco_intx->rank();
#endif
    moab::DebugOutput outputFormatter( std::cout, rank, 0 );
    outputFormatter.set_prefix( "[iMOAB_ComputeMeshIntersectionOnSphere]: " );

    // Next, compute intersections with MOAB.
    // for bilinear, this is an overkill
    MB_CHK_ERR( tdata.remapper->ComputeOverlapMesh( use_kdtree_search, false ) );

    // Mapping computation done
    if( validate )
    {
        IntxAreaUtils areaAdaptor( IntxAreaUtils::DEFAULT_AREA_METHOD );
        double local_areas[3] = { 0.0, 0.0, 0.0 }, global_areas[3] = { 0.0, 0.0, 0.0 };
        local_areas[0] = areaAdaptor.area_on_sphere( context.MBI, data_src.file_set, defaultradius /*radius_source*/ );
        local_areas[1] = areaAdaptor.area_on_sphere( context.MBI, data_tgt.file_set, defaultradius /*radius_target*/ );
        local_areas[2] = areaAdaptor.area_on_sphere( context.MBI, data_intx.file_set, defaultradius );
#ifdef MOAB_HAVE_MPI
        global_areas[0] = global_areas[1] = global_areas[2] = 0.0;
        MPI_Reduce( &local_areas[0], &global_areas[0], 3, MPI_DOUBLE, MPI_SUM, 0, pco_intx->comm() );
#else
        global_areas[0] = local_areas[0];
        global_areas[1] = local_areas[1];
        global_areas[2] = local_areas[2];
#endif

        if( rank == 0 )
        {
            outputFormatter.printf( 0,
                                    "initial area: source mesh = %12.14f, target mesh = %12.14f, "
                                    "overlap mesh = %12.14f\n",
                                    global_areas[0], global_areas[1], global_areas[2] );
            outputFormatter.printf( 0, " relative error w.r.t source = %12.14e, and target = %12.14e\n",
                                    fabs( global_areas[0] - global_areas[2] ) / global_areas[0],
                                    fabs( global_areas[1] - global_areas[2] ) / global_areas[1] );
        }
    }

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_ComputePointDoFIntersection( iMOAB_AppID pid_src, iMOAB_AppID pid_tgt, iMOAB_AppID pid_intx )
{
    ErrCode ierr;

    double radius_source = 1.0;
    double radius_target = 1.0;
    const double epsrel  = ReferenceTolerance;  // ReferenceTolerance is defined in Defines.h in tempestremap source ;
    const double boxeps  = 1.e-8;

    // Get the source and target data and pcomm objects
    appData& data_src  = context.appDatas[*pid_src];
    appData& data_tgt  = context.appDatas[*pid_tgt];
    appData& data_intx = context.appDatas[*pid_intx];
#ifdef MOAB_HAVE_MPI
    ParallelComm* pco_intx = context.appDatas[*pid_intx].pcomm;
#endif

    // Mesh intersection has already been computed; Return early.
    TempestMapAppData& tdata = data_intx.tempestData;
    if( tdata.remapper != nullptr ) return moab::MB_SUCCESS;  // nothing to do

#ifdef MOAB_HAVE_MPI
    if( pco_intx )
    {
        MB_CHK_ERR( pco_intx->check_all_shared_handles() );
    }
#endif

    ierr = iMOAB_UpdateMeshInfo( pid_src );MB_CHK_ERR( ierr );
    ierr = iMOAB_UpdateMeshInfo( pid_tgt );MB_CHK_ERR( ierr );

    // Rescale the radius of both to compute the intersection
    ComputeSphereRadius( pid_src, &radius_source );
    ComputeSphereRadius( pid_tgt, &radius_target );

    IntxAreaUtils areaAdaptor;
    // print verbosely about the problem setting
    {
        moab::Range rintxverts, rintxelems;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_src.file_set, 0, rintxverts ) );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_src.file_set, 2, rintxelems ) );
        MB_CHK_ERR( IntxUtils::fix_degenerate_quads( context.MBI, data_src.file_set ) );
        MB_CHK_ERR( areaAdaptor.positive_orientation( context.MBI, data_src.file_set, radius_source ) );
#ifdef VERBOSE
        std::cout << "The red set contains " << rintxverts.size() << " vertices and " << rintxelems.size()
                  << " elements \n";
#endif

        moab::Range bintxverts, bintxelems;
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_tgt.file_set, 0, bintxverts ) );
        MB_CHK_ERR( context.MBI->get_entities_by_dimension( data_tgt.file_set, 2, bintxelems ) );
        MB_CHK_ERR( IntxUtils::fix_degenerate_quads( context.MBI, data_tgt.file_set ) );
        MB_CHK_ERR( areaAdaptor.positive_orientation( context.MBI, data_tgt.file_set, radius_target ) );
#ifdef VERBOSE
        std::cout << "The blue set contains " << bintxverts.size() << " vertices and " << bintxelems.size()
                  << " elements \n";
#endif
    }

    data_intx.dimension = data_tgt.dimension;
    // set the context for the source and destination applications
    // set the context for the source and destination applications
    tdata.pid_src         = pid_src;
    tdata.pid_dest        = pid_tgt;
    data_intx.point_cloud = ( data_src.point_cloud || data_tgt.point_cloud );
    assert( data_intx.point_cloud == true );

    // Now allocate and initialize the remapper object
#ifdef MOAB_HAVE_MPI
    tdata.remapper = new moab::TempestRemapper( context.MBI, pco_intx );
#else
    tdata.remapper = new moab::TempestRemapper( context.MBI );
#endif
    tdata.remapper->meshValidate     = true;
    tdata.remapper->constructEdgeMap = true;

    // Do not create new filesets; Use the sets from our respective applications
    tdata.remapper->initialize( false );
    tdata.remapper->GetMeshSet( moab::Remapper::SourceMesh )  = data_src.file_set;
    tdata.remapper->GetMeshSet( moab::Remapper::TargetMesh )  = data_tgt.file_set;
    tdata.remapper->GetMeshSet( moab::Remapper::OverlapMesh ) = data_intx.file_set;

    /* Let make sure that the radius match for source and target meshes. If not, rescale now and
     * unscale later. */
    if( fabs( radius_source - radius_target ) > 1e-10 )
    { /* the radii are different */
        MB_CHK_ERR( IntxUtils::ScaleToRadius( context.MBI, data_src.file_set, 1.0 ) );
        MB_CHK_ERR( IntxUtils::ScaleToRadius( context.MBI, data_tgt.file_set, 1.0 ) );
    }

    MB_CHK_ERR( tdata.remapper->ConvertMeshToTempest( moab::Remapper::SourceMesh ) );
    MB_CHK_ERR( tdata.remapper->ConvertMeshToTempest( moab::Remapper::TargetMesh ) );

    // First, compute the covering source set.
    MB_CHK_ERR( tdata.remapper->ConstructCoveringSet( epsrel, 1.0, 1.0, boxeps, false ) );

#ifdef MOAB_HAVE_MPI
    /* VSM: This context should be set on the data_src but it would overwrite the source
       covering set context in case it is coupled to another APP as well.
       This needs a redesign. */
    // data_intx.covering_set = tdata.remapper->GetCoveringSet();
    // data_src.covering_set = tdata.remapper->GetCoveringSet();
#endif

    // Now let us re-convert the MOAB mesh back to Tempest representation
    // MB_CHK_ERR( tdata.remapper->ComputeGlobalLocalMaps() );

    return moab::MB_SUCCESS;
}

ErrCode iMOAB_ComputeScalarProjectionWeights(
    iMOAB_AppID pid_intx,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String disc_method_source,
    int* disc_order_source,
    const iMOAB_String disc_method_target,
    int* disc_order_target,
    const iMOAB_String fv_method,
    int* fNoBubble,
    int* fMonotoneTypeID,
    int* fVolumetric,
    int* fInverseDistanceMap,
    int* fNoConservation,
    int* fValidate,
    const iMOAB_String source_solution_tag_dof_name,
    const iMOAB_String target_solution_tag_dof_name )
{
    assert( disc_order_source && disc_order_target && *disc_order_source > 0 && *disc_order_target > 0 );
    assert( solution_weights_identifier && strlen( solution_weights_identifier ) );
    assert( disc_method_source && strlen( disc_method_source ) );
    assert( disc_method_target && strlen( disc_method_target ) );
    assert( source_solution_tag_dof_name && strlen( source_solution_tag_dof_name ) );
    assert( target_solution_tag_dof_name && strlen( target_solution_tag_dof_name ) );

    // Get the source and target data and pcomm objects
    appData& data_intx       = context.appDatas[*pid_intx];
    TempestMapAppData& tdata = data_intx.tempestData;

    // Setup computation of weights
    // Set the context for the remapping weights computation
    tdata.weightMaps[std::string( solution_weights_identifier )] = new moab::TempestOnlineMap( tdata.remapper );

    // Call to generate the remap weights with the tempest meshes
    moab::TempestOnlineMap* weightMap = tdata.weightMaps[std::string( solution_weights_identifier )];
    assert( weightMap != nullptr );

    GenerateOfflineMapAlgorithmOptions mapOptions;
    mapOptions.nPin           = *disc_order_source;
    mapOptions.nPout          = *disc_order_target;
    mapOptions.fSourceConcave = false;
    mapOptions.fTargetConcave = false;

    mapOptions.strMethod = "";
    if( fv_method ) mapOptions.strMethod += std::string( fv_method ) + ";";
    if( fMonotoneTypeID )
    {
        switch( *fMonotoneTypeID )
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
    }
    else
        mapOptions.fMonotone = false;

    mapOptions.fNoBubble       = ( fNoBubble ? *fNoBubble : false );
    mapOptions.fNoConservation = ( fNoConservation ? *fNoConservation > 0 : false );
    mapOptions.fNoCorrectAreas = false;
    // mapOptions.fNoCheck        = !( fValidate ? *fValidate : true );
    mapOptions.fNoCheck = true;
    if( fVolumetric && *fVolumetric ) mapOptions.strMethod += "volumetric;";
    if( fInverseDistanceMap && *fInverseDistanceMap ) mapOptions.strMethod += "invdist;";

    std::string metadataStr = mapOptions.strMethod + ";" + std::string( disc_method_source ) + ":" +
                              std::to_string( *disc_order_source ) + ":" + std::string( source_solution_tag_dof_name ) +
                              ";" + std::string( disc_method_target ) + ":" + std::to_string( *disc_order_target ) +
                              ":" + std::string( target_solution_tag_dof_name );

    data_intx.metadataMap[std::string( solution_weights_identifier )] = metadataStr;

    // Now let us compute the local-global mapping and store it in the context
    // We need this mapping when computing matvec products and to do reductions in parallel
    // Additionally, the call below will also compute weights with TempestRemap
    MB_CHK_ERR( weightMap->GenerateRemappingWeights(
        std::string( disc_method_source ),            // const std::string& strInputType
        std::string( disc_method_target ),            // const std::string& strOutputType
        mapOptions,                                   // GenerateOfflineMapAlgorithmOptions& mapOptions
        std::string( source_solution_tag_dof_name ),  // const std::string& srcDofTagName = "GLOBAL_ID"
        std::string( target_solution_tag_dof_name )   // const std::string& tgtDofTagName = "GLOBAL_ID"
        ) );

    // print some map statistics
    weightMap->PrintMapStatistics();

    if( fValidate && *fValidate )
    {
        const double dNormalTolerance = 1.0E-8;
        const double dStrictTolerance = 1.0E-12;
        weightMap->CheckMap( true, true, ( fMonotoneTypeID && *fMonotoneTypeID ), dNormalTolerance, dStrictTolerance );
    }

    return moab::MB_SUCCESS;
}

// Forward declaration of internal helper
static ErrCode ComputeRowBounds( iMOAB_AppID pid_intersection,
                                 const iMOAB_String solution_weights_identifier,
                                 const iMOAB_String source_solution_tag_name,
                                 const iMOAB_String lower_bound_tag_name,
                                 const iMOAB_String upper_bound_tag_name );

ErrCode iMOAB_ApplyScalarProjectionWeights(
    iMOAB_AppID pid_intersection,
    int* filter_type,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String source_solution_tag_name,
    const iMOAB_String target_solution_tag_name )
{
    assert( solution_weights_identifier && strlen( solution_weights_identifier ) );
    assert( source_solution_tag_name && strlen( source_solution_tag_name ) );
    assert( target_solution_tag_name && strlen( target_solution_tag_name ) );

    moab::ErrorCode rval;

    // Get the source and target data and pcomm objects
    appData& data_intx       = context.appDatas[*pid_intersection];
    TempestMapAppData& tdata = data_intx.tempestData;

    if( !tdata.weightMaps.count( std::string( solution_weights_identifier ) ) )  // key does not exist
        return moab::MB_INDEX_OUT_OF_RANGE;
    moab::TempestOnlineMap* weightMap = tdata.weightMaps[std::string( solution_weights_identifier )];

    // we assume that there are separators ":" between the tag names
    std::vector< std::string > srcNames;
    std::vector< std::string > tgtNames;
    std::vector< Tag > srcTagHandles;
    std::vector< Tag > tgtTagHandles;
    std::string separator( ":" );
    std::string src_name( source_solution_tag_name );
    std::string tgt_name( target_solution_tag_name );
    split_tag_names( src_name, separator, srcNames );
    split_tag_names( tgt_name, separator, tgtNames );
    if( srcNames.size() != tgtNames.size() )
    {
        std::cout << " error in parsing source and target tag names. \n";
        return moab::MB_FAILURE;
    }

    // get both source and target tag handles
    for( size_t i = 0; i < srcNames.size(); i++ )
    {
        Tag tagHandle;
        // get source tag handles and arrange in order
        rval = context.MBI->tag_get_handle( srcNames[i].c_str(), tagHandle );
        if( MB_SUCCESS != rval || nullptr == tagHandle )
        {
            return moab::MB_TAG_NOT_FOUND;
        }
        srcTagHandles.push_back( tagHandle );

        // get corresponding target handles and arrange in the same order
        rval = context.MBI->tag_get_handle( tgtNames[i].c_str(), tagHandle );
        if( MB_SUCCESS != rval || nullptr == tagHandle )
        {
            return moab::MB_TAG_NOT_FOUND;
        }
        tgtTagHandles.push_back( tagHandle );
    }

#ifdef VERBOSE
    // Now get reference to the remapper object
    moab::TempestRemapper* remapper = tdata.remapper;
    std::vector< double > solSTagVals;
    std::vector< double > solTTagVals;

    moab::Range sents, tents;
    if( data_intx.point_cloud )
    {
        appData& data_src = context.appDatas[*tdata.pid_src];
        appData& data_tgt = context.appDatas[*tdata.pid_dest];
        if( data_src.point_cloud )
        {
            moab::Range& covSrcEnts = remapper->GetMeshVertices( moab::Remapper::CoveringMesh );
            solSTagVals.resize( covSrcEnts.size(), 0. );
            sents = covSrcEnts;
        }
        else
        {
            moab::Range& covSrcEnts = remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
            solSTagVals.resize(
                covSrcEnts.size() * weightMap->GetSourceNDofsPerElement() * weightMap->GetSourceNDofsPerElement(), 0. );
            sents = covSrcEnts;
        }
        if( data_tgt.point_cloud )
        {
            moab::Range& tgtEnts = remapper->GetMeshVertices( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size(), 0. );
            tents = tgtEnts;
        }
        else
        {
            moab::Range& tgtEnts = remapper->GetMeshEntities( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size() * weightMap->GetDestinationNDofsPerElement() *
                                    weightMap->GetDestinationNDofsPerElement(),
                                0. );
            tents = tgtEnts;
        }
    }
    else
    {
        moab::Range& covSrcEnts = remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
        moab::Range& tgtEnts    = remapper->GetMeshEntities( moab::Remapper::TargetMesh );
        solSTagVals.resize(
            covSrcEnts.size() * weightMap->GetSourceNDofsPerElement() * weightMap->GetSourceNDofsPerElement(), -1.0 );
        solTTagVals.resize( tgtEnts.size() * weightMap->GetDestinationNDofsPerElement() *
                                weightMap->GetDestinationNDofsPerElement(),
                            0. );

        sents = covSrcEnts;
        tents = tgtEnts;
    }
#endif

    moab::TempestOnlineMap::CAASType caasType = moab::TempestOnlineMap::CAAS_NONE;
    if( filter_type )
    {
        switch( *filter_type )
        {
            case 1:
                caasType = moab::TempestOnlineMap::CAAS_GLOBAL;
                break;
            case 2:
                caasType = moab::TempestOnlineMap::CAAS_LOCAL;
                break;
            case 3:
                caasType = moab::TempestOnlineMap::CAAS_LOCAL_ADJACENT;
                break;
            default:
                caasType = moab::TempestOnlineMap::CAAS_NONE;
        }
    }

    // Look up the low-order weight map
    std::string lo_weights_identifier =
        ( strlen( solution_weights_identifier ) > 3 ? std::string( solution_weights_identifier + 3 ) : "" );
    moab::TempestOnlineMap* loWeightMap = nullptr;
    bool useDualMapBounds               = false;
    if( lo_weights_identifier.size() && tdata.weightMaps.count( lo_weights_identifier ) )
    {
        // store the reference to the weight map
        loWeightMap = tdata.weightMaps[lo_weights_identifier];
        // Dual-map nonlinear remapping: use low-order map stencil for CAAS bounds
        useDualMapBounds = ( loWeightMap && caasType != moab::TempestOnlineMap::CAAS_NONE );
    }

    if( useDualMapBounds )
    {
        for( size_t i = 0; i < srcTagHandles.size(); i++ )
        {
            // Apply high-order projection with dual-map CAAS bounds from low-order map
            // If caasType == CAAS_NONE, this just returns the low-order projection back
            MB_CHK_ERR(
                weightMap->ApplyWeightsWithDualMap( srcTagHandles[i], tgtTagHandles[i], loWeightMap, caasType ) );
        }

        // Store diagnostic per-row bounds from the HIGH-ORDER stencil on
        // target entities. ApplyWeightsWithDualMap uses the high-order
        // stencil to define [lo,hi] (per Bradley et al. 2019); the test
        // bounds must match for verification to be meaningful.
        for( size_t i = 0; i < srcNames.size(); i++ )
        {
            std::string loBoundName = srcNames[i] + "_DualMapLoBound";
            std::string hiBoundName = srcNames[i] + "_DualMapHiBound";
            MB_CHK_ERR( ComputeRowBounds( pid_intersection, solution_weights_identifier, srcNames[i].c_str(),
                                          loBoundName.c_str(), hiBoundName.c_str() ) );
        }
    }
    else
    {
        for( size_t i = 0; i < srcTagHandles.size(); i++ )
        {
            // Standard: apply projection with optional overlap-based CAAS
            MB_CHK_ERR( weightMap->ApplyWeights( srcTagHandles[i], tgtTagHandles[i], false, caasType ) );
        }
    }

// #define VERBOSE
#ifdef VERBOSE
    ParallelComm* pco_intx = context.appDatas[*pid_intersection].pcomm;

    int ivar = 0;
    {
        Tag ssolnTag = srcTagHandles[ivar];
        std::stringstream sstr;
        sstr << "covsrcTagData_" << *pid_intersection << "_" << ivar << "_" << pco_intx->rank() << ".txt";
        std::ofstream output_file( sstr.str().c_str() );
        for( unsigned i = 0; i < sents.size(); ++i )
        {
            EntityHandle elem = sents[i];
            std::vector< double > locsolSTagVals( 16 );
            MB_CHK_ERR( context.MBI->tag_get_data( ssolnTag, &elem, 1, &locsolSTagVals[0] ) );
            output_file << "\n" << remapper->GetGlobalID( Remapper::CoveringMesh, i ) << "-- \n\t";
            for( unsigned j = 0; j < 16; ++j )
                output_file << locsolSTagVals[j] << " ";
        }
        output_file.flush();  // required here
        output_file.close();
    }
    {
        std::stringstream sstr;
        sstr << "outputSrcDest_" << *pid_intersection << "_" << ivar << "_" << pco_intx->rank() << ".h5m";
        EntityHandle sets[2] = { context.appDatas[*tdata.pid_src].file_set,
                                 context.appDatas[*tdata.pid_dest].file_set };
        MB_CHK_ERR( context.MBI->write_file( sstr.str().c_str(), nullptr, "", sets, 2 ) );
    }
    {
        std::stringstream sstr;
        sstr << "outputCovSrcDest_" << *pid_intersection << "_" << ivar << "_" << pco_intx->rank() << ".h5m";
        // EntityHandle sets[2] = {data_intx.file_set, data_intx.covering_set};
        EntityHandle covering_set = remapper->GetCoveringSet();
        EntityHandle sets[2]      = { covering_set, context.appDatas[*tdata.pid_dest].file_set };
        MB_CHK_ERR( context.MBI->write_file( sstr.str().c_str(), nullptr, "", sets, 2 ) );
    }
    {
        std::stringstream sstr;
        sstr << "covMesh_" << *pid_intersection << "_" << pco_intx->rank() << ".vtk";
        // EntityHandle sets[2] = {data_intx.file_set, data_intx.covering_set};
        EntityHandle covering_set = remapper->GetCoveringSet();
        MB_CHK_ERR( context.MBI->write_file( sstr.str().c_str(), nullptr, "", &covering_set, 1 ) );
    }
    {
        std::stringstream sstr;
        sstr << "tgtMesh_" << *pid_intersection << "_" << pco_intx->rank() << ".vtk";
        // EntityHandle sets[2] = {data_intx.file_set, data_intx.covering_set};
        EntityHandle target_set = remapper->GetMeshSet( Remapper::TargetMesh );
        MB_CHK_ERR( context.MBI->write_file( sstr.str().c_str(), nullptr, "", &target_set, 1 ) );
    }
    {
        std::stringstream sstr;
        sstr << "colvector_" << *pid_intersection << "_" << ivar << "_" << pco_intx->rank() << ".txt";
        std::ofstream output_file( sstr.str().c_str() );
        for( unsigned i = 0; i < solSTagVals.size(); ++i )
            output_file << i << " " << weightMap->col_dofmap[i] << " " << weightMap->col_gdofmap[i] << " "
                        << solSTagVals[i] << "\n";
        output_file.flush();  // required here
        output_file.close();
    }
#endif
    // #undef VERBOSE

    return moab::MB_SUCCESS;
}

static ErrCode ComputeRowBounds( iMOAB_AppID pid_intersection,
                                 const iMOAB_String solution_weights_identifier,
                                 const iMOAB_String source_solution_tag_name,
                                 const iMOAB_String lower_bound_tag_name,
                                 const iMOAB_String upper_bound_tag_name )
{
    assert( solution_weights_identifier && strlen( solution_weights_identifier ) );
    assert( source_solution_tag_name && strlen( source_solution_tag_name ) );
    assert( lower_bound_tag_name && strlen( lower_bound_tag_name ) );
    assert( upper_bound_tag_name && strlen( upper_bound_tag_name ) );

    appData& data_intx       = context.appDatas[*pid_intersection];
    TempestMapAppData& tdata = data_intx.tempestData;

    if( !tdata.weightMaps.count( std::string( solution_weights_identifier ) ) ) return moab::MB_INDEX_OUT_OF_RANGE;
    moab::TempestOnlineMap* weightMap = tdata.weightMaps[std::string( solution_weights_identifier )];

    // Get entity ranges
    moab::TempestRemapper* remapper = tdata.remapper;
    moab::Range covSrcEnts          = remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
    moab::Range tgtEnts             = remapper->GetMeshEntities( moab::Remapper::TargetMesh );

    int srcNDof = weightMap->GetSourceNDofsPerElement();
    int tgtNDof = weightMap->GetDestinationNDofsPerElement();

    // Read source field values from coverage mesh
    Tag srcTag;
    moab::ErrorCode rval = context.MBI->tag_get_handle( source_solution_tag_name, srcTag );
    if( rval != moab::MB_SUCCESS ) return moab::MB_TAG_NOT_FOUND;

    std::vector< double > srcVals( covSrcEnts.size() * srcNDof * srcNDof, 0.0 );
    MB_CHK_ERR( context.MBI->tag_get_data( srcTag, covSrcEnts, &srcVals[0] ) );

    // Get or create lower/upper bound tags on target entities
    size_t nTargetDofs = tgtEnts.size() * tgtNDof * tgtNDof;
    Tag loTag, hiTag;
    MB_CHK_SET_ERR( context.MBI->tag_get_handle( lower_bound_tag_name, tgtNDof * tgtNDof, MB_TYPE_DOUBLE, loTag,
                                                 MB_TAG_DENSE | MB_TAG_CREAT ),
                    "Failed to get or create lower bound tag on target entities" );
    MB_CHK_SET_ERR( context.MBI->tag_get_handle( upper_bound_tag_name, tgtNDof * tgtNDof, MB_TYPE_DOUBLE, hiTag,
                                                 MB_TAG_DENSE | MB_TAG_CREAT ),
                    "Failed to get or create upper bound tag on target entities" );

    // Compute per-row bounds from weight matrix stencil.
    // Matrix column indices are NOT the same as srcVals indices; we must map
    // matrix-col -> source-vector-index via the inverse of col_dtoc_dofmap.
    auto& W                            = weightMap->GetWeightMatrix();
    const std::vector< int >& col_dtoc = weightMap->GetColDofMap();
    const std::vector< int >& row_dtoc = weightMap->GetRowDofMap();
    int maxMatCol                      = -1;
    for( size_t k = 0; k < srcVals.size() && k < col_dtoc.size(); k++ )
        if( col_dtoc[k] > maxMatCol ) maxMatCol = col_dtoc[k];
    std::vector< int > col_inv( maxMatCol + 1, -1 );
    for( size_t k = 0; k < srcVals.size() && k < col_dtoc.size(); k++ )
        if( col_dtoc[k] >= 0 ) col_inv[col_dtoc[k]] = (int)k;

    std::vector< double > loBound( nTargetDofs, 1e308 );
    std::vector< double > hiBound( nTargetDofs, -1e308 );

    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        int r = ( i < row_dtoc.size() ) ? row_dtoc[i] : -1;
        if( r < 0 || r >= W.outerSize() )
        {
            loBound[i] = 0.0;
            hiBound[i] = 0.0;
            continue;
        }
        for( moab::TempestOnlineMap::WeightMatrix::InnerIterator it( W, r ); it; ++it )
        {
            int mc = (int)it.col();
            if( mc < 0 || mc > maxMatCol ) continue;
            int srcIdx = col_inv[mc];
            if( srcIdx < 0 || srcIdx >= (int)srcVals.size() ) continue;
            double v = srcVals[srcIdx];
            if( v < loBound[i] ) loBound[i] = v;
            if( v > hiBound[i] ) hiBound[i] = v;
        }
        if( loBound[i] > hiBound[i] )
        {
            loBound[i] = 0.0;
            hiBound[i] = 0.0;
        }
    }

    // Write bounds to tags
    MB_CHK_SET_ERR( context.MBI->tag_set_data( loTag, tgtEnts, &loBound[0] ),
                    "Failed to set lower bound tag data on target entities" );
    MB_CHK_SET_ERR( context.MBI->tag_set_data( hiTag, tgtEnts, &hiBound[0] ),
                    "Failed to set upper bound tag data on target entities" );

    return moab::MB_SUCCESS;
}

#endif  // MOAB_HAVE_TEMPESTREMAP

#ifdef __cplusplus
}
#endif
