/*
 * This imoab_read_map test will simulate coupling between 2 components
 * 2 meshes will be loaded from 2 files (src, tgt), and one map file
 * after the map is read, in parallel, on coupler pes, with distributed rows, the
 * coupler meshes for source and target will be generated, in a migration step,
 * in which we will migrate from target pes according to row ids, to coupler target mesh,
 *  and from source to coverage mesh mesh on coupler. During this migration, par comm graphs
 *  will be established between source and coupler and target and coupler, which will assist
 *  in field transfer from source to target, through coupler
 *
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
#error mbtempest tool requires MPI configuration
#endif

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>
#include <sstream>

#include "imoab_coupler_utils.hpp"

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

#define COMPUTE_FILE_MAP
// #define COMPUTE_ONLINE_MAP

#if( !defined( COMPUTE_FILE_MAP ) && !defined( COMPUTE_ONLINE_MAP ) )
#error Enable either file-based map (COMPUTE_FILE_MAP) and/or online (COMPUTE_ONLINE_MAP) for coupling
#endif

int main( int argc, char* argv[] )
{
    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    // component ids are unique over all pes, and established in advance;
    int rankInAtmComm = -1, rankInOcnComm = -1, rankInCouComm = -1;

    std::string atmFilename = TestDir + "unittest/srcWithSolnTag.h5m";
    std::string ocnFilename = TestDir + "unittest/outTri15_8.h5m";
    std::string mapFilename = TestDir + "unittest/mapNE20_FV15.nc";  // this is a netcdf file!
    std::string baseline    = TestDir + "unittest/baseline2.txt";

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1, endG2 = numProcesses - 1;

    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout

    // Default: load ATM on 2 proc, ocean on 2,
    // Load map on 2 tasks also, in parallel, distributed by rows (which is very bad actually for ocean mesh, because
    // probably all source cells will be involved in coverage mesh on both tasks

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "ATM mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &ocnFilename );
    opts.addOpt< std::string >( "map_file,w", "map file from source to target", &mapFilename );

    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );

    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    int number_iterations = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &number_iterations );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 1", &no_regression_test );
    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " ATM file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n OCN file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2
                  << "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n";
        if( !no_regression_test )
        {
            std::cout << " check projection against baseline: " << baseline << "\n";
        }
    }

    // Load files on 3 different communicators, groups
    // First groups has task 0, second group tasks 0 and 1
    // Coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    CHECKIERR( create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm ),
               "Cannot create ATM MPI group and communicator" )

    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    CHECKIERR( create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm ),
               "Cannot create OCN MPI group and communicator" )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    CHECKIERR( create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm ),
               "Cannot create cpl MPI group and communicator" )

    // atm_coupler
    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    CHECKIERR( create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm ),
               "Cannot create joint ATM coupler communicator" )

    // ocn_coupler
    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    CHECKIERR( create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm ),
               "Cannot create joint OCN coupler communicator" )

    CHECKIERR( iMOAB_Initialize( argc, argv ),  // not really needed anything from argc, argv, yet; maybe we should
               "Cannot initialize iMOAB" )

    // -1 means it is not initialized
    int cmpAtmAppID = -1, cmpOcnAppID = -1;
    int cplOcnAppID = -1, cplAtmOcnFileAppID = -1, cplAtmOcnMemAppID = -1, cplAtmAppID = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;  // ATM on component PEs
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;  // OCN on component PEs
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // ATM on coupler PEs
    iMOAB_AppID cplOcnPID = &cplOcnAppID;  // OCN on coupler PEs
#ifdef COMPUTE_FILE_MAP
    iMOAB_AppID cplAtmOcnFilePID = &cplAtmOcnFileAppID;  // intx ATM - OCN on coupler PEs (file workflow)
#endif
#ifdef COMPUTE_ONLINE_MAP
    iMOAB_AppID cplAtmOcnMemPID  = &cplAtmOcnMemAppID;   // intx ATM - OCN on coupler PEs (memory workflow)
#endif

    // constant identifiers
    int cmpatm = 1, cmpocn = 2, cplatm = 3, cplocn = 5, atmocnfid = 7, atmocnmid = 8;

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        // ATM on coupler pes
        CHECKIERR( iMOAB_RegisterApplication( "CPLATMFILE", &couComm, &cplatm, cplAtmPID ),
                   "Cannot register ATM over coupler PEs" )
        // OCN on coupler pes
        CHECKIERR( iMOAB_RegisterApplication( "CPLOCNFILE", &couComm, &cplocn, cplOcnPID ),
                   "Cannot register OCN over coupler PEs" )

#ifdef COMPUTE_FILE_MAP
        // now load map between OCNx and ATMx on coupler PEs
        CHECKIERR( iMOAB_RegisterApplication( "ATMOCNFILE", &couComm, &atmocnfid, cplAtmOcnFilePID ),
                   "Cannot register ocn_atm map instance over coupler pes" )
#endif

#ifdef COMPUTE_ONLINE_MAP
        // now create an app to compute map  between OCNx and ATMx on coupler PEs
        CHECKIERR( iMOAB_RegisterApplication( "ATMOCNMEM", &couComm, &atmocnmid, cplAtmOcnMemPID ),
                   "Cannot register ocn_atm map instance over coupler pes" )
#endif
    }

    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        CHECKIERR( iMOAB_RegisterApplication( "ATMCMP", &atmComm, &cmpatm, cmpAtmPID ), "Cannot register ATM App" )
        CHECKIERR( iMOAB_LoadMesh( cmpAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay ), "Cannot load ATM mesh" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        CHECKIERR( iMOAB_RegisterApplication( "OCNCMP", &ocnComm, &cmpocn, cmpOcnPID ), "Cannot register OCN App" )
        CHECKIERR( iMOAB_LoadMesh( cmpOcnPID, ocnFilename.c_str(), readopts.c_str(), &nghlay ), "Cannot load OCN mesh" )
    }

    // migrate mesh from component to coupler

    // --------- ATM and OCN mesh migration ---------
    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;  // use the RCB/geometric partitioner
#endif
    if( atmComm != MPI_COMM_NULL )
    {
        // then send mesh to second coupler pes
        // send to  coupler pes
        CHECKIERR( iMOAB_SendMesh( cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm, &repartitioner_scheme ),
                   "cannot send atmosphere elements to coupler" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        // then send mesh to second coupler pes
        // send to  coupler pes
        CHECKIERR( iMOAB_SendMesh( cmpOcnPID, &ocnCouComm, &couPEGroup, &cplocn, &repartitioner_scheme ),
                   "cannot send ocean elements to coupler" )
    }

    // now, receive mesh, on coupler communicator; first mesh 1, atm
    if( couComm != MPI_COMM_NULL )
    {
        // receive from atmosphere component
        CHECKIERR( iMOAB_ReceiveMesh( cplAtmPID, &atmCouComm, &atmPEGroup, &cmpatm ),
                   "cannot receive atmosphere elements on coupler app" )

        // receive from ocean component
        CHECKIERR( iMOAB_ReceiveMesh( cplOcnPID, &ocnCouComm, &ocnPEGroup, &cmpocn ),
                   "cannot receive ocean elements on coupler app" )
    }

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ), "cannot free buffers used to send atmosphere mesh" )
    }

    // we can now free the sender buffers
    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpOcnPID, &cplocn ), "cannot free buffers used to send ocean mesh" )
    }

    // this model (recMeshOcn.h5m) has mixed meshes in it, we need to repair the comm graph
    // first delete the one created with migration, then compute a new one
    if( ocnCouComm != MPI_COMM_NULL )
    {
        int type = 3;  // type: 1 - SE, 2 - Vertex (point cloud), 3 - Element (FV scalars)
        CHECKIERR( iMOAB_ComputeCommGraph( cmpOcnPID, cplOcnPID, &ocnCouComm, &ocnPEGroup, &couPEGroup, &type, &type,
                                           &cmpocn, &cplocn ),
                   "cannot compute graph between ocn comp and ocn migrated to coupler" )
    }

    // write only for n==1 case
    if( couComm != MPI_COMM_NULL && 1 == number_iterations )
    {
        char outputFileATM[] = "recvAtmMem.h5m";
        CHECKIERR( iMOAB_WriteMesh( cplAtmPID, outputFileATM, fileWriteOptions ),
                   "cannot write second atm mesh after receiving" )
        char outputFileOCN[] = "recvOcnMem.h5m";
        CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileOCN, fileWriteOptions ),
                   "cannot write second atm mesh after receiving" )
    }
    // --------- ATM and OCN mesh migration ---------

    if( couComm != MPI_COMM_NULL )
    {
#ifdef COMPUTE_FILE_MAP
        PUSH_TIMER( "Compute ATM source coverage mesh for OCN (file-based)" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplAtmPID, cplOcnPID, cplAtmOcnFilePID ),
                   "cannot compute source ATM coverage mesh for OCN" )
        POP_TIMER( couComm, rankInCouComm )
#endif
#ifdef COMPUTE_ONLINE_MAP
        PUSH_TIMER( "Compute ATM source coverage mesh for OCN (in-memory)" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplAtmPID, cplOcnPID, cplAtmOcnMemPID ),
                   "cannot compute source ATM coverage mesh for OCN" )
        POP_TIMER( couComm, rankInCouComm )
#endif
    }

#ifdef COMPUTE_FILE_MAP
    if( atmCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from atm comp to coverage mesh;
        // it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        // results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from atm pes to coupler pes will use the new par comm
        // graph, that has more precise info about what to send for ocean cover ; every time, we
        // will use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute OCN coverage graph for ATM mesh" )
        CHECKIERR( iMOAB_CoverageGraphStandalone( &atmCouComm, cmpAtmPID, cplAtmPID, cplAtmOcnFilePID, &cmpatm, &cplatm,
                                                  &cplocn ),
                   "cannot recompute direct coverage graph for ocean" )
        POP_TIMER( atmCouComm, rankInAtmComm )  // hijack this rank
    }
#endif
    // --------- Compute ATM coverage mesh ---------

#ifdef COMPUTE_FILE_MAP
    const std::string map_from_file_identifier = "map-from-file";
    if( couComm != MPI_COMM_NULL )
    {
        // int dummyCpl     = -1;
        // int dummy_rowcol = 0;
        // int dummyType    = 1;
        // CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplAtmOcnFilePID, &dummyCpl, &dummy_rowcol, &dummyType,
        //                                              map_from_file_identifier.c_str(), mapFilename.c_str() ),
        //            "failed to load map file from disk" );
        int col_or_row   = 0; // row-based or target-based
        int disc_type    = 3; // element-based FV
        CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplAtmPID, cplOcnPID, cplAtmOcnFilePID, &col_or_row, &disc_type,
                                                     map_from_file_identifier.c_str(), mapFilename.c_str() ),
                   "failed to load map file from disk" );
    }

//     if( atmCouComm != MPI_COMM_NULL )
//     {
//         int type      = 1;  // quads in source set
//         int direction = 1;  // from source to coupler; will create a mesh on cplAtmPID
//         // because it is like "coverage", context will be cplocn
//         CHECKIERR( iMOAB_MigrateMapMesh( cmpAtmPID, cplAtmOcnFilePID, cplAtmPID, &atmCouComm, &atmPEGroup,
//                                          &couPEGroup, &type, &cmpatm, &cplatm, &direction ),
//                    "failed to migrate mesh for ATM on coupler" );
// #ifdef VERBOSE
//         if( *cplAtmPID >= 0 )
//         {
//             char prefix[] = "atmcov";
//             ierr          = iMOAB_WriteLocalMesh( cplAtmPID, prefix );
//             , "failed to write local mesh" );
//         }
// #endif
//     }

//     if( ocnCouComm != MPI_COMM_NULL )
//     {
//         int type      = 3;  // cells with GLOBAL_ID in ocean / target set
//         int direction = 2;  // from coupler to target; will create a mesh on cplOcnPID
//         // it will be like initial migrate cmpocn <-> cplocn
//         CHECKIERR( iMOAB_MigrateMapMesh( cmpOcnPID, cplAtmOcnFilePID, cplOcnPID, &ocnCouComm, &ocnPEGroup,
//                                          &couPEGroup, &type, &cmpocn, &cplocn, &direction ),
//                    "failed to migrate mesh for OCN on coupler" );
// #ifdef VERBOSE
//         if( *cplOcnPID >= 0 )
//         {
//             char prefix[] = "ocntgt";
//             CHECKIERR( iMOAB_WriteLocalMesh( cplOcnPID, prefix ), "failed to write local ocean mesh" );
//             char outputFileRec[] = "CoupOcn.h5m";
//             CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileRec, fileWriteOptions ),
//                        "failed to write ocean global mesh file" );
//         }
// #endif
//     }
#endif

    int tagTypes = DENSE_DOUBLE;
    int tagIndex[3]; /* OCN, ATM-File, ATM-Mem */
    int disc_orders[2] = { 1, 1 };
    int atmCompNDoFs = disc_orders[0] * disc_orders[0] /* SE */, ocnCompNDoFs = disc_orders[1] * disc_orders[1] /*FV*/;

#ifdef COMPUTE_ONLINE_MAP
    const std::string map_from_mem_identifier = "map-computed-online";
    if( couComm != MPI_COMM_NULL )
    {
        // PUSH_TIMER( "Compute ATM source coverage mesh for OCN" )
        // // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // // ocean partitions check if intx valid, write some h5m intx file
        // CHECKIERR( iMOAB_ComputeCoverageMesh( cplAtmPID, cplOcnPID, cplAtmOcnMemPID ),
        //            "cannot compute source ATM coverage mesh for OCN" )
        // POP_TIMER( couComm, rankInCouComm )

        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplAtmOcnMemPID ),
                   "cannot compute intersection for ATM/OCN" )
        POP_TIMER( couComm, rankInCouComm )
    }

    const std::string disc_methods[2]  = { "fv", "fv" };
    const std::string dof_tag_names[2] = { "GLOBAL_ID", "GLOBAL_ID" };
    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap" )
        CHECKIERR(  iMOAB_ComputeScalarProjectionWeights( cplAtmOcnMemPID, map_from_mem_identifier.c_str(),
                                                     disc_methods[0].c_str(), &disc_orders[0], disc_methods[1].c_str(),
                                                     &disc_orders[1], nullptr, &fNoBubble, &fMonotoneTypeID,
                                                     &fVolumetric, &fInverseDistanceMap, &fNoConserve, &fValidate,
                                                     dof_tag_names[0].c_str(), dof_tag_names[1].c_str() ),
         "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )

        {
            const std::string atmocn_map_file_name = "atm_ocn_map_computed.nc";
            CHECKIERR( iMOAB_WriteMappingWeightsToFile( cplAtmOcnMemPID, map_from_mem_identifier.c_str(),
                                                        atmocn_map_file_name.c_str() ),
                       "failed to write map file to disk" );
        }
    }

    // now create the linkage between the ATM component and the coverage source mesh in the coupler
    if( atmCouComm != MPI_COMM_NULL )
    {
        /*
         * This new graph will be used for sending data from ATM-cmp (cmpAtmPID) to ATM-coverage-cpl (cplAtmPID) application
         * The operation involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
         * results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
         * after this, the sending of tags from atm pes to coupler pes will use the new par comm
         * graph, that has more precise info about what to send for ocean cover ; every time, we
         * will use the element global id, which should uniquely identify the element
         */
        PUSH_TIMER( "Compute OCN coverage graph for ATM mesh" )
        CHECKIERR( iMOAB_CoverageGraph( &atmCouComm, cmpAtmPID, cplAtmPID, cplAtmOcnMemPID, &cmpatm, &cplatm, &cplocn ),
                   "cannot recompute direct coverage graph for ocean" )
        POP_TIMER( atmCouComm, rankInAtmComm )  // hijack this rank
    }

#endif

    int filter_type    = 0;
    const char* bottomFields = "AnalyticalSolnSrcExact";
#ifdef COMPUTE_FILE_MAP
    const char* bottomProjectedFieldsF = "Target_projF";
#endif
#ifdef COMPUTE_ONLINE_MAP
    const char* bottomProjectedFieldsM = "Target_projM";
#endif

    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, bottomFields, &tagTypes, &atmCompNDoFs, &tagIndex[0] ),
                   "failed to define the field tags AnalyticalSolnSrcExact" );

#ifdef COMPUTE_FILE_MAP
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFieldsF, &tagTypes, &ocnCompNDoFs, &tagIndex[1] ),
                   "failed to define the field tags bottomProjectedFieldsF" );
#endif
#ifdef COMPUTE_ONLINE_MAP
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFieldsM, &tagTypes, &ocnCompNDoFs, &tagIndex[2] ),
                   "failed to define the field tags bottomProjectedFieldsM" );
#endif
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
#ifdef COMPUTE_FILE_MAP
        CHECKIERR( iMOAB_DefineTagStorage( cmpOcnPID, bottomProjectedFieldsF, &tagTypes, &ocnCompNDoFs, &tagIndexIn2 ),
                   "failed to define the field tag for receiving back the tags "
                   "bottomProjectedFieldsF on OCN pes" );
#endif
#ifdef COMPUTE_ONLINE_MAP
        CHECKIERR( iMOAB_DefineTagStorage( cmpOcnPID, bottomProjectedFieldsM, &tagTypes, &ocnCompNDoFs, &tagIndexIn2 ),
                   "failed to define the field tag for receiving back the tags "
                   "bottomProjectedFieldsM on OCN pes" );
#endif
    }

    // start a virtual loop for number of iterations
    for( int iters = 0; iters < number_iterations; iters++ )
    {
        PUSH_TIMER( "Send/receive data from ATM component to coupler in OCN context" )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            CHECKIERR( iMOAB_SendElementTag( cmpAtmPID, bottomFields, &atmCouComm, &cplocn ), "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on ATM on coupler pes, that was redistributed according to coverage
            CHECKIERR( iMOAB_ReceiveElementTag( cplAtmPID, bottomFields, &atmCouComm, &cplocn ),
                       "cannot receive tag values" )
        }

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ),
                       "cannot free buffers used to resend ATM tag towards the coverage mesh" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

// #ifdef VERBOSE
        if( *cplAtmPID >= 0 && number_iterations == 1 )
        {
            char prefix[] = "atmcov_withdata";
            CHECKIERR( iMOAB_WriteLocalMesh( cplAtmPID, prefix ), "failed to write local ATM cov mesh with data" );
        }

        if( couComm != MPI_COMM_NULL && 1 == number_iterations )
        {
            // write only for n==1 case
            char outputFileRecvd[] = "cplAtmFile.h5m";
            CHECKIERR( iMOAB_WriteMesh( cplAtmPID, outputFileRecvd, fileWriteOptions ),
                       "could not write cplAtmFile.h5m to disk" )
        }
// #endif

        std::cout << "Finished sending/receiving the data from ATM to cpl-ATM" << std::endl;

        if( couComm != MPI_COMM_NULL )
        {
#ifdef COMPUTE_FILE_MAP
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply from file scalar projection weights" )
            CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplAtmOcnFilePID, &filter_type,
                                                           map_from_file_identifier.c_str(), bottomFields,
                                                           bottomProjectedFieldsF ),
                       "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
#endif

            std::cout << "Applying projection of the data from ATM to OCN on coupler is done.." << std::endl;

#ifdef COMPUTE_ONLINE_MAP
            PUSH_TIMER( "Apply in-memory scalar projection weights" )
            CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplAtmOcnMemPID, &filter_type,
                                                           map_from_mem_identifier.c_str(), bottomFields,
                                                           bottomProjectedFieldsM ),
                       "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
#endif
        }

        // send the projected tag back to ocean pes, with send/receive tag

        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm ocn_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            // need to use ocean comp id for context
#ifdef COMPUTE_FILE_MAP
            CHECKIERR( iMOAB_SendElementTag( cplOcnPID, bottomProjectedFieldsF, &ocnCouComm, &cmpocn ),
                       "cannot send tag values back to ocean pes" )
            {
                // write only for n==1 case
                char outputFileRecvd[] = "cplProjectedOCNFileF.h5m";
                CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileRecvd, fileWriteOptions ),
                           "could not write cplProjectedOCNFile.h5m to disk" )
            }
#endif
#ifdef COMPUTE_ONLINE_MAP
            CHECKIERR( iMOAB_SendElementTag( cplOcnPID, bottomProjectedFieldsM, &ocnCouComm, &cmpocn ),
                       "cannot send tag values back to ocean pes" )
            {
                // write only for n==1 case
                char outputFileRecvd[] = "cplProjectedOCNFileM.h5m";
                CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileRecvd, fileWriteOptions ),
                           "could not write cplProjectedOCNFile.h5m to disk" )
            }
#endif
        }

        // receive on component 2, ocean
        if( ocnComm != MPI_COMM_NULL )
        {
#ifdef COMPUTE_FILE_MAP
            CHECKIERR( iMOAB_ReceiveElementTag( cmpOcnPID, bottomProjectedFieldsF, &ocnCouComm, &cplocn ),
                       "cannot receive tag values from ocean mesh on coupler pes" )
#endif
#ifdef COMPUTE_ONLINE_MAP
            CHECKIERR( iMOAB_ReceiveElementTag( cmpOcnPID, bottomProjectedFieldsM, &ocnCouComm, &cplocn ),
                       "cannot receive tag values from ocean mesh on coupler pes" )
#endif
        }

        if( couComm != MPI_COMM_NULL )
        {
#ifdef COMPUTE_FILE_MAP
            CHECKIERR( iMOAB_FreeSenderBuffers( cplOcnPID, &cmpocn ), "Freeing buffers failed" )
#endif
#ifdef COMPUTE_ONLINE_MAP
            CHECKIERR( iMOAB_FreeSenderBuffers( cplOcnPID, &cmpocn ), "Freeing buffers failed" )
#endif
        }

        if( ocnComm != MPI_COMM_NULL && 1 == number_iterations )  // write only for n==1 case
        {
// #ifdef VERBOSE
            char outputFileOcnFile[] = "OcnWithProjection.h5m";
            CHECKIERR( iMOAB_WriteMesh( cmpOcnPID, outputFileOcnFile, fileWriteOptions ),
                       "could not write OcnWithProjection.h5m to disk" )
            // #endif
            // test results only for number_iterations== 1
            if( !no_regression_test )
            {
                // get global id storage
                const std::string GidStr = "GLOBAL_ID";  // hard coded too
                int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
                // first get GlobalIds from OCN, and fields:
                int nverts[3], nelem[3];
                std::vector< int > gidElems;
                std::vector< double > tempElems;
                int err_code = 1, ent_type = 1;

                CHECKIERR( iMOAB_DefineTagStorage( cmpOcnPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd ),
                           "failed to define global id tag" );

#ifdef COMPUTE_FILE_MAP
                CHECKIERR( iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 ), "failed to get OCN mesh info" );
                gidElems.resize( nelem[2] );
                tempElems.resize( nelem[2] );

                CHECKIERR( iMOAB_GetIntTagStorage( cmpOcnPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] ),
                           "failed to get global ids" );
                CHECKIERR( iMOAB_GetDoubleTagStorage( cmpOcnPID, bottomProjectedFieldsF, &nelem[2], &ent_type, &tempElems[0] ),
                           "failed to get bottomProjectedFieldsF field" );

                // check against the baseline
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn (file-based map projection) on ocean task "
                              << rankInOcnComm << "\n";
#endif

#ifdef COMPUTE_ONLINE_MAP
                CHECKIERR( iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 ), "failed to get OCN mesh info" );
                gidElems.resize( nelem[2] );
                tempElems.resize( nelem[2] );

                CHECKIERR( iMOAB_GetIntTagStorage( cmpOcnPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] ),
                           "failed to get global ids" );
                CHECKIERR( iMOAB_GetDoubleTagStorage( cmpOcnPID, bottomProjectedFieldsM, &nelem[2], &ent_type, &tempElems[0] ),
                           "failed to get bottomProjectedFieldsM field" );
                // check against the baseline
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn (in-memory map projection) on ocean task "
                              << rankInOcnComm << "\n";
#endif
            }
        }

    }  // end loop iterations

#ifdef COMPUTE_ONLINE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmOcnMemPID ), "cannot deregister app intx AO" )
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmPID ), "cannot deregister app ATMX" )
        CHECKIERR( iMOAB_DeregisterApplication( cplOcnPID ), "cannot deregister app OCNX" )
    }
#endif

#ifdef COMPUTE_FILE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmOcnFilePID ), "cannot deregister app intx AO" )
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmPID ), "cannot deregister app ATMX" )
        CHECKIERR( iMOAB_DeregisterApplication( cplOcnPID ), "cannot deregister app OCNX" )
    }
#endif

    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpOcnPID ), "cannot deregister app OCN" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpAtmPID ), "cannot deregister app ATM" )
    }

    //#endif
    CHECKIERR( iMOAB_Finalize(), "did not finalize iMOAB" )

    // free ATM coupler group and comm
    if( MPI_COMM_NULL != atmCouComm ) MPI_Comm_free( &atmCouComm );
    MPI_Group_free( &joinAtmCouGroup );
    if( MPI_COMM_NULL != atmComm ) MPI_Comm_free( &atmComm );

    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    // free OCN - coupler group and comm
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );

    MPI_Group_free( &ocnPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
