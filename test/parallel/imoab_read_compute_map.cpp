/*
 * This imoab_read_compute_map test will simulate coupling between 2 components
 * 2 meshes will be loaded from 2 files (src, tgt), and one map file
 * We will compute the map FV-FV between the meshes and we will also read the map
 * we will compare both workflows against the baseline test
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
#include <iomanip>
#include <sstream>

#include "imoab_coupler_utils.hpp"

// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

#define COMPUTE_FILE_MAP
#define COMPUTE_TRANSPOSE_FILE_MAP
#define COMPUTE_ONLINE_MAP

#if ( !defined( COMPUTE_FILE_MAP ) && !defined( COMPUTE_TRANSPOSE_FILE_MAP ) && !defined( COMPUTE_ONLINE_MAP ) )
#error Enable either file-based map (COMPUTE_FILE_MAP/COMPUTE_TRANSPOSE_FILE_MAP) and/or online (COMPUTE_ONLINE_MAP) for coupling
#endif

#define VERBOSE
int main( int argc, char* argv[] )
{
    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    const iMOAB_String readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    constexpr bool skip_apply_a2o = false;

    // component ids are unique over all pes, and established in advance;
    int rankInAtmComm = -1, rankInOcnComm = -1, rankInCouComm = -1;

    std::string atmFilename      = TestDir + "unittest/srcWithSolnTag.h5m";
    std::string ocnFilename      = TestDir + "unittest/outTri15_8.h5m";
    std::string mapFilename      = TestDir + "unittest/mapNE20_FV15.nc";    // this is a netcdf file!
    std::string mapFilenameTrans = TestDir + "unittest/mapNE20_FV15_T.nc";  // this is a netcdf file!
    std::string baseline         = TestDir + "unittest/baseline2.txt";
    std::string baselineATM      = TestDir + "unittest/baseline_atm.txt";

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
    opts.addOpt< std::string >( "transpose_map_file,x", "transposed map file from target to source",
                                &mapFilenameTrans );

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
                  << "\n map file:" << mapFilename.c_str() << "\n     on tasks : " << startG4 << ":" << endG4 << "\n";
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

    // constant identifiers
    int cmpatm = 1, cmpocn = 2, cplatm = 3, cplocn = 5;
    // -1 means it is not initialized
    int cmpAtmAppID = -1, cmpOcnAppID = -1;
    int cplOcnAppID = -1, cplAtmAppID = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;  // ATM on component PEs
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;  // OCN on component PEs
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // ATM on coupler PEs
    iMOAB_AppID cplOcnPID = &cplOcnAppID;  // OCN on coupler PEs
#ifdef COMPUTE_FILE_MAP
    int atmocnfid                = 7;
    int cplAtmOcnFileAppID       = -1;
    iMOAB_AppID cplAtmOcnFilePID = &cplAtmOcnFileAppID;  // intx ATM - OCN on coupler PEs (file workflow)
#endif
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
    int ocnatmfid                = 8;
    int cplOcnAtmFileAppID       = -1;
    iMOAB_AppID cplOcnAtmFilePID = &cplOcnAtmFileAppID;  // intx ATM - OCN on coupler PEs (file workflow)
#endif
#ifdef COMPUTE_ONLINE_MAP
    int atmocnmid               = 8;
    int cplAtmOcnMemAppID       = -1;
    iMOAB_AppID cplAtmOcnMemPID = &cplAtmOcnMemAppID;  // intx ATM - OCN on coupler PEs (memory workflow)
#endif

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        // ATM on coupler pes
        CHECKIERR( iMOAB_RegisterApplication( "CPLATM", &couComm, &cplatm, cplAtmPID ),
                   "Cannot register ATM over coupler PEs" )
        // OCN on coupler pes
        CHECKIERR( iMOAB_RegisterApplication( "CPLOCN", &couComm, &cplocn, cplOcnPID ),
                   "Cannot register OCN over coupler PEs" )

#ifdef COMPUTE_FILE_MAP
        // now load map between OCNx and ATMx on coupler PEs
        CHECKIERR( iMOAB_RegisterApplication( "ATMOCNFILE", &couComm, &atmocnfid, cplAtmOcnFilePID ),
                   "Cannot register ocn_atm map instance over coupler pes" )
#endif
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
        // now load map between OCNx and ATMx on coupler PEs
        CHECKIERR( iMOAB_RegisterApplication( "OCNATMFILE", &couComm, &ocnatmfid, cplOcnAtmFilePID ),
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
        CHECKIERR( iMOAB_LoadMesh( cmpAtmPID, atmFilename.c_str(), readopts, &nghlay ), "Cannot load ATM mesh" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        CHECKIERR( iMOAB_RegisterApplication( "OCNCMP", &ocnComm, &cmpocn, cmpOcnPID ), "Cannot register OCN App" )
        CHECKIERR( iMOAB_LoadMesh( cmpOcnPID, ocnFilename.c_str(), readopts, &nghlay ), "Cannot load OCN mesh" )
    }

    // migrate mesh from component to coupler

    // --------- ATM and OCN mesh migration ---------
    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;  // use the RCB/geometric partitioner
#endif
    if( atmComm != MPI_COMM_NULL )
    {
        // then send mesh to second coupler pes send to  coupler pes
        CHECKIERR( iMOAB_SendMesh( cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm, &repartitioner_scheme ),
                   "cannot send atmosphere elements to coupler" )
    }

    // now, receive mesh, on coupler communicator; first mesh 1, atm
    if( couComm != MPI_COMM_NULL )
    {
        // receive from atmosphere component
        CHECKIERR( iMOAB_ReceiveMesh( cplAtmPID, &atmCouComm, &atmPEGroup, &cmpatm ),
                   "cannot receive atmosphere elements on coupler app" )
    }

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ), "cannot free buffers used to send atmosphere mesh" )
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
        // receive from ocean component
        CHECKIERR( iMOAB_ReceiveMesh( cplOcnPID, &ocnCouComm, &ocnPEGroup, &cmpocn ),
                   "cannot receive ocean elements on coupler app" )
    }

    // we can now free the sender buffers
    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpOcnPID, &cplocn ), "cannot free buffers used to send ocean mesh" )
    }

    // write only for n==1 case
    if( couComm != MPI_COMM_NULL && 1 == number_iterations )
    {
        const iMOAB_String outputFileATM = "recvAtmMem.h5m";
        CHECKIERR( iMOAB_WriteMesh( cplAtmPID, outputFileATM, fileWriteOptions ),
                   "cannot write second atm mesh after receiving" )
        char outputFileOCN[] = "recvOcnMem.h5m";
        CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileOCN, fileWriteOptions ),
                   "cannot write second atm mesh after receiving" )
    }

    // --------- ATM and OCN mesh migration ---------
    if( couComm != MPI_COMM_NULL )
    {
#if defined( COMPUTE_ONLINE_MAP )
        PUSH_TIMER( "Compute ATM source coverage mesh for OCN (in-memory)" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplAtmPID, cplOcnPID, cplAtmOcnMemPID ),
                   "cannot compute source ATM coverage mesh for OCN" )
        POP_TIMER( couComm, rankInCouComm )
#endif
#ifdef COMPUTE_FILE_MAP
        PUSH_TIMER( "Compute ATM source coverage mesh for OCN (file-based)" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplAtmPID, cplOcnPID, cplAtmOcnFilePID ),
                   "cannot compute source ATM coverage mesh for OCN" )
        POP_TIMER( couComm, rankInCouComm )
        if( rankInCouComm == 0 )
        {
            CHECKIERR( iMOAB_WriteMesh( cplAtmPID, "atm.h5m", "" ), "cannot write ATM coverage mesh for OCN" )
            CHECKIERR( iMOAB_WriteCoverageMesh( cplAtmOcnFilePID, "atm_ocn_coverage.h5m", "" ),
                       "cannot write ATM coverage mesh for OCN" )
        }
#endif
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
        PUSH_TIMER( "Compute OCN source coverage mesh for ATM (file-based)" )
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplOcnPID, cplAtmPID, cplOcnAtmFilePID ),
                   "cannot compute source OCN coverage mesh for ATM" )
        POP_TIMER( couComm, rankInCouComm )
        if( rankInCouComm == 0 )
        {
            CHECKIERR( iMOAB_WriteMesh( cplOcnPID, "ocn.h5m", "" ), "cannot write OCN coverage mesh for ATM" )
            CHECKIERR( iMOAB_WriteCoverageMesh( cplOcnAtmFilePID, "ocn_atm_coverage.h5m", "" ),
                       "cannot write OCN coverage mesh for ATM" )
        }
#endif
    }

    // --------- Load map from disk or compute it online ---------
    const iMOAB_String map_from_file_identifier[2] = { "atm-ocn-file-map", "ocn-atm-file-map" };
#ifdef COMPUTE_FILE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        int src_disc_type = 3;  // element-based FV
        int tgt_disc_type = 3;  // element-based FV
        CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplAtmPID, cplOcnPID, cplAtmOcnFilePID, &src_disc_type,
                                                     &tgt_disc_type, map_from_file_identifier[0], mapFilename.c_str() ),
                   "failed to load ATM-OCN map file from disk" );

        // {
        //     const iMOAB_String atmocn_map_file_name = "atm_ocn_map_loaded.nc";
        //     CHECKIERR( iMOAB_WriteMappingWeightsToFile( cplAtmOcnFilePID, map_from_file_identifier[0],
        //                                                 atmocn_map_file_name ),
        //                "failed to write map file to disk" );
        // }

        int meshtype = 3;
        PUSH_TIMER( "Compute ATM coverage graph for OCN mesh" )
        CHECKIERR( iMOAB_ComputeCommGraph( cplAtmPID, cplAtmOcnFilePID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                           &meshtype, &cplatm, &atmocnfid ),
                   "cannot recompute ATM source coverage graph for OCN" )
        POP_TIMER( couComm, rankInCouComm )  // hijack this rank
    }
#endif
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        int src_disc_type = 3;  // element-based FV
        int tgt_disc_type = 3;  // element-based FV
        CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplOcnPID, cplAtmPID, cplOcnAtmFilePID, &src_disc_type,
                                                     &tgt_disc_type, map_from_file_identifier[1],
                                                     mapFilenameTrans.c_str() ),
                   "failed to load OCN-ATM map file from disk" );

        // {
        //     const iMOAB_String atmocn_map_file_name = "ocn_atm_map_loaded.nc";
        //     CHECKIERR( iMOAB_WriteMappingWeightsToFile( cplOcnAtmFilePID, map_from_file_identifier[1],
        //                                                 atmocn_map_file_name ),
        //                "failed to write map file to disk" );
        // }

        int meshtype = 3;
        PUSH_TIMER( "Compute OCN coverage graph for ATM mesh" )
        CHECKIERR( iMOAB_ComputeCommGraph( cplOcnPID, cplOcnAtmFilePID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                           &meshtype, &cplocn, &ocnatmfid ),
                   "cannot recompute OCN source coverage graph for ATM" )
        POP_TIMER( couComm, rankInCouComm )  // hijack this rank
    }
#endif

    int tagTypes = DENSE_DOUBLE;
    int tagIndex[3]; /* OCN, ATM-File, ATM-Mem */
    int disc_orders[2] = { 1, 1 };
    int atmCompNDoFs   = disc_orders[0] * disc_orders[0] /* FV */,
        ocnCompNDoFs   = disc_orders[1] * disc_orders[1] /* FV */;

#ifdef COMPUTE_ONLINE_MAP
    const iMOAB_String map_from_mem_identifier = "map-computed-online";
    const iMOAB_String disc_methods[2]         = { "fv", "fv" };
    const iMOAB_String dof_tag_names[2]        = { "GLOBAL_ID", "GLOBAL_ID" };
    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;
    if( couComm != MPI_COMM_NULL )
    {
        // compute the mesh intersection between ATM and OCN
        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        CHECKIERR( iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplAtmOcnMemPID ),
                   "cannot compute intersection for ATM/OCN" )
        POP_TIMER( couComm, rankInCouComm )

        if( !skip_apply_a2o )
        {
            // next compute the FV-FV map for runtime field projections
            PUSH_TIMER( "Compute the projection weights with TempestRemap" )
            CHECKIERR( iMOAB_ComputeScalarProjectionWeights(
                           cplAtmOcnMemPID, map_from_mem_identifier, disc_methods[0], &disc_orders[0], disc_methods[1],
                           &disc_orders[1], nullptr, &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap,
                           &fNoConserve, &fValidate, dof_tag_names[0], dof_tag_names[1] ),
                       "cannot compute scalar projection weights" )
            POP_TIMER( couComm, rankInCouComm )

            {
                const iMOAB_String atmocn_map_file_name = "atm_ocn_map_computed.nc";
                CHECKIERR( iMOAB_WriteMappingWeightsToFile( cplAtmOcnMemPID, map_from_mem_identifier,
                                                            atmocn_map_file_name ),
                           "failed to write map file to disk" );
            }
        }
        int meshtype = 3;
        PUSH_TIMER( "Compute ATM coverage graph for OCN mesh, for compute graph" )
        CHECKIERR( iMOAB_ComputeCommGraph( cplAtmPID, cplAtmOcnMemPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                           &meshtype, &cplatm, &atmocnmid ),
                   "cannot recompute ATM source coverage graph for ocean" )
        POP_TIMER( couComm, rankInCouComm )  // hijack this rank
    }
#endif

    int filter_type                           = 0;
    const iMOAB_String bottomFields           = "AnalyticalSolnSrcExact";
    const iMOAB_String bottomProjectedFieldsF = "Target_projF";
    std::string allProjectedFields( bottomProjectedFieldsF );
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
    const iMOAB_String bottomProjectedFieldsS = "Source_projF";
#endif
#ifdef COMPUTE_ONLINE_MAP
    const iMOAB_String bottomProjectedFieldsM = "Target_projM";
    allProjectedFields                        = std::string( bottomProjectedFieldsM );
#endif
#if defined( COMPUTE_FILE_MAP ) && defined( COMPUTE_ONLINE_MAP )
    allProjectedFields = std::string( bottomProjectedFieldsM ) + ":" + std::string( bottomProjectedFieldsF );
#endif

    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, bottomFields, &tagTypes, &atmCompNDoFs, &tagIndex[0] ),
                   "failed to define the field tags AnalyticalSolnSrcExact" );

        // just to be sure it is set, to be visible by iMOAB app
        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, "aream", &tagTypes, &atmCompNDoFs, &tagIndex[0] ),
                           "failed to define the field tags aream" );

        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, allProjectedFields.c_str(), &tagTypes, &ocnCompNDoFs,
                                           &tagIndex[1] ),
                   "failed to define the field tags allProjectedFields" );
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, "aream", &tagTypes, &ocnCompNDoFs,
                                           &tagIndex[1] ),
                   "failed to define the field tag aream" );

#ifdef COMPUTE_TRANSPOSE_FILE_MAP
        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, bottomProjectedFieldsS, &tagTypes, &atmCompNDoFs, &tagIndex[0] ),
                   "failed to define the field tags AnalyticalSolnSrcExact" );
#endif
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        CHECKIERR( iMOAB_DefineTagStorage( cmpOcnPID, allProjectedFields.c_str(), &tagTypes, &ocnCompNDoFs,
                                           &tagIndexIn2 ),
                   "failed to define the field tag for receiving back the tags "
                   "allProjectedFields on OCN pes" );
    }

#ifdef COMPUTE_TRANSPOSE_FILE_MAP
    if( ocnComm != MPI_COMM_NULL )
    {
        int tagIndexIn2;
        CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, bottomProjectedFieldsS, &tagTypes, &atmCompNDoFs, &tagIndexIn2 ),
                   "failed to define the field tag for receiving back the tags "
                   "bottomProjectedFieldsS on ATM pes" );
    }
#endif
    // start a virtual loop for number of iterations
    for( int iters = 0; iters < number_iterations; iters++ )
    {
        PUSH_TIMER( "Send/receive data from ATM component to coupler " )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            CHECKIERR( iMOAB_SendElementTag( cmpAtmPID, bottomFields, &atmCouComm, &cplatm ), "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on ATM on coupler pes, that was redistributed according to coverage
            CHECKIERR( iMOAB_ReceiveElementTag( cplAtmPID, bottomFields, &atmCouComm, &cmpatm ),
                       "cannot receive tag values" )
        }

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ),
                       "cannot free buffers used to resend ATM tag towards the coverage mesh" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        // second hop, is from atm towards ocean, on coupler
        //  it should send from each part on coupler towards the coverage set that should form the
        // rings around target cells (ocean)
        // basically we should send to more cells than needed just for intersection
#if defined( COMPUTE_FILE_MAP )
        if( couComm != MPI_COMM_NULL )
        {
            // send using the par comm graph computed by iMOAB_ComputeCommGraph
            CHECKIERR( iMOAB_SendElementTag( cplAtmPID, bottomFields, &couComm, &atmocnfid ),
                       "cannot send tag values towards coverage mesh for bilinear map" )

            CHECKIERR( iMOAB_ReceiveElementTag( cplAtmOcnFilePID, bottomFields, &couComm, &cplatm ),
                       "cannot receive tag values for bilinear map" )

            CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &atmocnfid ), "cannot free buffers" )
        }
#endif

#if defined( COMPUTE_ONLINE_MAP )
        if( couComm != MPI_COMM_NULL )
        {
            // send using the par comm graph computed by iMOAB_ComputeCommGraph
            CHECKIERR( iMOAB_SendElementTag( cplAtmPID, bottomFields, &couComm, &atmocnmid ),
                       "cannot send tag values towards coverage mesh for bilinear map" )

            CHECKIERR( iMOAB_ReceiveElementTag( cplAtmOcnMemPID, bottomFields, &couComm, &cplatm ),
                       "cannot receive tag values for bilinear map" )

            CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &atmocnmid ), "cannot free buffers" )
        }
#endif

#ifdef VERBOSE
        if( couComm != MPI_COMM_NULL && 1 == number_iterations )
        {
            // write only for n==1 case
            char outputFileRecvd[] = "cplAtmFile.h5m";
            CHECKIERR( iMOAB_WriteMesh( cplAtmPID, outputFileRecvd, fileWriteOptions ),
                       "could not write cplAtmFile.h5m to disk" )
        }
#endif

        if( couComm != MPI_COMM_NULL )
        {
#ifdef COMPUTE_FILE_MAP
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply from file scalar projection weights" )
            CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplAtmOcnFilePID, &filter_type, map_from_file_identifier[0],
                                                           bottomFields, bottomProjectedFieldsF ),
                       "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
#endif  // COMPUTE_FILE_MAP

// #undef COMPUTE_FILE_MAP
#ifdef COMPUTE_TRANSPOSE_FILE_MAP

            // send using the par comm graph computed by iMOAB_ComputeCommGraph
            CHECKIERR( iMOAB_SendElementTag( cplOcnPID, bottomProjectedFieldsF, &couComm, &ocnatmfid ),
                       "cannot send tag values towards OCN coverage mesh" )

            CHECKIERR( iMOAB_ReceiveElementTag( cplOcnAtmFilePID, bottomProjectedFieldsF, &couComm, &cplocn ),
                       "cannot receive tag values in OCN coverage from OCN coupler instance" )

            CHECKIERR( iMOAB_FreeSenderBuffers( cplOcnPID, &ocnatmfid ), "cannot free buffers" )

            PUSH_TIMER( "Apply from file scalar projection weights" )
            CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplOcnAtmFilePID, &filter_type, map_from_file_identifier[1],
                                                           bottomProjectedFieldsF, bottomProjectedFieldsS ),
                       "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
// #ifdef VERBOSE
//             char outputFileAtmFile[] = "AtmWithProjection.h5m";
//             CHECKIERR( iMOAB_WriteMesh( cplAtmPID, outputFileAtmFile, fileWriteOptions ),
//                        "could not write AtmWithProjection.h5m to disk" )
// #endif  // VERBOSE
#endif  // COMPUTE_TRANSPOSE_FILE_MAP

#ifdef COMPUTE_ONLINE_MAP
            PUSH_TIMER( "Apply in-memory scalar projection weights" )
            CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplAtmOcnMemPID, &filter_type, map_from_mem_identifier,
                                                           bottomFields, bottomProjectedFieldsM ),
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
            CHECKIERR( iMOAB_SendElementTag( cplOcnPID, allProjectedFields.c_str(), &ocnCouComm, &cmpocn ),
                       "cannot send tag values back to ocean pes" )
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
            // need to use ocean comp id for context
            CHECKIERR( iMOAB_SendElementTag( cplAtmPID, bottomProjectedFieldsS, &atmCouComm, &cmpatm ),
                       "cannot send tag values back to ocean pes" )
#endif
#ifdef VERBOSE
            {
                // write only for n==1 case
                char outputFileRecvd[] = "cplProjectedOCNFileMF.h5m";
                CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileRecvd, fileWriteOptions ),
                           "could not write cplProjectedOCNFile.h5m to disk" )
            }
#endif
        }

        // receive on component OCN
        if( ocnComm != MPI_COMM_NULL )
        {
            CHECKIERR( iMOAB_ReceiveElementTag( cmpOcnPID, allProjectedFields.c_str(), &ocnCouComm, &cplocn ),
                       "cannot receive tag values from OCN mesh on coupler pes" )
        }

#ifdef COMPUTE_TRANSPOSE_FILE_MAP
        // receive on component ATM
        if( atmComm != MPI_COMM_NULL )
        {
            CHECKIERR( iMOAB_ReceiveElementTag( cmpAtmPID, bottomProjectedFieldsS, &atmCouComm, &cplatm ),
                       "cannot receive tag values from ATM mesh on coupler pes" )
        }
#endif

        if( couComm != MPI_COMM_NULL )
        {
            CHECKIERR( iMOAB_FreeSenderBuffers( cplOcnPID, &cmpocn ), "Freeing buffers failed" )
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
            CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &cmpatm ), "Freeing buffers failed" )
#endif
        }

        if( ocnComm != MPI_COMM_NULL && 1 == number_iterations )  // write only for n==1 case
        {
#ifdef VERBOSE
            char outputFileOcnFile[] = "OcnWithProjection.h5m";
            CHECKIERR( iMOAB_WriteMesh( cmpOcnPID, outputFileOcnFile, fileWriteOptions ),
                       "could not write OcnWithProjection.h5m to disk" )
#endif
            // test results only for number_iterations== 1
            if( !no_regression_test )
            {
                // get global id storage
                const char* gidStr = "GLOBAL_ID";  // hard coded too
                int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
                // first get GlobalIds from OCN, and fields:
                int nverts[3], nelem[3];
                std::vector< int > gidElems;
                std::vector< double > tempElems;
                int err_code = 1, ent_type = 1;

                CHECKIERR( iMOAB_DefineTagStorage( cmpOcnPID, gidStr, &tag_type, &ncomp, &tagInd ),
                           "failed to define global id tag" );

#ifdef COMPUTE_ONLINE_MAP
                CHECKIERR( iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 ), "failed to get OCN mesh info" );
                gidElems.resize( nelem[2] );
                tempElems.resize( nelem[2] );

                CHECKIERR( iMOAB_GetIntTagStorage( cmpOcnPID, gidStr, &nelem[2], &ent_type, gidElems.data() ),
                           "failed to get global ids" );
                CHECKIERR( iMOAB_GetDoubleTagStorage( cmpOcnPID, bottomProjectedFieldsM, &nelem[2], &ent_type,
                                                      tempElems.data() ),
                           "failed to get bottomProjectedFieldsM field" );
                // check against the baseline
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn (in-memory map projection) on ocean task "
                              << rankInOcnComm << "\n";
#endif

#ifdef COMPUTE_TRANSPOSE_FILE_MAP
                CHECKIERR( iMOAB_GetMeshInfo( cmpAtmPID, nverts, nelem, 0, 0, 0 ), "failed to get OCN mesh info" );
                gidElems.resize( nelem[2] );
                tempElems.resize( nelem[2] );

                // we should not have to define global id tag, it is always defined
                CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, gidStr, &tag_type, &ncomp, &tagInd ),
                           "Failed to define GLOBAL_ID tag" );
                CHECKIERR( iMOAB_GetIntTagStorage( cmpAtmPID, gidStr, &nelem[2], &ent_type, gidElems.data() ),
                           "failed to get global ids" );
                CHECKIERR( iMOAB_GetDoubleTagStorage( cmpAtmPID, bottomProjectedFieldsS, &nelem[2], &ent_type,
                                                      tempElems.data() ),
                           "failed to get bottomProjectedFieldsS field" );

#ifdef VERBOSE
                char outputFileAtmFile[] = "AtmWithProjection.h5m";
                CHECKIERR( iMOAB_WriteMesh( cmpAtmPID, outputFileAtmFile, fileWriteOptions ),
                           "could not write AtmWithProjection.h5m to disk" )
#endif  // VERBOSE

                constexpr bool gen_baseline = false;
                if( gen_baseline )
                {
                    // int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
                    // // we should not have to define global id tag, it is always defined
                    // // CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, gidStr, &tag_type, &ncomp, &tagInd ), "Failed to define GLOBAL_ID tag" );

                    std::fstream fs;
                    fs.open( baselineATM, std::fstream::out );
                    fs << std::setprecision( 15 );  // maximum precision for doubles
                    for( int i = 0; i < nelem[2]; i++ )
                        fs << gidElems[i] << " " << tempElems[i] << "\n";
                    fs.close();
                }
                else
                {
                    // check against the baseline
                    check_baseline_file( baselineATM, gidElems, tempElems, 1.e-9, err_code );
                    if( 0 == err_code )
                        std::cout << " passed baseline test ocn2atm (file-based map projection) on ocean task "
                                  << rankInOcnComm << "\n";
                }
#endif

#ifdef COMPUTE_FILE_MAP
                CHECKIERR( iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 ), "failed to get OCN mesh info" );
                gidElems.resize( nelem[2] );
                tempElems.resize( nelem[2] );

                CHECKIERR( iMOAB_GetIntTagStorage( cmpOcnPID, gidStr, &nelem[2], &ent_type, gidElems.data() ),
                           "failed to get global ids" );
                CHECKIERR( iMOAB_GetDoubleTagStorage( cmpOcnPID, bottomProjectedFieldsF, &nelem[2], &ent_type,
                                                      tempElems.data() ),
                           "failed to get bottomProjectedFieldsF field" );

                // check against the baseline
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn (file-based map projection) on ocean task "
                              << rankInOcnComm << "\n";
#endif
            }
        }

    }  // end loop iterations

#ifdef COMPUTE_ONLINE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmOcnMemPID ), "cannot deregister app intx AO" )
    }
#endif

#ifdef COMPUTE_FILE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmOcnFilePID ), "cannot deregister app intx ATM-OCN" )
    }
#endif
#ifdef COMPUTE_TRANSPOSE_FILE_MAP
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplOcnAtmFilePID ), "cannot deregister app intx OCN-ATM" )
    }
#endif

    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmPID ), "cannot deregister app ATMX" )
        CHECKIERR( iMOAB_DeregisterApplication( cplOcnPID ), "cannot deregister app OCNX" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpOcnPID ), "cannot deregister app OCN" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpAtmPID ), "cannot deregister app ATM" )
    }

    // finalize iMOAB - all resources are now free
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

    // free all MPI resources
    MPI_Finalize();

    return 0;
}
