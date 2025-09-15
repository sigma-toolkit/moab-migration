/*
 * This imoab_coupler test will simulate coupling between 3 components
 * 3 meshes will be loaded from 3 files (fv, ocean, lnd), and they will be migrated to
 * all processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 * currently, the fv will send some data to be projected to ocean and land components
 *
 * first, intersect fv and se, and recompute comm graph 1 between fv and fv_cx, for se intx
 * second, intersect fv and lnd, and recompute comm graph 2 between fv and fv_cx for lnd intx

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

using namespace moab;

//#define GRAPH_INFO
// #define VERBOSE
#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
    std::string readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" );

    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;  // use the graph partitioner in that caseS
#endif

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    // std::string fvFilename = TestDir + "unittest/wholeRof_06.h5m";
    std::string fvFilename = TestDir + "../usgs-rawdata.nc";
    // on a regular case,  5 FV, 6 CPLFV (FVX), 17 SE     , 18 CPLSE (SEX)  ;
    // intx fv/se is not in e3sm yet, give a number
    //   6 * 100+ 18 = 618 : fvseid
    // 9 LND, 10 CPLLND
    //   6 * 100 + 10 = 610  fvlndid:
    // cmpfv is for fv on fv pes
    // cmpse is for ocean, on ocean pe
    // cplfv is for fv on coupler pes
    // cplse is for ocean on coupelr pes
    // fvseid is for intx fv / se on coupler pes
    //
    int rankInFVComm = -1;
    int cmpfv        = 5,
        cplfv        = 6;  // component ids are unique over all pes, and established in advance;
    std::string seFilename = TestDir + "unittest/wholeATM_T.h5m";
    int rankInSEComm       = -1;
    int cmpse = 17, cplse = 18,
        fvseid = 618;  // component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startFV = 0, startSE = 0, endFV = numProcesses - 1, endSE = numProcesses - 1;
    int startCPL = startFV, endCPL = endFV;  // these are for coupler layout
    int context_id = -1;                   // used now for freeing buffers

    // default: load fv on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "fvosphere,t", "fv mesh filename (source)", &fvFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &seFilename );
    opts.addOpt< int >( "startFV,a", "start task for fvosphere layout", &startFV );
    opts.addOpt< int >( "endFV,b", "end task for fvosphere layout", &endFV );
    opts.addOpt< int >( "startSE,c", "start task for ocean layout", &startSE );
    opts.addOpt< int >( "endSE,d", "end task for ocean layout", &endSE );
    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startCPL );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endCPL );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    int n = 1;  // number of send/receive / project / send back cycles
    opts.addOpt< int >( "iterations,n", "number of iterations for coupler", &n );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " fv file: " << fvFilename << "\n   on tasks : " << startFV << ":" << endFV <<
            "\n se file: " << seFilename << "\n     on tasks : " << startSE << ":" << endSE <<
            "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
    }

    // load files on 3 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group fvPEGroup;
    MPI_Comm fvComm;
    ierr = create_group_and_comm( startFV, endFV, jgroup, &fvPEGroup, &fvComm );
    CHECKIERR( ierr, "Cannot create fv MPI group and communicator " )

    MPI_Group sePEGroup;
    MPI_Comm seComm;
    ierr = create_group_and_comm( startSE, endSE, jgroup, &sePEGroup, &seComm );
    CHECKIERR( ierr, "Cannot create se MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startCPL, endCPL, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // fv_coupler
    MPI_Group joinFVCouGroup;
    MPI_Comm fvCouComm;
    ierr = create_joint_comm_group( fvPEGroup, couPEGroup, &joinFVCouGroup, &fvCouComm );
    CHECKIERR( ierr, "Cannot create joint fv cou communicator" )

    // se_coupler
    MPI_Group joinSECouGroup;
    MPI_Comm seCouComm;
    ierr = create_joint_comm_group( sePEGroup, couPEGroup, &joinSECouGroup, &seCouComm );
    CHECKIERR( ierr, "Cannot create joint se cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpFVAppID       = -1;
    iMOAB_AppID cmpFVPID = &cmpFVAppID;  // fv
    int cplFVAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplFVPID = &cplFVAppID;  // fv on coupler PEs
    int cmpSEAppID       = -1;
    iMOAB_AppID cmpSEPID = &cmpSEAppID;        // se
    int cplSEAppID = -1, cplFVSEAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplSEPID    = &cplSEAppID;     // se on coupler PEs
    iMOAB_AppID cplFVSEPID = &cplFVSEAppID;  // intx fv -se on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "FVX", &couComm, &cplfv,
                                          cplFVPID );  // fv on coupler pes
        CHECKIERR( ierr, "Cannot register FV over coupler PEs" )
        ierr = iMOAB_RegisterApplication( "SEX", &couComm, &cplse,
                                          cplSEPID );  // se on coupler pes
        CHECKIERR( ierr, "Cannot register SE over coupler PEs" )
    }

    if( fvComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( fvComm, &rankInFVComm );
        ierr = iMOAB_RegisterApplication( "FV1", &fvComm, &cmpfv, cmpFVPID );
        CHECKIERR( ierr, "Cannot register FV App" )
    }

    if( seComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( seComm, &rankInSEComm );
        ierr = iMOAB_RegisterApplication( "SE1", &seComm, &cmpse, cmpSEPID );
        CHECKIERR( ierr, "Cannot register SE App" )
    }

    // fv
    ierr =
        setup_component_coupler_meshes( cmpFVPID, cmpfv, cplFVPID, cplfv, &fvComm, &fvPEGroup, &couComm,
                                        &couPEGroup, &fvCouComm, fvFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate fv mesh" )

    int global_fv[2] = { 0, 0 };
    if( fvComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_GetGlobalInfo( cmpFVPID, &( global_fv[0] ), &( global_fv[1] ) );
        CHECKIERR( ierr, "cannot get fv global info" )
        if( !rankInFVComm )
            std::cout << " fv mesh: vertices: " << global_fv[0] << " cells: " << global_fv[1] << "\n";
    }
    // broadcast to all tasks
    int global_fv_red[2] = { 0, 0 };
    MPI_Allreduce( global_fv, global_fv_red, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    if( couComm != MPI_COMM_NULL )
    {
        int global_fv_cou[2];
        ierr = iMOAB_GetGlobalInfo( cplFVPID, &( global_fv_cou[0] ), &( global_fv_cou[1] ) );
        CHECKIERR( ierr, "cannot get fv global info on coupler" )
        CHECK_ARRAYS_EQUAL( global_fv_red, 2, global_fv_cou, 2 );
        if( !rankInCouComm )
            std::cout << " fv mesh on coupler: vertices: " << global_fv[0] << " cells: " << global_fv[1] << "\n";
    }

#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL && 1 == n )
    {  // write only for n==1 case
        char outputFileTgt3[] = "recvFV.h5m";
        ierr = iMOAB_WriteMesh( cplFVPID, outputFileTgt3, fileWriteOptions );
        CHECKIERR( ierr, "cannot write fv mesh after receiving" )
    }
#endif
#ifdef GRAPH_INFO
    if( fvComm != MPI_COMM_NULL )
    {
        int is_sender = 1;
        int context   = -1;
        iMOAB_DumpCommGraph( cmpFVPID, &context, &is_sender, "FVMigS" );
    }
    if( couComm != MPI_COMM_NULL )
    {
        int is_sender = 0;
        int context   = -1;
        iMOAB_DumpCommGraph( cplFVPID, &context, &is_sender, "FVMigR" );
    }
#endif
    MPI_Barrier( MPI_COMM_WORLD );

    // ocean
    ierr =
        setup_component_coupler_meshes( cmpSEPID, cmpse, cplSEPID, cplse, &seComm, &sePEGroup, &couComm,
                                        &couPEGroup, &seCouComm, seFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate se mesh" )

    MPI_Barrier( MPI_COMM_WORLD );

#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL && 1 == n )
    {  // write only for n==1 case
        char outputFileTgt3[] = "recvSE.h5m";
        ierr                  = iMOAB_WriteMesh( cplSEPID, outputFileTgt3, fileWriteOptions );
        CHECKIERR( ierr, "cannot write se mesh after receiving" )
    }
#endif

    // this model (recMeshSE.h5m) has mixed meshes in it, we need to repair the comm graph
    // first delete the one created with migration, then compute a new one
    if( seCouComm != MPI_COMM_NULL )
    {
        int srctype = 2;  // type: 1 - SE, 2 - Vertex (point cloud), 3 - Element (FV scalars)
        int dsttype = 1;  // type: 1 - SE, 2 - Vertex (point cloud), 3 - Element (FV scalars)
        CHECKIERR( iMOAB_ComputeCommGraph( cmpSEPID, cplSEPID, &seCouComm, &sePEGroup, &couPEGroup, &srctype, &dsttype,
                                           &cmpse, &cplse ),
                   "cannot compute graph between se comp and se migrated to coupler" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between SEx and FVx on coupler PEs
        ierr = iMOAB_RegisterApplication( "FVSE", &couComm, &fvseid, cplFVSEPID );
        CHECKIERR( ierr, "Cannot register se_fv intx over coupler pes " )
    }

    int disc_orders[2]                        = { 1, 4 };
    const iMOAB_String weights_identifiers[1] = { "fvse-averaged" };
    const iMOAB_String disc_methods[2]        = { "pcloud", "cgll" };
    const iMOAB_String dof_tag_names[2]       = { "GLOBAL_ID", "GLOBAL_DOFS" };

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute coupler mesh coverage" )
        CHECKIERR( iMOAB_ComputeCoverageMesh( cplFVPID, cplSEPID, cplFVSEPID ), "cannot compute intersection" )
        // ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplFVPID, cplSEPID, cplFVSEPID );
        // coverage mesh was computed here, for cplFVPID, fv on coupler pes
        // basically, fv was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        POP_TIMER( couComm, rankInCouComm )
    }

    if( fvCouComm != MPI_COMM_NULL )
    {
        // the new graph will be for sending data from fv comp to coverage mesh;
        // it involves initial fv app; cmpFVPID; also migrate fv mesh on coupler pes, cplFVPID
        // results are in cplFVSEPID, intx mesh; remapper also has some info about coverage mesh
        // after this, the sending of tags from fv pes to coupler pes will use the new par comm
        // graph, that has more precise info about what to send for ocean cover ; every time, we
        // will
        //  use the element global id, which should uniquely identify the element
        PUSH_TIMER( "Compute SE coverage graph for FV mesh" )
        ierr = iMOAB_CoverageGraph( &fvCouComm, cmpFVPID, cplFVPID, cplFVSEPID, &cmpfv, &cplfv,
                                    &cplse );  // it happens over joint communicator
        CHECKIERR( ierr, "cannot recompute direct coverage graph for ocean" )
        POP_TIMER( fvCouComm, rankInFVComm )  // hijack this rank
    }

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;

#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL && 1 == n )
    {                                    // write only for n==1 case
        char serialWriteOptions[] = "";  // for writing in serial
        std::stringstream outf;
        outf << "intxFVSE_" << rankInCouComm << ".h5m";
        std::string intxfile = outf.str();  // write in serial the intx file, for debugging
        ierr                 = iMOAB_WriteMesh( cplFVSEPID, intxfile.c_str(), serialWriteOptions );
        CHECKIERR( ierr, "cannot write intx file result" )
    }
#endif

    if( couComm != MPI_COMM_NULL )
    {
        const char* fvmethod = "fvse-averaged";
        PUSH_TIMER( "Compute the projection weights with TempestRemap" )
        ierr = iMOAB_ComputeScalarProjectionWeights( cplFVSEPID, weights_identifiers[0], disc_methods[0],
                                                     &disc_orders[0], disc_methods[1], &disc_orders[1], fvmethod,
                                                     &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap,
                                                     &fNoConserve, &fValidate, dof_tag_names[0], dof_tag_names[1] );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )

        // Let us now write the map file to disk and then read it back to test the I/O API in iMOAB
#ifdef MOAB_HAVE_NETCDF
        if (false)
        {
            const iMOAB_String fvse_map_file_name = "fv_se_map.nc";
            ierr = iMOAB_WriteMapFile( cplFVSEPID, weights_identifiers[0], fvse_map_file_name );
            CHECKIERR( ierr, "failed to write map file to disk" );

            // const iMOAB_String intx_from_file_identifier = "fvse-map-from-file";
            // int src_disc_type                            = 1;  // element-based SE-4
            // int tgt_disc_type                            = 3;  // element-based FV
            // int arearead = 1; // read only area_a (fvosphere)
            // CHECKIERR( iMOAB_LoadMapFile( cplFVPID, cplSEPID, cplFVSEPID, &src_disc_type,
            //                                              &tgt_disc_type, &arearead, intx_from_file_identifier,
            //                                              fvse_map_file_name ),
            //            "failed to load map file from disk" );
        }
#endif
    }

    int tagIndex[3];
    int tagTypes[3]  = { DENSE_DOUBLE, DENSE_DOUBLE, DENSE_DOUBLE };
    int fvCompNDoFs = disc_orders[0] * disc_orders[0] /*FV*/, seCompNDoFs = disc_orders[1] * disc_orders[1] /*SE*/;
    int filter_type = 0;

    // const char* bottomFields          = "Firr_rofi:Flrr_volr:area";
    // const char* bottomProjectedFields = "Firr_rofi_proj:Flrr_volr_proj:area_proj";
    const char* bottomFields          = "landfract:htopo";
    const char* bottomProjectedFields = "landfract_proj:htopo_proj";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplFVPID, bottomFields, &tagTypes[0], &fvCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag Firr_rofi" );
        ierr = iMOAB_DefineTagStorage( cplSEPID, bottomProjectedFields, &tagTypes[1], &seCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag Firr_rofi_proj" );
    }

    // need to make sure that the coverage mesh (created during intx method) received the tag that
    // need to be projected to target so far, the coverage mesh has only the ids and global dofs;
    // need to change the migrate method to accommodate any GLL tag
    // now send a tag from original fvosphere (cmpFVPID) towards migrated coverage mesh
    // (cplFVPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        // if( cplFVAppID >= 0 )
        // {
        //     int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
        //     /*
        //      * Each process in the communicator will have access to a local mesh instance, which
        //      * will contain the original cells in the local partition and ghost entities. Number of
        //      * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
        //      * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
        //      * numbers.
        //      */
        //     ierr = iMOAB_GetMeshInfo( cplFVPID, nverts, nelem, nblocks, nsbc, ndbc );
        //     CHECKIERR( ierr, "failed to get num primary elems" );
        //     int numAllElem = nverts[2];
        //     std::vector< double > vals;
        //     int storLeng = fvCompNDoFs * numAllElem * 3;  // 3 tags
        //     int eetype   = 0; // vertices

        //     vals.resize( storLeng );
        //     for( int k = 0; k < storLeng; k++ )
        //         vals[k] = 0.;

        //     ierr = iMOAB_SetDoubleTagStorage( cplFVPID, bottomFields, &storLeng, &eetype, &vals[0] );
        //     CHECKIERR( ierr, "cannot make tag nul" )
        //     // set the tag to 0
        // }
    }

    // start a virtual loop for number of iterations
    for( int iters = 0; iters < n; iters++ )
    {
        PUSH_TIMER( "Send/receive data from fv component to coupler in se context" )
        if( fvComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpFVPID, bottomFields, &fvCouComm, &cplse );
            CHECKIERR( ierr, "cannot send tag values" )
#ifdef GRAPH_INFO
            int is_sender = 1;
            int context   = cplse;
            iMOAB_DumpCommGraph( cmpFVPID, &context, &is_sender, "FVCovSES" );
#endif
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on fv on coupler pes, that was redistributed according to coverage
            ierr = iMOAB_ReceiveElementTag( cplFVPID, bottomFields, &fvCouComm, &cplse );
            CHECKIERR( ierr, "cannot receive tag values" )
#ifdef GRAPH_INFO
            int is_sender = 0;
            int context   = cplse;  // the same context, cplse
            iMOAB_DumpCommGraph( cmpFVPID, &context, &is_sender, "FVCovSER" );
#endif
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

        // we can now free the sender buffers
        if( fvComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpFVPID, &cplse );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend fv tag towards the coverage mesh" )
        }
#ifdef VERBOSE
        if( couComm != MPI_COMM_NULL && 1 == n )
        {
            // write only for n==1 case
            char outputFileRecvd[] = "recvFVCoupSE.h5m";
            ierr                   = iMOAB_WriteMesh( cplFVPID, outputFileRecvd, fileWriteOptions );
            CHECKIERR( ierr, "could not write recvFVCoupSE.h5m to disk" )
        }
#endif

        if( couComm != MPI_COMM_NULL )
        {
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplFVSEPID, &filter_type, weights_identifiers[0], bottomFields,
                                                       bottomProjectedFields );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )
            if( 1 == n )  // write only for n==1 case
            {
                char outputFileTgt[] = "fSEOnCpl.h5m";
                ierr                 = iMOAB_WriteMesh( cplSEPID, outputFileTgt, fileWriteOptions );
                CHECKIERR( ierr, "could not write fSEOnCpl.h5m to disk" )
            }
        }

        // send the projected tag back to ocean pes, with send/receive tag
        if( seComm != MPI_COMM_NULL )
        {
            int tagIndexIn2;
            ierr =
                iMOAB_DefineTagStorage( cmpSEPID, bottomProjectedFields, &tagTypes[1], &seCompNDoFs, &tagIndexIn2 );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tags "
                             "Firr_rofi_proj, Flrr_volr_proj, area_proj on se pes" );
        }
        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm se_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            // need to use ocean comp id for context
            context_id = cmpse;  // id for ocean on comp
            ierr       = iMOAB_SendElementTag( cplSEPID, bottomProjectedFields, &seCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }

        // receive on component 2, ocean
        if( seComm != MPI_COMM_NULL )
        {
            context_id = cplse;  // id for ocean on coupler
            ierr       = iMOAB_ReceiveElementTag( cmpSEPID, bottomProjectedFields, &seCouComm, &context_id );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpse;
            ierr       = iMOAB_FreeSenderBuffers( cplSEPID, &context_id );
            CHECKIERR( ierr, "cannot free send/receive buffers for SE context" )
        }
        if( seComm != MPI_COMM_NULL && 1 == n )  // write only for n==1 case
        {
            char outputFileSE[] = "SEWithProj.h5m";
            ierr                 = iMOAB_WriteMesh( cmpSEPID, outputFileSE, fileWriteOptions );
            CHECKIERR( ierr, "could not write SEWithProj.h5m to disk" )
        }


    }  // end loop iterations n

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplFVSEPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
    if( seComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpSEPID );
        CHECKIERR( ierr, "cannot deregister app SE1" )
    }

    if( fvComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpFVPID );
        CHECKIERR( ierr, "cannot deregister app FV1" )
    }


    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplSEPID );
        CHECKIERR( ierr, "cannot deregister app SEX" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplFVPID );
        CHECKIERR( ierr, "cannot deregister app FVX" )
    }

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free fv coupler group and comm
    if( MPI_COMM_NULL != fvCouComm ) MPI_Comm_free( &fvCouComm );
    MPI_Group_free( &joinFVCouGroup );
    if( MPI_COMM_NULL != fvComm ) MPI_Comm_free( &fvComm );

    if( MPI_COMM_NULL != seComm ) MPI_Comm_free( &seComm );
    // free se - coupler group and comm
    if( MPI_COMM_NULL != seCouComm ) MPI_Comm_free( &seCouComm );
    MPI_Group_free( &joinSECouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &fvPEGroup );
    MPI_Group_free( &sePEGroup );
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
