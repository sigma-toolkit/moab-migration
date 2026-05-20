/*
 * This imoab_coupler_bilin test will simulate coupling between atm and ocean using bilinear map
 * 2 meshes will be loaded from 2 files (atm, ocean), and they will be migrated to
 * coupler processors (coupler pes); then, intx will be performed between migrated meshes
 * and weights will be generated, such that a field from one component will be transferred to
 * the other component
 * currently, the atm will send some data to be projected to ocean comp
 *
 * first, intersect atm and ocn, compute bilinear map
 * use the current 2 hop strategy; first send from atm comp to atm coupler, then,
 * in a second hop, send from atm coupler to coverage for ocn intx
 *
 * project then 2 fields to the ocean (Sa_pbot and Sa_dens)
 * we also write the map in parallel
 *
 * when we run on 1 process and 2 processes, we get a different map and
 * different field projected
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
//#include <iomanip>
#include "imoab_coupler_utils.hpp"

using namespace moab;

//#define GRAPH_INFO

#ifndef MOAB_HAVE_TEMPESTREMAP
#error The climate coupler test example requires MOAB configuration with TempestRemap
#endif

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

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

    std::string atmFilename = TestDir + "unittest/atm_c2x.h5m";
    // on a regular case,  5 ATM, 6 CPLATM (ATMX), 17 OCN     , 18 CPLOCN (OCNX)  ;
    // intx atm/ocn is not in e3sm yet, give a number
    //   6 * 100+ 18 = 618 : atmocnid
    // 9 LND, 10 CPLLND
    //   6 * 100 + 10 = 610  atmlndid:
    // cmpatm is for atm on atm pes
    // cmpocn is for ocean, on ocean pe
    // cplatm is for atm on coupler pes
    // cplocn is for ocean on coupelr pes
    // atmocnid is for intx atm / ocn on coupler pes
    //
    int rankInAtmComm       = -1;
    int cmpatm              = 5,
        cplatm              = 6;  // component ids are unique over all pes, and established in advance;
    std::string ocnFilename = TestDir + "unittest/wholeOcn.h5m";
    int rankInOcnComm       = -1;
    int cmpocn = 17, cplocn = 18,
        atmocnid = 618;  // component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0, startG2 = 0, endG1 = numProcesses - 1, endG2 = numProcesses - 1;
    // Support launch of imoab_coupler test on any combo of 2*x processes
    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id;                        // used now for freeing buffers

    // default: load atm on 2 proc, ocean on 2, land on 2; migrate to 2 procs, then compute intx
    // later, we need to compute weight matrix with tempestremap

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "ocean mesh filename (target)", &ocnFilename );
    std::string baseline = TestDir + "unittest/baseline3.txt";
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.addOpt< int >( "partitioning,p", "partitioning option for migration", &repartitioner_scheme );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 3", &no_regression_test );

    std::string digestPrefix;
    opts.addOpt< std::string >( "digest_prefix", "prefix for BfB digest output files", &digestPrefix );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n ocn file: " << ocnFilename << "\n     on tasks : " << startG2 << ":" << endG2
                  << "\n  partitioning (0 trivial, 1 graph, 2 geometry) " << repartitioner_scheme << "\n  ";
    }

    // load files on 3 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm );
    CHECKIERR( ierr, "Cannot create atm MPI group and communicator " )

    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm );
    CHECKIERR( ierr, "Cannot create ocn MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // atm_coupler
    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    ierr = create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm );
    CHECKIERR( ierr, "Cannot create joint atm cou communicator" )

    // ocn_coupler
    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    ierr = create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm );
    CHECKIERR( ierr, "Cannot create joint ocn cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpAtmAppID       = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;  // atm
    int cplAtmAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // atm on coupler PEs
    int cmpOcnAppID       = -1;
    iMOAB_AppID cmpOcnPID = &cmpOcnAppID;        // ocn
    int cplOcnAppID = -1, cplAtmOcnAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplAtmOcnPID = &cplAtmOcnAppID;  // intx atm -ocn on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMCPL", &couComm, &cplatm,
                                          cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )
        ierr = iMOAB_RegisterApplication( "OCNCPL", &couComm, &cplocn,
                                          cplOcnPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        ierr = iMOAB_RegisterApplication( "ATMCMP", &atmComm, &cmpatm, cmpAtmPID );
        CHECKIERR( ierr, "Cannot register ATM App" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        ierr = iMOAB_RegisterApplication( "OCNCMP", &ocnComm, &cmpocn, cmpOcnPID );
        CHECKIERR( ierr, "Cannot register OCN App" )
    }

    // atm
    ierr =
        setup_component_coupler_meshes( cmpAtmPID, cmpatm, cplAtmPID, cplatm, &atmComm, &atmPEGroup, &couComm,
                                        &couPEGroup, &atmCouComm, atmFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate atm mesh" )

    if( atmCouComm != MPI_COMM_NULL )
    {
        // We just need to create a comm graph to internally transfer data from coupler atm to coupler ocean
        // ierr = iMOAB_CoverageGraph( &couComm, cplAtm2PID, cplAtm2OcnPID, cplAtm2OcnPID, &cplatm2, &atm2ocnid,
        // &cplocn );  // it happens over joint communicator
        int type = 3;
        CHECKIERR( iMOAB_ComputeCommGraph( cmpAtmPID, cplAtmPID, &atmCouComm, &atmPEGroup, &couPEGroup, &type, &type,
                                           &cmpatm, &cplatm ),
                   "cannot compute graph between ATM comp and ATM migrated to coupler" )
    }

#ifdef VERBOSE
    if( couComm != MPI_COMM_NULL )
    {  // write only for n==1 case
        char outputFileTgt3[] = "recvAtmx.h5m";
        ierr                  = iMOAB_WriteMesh( cplAtmPID, outputFileTgt3, fileWriteOptions );
        CHECKIERR( ierr, "cannot write atm mesh after receiving" )
    }
#endif

    // ocean
    ierr =
        setup_component_coupler_meshes( cmpOcnPID, cmpocn, cplOcnPID, cplocn, &ocnComm, &ocnPEGroup, &couComm,
                                        &couPEGroup, &ocnCouComm, ocnFilename, readopts, nghlay, repartitioner_scheme );
    CHECKIERR( ierr, "Cannot load and migrate ocn mesh" )

    if( couComm != MPI_COMM_NULL )
    {
        // now compute intersection between OCNx and ATMx on coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMOCNCPL", &couComm, &atmocnid, cplAtmOcnPID );
        CHECKIERR( ierr, "Cannot register ocn_atm intx over coupler pes " )
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

    int disc_orders[1]                       = { 1 };
    const std::string weights_identifiers[1] = { "bilinear" };
    const std::string disc_methods[1]        = { "fv" };
    const std::string dof_tag_names[1]       = { "GLOBAL_ID" };
    const std::string method                 = "bilin";

    if( couComm != MPI_COMM_NULL )
    {
        // set the ghost layers on the coupler for the source mesh
        nghlay         = 1;  // number of ghost layers for source and target
        int nghlay_tgt = 0;
        ierr           = iMOAB_SetMapGhostLayers( cplAtmOcnPID, &nghlay, &nghlay_tgt );
        CHECKIERR( ierr, "cannot set ghost layers" )

        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplAtmOcnPID );
        // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        // basically, atm was redistributed according to target (ocean) partition, to "cover" the
        // ocean partitions check if intx valid, write some h5m intx file
        CHECKIERR( ierr, "cannot compute intersection" )
        POP_TIMER( couComm, rankInCouComm )
    }

    if( couComm != MPI_COMM_NULL )
    {
        // We just need to create a comm graph to internally transfer data from coupler atm to coupler ocean
        // ierr = iMOAB_CoverageGraph( &couComm, cplAtm2PID, cplAtm2OcnPID, cplAtm2OcnPID, &cplatm2, &atm2ocnid,
        // &cplocn );  // it happens over joint communicator
        int type = 3;
        ierr     = iMOAB_ComputeCommGraph( cplAtmPID, cplAtmOcnPID, &couComm, &couPEGroup, &couPEGroup, &type, &type,
                                           &cplatm, &atmocnid );
        CHECKIERR( ierr, "cannot recompute direct coverage graph for ocean from atm" )
    }

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;

    if( couComm != MPI_COMM_NULL )
    {
        PUSH_TIMER( "Compute the projection weights with TempestRemap" )
        ierr = iMOAB_ComputeScalarProjectionWeights( cplAtmOcnPID, weights_identifiers[0].c_str(),
                                                     disc_methods[0].c_str(), &disc_orders[0], disc_methods[0].c_str(),
                                                     &disc_orders[0], method.c_str(), &fNoBubble, &fMonotoneTypeID,
                                                     &fVolumetric, &fInverseDistanceMap, &fNoConserve, &fValidate,
                                                     dof_tag_names[0].c_str(), dof_tag_names[0].c_str() );
        CHECKIERR( ierr, "cannot compute scalar projection weights" )
        POP_TIMER( couComm, rankInCouComm )

        // Let us now write the map file to disk and then read it back to test the I/O API in iMOAB
#ifdef MOAB_HAVE_PNETCDF
        {
            std::stringstream outf;
            outf << "atm_ocn_bilin_map_p" << endG4 - startG4 + 1 << ".nc";  // number of tasks on coupler
            std::string mapfile = outf.str();  // write in parallel the map file, for debugging
            ierr = iMOAB_WriteMapFile( cplAtmOcnPID, weights_identifiers[0].c_str(), outf.str().c_str() );
            CHECKIERR( ierr, "failed to write map file to disk" );
        }
#endif
    }

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = disc_orders[0] * disc_orders[0], ocnCompNDoFs = 1 /*FV*/;
    int filter_type = 0;

    const char* bottomFields          = "Sa_dens:Sa_pbot";
    const char* bottomProjectedFields = "Sa_dens:Sa_pbot";

    if( ocnComm != MPI_COMM_NULL )
    {
        context_id = cplocn;  // id for ocean on coupler
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomProjectedFields, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag Sa_dens:Sa_pbot" );
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomFields, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag Sa_dens:Sa_pbot" );
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomProjectedFields, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag Sa_dens:Sa_pbot" );
    }

    // need to make sure that the coverage mesh (created during intx method) received the tag that
    // need to be projected to target so far, the coverage mesh has only the ids and global dofs;
    // need to change the migrate method to accommodate any GLL tag
    // now send a tag from original atmosphere (cmpAtmPID) towards migrated coverage mesh
    // (cplAtmPID), using the new coverage graph communicator

    // make the tag 0, to check we are actually sending needed data
    {
        if( cplAtmAppID >= 0 )
        {
            int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
            /*
             * Each process in the communicator will have access to a local mesh instance, which
             * will contain the original cells in the local partition and ghost entities. Number of
             * vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
             * conditions will be returned in numProcesses 3 arrays, for local, ghost and total
             * numbers.
             */
            ierr = iMOAB_GetMeshInfo( cplAtmPID, nverts, nelem, nblocks, nsbc, ndbc );
            CHECKIERR( ierr, "failed to get num primary elems" );
            int numAllElem = nelem[2];
            std::vector< double > vals;
            int storLeng = atmCompNDoFs * numAllElem * 3;  // 2 tags
            int eetype   = 1;

            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;

            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomFields, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag nul" )
            // set the tag to 0
        }
    }

    // first hop
    PUSH_TIMER( "Send/receive data from atm component to coupler in atm context" )
    if( atmComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cmpAtmPID, bottomFields, &atmCouComm, &cplatm );
        CHECKIERR( ierr, "cannot send tag values" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on atm on coupler pes
        ierr = iMOAB_ReceiveElementTag( cplAtmPID, bottomFields, &atmCouComm, &cmpatm );
        CHECKIERR( ierr, "cannot receive tag values" )
    }

    // we can now free the sender buffers
    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm );  // context is for ocean
        CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
    }
    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )

    // second hop, is from atm towards ocean, on coupler
    //  it should send from each part on coupler towards the coverage set that should form the
    // rings around target cells (ocean)
    // basically we should send to more cells than needed just for intersection
    //TODO
    if( couComm != MPI_COMM_NULL )
    {
        // send using the par comm graph computed by iMOAB_ComputeCommGraph
        ierr = iMOAB_SendElementTag( cplAtmPID, bottomFields, &couComm, &atmocnid );
        CHECKIERR( ierr, "cannot send tag values towards coverage mesh for bilinear map" )

        ierr = iMOAB_ReceiveElementTag( cplAtmOcnPID, bottomFields, &couComm, &cplatm );
        CHECKIERR( ierr, "cannot receive tag values for bilinear map" )

        ierr = iMOAB_FreeSenderBuffers( cplAtmPID, &atmocnid );
        CHECKIERR( ierr, "cannot free buffers" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
		   on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, &filter_type, weights_identifiers[0].c_str(),
                                                   bottomFields, bottomProjectedFields , nullptr);
        CHECKIERR( ierr, "failed to compute projection weight application" );
        POP_TIMER( couComm, rankInCouComm )
        {
            char outputFileTgt[] = "fOcnBilinOnCpl.h5m";
            ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileTgt, fileWriteOptions );
            CHECKIERR( ierr, "could not write fOcnBilinOnCpl.h5m to disk" )
        }
        if( !digestPrefix.empty() )
        {
            int couSize;
            MPI_Comm_size( couComm, &couSize );
            std::ostringstream oss;
            oss << digestPrefix << "_ocn_" << couSize << ".txt";
            ierr = gather_and_write_proj_tag( couComm, rankInCouComm, cplOcnPID, "Sa_pbot", oss.str() );
            if( ierr )
                std::cerr << "WARNING: could not write digest " << oss.str() << "\n";
        }
    }

    if( couComm != MPI_COMM_NULL )
    {
        // need to use ocean comp id for context
        context_id = cmpocn;  // id for ocean on comp
        ierr       = iMOAB_SendElementTag( cplOcnPID, bottomProjectedFields, &ocnCouComm, &context_id );
        CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
    }

    // receive on component 2, ocean
    if( ocnComm != MPI_COMM_NULL )
    {
        context_id = cplocn;  // id for ocean on coupler
        ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, bottomProjectedFields, &ocnCouComm, &context_id );
        CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        context_id = cmpocn;
        ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
        CHECKIERR( ierr, "cannot free send/receive buffers for OCN context" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        char outputFileOcn[] = "OcnWithProjBilin.h5m";
        ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions );
        CHECKIERR( ierr, "could not write OcnWithProjBilin.h5m to disk" )
    }
    // do a check agains a baseline test
    if( !no_regression_test && ( ocnComm != MPI_COMM_NULL ) )
    {
        // the same as remap test
        // get temp field on ocean, from conservative, the global ids, and check to the baseline file
        // first get GlobalIds from ocn, and fields:
        int nverts[3], nelem[3];
        ierr = iMOAB_GetMeshInfo( cmpOcnPID, nverts, nelem, 0, 0, 0 );
        CHECKIERR( ierr, "failed to get ocn mesh info" );
        std::vector< int > gidElems;
        gidElems.resize( nelem[2] );
        std::vector< double > tempElems;
        tempElems.resize( nelem[2] );
        // get global id storage
        const std::string GidStr = "GLOBAL_ID";  // hard coded too
        int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
        ierr = iMOAB_DefineTagStorage( cmpOcnPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd );
        CHECKIERR( ierr, "failed to define global id tag" );

        int ent_type = 1;
        ierr         = iMOAB_GetIntTagStorage( cmpOcnPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] );
        CHECKIERR( ierr, "failed to get global ids" );
        ierr = iMOAB_GetDoubleTagStorage( cmpOcnPID, "Sa_pbot", &nelem[2], &ent_type, &tempElems[0] );
        CHECKIERR( ierr, "failed to get temperature field" );
        //        {
        //            // write baseline file
        //            std::fstream fs;
        //            fs.open( "baseline3.txt", std::fstream::out );
        //            fs << std::setprecision( 15 );  // maximum precision for doubles
        //            for( size_t i = 0; i < tempElems.size(); i++ )
        //                fs << gidElems[i] << " " << tempElems[i] << "\n";
        //            fs.close();
        //        }
        int err_code = 1;
        check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
        if( 0 == err_code ) std::cout << " passed baseline test atm2ocn on ocean task " << rankInOcnComm << "\n";
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmOcnPID );
        CHECKIERR( ierr, "cannot deregister app intx AO" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCN1" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATM1" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplOcnPID );
        CHECKIERR( ierr, "cannot deregister app OCNX" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmPID );
        CHECKIERR( ierr, "cannot deregister app ATMX" )
    }

    //#endif
    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free atm coupler group and comm
    if( MPI_COMM_NULL != atmCouComm ) MPI_Comm_free( &atmCouComm );
    MPI_Group_free( &joinAtmCouGroup );
    if( MPI_COMM_NULL != atmComm ) MPI_Comm_free( &atmComm );

    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    // free ocn - coupler group and comm
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    // free up MPI resources
    MPI_Group_free( &atmPEGroup );
    MPI_Group_free( &ocnPEGroup );
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
