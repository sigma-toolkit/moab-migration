/*
 * This imoab_map_r2l test will simulate coupling between land and river
 * 2 meshes will be loaded from 2 files (land domain file, source,
 * and target scrip file, tgt), and one map file read from disk
 * the migrate map mesh will be used to generate coverage set over target
 *  and will help for projection application
 */

#include "moab/Core.hpp"

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

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;

    // Timer data
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    int cplrof        = 22,
        cpllnd        = 10;  // component ids are unique over all pes, and established in advance;

    std::string rof_data = TestDir + "unittest/rof_comp_p32.h5m";
    std::string rof_mesh = TestDir + "unittest/SCRIPgrid_2x2_nomask_c210211.nc";
    std::string lndFilename = TestDir + "unittest/domain.lnd.ne4pg2_oQU480.200527.nc";
#ifdef MOAB_HAVE_ZOLTAN
    std::string readopts_lnd( "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION" );
    std::string readopts_rof( "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN" );
#else
    std::string readopts_lnd( "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=" );
    std::string readopts_rof( "PARALLEL=READ_PART;");
#endif

    std::string mapFilename = TestDir + "unittest/map_r2_to_ne4pg2_mono.210211.nc";  // this is a netcdf file!

    char field[] = "Forr_rofl";  // this is a tag name, on the exported rof file
    // this will be projected and generate a baseline after sending it to coupler

    std::string baseline = TestDir + "unittest/baseline_lnd.txt";
    int cmprof = 21,
        roflndid = 2210;  // 100*22 + 10 component ids are unique over all pes, and established in advance;

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0,  endG1 = numProcesses - 1;

    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout

    // default: load rof / source on 2 proc, land / target on 2,
    // load map on 2 also, in parallel, distributed by rows

    ProgOptions opts;
    opts.addOpt< std::string >( "rof,s", "rof mesh scrip filename (source)", &rof_mesh );
    opts.addOpt< std::string >( "rofdata,d", "rof data filename (source)", &rof_data );
    opts.addOpt< std::string >( "lnd,t", "land domain mesh filename (target)", &lndFilename );
    opts.addOpt< std::string >( "map_file,w", "map file from source to target", &mapFilename );

    opts.addOpt< int >( "startAtm,a", "start task for rof layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for rof layout", &endG1 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    int disc_orders[2] = { 1, 1 };  // 1 is for FV

    std::string fieldstr;
    opts.addOpt< std::string >( "field,f", "field to project using the map ", &fieldstr );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 1", &no_regression_test );
    opts.addOpt< std::string >( "newbaseline,n", "baseline to use for test ", &baseline );
    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " rof_data file: " << rof_data << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n rof_mesh scrip file on coupler: " << rof_mesh << "\n   on tasks : " << startG4 << ":" << endG4
                  << "\n lnd domain file on coupler " << lndFilename << "\n     on tasks : " << startG4 << ":" << endG4
                  << "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n"
                  << " baseline: " << baseline << "\n";
        if( !no_regression_test )
        {
            std::cout << " check projection against baseline: " << baseline << "\n";
        }
    }

    // load files on 2 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    MPI_Group rofPEGroup;
    MPI_Comm rofComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &rofPEGroup, &rofComm );
    CHECKIERR( ierr, "Cannot create atm MPI group and communicator " )

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )

    // rof comp + coupler communicator and group
    MPI_Group joinRofCouGroup;
    MPI_Comm rofCouComm;
    ierr = create_joint_comm_group( rofPEGroup, couPEGroup, &joinRofCouGroup, &rofCouComm );
    CHECKIERR( ierr, "Cannot create joint rof cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpRofAppID       = -1;
    iMOAB_AppID cmpRofPID = &cmpRofAppID;  // rof
    int cplRofAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplRofPID = &cplRofAppID;  // rof on coupler PEs

    int cplLndAppID = -1;
    iMOAB_AppID cplLndPID    = &cplLndAppID;     // lnd on coupler PEs
    int cplRofLndAppID = -1;
    iMOAB_AppID cplRofLndPID = &cplRofLndAppID;  // map rof -lnd on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ROFX", &couComm, &cplrof,
                cplRofPID );  // rof on coupler pes
        CHECKIERR( ierr, "Cannot register ROF over coupler PEs" )

        ierr = iMOAB_RegisterApplication( "LNDX", &couComm, &cpllnd,
                                          cplLndPID );  // lnd on coupler pes
        CHECKIERR( ierr, "Cannot register LND over coupler PEs" )
    }

    int rankInRofComm = -1;
    if( rofComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( rofComm, &rankInRofComm );
        ierr = iMOAB_RegisterApplication( "ROF1", &rofComm, &cmprof, cmpRofPID );
        CHECKIERR( ierr, "Cannot register ROF cmp App" )
    }

    MPI_Barrier( MPI_COMM_WORLD );


    if( couComm != MPI_COMM_NULL )
    {
        //
        ierr = iMOAB_RegisterApplication( "ROFLNDMAP", &couComm, &roflndid, cplRofLndPID );
        CHECKIERR( ierr, "Cannot register rof2lnd map instance over coupler pes " )
    }

    // load rof data on component rof; fake a time step export
    // rof data is on a point cloud mesh
    if (rofComm != MPI_COMM_NULL)
    {
        std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");// load a point cloud with rof data
        ierr = iMOAB_LoadMesh( cmpRofPID, rof_data.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load component data file" )
    }

    // load rof mesh, lnd mesh on coupler
    if (couComm != MPI_COMM_NULL)
    {
        ierr = iMOAB_LoadMesh( cplRofPID, rof_mesh.c_str(), readopts_rof.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load scrip file for rof " )

        ierr = iMOAB_WriteMesh( cplRofPID, "RofCpl1.h5m", fileWriteOptions );
        CHECKIERR( ierr, "Cannot write rof file on cpl " )
        ierr = iMOAB_LoadMesh( cplLndPID, lndFilename.c_str(), readopts_lnd.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load domain file for lnd " )
    }

    if( rofCouComm != MPI_COMM_NULL )
    {
        int type1 = 2; // 2: Vertex (point cloud); rof is point cloud on component side
        int type2 = 3; // 3: fv on coupler
        CHECKIERR( iMOAB_ComputeCommGraph( cmpRofPID, cplRofPID, &rofCouComm, &rofPEGroup, &couPEGroup, &type1, &type2,
                                           &cmprof, &cplrof ),
                   "cannot compute graph between rof on comp and rof on coupler" )
    }
    const std::string intx_from_file_identifier = "map-from-file";

    if( couComm != MPI_COMM_NULL )
    {
        int src_disc_type = 3;  // element-based FV
        int tgt_disc_type = 3;  // element-based FV
        int arearead = 0; // no need for aream
        CHECKIERR( iMOAB_LoadFromMappingFile( cplRofPID, cplLndPID, cplRofLndPID, &src_disc_type, &tgt_disc_type,
                                              &arearead, intx_from_file_identifier.c_str(), mapFilename.c_str() ),
                   "failed to load map file from disk" );
        int type      = 3;  // FV
        // because it is like "coverage", context will be roflndid
        ierr = iMOAB_MigrateMapMesh( cplRofPID, cplRofLndPID, &couComm, &couPEGroup, &couPEGroup, &type,
                                     &cplrof, &roflndid);
        CHECKIERR( ierr, "failed to migrate mesh for rof on coupler" );
    }
    MPI_Barrier( MPI_COMM_WORLD );

    int tagIndex;
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int compOrder = disc_orders[0] * disc_orders[0] /*FV*/;
    int filter_type = 0;

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplRofPID, field, &tagTypes[0], &compOrder, &tagIndex );
        CHECKIERR( ierr, "failed to define the field tag" );

        ierr = iMOAB_DefineTagStorage( cplLndPID, field, &tagTypes[1], &compOrder, &tagIndex );
        CHECKIERR( ierr, "failed to define the field tag on projection" );
    }

    if( rofComm != MPI_COMM_NULL  )  // we are on source rof pes
    {
        ierr = iMOAB_DefineTagStorage( cmpRofPID, field, &tagTypes[0], &compOrder, &tagIndex);
        CHECKIERR( ierr, "failed to define the field tag" );
    }

    // first hop
    if( rofComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cmpRofPID, field, &rofCouComm, &cplrof );
        CHECKIERR( ierr, "cannot send tag values" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on atm on coupler pes
        ierr = iMOAB_ReceiveElementTag( cplRofPID, field, &rofCouComm, &cmprof );
        CHECKIERR( ierr, "cannot receive tag values" )
    }

    // we can now free the sender buffers
    if( rofComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cmpRofPID, &cplrof );
        CHECKIERR( ierr, "cannot free buffers used to send rof towards coupler" )
    }
    if( couComm != MPI_COMM_NULL ){
        ierr = iMOAB_WriteMesh( cplRofPID, "RofCpl2.h5m", fileWriteOptions );
        CHECKIERR( ierr, "cannot write rof on coupler" )
    }
    // start the second hop, from rof cpl to rof coverage for lnd
    // the data is now on cpl Rof, need to be sent to rof coverage over lnd
    PUSH_TIMER( "Send/receive data from rof cpl to coverage in lnd context" )
    if( couComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cplRofPID, field, &couComm, &roflndid );
        CHECKIERR( ierr, "cannot send tag values" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on atm on coupler pes, that was redistributed according to coverage
        // the trick is we use the map imoab app
        ierr = iMOAB_ReceiveElementTag( cplRofLndPID, field, &couComm, &cplrof );
        CHECKIERR( ierr, "cannot receive tag values" )
    }


    // we can now free the sender buffers
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cplRofPID, &roflndid );  // context is for ocean
        CHECKIERR( ierr, "cannot free buffers " )
    }
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_WriteCoverageMesh( cplRofLndPID, "rof_cover_lnd");
        CHECKIERR( ierr, "cannot write coverage mesh" )
    }

    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )


    if( couComm != MPI_COMM_NULL )
    {
        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
           on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        ierr = iMOAB_ApplyScalarProjectionWeights( cplRofLndPID, &filter_type, intx_from_file_identifier.c_str(),
                                                   field, field );
        CHECKIERR( ierr, "failed to compute projection weight application" );
        POP_TIMER( couComm, rankInCouComm )

        {
            int numTasksCpl=endG4-startG4+1;
            std::ostringstream outfile;
            outfile << "fLndOnCpl_" << numTasksCpl << ".h5m";
            ierr = iMOAB_WriteMesh( cplLndPID, outfile.str().c_str(), fileWriteOptions );
            CHECKIERR( ierr, "could not write fLndOnCpl5.h5m to disk" )
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        if( !no_regression_test )
        {
            // the same as remap test
            // get rofl field on land, the global ids, and dump to the baseline file
            // first get GlobalIds from lnd, and fields:
            int nverts[3], nelem[3];
            ierr = iMOAB_GetMeshInfo( cplLndPID, nverts, nelem, 0, 0, 0 );
            CHECKIERR( ierr, "failed to get lnd mesh info" );
            std::vector< int > gidElems;
            gidElems.resize( nelem[2] );
            std::vector< double > tempElems;
            tempElems.resize( nelem[2] );
            // get global id storage
            const std::string GidStr = "GLOBAL_ID";  // hard coded too
            int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
            ierr = iMOAB_DefineTagStorage( cplLndPID, GidStr.c_str(), &tag_type, &ncomp, &tagInd );
            CHECKIERR( ierr, "failed to define global id tag" );

            int ent_type = 1;
            ierr         = iMOAB_GetIntTagStorage( cplLndPID, GidStr.c_str(), &nelem[2], &ent_type, &gidElems[0] );
            CHECKIERR( ierr, "failed to get global ids" );
            ierr = iMOAB_GetDoubleTagStorage( cplLndPID, field, &nelem[2], &ent_type,
                                              &tempElems[0] );
            CHECKIERR( ierr, "failed to get temperature field" );
            int err_code = 1;
            check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
            if( 0 == err_code )
                std::cout << " passed baseline test atm2ocn on ocean task " << rankInGlobalComm << "\n";
        }
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplRofLndPID );
        CHECKIERR( ierr, "cannot deregister app intx RL" )
    }
    if( rofComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cmpRofPID );
        CHECKIERR( ierr, "cannot deregister app " )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplRofPID );
        CHECKIERR( ierr, "cannot deregister app " )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplLndPID );
        CHECKIERR( ierr, "cannot deregister app " )
    }

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    // free rof coupler group and comm
    if( MPI_COMM_NULL != rofCouComm ) MPI_Comm_free( &rofCouComm );
    MPI_Group_free( &joinRofCouGroup );
    if( MPI_COMM_NULL != rofComm ) MPI_Comm_free( &rofComm );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &rofPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
