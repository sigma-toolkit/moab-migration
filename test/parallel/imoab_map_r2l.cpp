/*
 * This imoab_map_r2l test will simulate coupling between land and river
 * 2 meshes will be loaded from 2 files (land domain file, source,
 * and target scrip file, tgt), and one map file read from disk
 * the 
 * the migrate map mesh will be used to generate coverage set over target
 *  and will help for projection application
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

    std::string lndFilename = TestDir + "unittest/domain.lnd.ne4pg2_oQU480.200527.nc";
    std::string readopts_lnd( "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;REPARTITION" );

    int cplrof        = 22,
        cpllnd        = 10;  // component ids are unique over all pes, and established in advance;

    std::string rof_data = TestDir + "unittest/io/rof_comp_p32.h5m";

    std::string rof_mesh = TestDir + "unittest/io/SCRIPgrid_2x2_nomask_c210211.nc";
    std::string readopts_rof( "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN" );

    std::string mapFilename = TestDir + "unittest/map_r2_to_ne4pg2_mono.210211.nc";  // this is a netcdf file!

    std::string field_source = "Forr_rofl";  // this is a tag name, on the exported rof file
    // this will be projected and generate a baseline after sending it to coupler

    std::string baseline = TestDir + "unittest/baseline_lnd.txt";
    int rankInOcnComm    = -1;
    int cmprof = 21,
        roflndid = 1022;  // 100*10 + 22 component ids are unique over all pes, and established in advance;

    // we should modify the MigrateMapMesh to work with source coverage directly, like an intersection app

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;
    int startG1 = 0,  endG1 = numProcesses - 1;

    int startG4 = startG1, endG4 = endG1;  // these are for coupler layout
    int context_id;                        // used now for freeing buffers

    int repartitioner_scheme = 0;

    // default: load rof / source on 2 proc, land / target on 2,
    // load map on 2 also, in parallel, distributed by rows
    // probably all source cells will be involved in coverage mesh on both tasks

    ProgOptions opts;
    opts.addOpt< std::string >( "rof,s", "rof mesh scrip filename (source)", &rof_mesh );
    opts.addOpt< std::string >( "rofdata,d", "rof data filename (source)", &rof_data );
    opts.addOpt< std::string >( "lnd,t", "land domain mesh filename (target)", &lndFilename );
    opts.addOpt< std::string >( "map_file,w", "map file from source to target", &mapFilename );

    opts.addOpt< int >( "startAtm,a", "start task for rof layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for rof layout", &endG1 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    int types[2]       = { 3, 3 };  // type of source and target;  1 = SE, 2,= PC, 3 = FV
    int disc_orders[2] = { 1, 1 };  // 1 is for FV and PC; 4 could be for SE

    opts.addOpt< std::string >( "field,f", "field to project using the map ", &field_source );

    bool no_regression_test = false;
    opts.addOpt< void >( "no_regression,r", "do not do regression test against baseline 1", &no_regression_test );
    opts.addOpt< std::string >( "newbaseline,n", "baseline to use for test ", &baseline );
    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << " rof_data file: " << rof_data << "\n   on tasks : " << startG1 << ":" << endG1
                  << " rof_mesh scrip file on coupler: " << rof_mesh << "\n   on tasks : " << startG4 << ":" << endG4
                  << "\n lnd domain file on coupler " << lndFilename << "\n     on tasks : " << startG4 << ":" << endG4
                  << "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n" <<
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

    // rof_coupler
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
    if (rofComm != MPI_COMM_NULL)
    {
        std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");// load a point cloud with rof data
        ierr = iMOAB_LoadMesh( cmpRofPID, rof_data.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load component data file" )
    }

    // load rof mesh, lnd mesh on coupler
    if (couComm != MPI_COMM_NULL)
    {

    }

    const std::string intx_from_file_identifier = "map-from-file";

    if( couComm != MPI_COMM_NULL )
    {
        int src_disc_type = 3;  // element-based FV
        int tgt_disc_type = 3;  // element-based FV
        CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplRofPID, cplLndPID, cplRofLndPID, &src_disc_type, &tgt_disc_type,
                                                     intx_from_file_identifier.c_str(), mapFilename.c_str() ),
                   "failed to load map file from disk" );
        int type      = types[0];  // FV
        // because it is like "coverage", context will be atmocnid
        ierr = iMOAB_MigrateMapMesh( cplAtmPID, cplAtmOcnPID, &couComm, &couPEGroup, &couPEGroup, &type,
                                     &cplatm, &atmocnid);
        CHECKIERR( ierr, "failed to migrate mesh for atm on coupler" );
#ifdef VERBOSE
        if( *cplAtmPID >= 0 )
        {
            char prefix[] = "atmcov";
            ierr          = iMOAB_WriteLocalMesh( cplAtmPID, prefix );
            CHECKIERR( ierr, "failed to write local mesh" );
        }
#endif
    }
    MPI_Barrier( MPI_COMM_WORLD );

    int tagIndex[2];
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int atmCompNDoFs = disc_orders[0] * disc_orders[0], ocnCompNDoFs = disc_orders[1] * disc_orders[1] /*FV*/;
    int filter_type = 0;

    const char* bottomTempField          = "AnalyticalSolnSrcExact";
    const char* bottomTempProjectedField = "Target_proj";

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag AnalyticalSolnSrcExact" );

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs, &tagIndex[1] );
        CHECKIERR( ierr, "failed to define the field tag Target_proj" );
    }

    if( analytic_field && ( atmComm != MPI_COMM_NULL ) )  // we are on source /atm  pes
    {
        // cmpOcnPID, "T_proj;u_proj;v_proj;"
        ierr = iMOAB_DefineTagStorage( cmpAtmPID, bottomTempField, &tagTypes[0], &atmCompNDoFs, &tagIndex[0] );
        CHECKIERR( ierr, "failed to define the field tag AnalyticalSolnSrcExact" );

        int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];
        /*
         * Each process in the communicator will have access to a local mesh instance, which will contain the
         * original cells in the local partition and ghost entities. Number of vertices, primary cells, visible
         * blocks, number of sidesets and nodesets boundary conditions will be returned in numProcesses 3 arrays,
         * for local, ghost and total numbers.
         */
        ierr = iMOAB_GetMeshInfo( cmpAtmPID, nverts, nelem, nblocks, nsbc, ndbc );
        CHECKIERR( ierr, "failed to get num primary elems" );
        int numAllElem = nelem[2];
        int eetype     = 1;

        if( types[0] == 2 )  // point cloud
        {
            numAllElem = nverts[2];
            eetype     = 0;
        }
        std::vector< double > vals;
        int storLeng = atmCompNDoFs * numAllElem;
        vals.resize( storLeng );
        for( int k = 0; k < storLeng; k++ )
            vals[k] = k;

        ierr = iMOAB_SetDoubleTagStorage( cmpAtmPID, bottomTempField, &storLeng, &eetype, &vals[0] );
        CHECKIERR( ierr, "cannot make analytical tag" )
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
            int eetype     = 1;
            if( types[0] == 2 )  // Point cloud
            {
                eetype     = 0;  // vertices
                numAllElem = nverts[2];
            }
            std::vector< double > vals;
            int storLeng = atmCompNDoFs * numAllElem;

            vals.resize( storLeng );
            for( int k = 0; k < storLeng; k++ )
                vals[k] = 0.;

            ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomTempField, &storLeng, &eetype, &vals[0] );
            CHECKIERR( ierr, "cannot make tag nul" )

            // set the tag to 0
        }
    }

    const char* concat_fieldname  = field_source.c_str();
    const char* concat_fieldnameT = "Target_proj";

    {
        // first hop
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cmpAtmPID, concat_fieldname, &atmCouComm, &cplatm );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on atm on coupler pes
            ierr = iMOAB_ReceiveElementTag( cplAtmPID, concat_fieldname, &atmCouComm, &cmpatm );
            CHECKIERR( ierr, "cannot receive tag values" )
        }

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }

        // start the second hop, from atm cpl to atm coverage for ocn
        // the data is now on cpl Atm, need to be sent to atm coverage over ocean
        PUSH_TIMER( "Send/receive data from atm cpl to coverage in ocn context" )
        if( atmComm != MPI_COMM_NULL )
        {
            // as always, use nonblocking sends
            // this is for projection to ocean:
            ierr = iMOAB_SendElementTag( cplAtmPID, concat_fieldname, &couComm, &atmocnid );
            CHECKIERR( ierr, "cannot send tag values" )
        }
        if( couComm != MPI_COMM_NULL )
        {
            // receive on atm on coupler pes, that was redistributed according to coverage
            // the trick is we use the map imoab app
            ierr = iMOAB_ReceiveElementTag( cplAtmOcnPID, concat_fieldname, &couComm, &cplatm );
            CHECKIERR( ierr, "cannot receive tag values" )
        }

        // we can now free the sender buffers
        if( atmComm != MPI_COMM_NULL )
        {
            ierr = iMOAB_FreeSenderBuffers( cplAtmPID, &atmocnid );  // context is for ocean
            CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )
        }
        POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )


        if( couComm != MPI_COMM_NULL )
        {
            /* We have the remapping weights now. Let us apply the weights onto the tag we defined
               on the source mesh and get the projection on the target mesh */
            PUSH_TIMER( "Apply Scalar projection weights" )
            ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, &filter_type, intx_from_file_identifier.c_str(),
                                                       concat_fieldname, concat_fieldnameT );
            CHECKIERR( ierr, "failed to compute projection weight application" );
            POP_TIMER( couComm, rankInCouComm )

            {
                char outputFileTgt[] = "fOcnOnCpl5.h5m";
                ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileTgt, fileWriteOptions );
                CHECKIERR( ierr, "could not write fOcnOnCpl5.h5m to disk" )
            }
        }

        // send the projected tag back to ocean pes, with send/receive tag
        if( ocnComm != MPI_COMM_NULL )
        {
            int tagIndexIn2;
            ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomTempProjectedField, &tagTypes[1], &ocnCompNDoFs,
                                           &tagIndexIn2 );
            CHECKIERR( ierr, "failed to define the field tag for receiving back the tag "
                             "Target_proj on ocn pes" );
        }
        // send the tag to ocean pes, from ocean mesh on coupler pes
        //   from couComm, using common joint comm ocn_coupler
        // as always, use nonblocking sends
        // original graph (context is -1_
        if( couComm != MPI_COMM_NULL )
        {
            // need to use ocean comp id for context
            context_id = cmpocn;  // id for ocean on comp
            ierr       = iMOAB_SendElementTag( cplOcnPID, "Target_proj", &ocnCouComm, &context_id );
            CHECKIERR( ierr, "cannot send tag values back to ocean pes" )
        }

        // receive on component 2, ocean
        if( ocnComm != MPI_COMM_NULL )
        {
            context_id = cplocn;  // id for ocean on coupler
            ierr       = iMOAB_ReceiveElementTag( cmpOcnPID, "Target_proj", &ocnCouComm, &context_id );
            CHECKIERR( ierr, "cannot receive tag values from ocean mesh on coupler pes" )
        }

        if( couComm != MPI_COMM_NULL )
        {
            context_id = cmpocn;
            ierr       = iMOAB_FreeSenderBuffers( cplOcnPID, &context_id );
            CHECKIERR( ierr, "cannot free buffers for Target_proj tag migration " )
        }
        MPI_Barrier( MPI_COMM_WORLD );

        if( ocnComm != MPI_COMM_NULL )
        {
#ifdef VERBOSE
            char outputFileOcn[] = "OcnWithProj6.h5m";
            ierr                 = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions );
            CHECKIERR( ierr, "could not write OcnWithProj6.h5m to disk" )
#endif
            // test results only for n == 1, for bottomTempProjectedField
            if( !no_regression_test )
            {
                // the same as remap test
                // get temp field on ocean, from conservative, the global ids, and dump to the baseline file
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
                ierr = iMOAB_GetDoubleTagStorage( cmpOcnPID, bottomTempProjectedField, &nelem[2], &ent_type,
                                                  &tempElems[0] );
                CHECKIERR( ierr, "failed to get temperature field" );
                int err_code = 1;
                check_baseline_file( baseline, gidElems, tempElems, 1.e-9, err_code );
                if( 0 == err_code )
                    std::cout << " passed baseline test atm2ocn on ocean task " << rankInOcnComm << "\n";
            }
        }

    }  // end loop iterations n

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

    MPI_Group_free( &atmPEGroup );

    MPI_Group_free( &ocnPEGroup );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}
