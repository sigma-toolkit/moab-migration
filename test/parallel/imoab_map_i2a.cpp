/*
 * This imoab_map_i2a test will simulate coupling between ice and atm
 * 2 meshes will be loaded from 2 files (ice on coupler, with ofrac), source,
 * and target atm file (recMeshAtm_PG2.h5m) file, and one map file read from disk
 * the migrate map mesh will be used to generate coverage set over target
 *  and will help for projection application of ofrac; fraction of ocean over atm
 *  We have some problems for this using intel compiler, in e3sm
 *  maybe we can replicate the differences we are seeing
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
//#include <iomanip>

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

    int cplIce        = 14,
        cplAtm        =  6;  // component ids are unique over all pes, and established in advance;

    std::string ice_mesh = "iceCplInit1Fr_P32.h5m";
    std::string atmFilename = "recMeshAtmPG_P32.h5m";

     std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );


    std::string mapFilename = "map_oQU480_to_ne4pg2_mono.200527.nc";  // this is a netcdf file!

    char field[] = "ofrac";  // this is a tag name, on the exported ice file
    // this will be projected and generate a baseline after sending it to coupler

    int rankInCouComm = -1;

    int nghlay = 0;  // number of ghost layers for loading the file
    std::vector< int > groupTasks;

    int startG4 = 0, endG4 = numProcesses - 1;  // these are for coupler layout

    // default: load ice / source on 2 proc, land / target on 2,
    // load map on 2 also, in parallel, distributed by rows

    ProgOptions opts;
    opts.addOpt< std::string >( "ice,s", "ice filename (source)", &ice_mesh );
    opts.addOpt< std::string >( "atm,t", "atm filename (target)", &atmFilename );
    opts.addOpt< std::string >( "map_file,w", "map file from source to target", &mapFilename );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    int disc_orders[2] = { 1, 1 };  // 1 is for FV

    std::string fieldstr;
    opts.addOpt< std::string >( "field,f", "field to project using the map ", &fieldstr );

    opts.parseCommandLine( argc, argv );

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";

    if( !rankInGlobalComm )
    {
        std::cout << "\n ice mesh file on coupler: " << ice_mesh << "\n   on tasks : " << startG4 << ":" << endG4
                  << "\n atm file on coupler " << atmFilename << "\n     on tasks : " << startG4 << ":" << endG4
                  << "\n map file:" << mapFilename << "\n     on tasks : " << startG4 << ":" << endG4 << "\n";
    }

    // load files on 2 different communicators, groups
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)
    // first groups has task 0, second group tasks 0 and 1
    // coupler will be on joint tasks, will be on a third group (0 and 1, again)

    // we will always have a coupler
    MPI_Group couPEGroup;
    MPI_Comm couComm;
    ierr = create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm );
    CHECKIERR( ierr, "Cannot create cpl MPI group and communicator " )


    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cplIceAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplIcePID = &cplIceAppID;  // ice on coupler PEs

    int cplAtmAppID = -1;
    iMOAB_AppID cplAtmPID    = &cplAtmAppID;     // Atm on coupler PEs
    int cplIceAtmAppID = -1;
    iMOAB_AppID cplIceAtmPID = &cplIceAtmAppID;  // map Ice -Atm on coupler PEs

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "IceX", &couComm, &cplIce,
                cplIcePID );  // ice on coupler pes
        CHECKIERR( ierr, "Cannot register ice over coupler PEs" )

        ierr = iMOAB_RegisterApplication( "AtmX", &couComm, &cplAtm,
                                          cplAtmPID );  // Atm on coupler pes
        CHECKIERR( ierr, "Cannot register Atm over coupler PEs" )
    }

    MPI_Barrier( MPI_COMM_WORLD );

    int IceAtmid = 100*14+6;
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "iceAtmMAP", &couComm, &IceAtmid, cplIceAtmPID );
        CHECKIERR( ierr, "Cannot register Ice2Atm map instance over coupler pes " )
    }


    // load Ice mesh, Atm mesh on coupler
    if (couComm != MPI_COMM_NULL)
    {
        ierr = iMOAB_LoadMesh( cplIcePID, ice_mesh.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load mesh file for Ice " )

        ierr = iMOAB_WriteMesh( cplIcePID, "IceCpl1.h5m", fileWriteOptions );
        CHECKIERR( ierr, "Cannot write Ice file on cpl " )
        ierr = iMOAB_LoadMesh( cplAtmPID, atmFilename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load domain file for lnd " )
    }

    const std::string intx_from_file_identifier = "map-from-file";

    if( couComm != MPI_COMM_NULL )
    {
        int src_disc_type = 3;  // element-based FV
        int tgt_disc_type = 3;  // element-based FV
        CHECKIERR( iMOAB_LoadMappingWeightsFromFile( cplIcePID, cplAtmPID, cplIceAtmPID, &src_disc_type, &tgt_disc_type,
                                                     intx_from_file_identifier.c_str(), mapFilename.c_str() ),
                   "failed to load map file from disk" );
        int type      = 3;  // FV
        // because it is like "coverage", context will be IceAtmid
        ierr = iMOAB_MigrateMapMesh( cplIcePID, cplIceAtmPID, &couComm, &couPEGroup, &couPEGroup, &type,
                                     &cplIce, &IceAtmid);
        CHECKIERR( ierr, "failed to migrate mesh for Ice on coupler" );
    }
    MPI_Barrier( MPI_COMM_WORLD );

    int tagIndex;
    int tagTypes[2]  = { DENSE_DOUBLE, DENSE_DOUBLE };
    int compOrder = disc_orders[0] * disc_orders[0] /*FV*/;
    int filter_type = 0;

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DefineTagStorage( cplIcePID, field, &tagTypes[0], &compOrder, &tagIndex );
        CHECKIERR( ierr, "failed to define the field tag" );

        ierr = iMOAB_DefineTagStorage( cplAtmPID, field, &tagTypes[1], &compOrder, &tagIndex );
        CHECKIERR( ierr, "failed to define the field tag on projection" );
    }

    // start the second hop, from Ice cpl to Ice coverage for Atm
    // the data is now on cpl Ice, need to be sent to Ice coverage over Atm
    PUSH_TIMER( "Send/receive data from Ice cpl to coverage in Atm context" )
    if( couComm != MPI_COMM_NULL )
    {
        // as always, use nonblocking sends
        // this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cplIcePID, field, &couComm, &IceAtmid );
        CHECKIERR( ierr, "cannot send tag values" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        // receive on atm on coupler pes, that was redistributed according to coverage
        // the trick is we use the map imoab app
        ierr = iMOAB_ReceiveElementTag( cplIceAtmPID, field, &couComm, &cplIce );
        CHECKIERR( ierr, "cannot receive tag values" )
    }


    // we can now free the sender buffers
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_FreeSenderBuffers( cplIcePID, &IceAtmid );  // context is for ocean
        CHECKIERR( ierr, "cannot free buffers " )
    }
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_WriteCoverageMesh( cplIceAtmPID, "Ice_cover_Atm");
        CHECKIERR( ierr, "cannot write coverage mesh" )
    }

    POP_TIMER( MPI_COMM_WORLD, rankInGlobalComm )


    if( couComm != MPI_COMM_NULL )
    {
        /* We have the remapping weights now. Let us apply the weights onto the tag we defined
           on the source mesh and get the projection on the target mesh */
        PUSH_TIMER( "Apply Scalar projection weights" )
        ierr = iMOAB_ApplyScalarProjectionWeights( cplIceAtmPID, &filter_type, intx_from_file_identifier.c_str(),
                                                   field, field );
        CHECKIERR( ierr, "failed to compute projection weight application" );
        POP_TIMER( couComm, rankInCouComm )

        {
            int numTasksCpl=endG4-startG4+1;
            std::ostringstream outfile;
            outfile << "fAtmOnCpl_" << numTasksCpl << ".h5m";
            ierr = iMOAB_WriteMesh( cplAtmPID, outfile.str().c_str(), fileWriteOptions );
            CHECKIERR( ierr, "could not write fAtmOnCpl5.h5m to disk" )
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplIceAtmPID );
        CHECKIERR( ierr, "cannot deregister app intx RL" )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplIcePID );
        CHECKIERR( ierr, "cannot deregister app " )
    }

    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplAtmPID );
        CHECKIERR( ierr, "cannot deregister app " )
    }

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();

    return 0;
}





