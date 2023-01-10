/*
 * This test will load a file and duplicate it
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

using namespace moab;

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;

    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    std::string rofInp = TestDir + "unittest/recMeshOcn.h5m";
    std::string filename( "outmesh.h5m" );

    int cmpRof = 21, cplRof = 22;
    int nghlay = 0;  // no ghost layers

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );
    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    int startG1 = 0, startG4 = 0;
    int endG1, endG4;
    endG1 = endG4 = numProcesses - 1;

    ProgOptions opts;

    opts.addOpt< std::string >( "file,f", " imoab mesh file", &rofInp );

    opts.addOpt< std::string >( "outfile,o", "output mesh file", &filename );

    opts.addOpt< int >( "startAtm,a", "start task for input layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for input layout", &endG1 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.parseCommandLine( argc, argv );

    if( !rankInGlobalComm )
    {
        std::cout << " input file: " << rofInp << "\n   on tasks : " << startG1 << ":" << endG1
                  << "\n coupler    on tasks : " << startG4 << ":" << endG4 << "\n";
    }

    // load files on 2 different communicators, groups
    // coupler will be on group 4

    MPI_Group rofPEGroup;
    MPI_Comm rofComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &rofPEGroup, &rofComm );
    CHECKIERR( ierr, "Cannot create rof MPI group and communicator " )

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

    int cmpRofID       = -1;
    iMOAB_AppID rofPID = &cmpRofID;
    if( rofComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "ROF", &rofComm, &cmpRof, rofPID );
        CHECKIERR( ierr, "Cannot register Rof App" )
    }

    int cplRofAppID       = -1;
    iMOAB_AppID cplRofPID = &cplRofAppID;

    int rankInCouComm = -1;
    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs

        ierr = iMOAB_RegisterApplication( "ROFX", &couComm, &cplRof,
                                          cplRofPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register ROFX over coupler PEs" )
    }

    // load atm mesh and migrate, not used actually
    int repartitioner_scheme = 2;  // zoltan is used
    if( rofComm != MPI_COMM_NULL )
    {
        ierr =
            setup_component_coupler_meshes( rofPID, cmpRof, cplRofPID, cplRof, &rofComm, &rofPEGroup, &couComm,
                                            &couPEGroup, &rofCouComm, rofInp, readopts, nghlay, repartitioner_scheme );
        CHECKIERR( ierr, "Cannot load and migrate rof mesh " )
    }

    int cplCopyAppID       = -1;
    iMOAB_AppID cplCopyPID = &cplCopyAppID;
    int copyId             = 100 + cplRof;
    // make a copy on coupler comm
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "COPY", &couComm, &copyId,
                                          cplCopyPID );  // copy on coupler pes
        CHECKIERR( ierr, "Cannot register COPY app over coupler PEs" )

        ierr = iMOAB_DuplicateAppMesh( cplRofPID, cplCopyPID );
        CHECKIERR( ierr, "Cannot duplicate mesh over coupler PEs" )

        char fileWriteOptions[] = "PARALLEL=WRITE_PART";
        ierr                    = iMOAB_WriteMesh( cplCopyPID, filename.c_str(), fileWriteOptions );
        CHECKIERR( ierr, "cannot write duplicated mesh" )
    }

    // we could deregister cplLndAtmPID
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplCopyPID );
        CHECKIERR( ierr, "cannot deregister copy app" )
    }

    // we could deregister cplRofPID
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( cplRofPID );
        CHECKIERR( ierr, "cannot deregister coupler app" )
    }

    if( rofComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( rofPID );
        CHECKIERR( ierr, "cannot deregister app rofPID " )
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
