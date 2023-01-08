/*
 * This test will load a scrip file in parallel and write it back
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
    std::string readopts( "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN" );
    std::string readopts2( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
    std::string readopts3("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION");
    std::string filename=TestDir + "unittest/SCRIPgrid_2x2_nomask_c210211.nc";
    std::string atmFilename=TestDir + "unittest/wholeATM_T.h5m";
    std::string rofInp=TestDir + "unittest/wholeRof_06.h5m";
    int cmpAtm = 5, cmpRof = 21, cplRof=22;
    int cplatm        = 6;  // component ids are unique over all pes, and established in advance;
    int nghlay = 0;// no ghost layers


    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );
    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );  // all processes in jgroup

    int startG1 = 0, startG2 = 0, startG4 = 0;
    int endG1, endG2, endG4;
    endG1 = endG2 = endG4 = numProcesses-1;

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "atm mesh filename ", &atmFilename );

    opts.addOpt< std::string >( "mosart,m", " mosart with data", &rofInp );


    opts.addOpt< std::string >( "scrip,s", "scrip mesh file", &filename );

    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );

    opts.addOpt< int >( "startOcn,c", "start task for mosart layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for mosart layout", &endG2 );

    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );

    opts.parseCommandLine( argc, argv );

    if( !rankInGlobalComm )
    {
        std::cout << " atm file: " << atmFilename << "\n   on tasks : " << startG1 << ":" << endG1 <<
            "\n mosart input file file: " << rofInp << "\n     on tasks : " << startG2 << ":" << endG2 <<
            "\n scrip file on coupler: " << filename <<
            "\n coupler    on tasks : " << startG4 << ":" << endG4 << "\n";
    }

    // load files on 2 different communicators, groups
    // coupler will be on group 4
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    ierr = create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm );
    CHECKIERR( ierr, "Cannot create atm MPI group and communicator " )

    MPI_Group rofPEGroup;
    MPI_Comm rofComm;
    ierr = create_group_and_comm( startG2, endG2, jgroup, &rofPEGroup, &rofComm );
    CHECKIERR( ierr, "Cannot create rof MPI group and communicator " )


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

    // rof_coupler
    MPI_Group joinRofCouGroup;
    MPI_Comm rofCouComm;
    ierr = create_joint_comm_group( rofPEGroup, couPEGroup, &joinRofCouGroup, &rofCouComm );
    CHECKIERR( ierr, "Cannot create joint rof cou communicator" )

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    int cmpRofID       = -1;
    iMOAB_AppID rofPID = &cmpRofID;
    if (rofComm != MPI_COMM_NULL ) {
        ierr = iMOAB_RegisterApplication( "ROF", &rofComm, &cmpRof, rofPID );
        CHECKIERR( ierr, "Cannot register Rof App" )
    }

    int cmpAtmAppID       = -1;
    iMOAB_AppID cmpAtmPID = &cmpAtmAppID;
    if (atmComm != MPI_COMM_NULL) {
        ierr = iMOAB_RegisterApplication( "ATM", &atmComm, &cmpAtm, cmpAtmPID );
        CHECKIERR( ierr, "Cannot register Atm App" )
    }
    int cplAtmAppID       = -1;
    iMOAB_AppID cplAtmPID = &cplAtmAppID;

    int cplRofAppID = -1;
    iMOAB_AppID  cplRofPID = &cplRofAppID;

    int rankInCouComm = -1;
    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        // Register all the applications on the coupler PEs
        ierr = iMOAB_RegisterApplication( "ATMX", &couComm, &cplatm,
                                          cplAtmPID );  // atm on coupler pes
        CHECKIERR( ierr, "Cannot register ATM over coupler PEs" )

        ierr = iMOAB_RegisterApplication( "ROFX", &couComm, &cplRof,
                                          cplRofPID );  // ocn on coupler pes
        CHECKIERR( ierr, "Cannot register ROFX over coupler PEs" )
    }

    // load atm mesh and migrate, not used actually
    int repartitioner_scheme = 2; // zoltan is used
    /*ierr =
        setup_component_coupler_meshes( cmpAtmPID, cmpAtm, cplAtmPID, cplatm, &atmComm, &atmPEGroup, &couComm,
                                        &couPEGroup, &atmCouComm, atmFilename, readopts, nghlay, repartitioner_scheme );*/

    if( cmpRofID >= 0 ) {
        // load  rof mesh with data on it
        ierr = iMOAB_LoadMesh( rofPID, rofInp.c_str(), readopts3.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load mosart data mesh" )
    }
    // load rof scrip file on coupler only
    if( couComm != MPI_COMM_NULL )
    {
        ierr = iMOAB_LoadMesh( cplRofPID, filename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load atm component mesh" )
        // test what we read from scrip file
        char outputFileTgt[] = "readCplRof.h5m";
        char fileWriteOptions[] = "PARALLEL=WRITE_PART;DEBUG_IO=2";
        ierr                  = iMOAB_WriteMesh( cplRofPID, outputFileTgt, fileWriteOptions);
        CHECKIERR( ierr, "cannot write Rof mesh after receiving" )
    }
    // compute comm graph between coupler and wholeRof

    MPI_Finalize();

    return 0;
}
