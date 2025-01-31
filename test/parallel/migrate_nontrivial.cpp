/*
 * migrate_nontrivial.cpp
 *
 *  migrate_nontrivial will contain tests for migrating meshes in parallel environments, with iMOAB
 * api, using nontrivial partitions; for starters, we will use zoltan to compute partition in
 * parallel, and then use the result to migrate the mesh trivial partition gave terrible results for
 * mpas type meshes, that were numbered like a fractal; the resulting migrations were like Swiss
 * cheese, full of holes a mesh is read on senders tasks we will use graph like methods or geometry
 * methods, and will modify the ZoltanPartitioner to add needed methods
 *
 *  mesh will be sent to receivers tasks, with nonblocking MPI_Isend calls, and then received
 *    with blocking MPI_Recv calls;
 *
 *  we will not modify the GLOBAL_ID tag, we assume it was set correctly before we started
 */

#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/ProgOptions.hpp"

#define RUN_TEST_ARG2( A, B ) run_test( &( A ), #A, B )

using namespace moab;

//#define GRAPH_INFO

#define CHECKRC( rc, message )            \
    if( 0 != ( rc ) )                     \
    {                                     \
        printf( "Error: %s\n", message ); \
        return MB_FAILURE;                \
    }

struct RunContext
{
    std::string filename;  // name of mesh file
    int startG1;           // starting process for group 1
    int startG2;           // starting process for group 2
    int endG1;             // ending process for group 1
    int endG2;             // ending process for group 2
    MPI_Comm jcomm;        // will be a copy of the global communicator
    MPI_Group jgroup;      // will be a copy of the global group
    int rank, size;        // current rank and total PE size
};

int is_any_proc_error( int is_my_error )
{
    int result = 0;
    int err    = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    return err || result;
}

int run_test( ErrorCode ( *func )( RunContext& ), const char* func_name, RunContext& context )
{
    ErrorCode result = ( *func )( context );
    int is_err       = is_any_proc_error( ( MB_SUCCESS != result ) );
    if( context.rank == 0 )
    {
        if( is_err )
            std::cout << func_name << " : FAILED!!" << std::endl;
        else
            std::cout << func_name << " : success" << std::endl;
    }

    return is_err;
}

ErrorCode migrate_graph( RunContext& );
ErrorCode migrate_geom( RunContext& );
ErrorCode migrate_trivial( RunContext& );

ErrorCode migrate_smart( RunContext& context, const char* outfile, int partMethod )
{
    int ierr;
    int compid1, compid2;           // component ids are unique over all pes, and established in advance;
    std::vector< int > groupTasks;  // at most 4 tasks

    // first create MPI groups
    MPI_Group group1, group2;
    groupTasks.resize( context.endG1 - context.startG1 + 1 );
    for( int i = context.startG1; i <= context.endG1; i++ )
        groupTasks[i - context.startG1] = i;

    ierr = MPI_Group_incl( context.jgroup, context.endG1 - context.startG1 + 1, &groupTasks[0], &group1 );
    CHECKRC( ierr, "can't create group1" )

    groupTasks.resize( context.endG2 - context.startG2 + 1 );
    for( int i = context.startG2; i <= context.endG2; i++ )
        groupTasks[i - context.startG2] = i;

    ierr = MPI_Group_incl( context.jgroup, context.endG2 - context.startG2 + 1, &groupTasks[0], &group2 );
    CHECKRC( ierr, "can't create group2" )

    // create 2 communicators, one for each group
    int tagcomm1 = 1, tagcomm2 = 2;
    int context_id = -1;  // plain migrate, default context
    MPI_Comm comm1, comm2;
    ierr = MPI_Comm_create_group( context.jcomm, group1, tagcomm1, &comm1 );
    CHECKRC( ierr, "can't create comm1" )

    ierr = MPI_Comm_create_group( context.jcomm, group2, tagcomm2, &comm2 );
    CHECKRC( ierr, "can't create comm2" )

    ierr = iMOAB_Initialize( 0, 0 );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKRC( ierr, "can't initialize iMOAB" )

    // give some dummy values to component ids, just to differentiate between them
    // the par comm graph is unique between components
    compid1 = 4;
    compid2 = 7;

    int appID1;
    iMOAB_AppID pid1 = &appID1;
    int appID2;
    iMOAB_AppID pid2 = &appID2;

    if( comm1 != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "APP1", &comm1, &compid1, pid1 );
        CHECKRC( ierr, "can't register app1 " )
    }
    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "APP2", &comm2, &compid2, pid2 );
        CHECKRC( ierr, "can't register app2 " )
    }

    if( comm1 != MPI_COMM_NULL )
    {

        std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );
        int nghlay = 0;  // number of ghost layers for loading the file

        ierr = iMOAB_LoadMesh( pid1, context.filename.c_str(), readopts.c_str(), &nghlay );
        CHECKRC( ierr, "can't load mesh " )
        ierr = iMOAB_SendMesh( pid1, &context.jcomm, &group2, &compid2, &partMethod );  // send to component 2
        CHECKRC( ierr, "cannot send elements" )
#ifdef GRAPH_INFO
        int is_sender = 1;
        int context   = compid2;
        iMOAB_DumpCommGraph( pid1, &context, &is_sender, "MigrateS" );
#endif
    }

    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_ReceiveMesh( pid2, &context.jcomm, &group1, &compid1 );  // receive from component 1
        CHECKRC( ierr, "cannot receive elements" )
        std::string wopts;
        wopts = "PARALLEL=WRITE_PART;";
        ierr  = iMOAB_WriteMesh( pid2, outfile, wopts.c_str() );
        CHECKRC( ierr, "cannot write received mesh" )
#ifdef GRAPH_INFO
        int is_sender = 0;
        int context   = compid1;
        iMOAB_DumpCommGraph( pid2, &context, &is_sender, "MigrateR" );
#endif
    }

    MPI_Barrier( context.jcomm );

    // we can now free the sender buffers
    context_id = compid2;  // even for default migrate, be more explicit
    if( comm1 != MPI_COMM_NULL ) ierr = iMOAB_FreeSenderBuffers( pid1, &context_id );

    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( pid2 );
        CHECKRC( ierr, "cannot deregister app 2 receiver" )
    }

    if( comm1 != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( pid1 );
        CHECKRC( ierr, "cannot deregister app 1 sender" )
    }

    ierr = iMOAB_Finalize();
    CHECKRC( ierr, "did not finalize iMOAB" )

    if( MPI_COMM_NULL != comm1 ) MPI_Comm_free( &comm1 );
    if( MPI_COMM_NULL != comm2 ) MPI_Comm_free( &comm2 );

    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    return MB_SUCCESS;
}

// migrate from 2 tasks to 3 tasks
ErrorCode migrate_graph( RunContext& context )
{
    return migrate_smart( context, "migrate_graph.h5m", 1 );
}

ErrorCode migrate_geom( RunContext& context )
{
    return migrate_smart( context, "migrate_geom.h5m", 2 );
}

ErrorCode migrate_trivial( RunContext& context )
{
    return migrate_smart( context, "migrate_trivial.h5m", 0 );
}

int main( int argc, char* argv[] )
{
    // some global variables, used by all tests
    ProgOptions opts;
    int typeTest = 2;  // default: geometric partitioner
    RunContext context;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &context.rank );
    MPI_Comm_size( MPI_COMM_WORLD, &context.size );

    MPI_Comm_dup( MPI_COMM_WORLD, &context.jcomm );
    MPI_Comm_group( context.jcomm, &context.jgroup );

    // set default run arguments
    context.filename = TestDir + "unittest/field1.h5m";
    context.startG1  = 0;
    context.startG2  = 0;
    context.endG1    = 0;
    context.endG2    = 1;

    opts.addOpt< std::string >( "file,f", "source file", &context.filename );
    opts.addOpt< int >( "startSender,a", "start task for source layout", &context.startG1 );
    opts.addOpt< int >( "endSender,b", "end task for source layout", &context.endG1 );
    opts.addOpt< int >( "startRecv,c", "start task for receiver layout", &context.startG2 );
    opts.addOpt< int >( "endRecv,d", "end task for receiver layout", &context.endG2 );

    opts.addOpt< int >( "typeTest,t", "test types (0 - trivial, 1 graph, 2 geom, 3 both  graph and geometry",
                        &typeTest );

    opts.parseCommandLine( argc, argv );

    if( context.rank == 0 )
    {
        std::cout << " input file : " << context.filename << "\n";
        std::cout << " sender   on tasks: " << context.startG1 << ":" << context.endG1 << "\n";
        std::cout << " receiver on tasks: " << context.startG2 << ":" << context.endG2 << "\n";
        std::cout << " type migrate: " << typeTest << " (0 - trivial, 1 graph , 2 geom, 3 both graph and geom  ) \n";
    }

    int num_errors = 0;

    if( 0 == typeTest ) num_errors += RUN_TEST_ARG2( migrate_trivial, context );
    if( 3 == typeTest || 1 == typeTest ) num_errors += RUN_TEST_ARG2( migrate_graph, context );
    if( 3 == typeTest || 2 == typeTest ) num_errors += RUN_TEST_ARG2( migrate_geom, context );

    if( context.rank == 0 )
    {
        if( !num_errors )
            std::cout << "All tests passed" << std::endl;
        else
            std::cout << num_errors << " TESTS FAILED!" << std::endl;
    }

    MPI_Group_free( &context.jgroup );
    MPI_Comm_free( &context.jcomm );
    MPI_Finalize();
    return num_errors;
}

#undef VERBOSE
