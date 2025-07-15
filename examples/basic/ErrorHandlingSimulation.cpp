/**
 * @file ErrorHandlingSimulation.cpp
 * @brief Example demonstrating MOAB's enhanced error handling in parallel
 *
 * This example shows how to:
 * - Initialize and finalize MOAB's error handler
 * - Simulate different types of errors in parallel execution
 * - Handle global fatal errors vs per-processor errors
 * - Use error propagation through function call hierarchy
 * - Demonstrate error reporting behavior across processors
 *
 * The example demonstrates four test cases:
 * - Test case 1: Global fatal error (MB_NOT_IMPLEMENTED) on all processors
 * - Test case 2: Per-processor error (MB_INDEX_OUT_OF_RANGE) on all processors
 * - Test case 3: Per-processor error (MB_TYPE_OUT_OF_RANGE) on non-root processors
 * - Test case 4: Different errors on specific processors (1 and 3)
 *
 * @author MOAB Development Team
 * @date 2024
 *

 * Description: This example simulates MOAB's enhanced error handling in parallel. \n
 * All of the errors are contrived, used for simulation purpose only. \n
 *
 * Note: We do not need a moab instance for this example
 *
 * <b>To run</b>: mpiexec -np 4 ./ErrorHandlingSimulation <test_case_num(1 to 4)> \n
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */

#include "moab/MOABConfig.h"
#include "moab/ErrorHandler.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>
#include <cstdlib>

using namespace moab;
using namespace std;

// Functions that create and handle contrived errors
// Call hierarchy: A calls B, and B calls C
ErrorCode FunctionC( int test_case_num, int rank )
{
    switch( test_case_num )
    {
        case 1:
            // Simulate a global fatal error MB_NOT_IMPLEMENTED on all processors
            // Note, it is printed by root processor 0 only
            MB_SET_GLB_ERR( MB_NOT_IMPLEMENTED, "A contrived global error MB_NOT_IMPLEMENTED" );
            break;
        case 2:
            // Simulate a per-processor relevant error MB_INDEX_OUT_OF_RANGE on all processors
            // Note, it is printed by all processors
            MB_SET_ERR( MB_INDEX_OUT_OF_RANGE, "A contrived error MB_INDEX_OUT_OF_RANGE on processor " << rank );
            break;
        case 3:
            // Simulate a per-processor relevant error MB_TYPE_OUT_OF_RANGE on all processors except
            // root Note, it is printed by all non-root processors
            if( 0 != rank )
                MB_SET_ERR( MB_TYPE_OUT_OF_RANGE, "A contrived error MB_TYPE_OUT_OF_RANGE on processor " << rank );
            break;
        case 4:
            // Simulate a per-processor relevant error MB_INDEX_OUT_OF_RANGE on processor 1
            // Note, it is printed by processor 1 only
            if( 1 == rank )
                MB_SET_ERR( MB_INDEX_OUT_OF_RANGE, "A contrived error MB_INDEX_OUT_OF_RANGE on processor 1" );

            // Simulate a per-processor relevant error MB_TYPE_OUT_OF_RANGE on processor 3
            // Note, it is printed by processor 3 only
            if( 3 == rank ) MB_SET_ERR( MB_TYPE_OUT_OF_RANGE, "A contrived error MB_TYPE_OUT_OF_RANGE on processor 3" );
            break;
        default:
            break;
    }

    return MB_SUCCESS;
}

ErrorCode FunctionB( int test_case_num, int rank )
{
    MB_CHK_ERR( FunctionC( test_case_num, rank ) );

    return MB_SUCCESS;
}

ErrorCode FunctionA( int test_case_num, int rank )
{
    MB_CHK_ERR( FunctionB( test_case_num, rank ) );

    return MB_SUCCESS;
}

int main( int argc, char** argv )
{
    if( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <test_case_num(1 to 4)>" << endl;
        return 0;
    }

#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    // Initialize error handler, required for this example (not using a moab instance)
    MBErrorHandler_Init();

    int test_case_num = atoi( argv[1] );
    int rank          = 0;
#ifdef MOAB_HAVE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif

    MB_CHK_ERR( FunctionA( test_case_num, rank ) );

    // Finalize error handler, required for this example (not using a moab instance)
    MBErrorHandler_Finalize();

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
