/**
 * @file TestErrorHandling.cpp
 * @brief Example demonstrating MOAB's trace back error handler in serial
 *
 * This example shows how to:
 * - Initialize and finalize MOAB's error handler
 * - Simulate different types of errors in serial execution
 * - Handle various MOAB error codes and error propagation
 * - Test error handling for file loading, tag creation, and tag iteration
 *
 * The example demonstrates four test cases:
 * - Test case 1: Error MB_NOT_IMPLEMENTED (unsupported variable on edges)
 * - Test case 2: Error MB_TYPE_OUT_OF_RANGE (invalid GATHER_SET option)
 * - Test case 3: Error MB_FAILURE (NOMESH option with NULL file set)
 * - Test case 4: Error MB_VARIABLE_DATA_LENGTH (variable-length tag iteration)
 *
 * @author MOAB Development Team
 * @date 2024
 *

 * Description: This example tests MOAB's trace back error handler in serial. \n
 *
 * <b>To run</b>: ./TestErrorHandling <test_case_num(1 to 4)> \n
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */

#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>
#include <memory>

using namespace moab;
using namespace std;

// In this test case, an error MB_NOT_IMPLEMENTED is returned by MOAB
ErrorCode TestErrorHandling_1()
{
    Core moab;
    Interface& mb = moab;

    // Load a CAM-FV file and read a variable on edges (not supported yet)
    string test_file = string( MESH_DIR ) + string( "/io/fv3x46x72.t.3.nc" );
    MB_CHK_ERR( mb.load_file( test_file.c_str(), NULL, "VARIABLE=US" ) );

    return MB_SUCCESS;
}

// In this test case, an error MB_TYPE_OUT_OF_RANGE is returned by MOAB
ErrorCode TestErrorHandling_2()
{
    Core moab;
    Interface& mb = moab;

    // Load a HOMME file with an invalid GATHER_SET option
    string test_file = string( MESH_DIR ) + string( "/io/homme3x3458.t.3.nc" );
    MB_CHK_ERR( mb.load_file( test_file.c_str(), NULL, "VARIABLE=T;GATHER_SET=0.1" ) );

    return MB_SUCCESS;
}

// In this test case, an error MB_FAILURE is returned by MOAB
ErrorCode TestErrorHandling_3()
{
    Core moab;
    Interface& mb = moab;

    // Load a CAM-FV file with NOMESH option and a NULL file set
    string test_file = string( MESH_DIR ) + string( "/io/fv3x46x72.t.3.nc" );
    MB_CHK_ERR( mb.load_file( test_file.c_str(), NULL, "NOMESH;VARIABLE=" ) );

    return MB_SUCCESS;
}

// In this test case, an error MB_VARIABLE_DATA_LENGTH is returned by MOAB
ErrorCode TestErrorHandling_4()
{
    Core moab;
    Interface& mb = moab;

    // Create 100 vertices
    const int NUM_VTX = 100;
    vector< double > coords( 3 * NUM_VTX );
    Range verts;
    ErrorCode rval = mb.create_vertices( &coords[0], NUM_VTX, verts );MB_CHK_SET_ERR( rval, "Failed to create vertices" );

    // Create a variable-length dense tag
    Tag tag;
    rval = mb.tag_get_handle( "var_len_den", 1, MB_TYPE_INTEGER, tag, MB_TAG_VARLEN | MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Failed to create a tag" );

    // Attempt to iterate over a variable-length tag, which will never be possible
    void* ptr = NULL;
    int count = 0;
    rval      = mb.tag_iterate( tag, verts.begin(), verts.end(), count, ptr );MB_CHK_SET_ERR( rval, "Failed to iterate over tag on " << NUM_VTX << " vertices" );

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

    // Initialize error handler, optional for this example (using moab instances)
    MBErrorHandler_Init();

    ErrorCode rval = MB_SUCCESS;

    int test_case_num = atoi( argv[1] );
    switch( test_case_num )
    {
        case 1:
            rval = TestErrorHandling_1();MB_CHK_ERR( rval );
            break;
        case 2:
            rval = TestErrorHandling_2();MB_CHK_ERR( rval );
            break;
        case 3:
            rval = TestErrorHandling_3();MB_CHK_ERR( rval );
            break;
        case 4:
            rval = TestErrorHandling_4();MB_CHK_ERR( rval );
            break;
        default:
            break;
    }

    // Finalize error handler, optional for this example (using moab instances)
    MBErrorHandler_Finalize();

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
