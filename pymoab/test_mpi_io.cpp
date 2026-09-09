#include <mpi.h>
#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    std::cout << "MPI I/O Test: Rank " << rank << " of " << size << std::endl;

    // Test MPI I/O functionality
    MPI_File file;
    MPI_Offset offset = rank * 1024;

    int ierr = MPI_File_open( MPI_COMM_WORLD, "test_mpiio.tmp", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file );

    if( ierr == MPI_SUCCESS )
    {
        std::cout << "Rank " << rank << ": MPI_File_open successful" << std::endl;
        MPI_File_close( &file );
    }
    else
    {
        std::cout << "Rank " << rank << ": MPI_File_open failed" << std::endl;
    }

    MPI_Finalize();
    return 0;
}