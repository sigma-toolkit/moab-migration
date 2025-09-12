/**
 * @file test_distributed_reduction.cpp
 * @brief Test suite for the distributed sparse matrix reduction functionality.
 *
 * This program tests the `determine_communication_pattern` and `perform_distributed_reduction` functions.
 * It sets up a parallel environment using MPI and performs the following steps:
 * 1. Initializes MPI and MOAB.
 * 2. Creates a distributed sparse matrix where different MPI ranks own overlapping rows, simulating a
 *    scenario where data is shared across process boundaries.
 * 3. Prints the initial state of the distributed matrix.
 * 4. Calls `determine_communication_pattern` to discover which ranks need to communicate based on the
 *    shared rows and prints the resulting communication pattern for verification.
 * 5. Calls `perform_distributed_reduction` to execute the in-place reduction of the shared rows.
 * 6. Prints the final state of the matrix to verify that:
 *    - Shared rows have been correctly summed across all participating ranks.
 *    - Non-shared (local) rows have been preserved without modification.
 * 7. Cleans up resources and finalizes MPI.
 */
#include "determine_communication_pattern.hpp"
#include "moab/Core.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <mpi.h>

using namespace moab;

/**
 * @brief A helper function to print the contents of a distributed sparse matrix in a synchronized manner.
 *
 * This function iterates through all MPI ranks, and each rank prints its local portion of the matrix to standard output.
 * An `MPI_Barrier` is used after each rank prints to ensure that the output from different ranks does not get jumbled,
 * making the overall matrix structure easy to inspect.
 *
 * @param matrix The local sparse matrix on the current rank.
 * @param rank The MPI rank of the current process.
 * @param size The total number of MPI ranks.
 */
void print_distributed_matrix(
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
    int rank,
    int size)
{
    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            std::cout << "Rank " << rank << " local matrix (" << matrix.rows() << "x"
                      << matrix.cols() << ", " << matrix.nonZeros() << " non-zeros):" << std::endl;
            for (int k = 0; k < matrix.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, k); it; ++it) {
                    if (fabs(it.value()) > 1e-12) {
                        std::cout << "  (" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize MOAB and the ParallelComm interface
    moab::Core* mb = new moab::Core();
    moab::ParallelComm* pcomm = new moab::ParallelComm(mb, MPI_COMM_WORLD);

    try {
        // Define the global dimensions of the test matrix.
        const int global_rows = 15;
        const int global_cols = 15;

        // --- Create a distributed test case with overlapping rows ---
        // Each process is assigned a block of rows. To create sharing, the last row of each rank
        // (except the last rank) is made to overlap with the first row of the next rank.
        int rows_per_proc = (global_rows + size - 1) / size;
        int my_first_row = rank * rows_per_proc;
        int my_last_row = std::min((rank + 1) * rows_per_proc, global_rows) - 1;

        // Introduce overlap: rank `r` shares its last row with rank `r+1`.
        if (rank < size - 1) my_last_row++;

        my_first_row = std::max(0, my_first_row);
        my_last_row = std::min(global_rows - 1, my_last_row);

        // -- Populate the local sparse matrix --
        // Each rank creates a tri-diagonal matrix for its assigned rows.
        // The values are simple: 1.0 on the main diagonal and 0.1 on the off-diagonals.
        // An extra entry is added to a shared row to make the test case more robust.
        Eigen::SparseMatrix<double, Eigen::RowMajor> local_matrix(global_rows, global_cols);
        typedef Eigen::Triplet< double > Triplet;
        std::vector< Triplet > tripletList;
        tripletList.reserve( 3*(my_last_row-my_first_row+1) );

        for (int i = my_first_row; i <= my_last_row; ++i) {
            std::vector<double> row(global_cols, 0.0);

            // Set diagonal and some off-diagonal elements
            for (int j = std::max(0, i-1); j <= std::min(global_cols-1, i+1); ++j) {
                row[j] = (i == j) ? 1.0 : 0.1;
                tripletList.push_back(Triplet(i,j,row[j]));
            }
            // Add an extra element to the last row of this rank to test non-symmetric contributions.
            if (i == my_last_row) {
                tripletList.push_back(Triplet(i,i-2,1.0));
            }
        }

        local_matrix.setFromTriplets( tripletList.begin(), tripletList.end() );
        local_matrix.makeCompressed();

        // Print the initial state of the matrix across all ranks.
        if (rank == 0) {
            std::cout << "Initial distributed matrix:" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        print_distributed_matrix(local_matrix, rank, size);

        // --- Determine the communication pattern based on row distribution ---
        CommunicationPattern comm_pattern;
        std::map<int, std::set<int>> all_row_sharers;
        ErrorCode rval = determine_communication_pattern(pcomm, local_matrix, comm_pattern, all_row_sharers);
        MB_CHK_SET_ERR(rval, "determine_communication_pattern failed");

        // Print the determined communication pattern for verification.
        for (int i = 0; i < size; ++i) {
            if (rank == i) {
                printf("Rank %d Communication Pattern:\n", rank);
                if (!comm_pattern.send_map.empty()) {
                    for (const auto& pair : comm_pattern.send_map) {
                        printf("  Send to %d, rows: ", pair.first);
                        for (int row : pair.second) {
                            printf("%d ", row);
                        }
                        printf("\n");
                    }
                }
                if (!comm_pattern.recv_map.empty()) {
                    for (const auto& pair : comm_pattern.recv_map) {
                        printf("  Receive from %d, rows: ", pair.first);
                        for (int row : pair.second) {
                            printf("%d ", row);
                        }
                        printf("\n");
                    }
                }
            }
            MPI_Barrier(pcomm->comm());
        }

        // --- Perform the distributed reduction in-place ---
        perform_distributed_reduction(pcomm, comm_pattern, all_row_sharers, local_matrix);
        MB_CHK_SET_ERR(rval, "perform_distributed_reduction failed");

        // Print the final state of the matrix after reduction.
        if (rank == 0) {
            std::cout << "\nAfter reduction:" << std::endl;
        }
        print_distributed_matrix(local_matrix, rank, size);

        // --- Cleanup ---
        delete pcomm;
        delete mb;

        MPI_Finalize();
        return 0;
    }
    catch (std::exception& e) {
        std::cerr << "Rank " << rank << " Error: " << e.what() << std::endl;
        delete pcomm;
        delete mb;
        MPI_Finalize();
        return 1;
    }
}
