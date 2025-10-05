/** @example HelloParMOAB.cpp
 * @brief Demonstrates parallel mesh operations with MOAB, including shared entity resolution and ghost exchange.
 *
 * This example shows how to:
 * - Load a mesh in parallel with specific read options
 * - Work with MPI communicators and process groups
 * - Identify and report shared entities across processors
 * - Exchange ghost elements between processors
 * - Gather and report statistics about the parallel decomposition
 *
 * @note This example requires MOAB to be compiled with MPI and HDF5 support.
 *
 * @par Usage:
 * @code
 * mpiexec -np <num_procs> HelloParMOAB [filename] [num_communicators]
 * @endcode
 *
 * @par Example:
 * @code
 * mpiexec -np 8 HelloParMOAB mesh.h5m 2
 * @endcode
 *
 * This will load the mesh on 8 processes, split into 2 communicator groups (4 processes each).
 *
 * @see Core, ParallelComm, Range, MBParallelConventions
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBParallelConventions.h"
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>

// Only include MPI-specific headers if MOAB was built with MPI support
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include <mpi.h>
#endif

// Using declarations for cleaner code
using moab::Core;
using moab::ErrorCode;
using moab::Range;
using moab::ParallelComm;

namespace {
    // Default mesh file path
    const char* const DEFAULT_MESH_FILE = "/64bricks_512hex_256part.h5m";

    // Default read options
    const char* const DEFAULT_READ_OPTS =
        "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

    // Print usage information
    void print_usage(const char* program_name) {
        std::cout << "Usage: mpiexec -np <num_procs> " << program_name
                  << " [meshfile] [num_communicators]\n"
                  << "\nOptions:"
                  << "\n  meshfile         - Path to the mesh file (default: " << MESH_DIR << DEFAULT_MESH_FILE << ")"
                  << "\n  num_communicators - Number of MPI communicator groups to create (default: 1)\n";
    }

    // Print entity statistics
    void print_entity_stats(const std::vector<int>& counts, int nprocs, const std::string& title = "") {
        if (!title.empty()) {
            std::cout << "\n" << title << ":\n";
        }

        for (int i = 0; i < nprocs; ++i) {
            std::cout << "  Shared, owned entities on proc " << i << ": "
                      << counts[4 * i] << " verts, "
                      << counts[4 * i + 1] << " edges, "
                      << counts[4 * i + 2] << " faces, "
                      << counts[4 * i + 3] << " elements\n";
        }
    }
} // namespace

int main(int argc, char** argv) {
    // Initialize MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(&argc, &argv);
    }

    int global_rank = 0, global_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    try {
        // Parse command line arguments
        std::string mesh_file = std::string(MESH_DIR) + DEFAULT_MESH_FILE;
        int num_comms = 1;

        if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
            if (global_rank == 0) print_usage(argv[0]);
            MPI_Finalize();
            return 0;
        }

        if (argc > 1) {
            mesh_file = argv[1];
        }

        if (argc > 2) {
            num_comms = std::atoi(argv[2]);
            if (num_comms < 1) {
                if (global_rank == 0) {
                    std::cerr << "Error: Number of communicators must be at least 1\n";
                    print_usage(argv[0]);
                }
                MPI_Finalize();
                return 1;
            }
        }

        // Create MOAB instance with smart pointer for automatic cleanup
        auto moab = std::make_unique<Core>();
        if (!moab) {
            std::cerr << "Error: Failed to create MOAB instance\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Split the communicator if requested
        MPI_Comm comm;
        int color = global_rank % num_comms;  // Determine color based on requested number of communicators

        if (num_comms > 1) {
            // Split the communicator into the specified number of groups
            MPI_Comm_split(MPI_COMM_WORLD, color, global_rank, &comm);
        } else {
            // Use the default communicator
            comm = MPI_COMM_WORLD;
        }

        // Get the local communicator information
        int local_rank = 0, local_size = 1;
        MPI_Comm_rank(comm, &local_rank);
        MPI_Comm_size(comm, &local_size);

        // Print process information
        if (global_rank == 0) {
            std::cout << "\n=== MOAB Parallel Hello World ===\n"
                      << "Global processes: " << global_size << "\n"
                      << "Number of communicator groups: " << num_comms << "\n"
                      << "Processes per group: ~" << (global_size + num_comms - 1) / num_comms << "\n"
                      << "Reading file: " << mesh_file << "\n"
                      << "Read options: " << DEFAULT_READ_OPTS << "\n\n";
        }

        // Read the mesh file in parallel
        if (global_rank == 0) {
            std::cout << "Loading mesh file..." << std::endl;
        }

        // Create root sets for each mesh.  Then pass these
        // to the load_file functions to be populated.
        moab::EntityHandle rootset, partnset;
        MB_CHK_SET_ERR( moab->create_meshset( moab::MESHSET_SET, rootset ), "Creating root set failed" );
        MB_CHK_SET_ERR( moab->create_meshset( moab::MESHSET_SET, partnset ), "Creating partition set failed" );

        // Create the parallel communicator object with the partition handle associated with MOAB
        auto pcomm = std::auto_ptr<moab::ParallelComm>(ParallelComm::get_pcomm( moab.get(), partnset, &comm ));

        MB_CHK_SET_ERR(
            moab->load_file(mesh_file.c_str(), &rootset, DEFAULT_READ_OPTS),
            "Failed to load mesh file: " << mesh_file
        );

        // Create parallel communication interface
        // auto pcomm = std::make_unique<ParallelComm>(moab.get(), comm);
        // auto pcomm = std::make_unique<ParallelComm>(moab.get(), comm);

        // Get shared entities and report statistics
        Range shared_sets, shared_ents, owned_entities;

        // Range src_elems;
        // MB_CHK_ERR( pcomm->get_part_entities( src_elems, 3 ) );
        // src_elems.print("Source elements: ");

        // Get sets shared with all other processors
        MB_CHK_SET_ERR(
            pcomm->get_shared_sets(shared_sets),
            "Failed to get shared sets"
        );

        shared_sets.print("Shared sets: ");

        // Get entities shared with all other processors
        for (int dim = 0; dim < 4; ++dim) {
            MB_CHK_SET_ERR(
                pcomm->get_part_entities(shared_ents, dim),
                "Failed to get shared entities"
            );
            shared_ents.print("Shared entities: ");
        }

        // MB_CHK_SET_ERR(
        //     pcomm->get_shared_entities(-1, shared_ents),
        //     "Failed to get shared entities"
        // );

        // Filter to get only locally owned shared entities
        MB_CHK_SET_ERR(
            pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities),
            "Failed to filter owned entities"
        );

        shared_ents.print("Owned entities: ");

        // Count owned entities by dimension
        unsigned int entity_counts[4] = {0};
        for (int dim = 0; dim < 4; ++dim) {
            entity_counts[dim] = static_cast<unsigned int>(owned_entities.num_of_dimension(dim));
        }

        // Gather statistics on process 0 of each communicator
        std::vector<int> recv_buffer(local_size * 4, 0);
        MPI_Gather(
            entity_counts, 4, MPI_INT,
            recv_buffer.data(), 4, MPI_INT,
            0, comm
        );

        // Print initial statistics
        if (local_rank == 0) {
            std::cout << "\n=== Initial Shared Entity Statistics (Group " << color << ") ===\n";
            print_entity_stats(std::vector<int>(entity_counts, entity_counts+4), local_size);
            print_entity_stats(recv_buffer, local_size);

        }

        // Exchange ghost elements (1 layer, using vertices as bridge)
        if (global_rank == 0) std::cout << "\nExchanging ghost elements..." << std::endl;
        MB_CHK_SET_ERR(
            pcomm->exchange_ghost_cells(
                3,      // ghost_dim: exchange 3D elements (hexahedra)
                0,      // bridge_dim: use vertices as bridge
                1,      // num_layers: exchange 1 layer of ghost elements
                0,      // addl_ents: no additional entities
                true    // store_remote_handles: store remote handles for ghost entities
            ),
            "Failed to exchange ghost cells"
        );

        // Update shared entity information after ghost exchange
        shared_ents.clear();
        owned_entities.clear();

        MB_CHK_SET_ERR(
            pcomm->get_shared_entities(-1, shared_ents),
            "Failed to get shared entities after ghost exchange"
        );

        MB_CHK_SET_ERR(
            pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities),
            "Failed to filter owned entities after ghost exchange"
        );

        // Recalculate entity counts
        for (int dim = 0; dim < 4; ++dim) {
            entity_counts[dim] = static_cast<unsigned int>(owned_entities.num_of_dimension(dim));
        }

        // Gather updated statistics
        MPI_Gather(
            entity_counts, 4, MPI_INT,
            recv_buffer.data(), 4, MPI_INT,
            0, comm
        );

        // Print final statistics
        if (local_rank == 0) {
            std::cout << "\n=== Final Shared Entity Statistics After Ghost Exchange (Group "
                      << color << ") ===\n";
            print_entity_stats(recv_buffer, local_size);
            std::cout << "\n=== Parallel Example Completed Successfully ===\n\n";
        }

        // Clean up MPI communicator if we created it
        if (num_comms > 1) {
            MPI_Comm_free(&comm);
        }

        // Clean up MOAB and MPI
        pcomm.reset();
        moab.reset();

        // Only finalize MPI if we initialized it
        if (!mpi_initialized) {
            MPI_Finalize();
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error on rank " << global_rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown error occurred on rank " << global_rank << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
