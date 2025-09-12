/**
 * @file determine_communication_pattern.cpp
 * @brief Implements functions for determining communication patterns and performing distributed reductions on sparse matrices.
 *
 * This file contains the core logic for identifying which MPI ranks need to communicate with each other
 * based on shared data (rows in a sparse matrix or Degrees of Freedom in a mesh). It also provides a function
 * to perform an efficient, in-place distributed reduction (summation) of shared data, ensuring that all
 * participating ranks have the final, consistent result while preserving non-shared local data.
 */
#include "determine_communication_pattern.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include <algorithm> // For std::find
#include <unordered_set>

#include <fstream>

/**
 * @brief Writes a log message to a rank-specific file.
 * @param rank The MPI rank of the process.
 * @param msg The message to log.
 */
void file_log(int rank, const std::string& msg) {
    std::ostringstream filename;
    filename << "rank_" << rank << ".log";
    std::ofstream log_file(filename.str(), std::ios_base::app);
    log_file << msg << std::endl;
}

#define DEBUG_LOG(msg) do { \
    std::ostringstream oss; \
    int rank; \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
    oss << msg; \
    file_log(rank, oss.str()); \
} while(0)

/**
 * @brief Clears all rank-specific log files. Should be called by rank 0.
 * @param rank The MPI rank of the calling process.
 * @param size The total number of MPI ranks.
 */
void clear_log_files(int rank, int size) {
    if (rank == 0) { // Only rank 0 clears the files
        for (int i = 0; i < size; ++i) {
            std::ostringstream filename;
            filename << "rank_" << i << ".log";
            std::ofstream log_file(filename.str(), std::ios_base::trunc);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD); // Ensure files are cleared before anyone writes
}

using namespace moab;

/**
 * @brief Determines the communication pattern based on the distribution of sparse matrix rows across MPI ranks.
 *
 * This function discovers which ranks share common global row indices in a distributed sparse matrix.
 * It constructs a symmetric communication pattern where if rank A sends data for a row to rank B, rank B
 * also sends its data for the same row to rank A. This pattern is essential for performing collective
 * operations like distributed reductions.
 *
 * The algorithm proceeds as follows:
 * 1. Each rank identifies its unique set of local row indices.
 * 2. `MPI_Allgather` is used to collect the number of rows from each rank.
 * 3. `MPI_Allgatherv` gathers all row indices from all ranks into a single buffer on every process.
 *    This gives every rank a complete picture of the row distribution.
 * 4. From this global information, each rank builds a map (`all_row_sharers`) that lists, for each
 *    global row index, the set of ranks that hold data for that row.
 * 5. Using the `all_row_sharers` map, each rank constructs its `CommunicationPattern`, populating the
 *    send and receive maps. If a rank shares a row with other ranks, it adds that row to its send/receive
 *    lists for each of those other ranks.
 * 6. The send/receive lists are sorted and made unique for consistency.
 *
 * @param pcomm A pointer to the MOAB ParallelComm object for MPI communication.
 * @param local_matrix_data The local sparse matrix data for the current rank.
 * @param comm_pattern [out] The communication pattern to be filled.
 * @param all_row_sharers [out] A map where the key is the global row index and the value is the set of ranks sharing that row.
 * @return moab::ErrorCode Returns MB_SUCCESS on success.
 */
ErrorCode determine_communication_pattern(
    ParallelComm* pcomm,
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& local_matrix_data,
    CommunicationPattern& comm_pattern,
    std::map<int, std::set<int>>& all_row_sharers)
{
    int rank = pcomm->rank();
    int size = pcomm->size();
    MPI_Comm comm = pcomm->comm();

    // Step 1: Each rank identifies its unique global row indices.
    std::set<int> local_rows_set;
    for (int k = 0; k < local_matrix_data.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(local_matrix_data, k); it; ++it) {
            local_rows_set.insert(it.row());
        }
    }
    std::vector<int> local_rows(local_rows_set.begin(), local_rows_set.end());

    // Step 2: Gather the number of rows from each rank.
    int local_row_count = local_rows.size();
    std::vector<int> all_row_counts(size);
    MPI_Allgather(&local_row_count, 1, MPI_INT, all_row_counts.data(), 1, MPI_INT, comm);

    // Step 3: Gather all row indices from all ranks (Allgatherv).
    std::vector<int> all_rows_buffer;
    std::vector<int> displacements(size);
    int total_rows_gathered = 0;
    for (int i = 0; i < size; ++i) {
        displacements[i] = total_rows_gathered;
        total_rows_gathered += all_row_counts[i];
    }
    all_rows_buffer.resize(total_rows_gathered);

    MPI_Allgatherv(local_rows.data(), local_row_count, MPI_INT,
                   all_rows_buffer.data(), all_row_counts.data(), displacements.data(),
                   MPI_INT, comm);

    // Step 4: Build the global row sharing map on every rank.
    all_row_sharers.clear();
    int current_pos = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < all_row_counts[i]; ++j) {
            int global_row = all_rows_buffer[current_pos++];
            all_row_sharers[global_row].insert(i);
        }
    }

    // Step 5: Build the symmetric communication pattern.
    comm_pattern.send_map.clear();
    comm_pattern.recv_map.clear();

    for (const auto& pair : all_row_sharers) {
        int global_row = pair.first;
        const std::set<int>& sharers = pair.second;

        if (sharers.size() > 1 && sharers.count(rank)) {
            // If this rank is a sharer, it needs to communicate with all other sharers.
            for (int other_rank : sharers) {
                if (other_rank != rank) {
                    comm_pattern.send_map[other_rank].push_back(global_row);
                    comm_pattern.recv_map[other_rank].push_back(global_row);
                }
            }
        }
    }

    // Ensure send/recv lists are sorted and unique, which is good practice.
    for (auto& pair : comm_pattern.send_map) {
        std::sort(pair.second.begin(), pair.second.end());
        pair.second.erase(std::unique(pair.second.begin(), pair.second.end()), pair.second.end());
    }
    for (auto& pair : comm_pattern.recv_map) {
        std::sort(pair.second.begin(), pair.second.end());
        pair.second.erase(std::unique(pair.second.begin(), pair.second.end()), pair.second.end());
    }

    return MB_SUCCESS;
}

/**
 * @brief Performs an in-place distributed reduction (summation) on shared rows of a sparse matrix.
 *
 * This function updates the `local_matrix_data` directly, summing the values of shared rows from
 * different MPI ranks without requiring a separate output matrix. This approach minimizes memory usage.
 * The process ensures that after the reduction, all ranks that share a particular row have the same
 * final, summed values for that row, while non-shared rows remain untouched.
 *
 * The reduction is performed using a two-stage communication strategy to avoid deadlocks and ensure correctness:
 *
 * Stage 1: Reduction to Owner
 * - For each shared row, a single "owner" rank is determined (here, the lowest-ranking sharer).
 * - All other ranks that share the row send their local data for that row to the owner.
 * - The owner rank receives the data from all other sharers and adds their contributions to its own local values.
 * - This stage uses non-blocking sends and receives (`MPI_Isend`/`MPI_Irecv`) with `MPI_Probe` to handle variable message sizes.
 *
 * Stage 2: Broadcast from Owner
 * - The owner rank, which now holds the final, fully reduced row, sends this complete row back to all the other ranks that share it.
 * - The non-owning ranks receive this final version.
 * - To update their local matrix, non-owning ranks first clear their existing row (`*= 0`) and then insert the final, correct values received from the owner.
 * - A different MPI tag is used for the broadcast to distinguish it from the reduction messages.
 *
 * @param pcomm A pointer to the MOAB ParallelComm object.
 * @param comm_pattern The pre-determined communication pattern indicating who to send to and receive from.
 * @param all_row_sharers A map detailing which ranks share each global row.
 * @param local_matrix_data [in, out] The local sparse matrix, which will be modified in-place.
 * @return moab::ErrorCode Returns MB_SUCCESS on success.
 */
ErrorCode perform_distributed_reduction(
    ParallelComm* pcomm,
    const CommunicationPattern& comm_pattern,
    const std::map<int, std::set<int>>& all_row_sharers,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& local_matrix_data)
{
    int rank = pcomm->rank();
    MPI_Comm comm = pcomm->comm();

    // Communication happens in two stages: reduction to an owner, then broadcast from the owner.

    // Determine owner for each shared row (simplistic: lowest rank is owner).
    std::map<int, int> row_owners;
    for (const auto& pair : all_row_sharers) {
        if (!pair.second.empty()) {
            row_owners[pair.first] = *pair.second.begin();
        }
    }

    // Stage 1: Reduction - Non-owners send their row data to the owner.
    std::vector<MPI_Request> reduce_requests;
    std::map<int, std::vector<double>> reduce_recv_buffers; // Owner receives into these

    for (const auto& pair : all_row_sharers) {
        int global_row = pair.first;
        const auto& sharers = pair.second;
        int owner_rank = row_owners[global_row];

        if (sharers.count(rank)) { // If this rank has the row
            if (rank == owner_rank) {
                // Owner posts receives for this row from all other sharers.
                for (int other_rank : sharers) {
                    if (other_rank != rank) {
                        // For simplicity, probe for size. In a real scenario, sizes would be exchanged first.
                        MPI_Status status;
                        MPI_Probe(other_rank, global_row, comm, &status);
                        int recv_size;
                        MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

                        reduce_recv_buffers[global_row].resize(recv_size);
                        MPI_Request req;
                        MPI_Irecv(reduce_recv_buffers[global_row].data(), recv_size, MPI_DOUBLE, other_rank, global_row, comm, &req);
                        reduce_requests.push_back(req);
                    }
                }
            } else {
                // Non-owner sends its data for this row to the owner.
                std::vector<double> send_buffer;
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(local_matrix_data, global_row); it; ++it) {
                    send_buffer.push_back(it.col());
                    send_buffer.push_back(it.value());
                }
                MPI_Request req;
                MPI_Isend(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, owner_rank, global_row, comm, &req);
                reduce_requests.push_back(req);
            }
        }
    }

    if (!reduce_requests.empty()) {
        MPI_Waitall(reduce_requests.size(), reduce_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Owners perform the reduction.
    if (reduce_recv_buffers.size() > 0) {
        for (auto const& [global_row, buffer] : reduce_recv_buffers) {
            for (size_t i = 0; i < buffer.size(); i += 2) {
                int col = buffer[i];
                double val = buffer[i+1];
                local_matrix_data.coeffRef(global_row, col) += val;
            }
        }
    }

    // Stage 2: Broadcast - Owners send the final reduced row to all non-owning sharers.
    std::vector<MPI_Request> bcast_requests;
    std::map<int, std::vector<double>> bcast_recv_buffers; // Non-owners receive into these

    for (const auto& pair : all_row_sharers) {
        int global_row = pair.first;
        const auto& sharers = pair.second;
        int owner_rank = row_owners[global_row];

        if (sharers.count(rank)) {
            if (rank == owner_rank) {
                // Owner sends the final row to other sharers.
                std::vector<double> send_buffer;
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(local_matrix_data, global_row); it; ++it) {
                    send_buffer.push_back(it.col());
                    send_buffer.push_back(it.value());
                }
                for (int other_rank : sharers) {
                    if (other_rank != rank) {
                        MPI_Request req;
                        MPI_Isend(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, other_rank, global_row + local_matrix_data.rows(), comm, &req);
                        bcast_requests.push_back(req);
                    }
                }
            } else {
                // Non-owner receives the final row from the owner.
                MPI_Status status;
                MPI_Probe(owner_rank, global_row + local_matrix_data.rows(), comm, &status);
                int recv_size;
                MPI_Get_count(&status, MPI_DOUBLE, &recv_size);
                bcast_recv_buffers[global_row].resize(recv_size);
                MPI_Request req;
                MPI_Irecv(bcast_recv_buffers[global_row].data(), recv_size, MPI_DOUBLE, owner_rank, global_row + local_matrix_data.rows(), comm, &req);
                bcast_requests.push_back(req);
            }
        }
    }

    if (!bcast_requests.empty()) {
        MPI_Waitall(bcast_requests.size(), bcast_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Non-owners update their matrix with the final, reduced row data.
    if (bcast_recv_buffers.size() > 0) {
        for (auto const& [global_row, buffer] : bcast_recv_buffers) {
            // Clear the old row before inserting new values.
            local_matrix_data.row(global_row) *= 0;
            for (size_t i = 0; i < buffer.size(); i += 2) {
                int col = buffer[i];
                double val = buffer[i+1];
                local_matrix_data.coeffRef(global_row, col) = val;
            }
        }
    }

    return MB_SUCCESS;
}

/**
 * @brief Determines a communication pattern based on shared Degrees of Freedom (DoFs) from MOAB mesh tags.
 *
 * This function inspects the mesh geometry and connectivity data stored in MOAB to figure out which
 * DoFs are shared across MPI partition boundaries. It is an alternative to `determine_communication_pattern`
 * when the sharing information is derived from the mesh topology rather than an already-existing matrix.
 *
 * The workflow is as follows:
 * 1. Identify geometric entities (e.g., edges) on the partition boundaries using `pcomm->get_shared_entities`.
 * 2. For each shared entity, retrieve the associated DoF IDs from the provided MOAB `dof_tag`.
 * 3. A preliminary, potentially non-symmetric map of shared DoFs is created (`rank_to_dofs`).
 * 4. A two-phase, non-blocking MPI exchange is performed to symmetrize this map:
 *    a. First, ranks exchange the *number* of DoFs they will send to each neighbor.
 *    b. Second, they exchange the actual DoF IDs.
 * 5. Each rank computes the intersection of the DoFs it sent and the DoFs it received from each neighbor.
 *    This intersection represents the truly shared DoFs, creating a symmetric pattern.
 * 6. The final, symmetric `CommunicationPattern` and the `all_dof_sharers` map are constructed from this intersection.
 *
 * @param pcomm A pointer to the MOAB ParallelComm object.
 * @param tag_name The name of the MOAB tag containing the DoF data (used for logging).
 * @param dof_tag The MOAB handle for the DoF tag.
 * @param comm_pattern [out] The resulting symmetric communication pattern.
 * @param all_dof_sharers [out] A map detailing which ranks share each DoF.
 * @return moab::ErrorCode Returns MB_SUCCESS on success.
 */
moab::ErrorCode determine_communication_pattern_from_tag(
    moab::ParallelComm* pcomm,
    const char* tag_name,
    moab::Tag dof_tag,
    CommunicationPattern& comm_pattern,
    std::map<int, std::set<int>>& all_dof_sharers)
{
    moab::ErrorCode rval;
    moab::Interface* mb = pcomm->get_moab();
    int rank = pcomm->rank();
    int size = pcomm->size();
    MPI_Comm comm = pcomm->comm();

    // Clear output parameters
    comm_pattern.send_map.clear();
    comm_pattern.recv_map.clear();
    all_dof_sharers.clear();

    std::cout << "Rank " << rank << ": Starting communication pattern detection for tag '" << tag_name << "'" << std::endl;

    // Get all elements
    moab::Range elems;
    rval = mb->get_entities_by_dimension(0, 2, elems); // Assuming 2D elements
    if (rval != MB_SUCCESS) return rval;
    std::cout << "Rank " << rank << ": Found " << elems.size() << " elements" << std::endl;

    // Get shared entities and their sharing ranks
    moab::Range shared_ents_edge;
    rval = pcomm->get_shared_entities(-1, shared_ents_edge, 1, true, false); // Check dimension 1 (edges)
    if (rval != MB_SUCCESS) {
        std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
        return rval;
    }
    std::cout << "Rank " << rank << ": Found " << shared_ents_edge.size() << " interface edges" << std::endl;
    moab::Range shared_ents;
    rval = mb->get_adjacencies(shared_ents_edge, 2, true, shared_ents, moab::Interface::UNION); // Check dimension 1 (edges)
    if (rval != MB_SUCCESS) {
        std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
        return rval;
    }

    // Get skin entities (boundary)
    moab::Skinner skinner(mb);
    moab::Range skin_ents;
    rval = skinner.find_skin(0, elems, false, skin_ents, nullptr, true, true);
    if (rval != MB_SUCCESS) {
        std::cerr << "Rank " << rank << ": Error finding skin" << std::endl;
        return rval;
    }
    std::cout << "Rank " << rank << ": Found " << skin_ents.size() << " skin entities" << std::endl;

    // Find shared skin entities (on partition boundaries)
    // moab::Range shared_skin = moab::intersect(skin_ents, shared_ents);
    moab::Range shared_skin = shared_ents;
    std::cout << "Rank " << rank << ": Found " << shared_skin.size()
              << " shared skin entities (partition boundaries)" << std::endl;

    // First, collect all DoFs on the interface for each adjacent rank
    std::map<int, std::set<int>> rank_to_dofs;  // Maps rank to set of DoFs on interface

    std::set< unsigned int > neighbor_procs;
    rval = pcomm->get_interface_procs(neighbor_procs);

    // Debug: Print neighbor procs and initialize rank_to_dofs
    std::cout << "Rank " << rank << ": Found " << neighbor_procs.size() << " neighbor processes" << std::endl;
    {
        std::ostringstream oss;
        oss << "Identified neighbor procs: ";
        for (int p : neighbor_procs) oss << p << " ";
        DEBUG_LOG(oss.str());
    }
    for (int np : neighbor_procs) {
        // std::cout << "Rank " << rank << ": Neighbor rank " << np << std::endl;
        // Initialize map entry for each neighbor
        rank_to_dofs[np] = std::set<int>();
    }

    // For each shared skin entity, get its adjacent elements and their DoFs
    std::cout << "Rank " << rank << ": Processing " << shared_skin.size() << " shared skin entities" << std::endl;
    for (moab::Range::iterator it = shared_skin.begin(); it != shared_skin.end(); ++it) {
        moab::EntityHandle elem = *it;

        // Get DoF numbers for this element (16 DoFs per element)
        std::vector<int> dof_numbers(16);
        rval = mb->tag_get_data(dof_tag, &elem, 1, dof_numbers.data());
        if (rval != MB_SUCCESS) {
            std::cerr << "Rank " << rank << ": Failed to get DOF data for element " << elem << std::endl;
            return rval;
        }

        // Verify DOF numbers are valid
        for (int i = 0; i < 16; ++i) {
            if (dof_numbers[i] < 0) {
                std::cerr << "Rank " << rank << ": Invalid DOF number " << dof_numbers[i]
                         << " at index " << i << " for element " << elem << std::endl;
                return MB_FAILURE;
            }
        }

        // For each element, get its DoF numbers
        moab::Range shared_edges;
        rval = mb->get_adjacencies(&elem, 1, 1, false, shared_edges, moab::Interface::UNION); // Check dimension 1 (edges)
        if (rval != MB_SUCCESS) {
            std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
            return rval;
        }
        // Get sharing ranks for this element
        std::set<int> sharing_ranks_set;
        rval = pcomm->get_sharing_data(shared_edges, sharing_ranks_set, moab::Interface::UNION);
        if (rval != MB_SUCCESS) {
            std::cerr << "Rank " << rank << ": Failed to get sharing data for element " << elem << std::endl;
            return rval;
        }

        // Add current rank to sharing set for this element
        sharing_ranks_set.insert(rank);

        // Debug: Print sharing information
        if (sharing_ranks_set.size() > 1) {  // Only print if shared with other ranks
            // std::cout << "Rank " << rank << ": Element " << elem << " is shared with ranks: ";
            // for (int r : sharing_ranks_set) std::cout << r << " ";
            // std::cout << "(DOFs: ";
            // for (int i = 0; i < 16; ++i) std::cout << dof_numbers[i] << " ";
            // std::cout << ")" << std::endl;

            // Add DOFs to rank_to_dofs for each sharing rank
            for (int other_rank : sharing_ranks_set) {
                if (other_rank != rank) {
                    for (int i = 0; i < 16; ++i) {
                        rank_to_dofs[other_rank].insert(dof_numbers[i]);
                    }
                }
            }
        }

        // Convert to vector and add current rank
        std::vector<int> sharing_ranks(sharing_ranks_set.begin(), sharing_ranks_set.end());

        // Add this processor to sharing ranks
        sharing_ranks.push_back(rank);
        std::sort(sharing_ranks.begin(), sharing_ranks.end());
        sharing_ranks.erase(std::unique(sharing_ranks.begin(), sharing_ranks.end()),
                            sharing_ranks.end());

        // For each sharing rank, add these DoFs to their interface set
        for (int other_rank : sharing_ranks) {
            if (other_rank != rank) {
                for (int j = 0; j < 16; ++j) {
                    rank_to_dofs[other_rank].insert(dof_numbers[j]);
                }
            }
        }
    }

    // Symmetrize the communication pattern in two phases to avoid deadlock.

    // Phase 1: Exchange the number of DOFs to be sent.
    std::map<int, int> incoming_sizes;
    std::vector<MPI_Request> size_requests;

    for (int neighbor : neighbor_procs) {
        if (neighbor == rank) continue;
        // Post receive for the size of the neighbor's DOF list.
        size_requests.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(&incoming_sizes[neighbor], 1, MPI_INT, neighbor, 0, comm, &size_requests.back());
    }

    for (int neighbor : neighbor_procs) {
        if (neighbor == rank) continue;
        // Post send for the size of our DOF list for that neighbor.
        int send_size = rank_to_dofs[neighbor].size();
        size_requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&send_size, 1, MPI_INT, neighbor, 0, comm, &size_requests.back());
    }

    // Wait for all size exchanges to complete.
    if (!size_requests.empty()) {
        MPI_Waitall(size_requests.size(), size_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Phase 2: Exchange the actual DOF data.
    std::map<int, std::vector<int>> received_dofs;
    std::vector<MPI_Request> data_requests;

    for (int neighbor : neighbor_procs) {
        if (neighbor == rank) continue;
        // Post receive for the actual DOF list.
        if (incoming_sizes[neighbor] > 0) {
            received_dofs[neighbor].resize(incoming_sizes[neighbor]);
            data_requests.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(received_dofs[neighbor].data(), incoming_sizes[neighbor], MPI_INT, neighbor, 1, comm, &data_requests.back());
        }
    }

    for (int neighbor : neighbor_procs) {
        if (neighbor == rank) continue;
        // Post send for our DOF list.
        if (!rank_to_dofs[neighbor].empty()) {
            std::vector<int> dofs_to_send(rank_to_dofs[neighbor].begin(), rank_to_dofs[neighbor].end());
            data_requests.push_back(MPI_REQUEST_NULL);
            MPI_Isend(dofs_to_send.data(), dofs_to_send.size(), MPI_INT, neighbor, 1, comm, &data_requests.back());
        }
    }

    // Wait for all data exchanges to complete.
    if (!data_requests.empty()) {
        MPI_Waitall(data_requests.size(), data_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Now, compute the intersection to get the symmetric communication pattern
    std::map<int, std::set<int>> final_rank_to_dofs;
    for (int neighbor : neighbor_procs) {
        if (neighbor == rank) continue;

        std::set<int> local_dofs = rank_to_dofs[neighbor];
        std::set<int> remote_dofs(received_dofs[neighbor].begin(), received_dofs[neighbor].end());

        std::set<int> intersection;
        std::set_intersection(local_dofs.begin(), local_dofs.end(),
                              remote_dofs.begin(), remote_dofs.end(),
                              std::inserter(intersection, intersection.begin()));

        final_rank_to_dofs[neighbor] = intersection;
    }
    rank_to_dofs = final_rank_to_dofs; // Replace with the symmetrized map

    // Finally, build the communication pattern and all_dof_sharers map from the
    // now-symmetric rank_to_dofs map.
    for (const auto& pair : rank_to_dofs) {
        int other_rank = pair.first;
        const std::set<int>& dofs = pair.second;

        if (!dofs.empty()) {
            // We send these DOFs to the other rank
            comm_pattern.send_map[other_rank].assign(dofs.begin(), dofs.end());

            // Symmetrically, we will receive the same DOFs from that rank
            comm_pattern.recv_map[other_rank].assign(dofs.begin(), dofs.end());
        }

        // Update the global list of sharers for each DOF
        for (int dof : dofs) {
            all_dof_sharers[dof].insert(rank);
            all_dof_sharers[dof].insert(other_rank);
        }
    }


    return MB_SUCCESS;
}