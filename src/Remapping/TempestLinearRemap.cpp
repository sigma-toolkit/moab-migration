///////////////////////////////////////////////////////////////////////////////
///
/// \file    TempestLinearRemap.cpp
/// \author  Vijay Mahadevan
/// \version Mar 08, 2017
///

#ifdef WIN32               /* windows */
#define _USE_MATH_DEFINES  // For M_PI
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include "Announce.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
#include "FiniteVolumeTools.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MathHelper.h"
#include "kdtree.h"
#include "SparseMatrix.h"
#include "OverlapMesh.h"
#include "MeshUtilitiesFuzzy.h"

#include "DebugOutput.hpp"
#include "moab/nanoflann.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/Skinner.hpp"
#include "moab/TupleList.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/MeshTopoUtil.hpp"

#pragma GCC diagnostic pop

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <unordered_set>

// #define VERBOSE
#define USE_ComputeAdjacencyRelations

/// <summary>
///     Face index and distance metric pair.
/// </summary>
typedef std::pair< int, int > FaceDistancePair;

/// <summary>
///     Vector storing adjacent Faces.
/// </summary>
typedef std::vector< FaceDistancePair > AdjacentFaceVector;

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::LinearRemapNN_MOAB( bool use_GID_matching, bool strict_check )
{
    /* m_mapRemap size = (m_nTotDofs_Dest X m_nTotDofs_SrcCov)  */

#ifdef VVERBOSE
    {
        std::ofstream output_file( "rowcolindices.txt", std::ios::out );
        output_file << m_nTotDofs_Dest << " " << m_nTotDofs_SrcCov << " " << row_gdofmap.size() << " "
                    << row_ldofmap.size() << " " << col_gdofmap.size() << " " << col_ldofmap.size() << "\n";
        output_file << "Rows \n";
        for( unsigned iv = 0; iv < row_gdofmap.size(); iv++ )
            output_file << row_gdofmap[iv] << " " << row_dofmap[iv] << "\n";
        output_file << "Cols \n";
        for( unsigned iv = 0; iv < col_gdofmap.size(); iv++ )
            output_file << col_gdofmap[iv] << " " << col_dofmap[iv] << "\n";
        output_file.flush();  // required here
        output_file.close();
    }
#endif

    if( use_GID_matching )
    {
        std::map< unsigned, unsigned > src_gl;
        for( unsigned it = 0; it < col_gdofmap.size(); ++it )
            src_gl[col_gdofmap[it]] = it;

        std::map< unsigned, unsigned >::iterator iter;
        for( unsigned it = 0; it < row_gdofmap.size(); ++it )
        {
            unsigned row = row_gdofmap[it];
            iter         = src_gl.find( row );
            if( strict_check && iter == src_gl.end() )
            {
                std::cout << "Searching for global target DOF " << row
                          << " but could not find correspondence in source mesh.\n";
                assert( false );
            }
            else if( iter == src_gl.end() )
            {
                continue;
            }
            else
            {
                unsigned icol = src_gl[row];
                unsigned irow = it;

                // Set the permutation matrix in local space
                m_mapRemap( irow, icol ) = 1.0;
            }
        }

        return moab::MB_SUCCESS;
    }
    else
    {
        /* Create a Kd-tree to perform local queries to find nearest neighbors */

        return moab::MB_FAILURE;
    }
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapFVtoFV_Tempest_MOAB( int nOrder )
{
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if( m_meshInputCov->faces.size() > 0 && m_meshInputCov->revnodearray.size() == 0 )
    {
        _EXCEPTIONT( "ReverseNodeArray has not been calculated for m_meshInputCov" );
    }

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( TriQuadRuleOrder );

    // Number of coefficients needed at this order
#ifdef RECTANGULAR_TRUNCATION
    int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION
    int nCoefficients = nOrder * ( nOrder + 1 ) / 2;
#endif

    // Number of faces you need
    const int nRequiredFaceSetSize = nCoefficients;

    // Fit weight exponent
    const int nFitWeightsExponent = nOrder + 2;

    // Announcements
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapFVtoFV_Tempest_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Volume to Finite Volume Projection\n" );
        dbgprint.printf( 0, "Triangular quadrature rule order %i\n", TriQuadRuleOrder );
        dbgprint.printf( 0, "Number of coefficients: %i\n", nCoefficients );
        dbgprint.printf( 0, "Required adjacency set size: %i\n", nRequiredFaceSetSize );
        dbgprint.printf( 0, "Fit weights exponent: %i\n", nFitWeightsExponent );
    }

    // Current overlap face
    int ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    DataArray2D< double > dIntArray;
    DataArray1D< double > dConstraint( nCoefficients );

    // Loop through all faces on m_meshInputCov
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        // Output every 1000 elements
#ifdef VERBOSE
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Find the set of Faces that overlap faceFirst
        int ixOverlapBegin    = ixOverlap;
        unsigned ixOverlapEnd = ixOverlapBegin;

        for( ; ixOverlapEnd < m_meshOverlap->faces.size(); ixOverlapEnd++ )
        {
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapEnd] != 0 ) break;
        }

        unsigned nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

        if( nOverlapFaces == 0 ) continue;

        // Build integration array
        BuildIntegrationArray( *m_meshInputCov, *m_meshOverlap, triquadrule, ixFirst, ixOverlapBegin, ixOverlapEnd,
                               nOrder, dIntArray );

        // Set of Faces to use in building the reconstruction and associated
        // distance metric.
        AdjacentFaceVector vecAdjFaces;

        GetAdjacentFaceVectorByEdge( *m_meshInputCov, ixFirst, nRequiredFaceSetSize, vecAdjFaces );

        // Number of adjacent Faces
        int nAdjFaces = vecAdjFaces.size();

        // Determine the conservative constraint equation
        double dFirstArea = m_meshInputCov->vecFaceArea[ixFirst];
        dConstraint.Zero();
        for( int p = 0; p < nCoefficients; p++ )
        {
            for( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                dConstraint[p] += dIntArray[p][j];
            }
            dConstraint[p] /= dFirstArea;
        }

        // Build the fit array from the integration operator
        DataArray2D< double > dFitArray;
        DataArray1D< double > dFitWeights;
        DataArray2D< double > dFitArrayPlus;

        BuildFitArray( *m_meshInputCov, triquadrule, ixFirst, vecAdjFaces, nOrder, nFitWeightsExponent, dConstraint,
                       dFitArray, dFitWeights );

        // Compute the inverse fit array
        bool fSuccess = InvertFitArray_Corrected( dConstraint, dFitArray, dFitWeights, dFitArrayPlus );

        // Multiply integration array and fit array
        DataArray2D< double > dComposedArray( nAdjFaces, nOverlapFaces );
        if( fSuccess )
        {
            // Multiply integration array and inverse fit array
            for( int i = 0; i < nAdjFaces; i++ )
            {
                for( size_t j = 0; j < nOverlapFaces; j++ )
                {
                    for( int k = 0; k < nCoefficients; k++ )
                    {
                        dComposedArray( i, j ) += dIntArray( k, j ) * dFitArrayPlus( i, k );
                    }
                }
            }

            // Unable to invert fit array, drop to 1st order.  In this case
            // dFitArrayPlus(0,0) = 1 and all other entries are zero.
        }
        else
        {
            dComposedArray.Zero();
            for( size_t j = 0; j < nOverlapFaces; j++ )
            {
                dComposedArray( 0, j ) += dIntArray( 0, j );
            }
        }

        // Put composed array into map
        for( unsigned i = 0; i < vecAdjFaces.size(); i++ )
        {
            for( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                int& ixFirstFaceLoc  = vecAdjFaces[i].first;
                int& ixSecondFaceLoc = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];
                // int ixFirstFaceGlob = m_remapper->GetGlobalID(moab::Remapper::SourceMesh,
                // ixFirstFaceLoc); int ixSecondFaceGlob =
                // m_remapper->GetGlobalID(moab::Remapper::TargetMesh, ixSecondFaceLoc);

                // signal to not participate, because it is a ghost target
                if( ixSecondFaceLoc < 0 ) continue;  // do not do anything

                m_mapRemap( ixSecondFaceLoc, ixFirstFaceLoc ) +=
                    dComposedArray[i][j] / m_meshOutput->vecFaceArea[ixSecondFaceLoc];
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::PrintMapStatistics()
{
    int nrows = m_weightMatrix.rows();      // Number of rows
    int ncols = m_weightMatrix.cols();      // Number of columns
    int NNZ   = m_weightMatrix.nonZeros();  // Number of non zero values
#ifdef MOAB_HAVE_MPI
    // find out min/max for NNZ, ncols, nrows
    // should work on std c++ 11
    int arr3[6] = { NNZ, nrows, ncols, -NNZ, -nrows, -ncols };
    int rarr3[6];
    MPI_Reduce( arr3, rarr3, 6, MPI_INT, MPI_MIN, 0, m_pcomm->comm() );

    int total[3];
    MPI_Reduce( arr3, total, 3, MPI_INT, MPI_SUM, 0, m_pcomm->comm() );
    if( !rank )
        std::cout << "-> Rows (min/max/sum): (" << rarr3[1] << " / " << -rarr3[4] << " / " << total[1] << "), "
                  << " Cols (min/max/sum): (" << rarr3[2] << " / " << -rarr3[5] << " / " << total[2] << "), "
                  << " NNZ (min/max/sum): (" << rarr3[0] << " / " << -rarr3[3] << " / " << total[0] << ")\n";
#else
    std::cout << "-> Rows: " << nrows << ", Cols: " << ncols << ", NNZ: " << NNZ << "\n";
#endif
}

///////////////////////////////////////////////////////////////////////////////

#ifdef MOAB_HAVE_EIGEN3

struct CommunicationPattern {
    std::map<int, std::vector<int>> send_map; // rank -> [global_rows]
    std::map<int, std::vector<int>> recv_map; // rank -> [global_rows]
};

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
static moab::ErrorCode determine_communication_pattern(
    moab::ParallelComm* pcomm,
    const std::vector<Eigen::Triplet<double>>& local_triplets,
    const std::vector<mbGIDType>& row_gdofmap,
    CommunicationPattern& comm_pattern,
    std::map<mbGIDType, std::set<int>>& all_row_sharers)
{
    int rank = pcomm->rank();
    int size = pcomm->size();
    MPI_Comm comm = pcomm->comm();

    // Step 1: Each rank identifies its unique global row DOF IDs from its triplets.
    std::set<int> global_rows_set;
    for (const auto& triplet : local_triplets) {
        // Convert local row index to global DOF ID using row_gdofmap
        int global_row_id = row_gdofmap[triplet.row()];
        global_rows_set.insert(global_row_id);
    }
    std::vector<int> local_rows(global_rows_set.begin(), global_rows_set.end());

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

    return moab::MB_SUCCESS;
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
 * @param local_triplets [in, out] The local sparse matrix, which will be modified in-place.
 * @return moab::ErrorCode Returns MB_SUCCESS on success.
 */
static moab::ErrorCode perform_distributed_reduction(
    moab::ParallelComm* pcomm,
    const std::map<mbGIDType, mbGIDType>& rowMap,
    const std::map<mbGIDType, mbGIDType>& colMap,
    const std::vector<mbGIDType>& row_gdofmap,
    const std::vector<mbGIDType>& col_gdofmap,
    const std::map<mbGIDType, std::set<int>>& all_row_sharers,
    std::vector<Eigen::Triplet<double>>& local_triplets)
{
    int rank = pcomm->rank();
    MPI_Comm comm = pcomm->comm();

    if (rank == 0) {
        std::cout << "TempestLinearRemap:: Performing distributed reduction on shared rows." << std::endl;
    }

    // Communication happens in two stages: reduction to an owner, then broadcast from the owner.

    // Determine owner for each shared row (simplistic: lowest rank is owner).
    std::map<int, int> row_owners;
    for (const auto& pair : all_row_sharers) {
        if (!pair.second.empty()) {
            row_owners[pair.first] = *pair.second.begin();
        }
    }

    return moab::MB_SUCCESS;

    // Calculate a globally consistent tag offset for the broadcast stage.
    // This prevents tag collisions between the reduction and broadcast stages.
    int local_max_row = 0;
    for (const auto& triplet : local_triplets) {
        int global_row = rowMap.at(triplet.row());
        if (global_row > local_max_row) {
            local_max_row = global_row;
        }
    }

    int global_max_row = 0;
    MPI_Allreduce(&local_max_row, &global_max_row, 1, MPI_INT, MPI_MAX, comm);
    int tag_offset = global_max_row + 1;

    // Stage 1: Reduction - Non-owners send their row data to the owner rank.
    std::vector<MPI_Request> reduce_requests;
    std::map<int, std::vector<char>> reduce_recv_buffers; // Owners receive into these

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
                        int recv_bytes;
                        MPI_Get_count(&status, MPI_BYTE, &recv_bytes);

                        reduce_recv_buffers[global_row].resize(recv_bytes);
                        MPI_Request req;
                        MPI_Irecv(reduce_recv_buffers[global_row].data(), recv_bytes, MPI_BYTE, other_rank, global_row, comm, &req);
                        reduce_requests.push_back(req);
                    }
                }
            } else {
                // Non-owner sends its data for this row to the owner.
                std::vector<Eigen::Triplet<double>> triplets_to_send;
                for (const auto& triplet : local_triplets) {
                    // Check if this triplet's global row DOF ID matches the shared row
                    if (row_gdofmap[triplet.row()] == global_row) {
                        triplets_to_send.emplace_back(row_gdofmap[triplet.row()], col_gdofmap[triplet.col()], triplet.value());
                    }
                }
                if (!triplets_to_send.empty()) {
                    // Convert to global DOF IDs before sending
                    std::vector<Eigen::Triplet<double>> triplets_to_send_global;
                    triplets_to_send_global.reserve(triplets_to_send.size());
                    for (const auto& triplet : triplets_to_send) {
                        triplets_to_send_global.emplace_back(row_gdofmap[triplet.row()], col_gdofmap[triplet.col()], triplet.value());
                    }

                    MPI_Request req;
                    MPI_Isend(triplets_to_send_global.data(), triplets_to_send_global.size() * sizeof(Eigen::Triplet<double>),
                              MPI_BYTE, owner_rank, global_row, comm, &req);
                    reduce_requests.push_back(req);
                }
            }
        }
    }

    if (!reduce_requests.empty()) {
        MPI_Waitall(reduce_requests.size(), reduce_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Owners perform the reduction.
    if (reduce_recv_buffers.size() > 0) {
        for (const auto& pair : reduce_recv_buffers) {
            // int global_row = pair.first;
            const std::vector<char>& recv_buffer = pair.second;
            const Eigen::Triplet<double>* received_triplets = reinterpret_cast<const Eigen::Triplet<double>*>(recv_buffer.data());
            int num_triplets = recv_buffer.size() / sizeof(Eigen::Triplet<double>);

            for (int i = 0; i < num_triplets; ++i) {
                const Eigen::Triplet<double>& global_triplet = received_triplets[i];
                auto row_it = rowMap.find(global_triplet.row());
                auto col_it = colMap.find(global_triplet.col());

                if (row_it != rowMap.end() && col_it != colMap.end()) {
                    int lRow = row_it->second;
                    int lCol = col_it->second;
                    double value = global_triplet.value();

                    // Find existing triplet and add to it, or create new one
                    bool found = false;
                    for (auto& local_triplet : local_triplets) {
                        if (local_triplet.row() == lRow && local_triplet.col() == lCol) {
                            const_cast<double&>(local_triplet.value()) += value;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        local_triplets.emplace_back(lRow, lCol, value);
                    }
                } else {
                    std::cout << "Warning: Rank " << rank << " received global DOF IDs ("
                              << global_triplet.row() << ", " << global_triplet.col()
                              << ") that are not found in local maps. Skipping." << std::endl;
                }
            }
        }
    }

    // Stage 2: Broadcast - Owners send the final reduced row to all non-owning sharers.
    std::vector<MPI_Request> bcast_requests;
    std::map<int, std::vector<Eigen::Triplet<double>>> bcast_recv_buffers; // Non-owners receive into these

    for (const auto& pair : all_row_sharers) {
        int global_row = pair.first;
        const auto& sharers = pair.second;
        int owner_rank = row_owners[global_row];

        if (sharers.count(rank)) {
            if (rank == owner_rank) {
                // Owner sends the final row to other sharers.
                std::vector<Eigen::Triplet<double>> send_buffer;
                for (const auto& triplet : local_triplets) {
                    // Check if this triplet's global row DOF ID matches the shared row
                    if (row_gdofmap[triplet.row()] == global_row) {
                        send_buffer.push_back(triplet);
                    }
                }
                // Convert to global DOF IDs before broadcasting
                std::vector<Eigen::Triplet<double>> triplets_for_bcast_global;
                triplets_for_bcast_global.reserve(send_buffer.size());
                for (const auto& triplet : send_buffer) {
                    triplets_for_bcast_global.emplace_back(row_gdofmap[triplet.row()], col_gdofmap[triplet.col()], triplet.value());
                }
                for (int other_rank : sharers) {
                    if (other_rank != rank) {
                        MPI_Request req;
                        MPI_Isend(triplets_for_bcast_global.data(), triplets_for_bcast_global.size() * sizeof(Eigen::Triplet<double>),
                                  MPI_BYTE, other_rank, global_row + tag_offset, comm, &req);
                        bcast_requests.push_back(req);
                    }
                }
            } else {
                // Non-owner receives the final row from the owner.
                MPI_Status status;
                MPI_Probe(owner_rank, global_row + tag_offset, comm, &status);
                int recv_bytes;
                MPI_Get_count(&status, MPI_BYTE, &recv_bytes);
                int recv_count = recv_bytes / sizeof(Eigen::Triplet<double>);
                bcast_recv_buffers[global_row].resize(recv_count);
                MPI_Request req;
                MPI_Irecv(bcast_recv_buffers[global_row].data(), recv_bytes, MPI_BYTE, owner_rank, global_row + tag_offset, comm, &req);
                bcast_requests.push_back(req);
            }
        }
    }

    if (!bcast_requests.empty()) {
        MPI_Waitall(bcast_requests.size(), bcast_requests.data(), MPI_STATUSES_IGNORE);
    }

    // Non-owners update their matrix with the final, reduced row data.
    if (!bcast_recv_buffers.empty()) {
        // First, create a set of all rows that this rank will receive from owners.
        std::set<int> shared_rows_to_update;
        for (const auto& pair : bcast_recv_buffers) {
            shared_rows_to_update.insert(pair.first);
        }

        // Remove all existing triplets for these shared rows (using global DOF IDs).
        local_triplets.erase(
            std::remove_if(local_triplets.begin(), local_triplets.end(),
                [&](const Eigen::Triplet<double>& t) {
                    int global_row_id = row_gdofmap[t.row()];
                    return shared_rows_to_update.count(global_row_id);
                }),
            local_triplets.end());

        // Now, add the final, authoritative triplets received from the owners.
        for (const auto& pair : bcast_recv_buffers) {
            const auto& received_triplets = pair.second;
            for (const auto& global_triplet : received_triplets) {
                auto row_it = rowMap.find(global_triplet.row());
                auto col_it = colMap.find(global_triplet.col());

                if (row_it != rowMap.end() && col_it != colMap.end()) {
                    local_triplets.emplace_back(row_it->second, col_it->second, global_triplet.value());
                } else {
                    std::cout << "Warning: Rank " << rank << " received global DOF IDs ("
                              << global_triplet.row() << ", " << global_triplet.col()
                              << ") that are not found in local maps. Skipping." << std::endl;
                }
            }
        }
    }

    return moab::MB_SUCCESS;
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
    std::map<mbGIDType, std::set<int>>& all_dof_sharers)
{
    moab::ErrorCode rval;
    moab::Interface* mb = pcomm->get_moab();
    int rank = pcomm->rank();
    //int size = pcomm->size();
    MPI_Comm comm = pcomm->comm();

    // Clear output parameters
    comm_pattern.send_map.clear();
    comm_pattern.recv_map.clear();
    all_dof_sharers.clear();

    std::cout << "Rank " << rank << ": Starting communication pattern detection for tag '" << tag_name << "'" << std::endl;

    // Get all elements
    moab::Range elems;
    rval = mb->get_entities_by_dimension(0, 2, elems); // Assuming 2D elements
    if (rval != moab::MB_SUCCESS) return rval;
    std::cout << "Rank " << rank << ": Found " << elems.size() << " elements" << std::endl;

    // Get shared entities and their sharing ranks
    moab::Range shared_ents_edge;
    rval = pcomm->get_shared_entities(-1, shared_ents_edge, 1, true, false); // Check dimension 1 (edges)
    if (rval != moab::MB_SUCCESS) {
        std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
        return rval;
    }
    std::cout << "Rank " << rank << ": Found " << shared_ents_edge.size() << " interface edges" << std::endl;
    moab::Range shared_ents;
    rval = mb->get_adjacencies(shared_ents_edge, 2, true, shared_ents, moab::Interface::UNION); // Check dimension 1 (edges)
    if (rval != moab::MB_SUCCESS) {
        std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
        return rval;
    }

    // Get skin entities (boundary)
    moab::Skinner skinner(mb);
    moab::Range skin_ents;
    rval = skinner.find_skin(0, elems, false, skin_ents, nullptr, true, true);
    if (rval != moab::MB_SUCCESS) {
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
        std::cout << oss.str() << std::endl;
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
        std::vector<mbGIDType> dof_numbers(16);
        rval = mb->tag_get_data(dof_tag, &elem, 1, dof_numbers.data());
        if (rval != moab::MB_SUCCESS) {
            std::cerr << "Rank " << rank << ": Failed to get DOF data for element " << elem << std::endl;
            return rval;
        }

        // Verify DOF numbers are valid
        for (int i = 0; i < 16; ++i) {
            if (dof_numbers[i] < 0) {
                std::cerr << "Rank " << rank << ": Invalid DOF number " << dof_numbers[i]
                         << " at index " << i << " for element " << elem << std::endl;
                return moab::MB_FAILURE;
            }
        }

        // For each element, get its DoF numbers
        moab::Range shared_edges;
        rval = mb->get_adjacencies(&elem, 1, 1, false, shared_edges, moab::Interface::UNION); // Check dimension 1 (edges)
        if (rval != moab::MB_SUCCESS) {
            std::cerr << "Rank " << rank << ": Error getting shared entities" << std::endl;
            return rval;
        }
        // Get sharing ranks for this element
        std::set<int> sharing_ranks_set;
        rval = pcomm->get_sharing_data(shared_edges, sharing_ranks_set, moab::Interface::UNION);
        if (rval != moab::MB_SUCCESS) {
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

    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::copy_tempest_sparsemat_to_eigen3(bool perform_reduction)
{
    // std::cout << "copy_tempest_sparsemat_to_eigen3\n";
    // std::cout << "m_weightMatrix.rows() = " << m_weightMatrix.rows() << "\n";
    // std::cout << "m_weightMatrix.cols() = " << m_weightMatrix.cols() << "\n";
#ifndef VERBOSE
#define VERBOSE_ACTIVATED
// #define VERBOSE
#endif

    /* Should the columns be the global size of the matrix ? */
    m_weightMatrix.resize( m_nTotDofs_Dest, m_nTotDofs_SrcCov );
    m_rowVector.resize( m_weightMatrix.rows() );
    m_colVector.resize( m_weightMatrix.cols() );

#ifdef VERBOSE
    int locrows = std::max( m_mapRemap.GetRows(), m_nTotDofs_Dest );
    int loccols = std::max( m_mapRemap.GetColumns(), m_nTotDofs_SrcCov );

    std::cout << m_weightMatrix.rows() << ", " << locrows << ", " << m_weightMatrix.cols() << ", " << loccols << "\n";
    // assert(m_weightMatrix.rows() == locrows && m_weightMatrix.cols() == loccols);
#endif

    // Extract triplets from m_mapRemap and populate local_triplets
    DataArray1D< int > lrows;
    DataArray1D< int > lcols;
    DataArray1D< double > lvals;
    m_mapRemap.GetEntries( lrows, lcols, lvals );
    size_t locvals = lvals.GetRows();

    // Clear and populate local_triplets
    this->local_triplets.clear();
    this->local_triplets.reserve( locvals );
    for( size_t iv = 0; iv < locvals; iv++ )
    {
        this->local_triplets.emplace_back( lrows[iv], lcols[iv], lvals[iv] );
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (perform_reduction && false) {
        // --- Determine the communication pattern based on row distribution ---
        CommunicationPattern comm_pattern;
        std::map<mbGIDType, std::set<int>> all_row_sharers;
        constexpr int method = 1;

        if (method == 0)
        {
            std::cout << rank << ": Determining communication pattern. " << std::endl;
            MB_CHK_SET_ERR( determine_communication_pattern(m_pcomm, this->local_triplets, this->row_gdofmap, comm_pattern, all_row_sharers),
                                "determine_communication_pattern failed");

            // --- Perform the distributed reduction in-place on the triplet list ---
            if (rank == 0) std::cout << "Performing distributed reduction on triplets." << std::endl;
            MB_CHK_SET_ERR( perform_distributed_reduction(m_pcomm, this->rowMap, this->colMap, this->row_gdofmap, this->col_gdofmap, all_row_sharers, this->local_triplets),
                                "perform_distributed_reduction failed" );
            if (rank == 0) std::cout << rank << ": Performing distributed reduction all done. Operator is synced!!" << std::endl;
        }
        else
        {
            MB_CHK_SET_ERR( determine_communication_pattern_from_tag(m_pcomm, "GLOBAL_DOFS", m_dofTagDest, comm_pattern, all_row_sharers),
                                "determine_communication_pattern_from_tag failed");

            // For now, just print some debug information about the triplets
            std::cout << rank << ": Found " << this->local_triplets.size() << " local triplets for potential reduction" << std::endl;
        }
    }

    m_weightMatrix.setFromTriplets( this->local_triplets.begin(), this->local_triplets.end() );
    m_weightMatrix.makeCompressed();

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "tempestmatrix.txt.0000" << rank;
    std::ofstream output_file( sstr.str(), std::ios::out );
    output_file << "0 " << locrows << " 0 " << loccols << "\n";
    for( const auto& triplet : this->local_triplets )
    {
        output_file << row_gdofmap[triplet.row()] << " " << col_gdofmap[triplet.col()] << " "
                    << triplet.value() << "\n";
    }
    output_file.flush();  // required here
    output_file.close();
#endif

#ifdef VERBOSE_ACTIVATED
#undef VERBOSE_ACTIVATED
#undef VERBOSE
#endif
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

template < typename T >
static std::vector< size_t > sort_indexes( const std::vector< T >& v )
{
    // initialize original index locations
    std::vector< size_t > idx( v.size() );
    std::iota( idx.begin(), idx.end(), 0 );

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort( idx.begin(), idx.end(), [&v]( size_t i1, size_t i2 ) { return fabs( v[i1] ) > fabs( v[i2] ); } );

    return idx;
}

double moab::TempestOnlineMap::QLTLimiter( int caasIteration,
                                           std::vector< double >& dataCorrectedField,
                                           std::vector< double >& dataLowerBound,
                                           std::vector< double >& dataUpperBound,
                                           std::vector< double >& dMassDefect )
{
    const size_t nrows = dataCorrectedField.size();
    double dMassL      = 0.0;
    double dMassU      = 0.0;
    std::vector< double > dataCorrection( nrows );
    double dMassDiffCum                       = 0.0;
    double dLMinusU                           = fabs( dataUpperBound[0] - dataLowerBound[0] );
    const DataArray1D< double >& dTargetAreas = this->m_remapper->m_target->vecFaceArea;

    // std::vector< size_t > sortedIdx = sort_indexes( dMassDefect );
    std::vector< std::unordered_set< int > > vecAdjTargetFaces( nrows );
    constexpr bool useMOABAdjacencies = true;
#ifdef USE_ComputeAdjacencyRelations
    if( useMOABAdjacencies )
        ComputeAdjacencyRelations( vecAdjTargetFaces, caasIteration, m_remapper->m_target_entities,
                                   useMOABAdjacencies );
    else
        ComputeAdjacencyRelations( vecAdjTargetFaces, caasIteration, m_remapper->m_target_entities, useMOABAdjacencies,
                                   this->m_remapper->m_target );
#else
    moab::MeshTopoUtil mtu( m_interface );
    ;
#endif

    for( size_t i = 0; i < nrows; i++ )
    {
        // size_t index = sortedIdx[i];
        size_t index          = i;
        dataCorrection[index] = fmax( dataLowerBound[index], fmin( dataUpperBound[index], 0.0 ) );
        // dMassDiff[index] = dMassDefect[index] - dTargetAreas[index] * dataCorrection[index];
        // dMassDiff[index] = dMassDefect[index];

        dMassL += dTargetAreas[index] * dataLowerBound[index];
        dMassU += dTargetAreas[index] * dataUpperBound[index];
        dLMinusU = fmax( dLMinusU, fabs( dataUpperBound[index] - dataLowerBound[index] ) );
        dMassDiffCum += dMassDefect[index] - dTargetAreas[index] * dataCorrection[index];

#ifndef USE_ComputeAdjacencyRelations
        vecAdjTargetFaces[index].insert( index );  // add self target face first
        {
            // Compute the adjacent faces to the target face
            if( useMOABAdjacencies )
            {
                moab::Range ents;
                // ents.insert( m_remapper->m_target_entities.index( m_remapper->m_target_entities[index] ) );
                ents.insert( m_remapper->m_target_entities[index] );
                moab::Range adjEnts;
                moab::ErrorCode rval = mtu.get_bridge_adjacencies( ents, 0, 2, adjEnts, caasIteration );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
                for( moab::Range::iterator it = adjEnts.begin(); it != adjEnts.end(); ++it )
                {
                    // int adjIndex = m_interface->id_from_handle(*it)-1;
                    int adjIndex = m_remapper->m_target_entities.index( *it );
                    // printf("rank: %d, Element %lu, entity: %lu, adjIndex %d\n", rank, index, *it, adjIndex);
                    if( adjIndex >= 0 ) vecAdjTargetFaces[index].insert( adjIndex );
                }
            }
            else
            {
                AdjacentFaceVector vecAdjFaces;
                GetAdjacentFaceVectorByEdge( *this->m_remapper->m_target, index,
                                             ( m_output_order + 1 ) * ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             //  ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             //  ( 4 ) * ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             vecAdjFaces );

                // Add the adjacent faces to the target face list
                for( auto adjFace : vecAdjFaces )
                    if( adjFace.first >= 0 )
                        vecAdjTargetFaces[index].insert( adjFace.first );  // map target face to source face
            }
        }
#endif
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localDefects( 5, 0.0 ), globalDefects( 5, 0.0 );
    localDefects[0] = dMassL;
    localDefects[1] = dMassU;
    localDefects[2] = dMassDiffCum;
    localDefects[3] = dLMinusU;
    // localDefects[4] = dMassCorrectU;

    MPI_Allreduce( localDefects.data(), globalDefects.data(), 4, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );

    dMassL       = globalDefects[0];
    dMassU       = globalDefects[1];
    dMassDiffCum = globalDefects[2];
    dLMinusU     = globalDefects[3];
    // dMassCorrectU = globalDefects[4];
#endif

    //If the upper and lower bounds are too close together, just clip
    if( fabs( dMassDiffCum ) < 1e-15 || dLMinusU < 1e-15 )
    {
        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
        return dMassDiffCum;
    }
    else
    {
        if( dMassL > dMassDiffCum )
        {
            Announce( "Lower bound mass exceeds target mass by %1.15e: CAAS will need another iteration",
                      dMassL - dMassDiffCum );
            dMassDiffCum = dMassL;
            // dMass -= dMassL;
        }
        else if( dMassU < dMassDiffCum )
        {
            Announce( "Target mass exceeds upper bound mass by %1.15e: CAAS will need another iteration",
                      dMassDiffCum - dMassU );
            dMassDiffCum = dMassU;
            // dMass -= dMassU;
        }

        // TODO: optimize away dataMassVec by a simple transient double within the loop
        // DataArray1D< double > dataMassVec( nrows );  //vector of mass redistribution
        for( size_t i = 0; i < nrows; i++ )
        {
            // size_t index   = sortedIdx[i];
            size_t index                               = i;
            const std::unordered_set< int >& neighbors = vecAdjTargetFaces[index];
            if( dMassDefect[index] > 0.0 )
            {
                double dMassCorrectU = 0.0;
                for( auto it : neighbors )
                    dMassCorrectU += dTargetAreas[it] * ( dataUpperBound[it] - dataCorrection[it] );

                // double dMassDiffCumOld = dMassDefect[index];
                for( auto it : neighbors )
                    dataCorrection[it] +=
                        dMassDefect[index] * ( dataUpperBound[it] - dataCorrection[it] ) / dMassCorrectU;
            }
            else
            {
                double dMassCorrectL = 0.0;
                for( auto it : neighbors )
                    dMassCorrectL += dTargetAreas[it] * ( dataCorrection[it] - dataLowerBound[it] );

                // double dMassDiffCumOld = dMassDefect[index];
                for( auto it : neighbors )
                    dataCorrection[it] +=
                        dMassDefect[index] * ( dataCorrection[it] - dataLowerBound[it] ) / dMassCorrectL;
            }
        }

        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
    }

    return dMassDiffCum;
}

void moab::TempestOnlineMap::CAASLimiter( std::vector< double >& dataCorrectedField,
                                          std::vector< double >& dataLowerBound,
                                          std::vector< double >& dataUpperBound,
                                          double& dMass )
{
    const size_t nrows = dataCorrectedField.size();
    double dMassL      = 0.0;
    double dMassU      = 0.0;
    std::vector< double > dataCorrection( nrows );
    const DataArray1D< double >& dTargetAreas = this->m_remapper->m_target->vecFaceArea;
    double dMassDiff                          = dMass;
    double dLMinusU                           = fabs( dataUpperBound[0] - dataLowerBound[0] );
    double dMassCorrectU                      = 0.0;
    double dMassCorrectL                      = 0.0;
    for( size_t i = 0; i < nrows; i++ )
    {
        dataCorrection[i] = fmax( dataLowerBound[i], fmin( dataUpperBound[i], 0.0 ) );
        dMassL += dTargetAreas[i] * dataLowerBound[i];
        dMassU += dTargetAreas[i] * dataUpperBound[i];
        dMassDiff -= dTargetAreas[i] * dataCorrection[i];
        dLMinusU = fmax( dLMinusU, fabs( dataUpperBound[i] - dataLowerBound[i] ) );
        dMassCorrectL += dTargetAreas[i] * ( dataCorrection[i] - dataLowerBound[i] );
        dMassCorrectU += dTargetAreas[i] * ( dataUpperBound[i] - dataCorrection[i] );
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localDefects( 5, 0.0 ), globalDefects( 5, 0.0 );
    localDefects[0] = dMassL;
    localDefects[1] = dMassU;
    localDefects[2] = dMassDiff;
    localDefects[3] = dMassCorrectL;
    localDefects[4] = dMassCorrectU;

    MPI_Allreduce( localDefects.data(), globalDefects.data(), 5, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );

    dMassL        = globalDefects[0];
    dMassU        = globalDefects[1];
    dMassDiff     = globalDefects[2];
    dMassCorrectL = globalDefects[3];
    dMassCorrectU = globalDefects[4];
#endif

    //If the upper and lower bounds are too close together, just clip
    if( fabs( dMassDiff ) < 1e-15 || fabs( dLMinusU ) < 1e-15 )
    {
        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
        return;
    }
    else
    {
        if( dMassL > dMassDiff )
        {
            Announce( "%d: Lower bound mass exceeds target mass by %1.15e: CAAS will need another iteration", rank,
                      dMassL - dMassDiff );
            dMassDiff = dMassL;
            dMass -= dMassL;
        }
        else if( dMassU < dMassDiff )
        {
            Announce( "%d: Target mass exceeds upper bound mass by %1.15e: CAAS will need another iteration", rank,
                      dMassDiff - dMassU );
            dMassDiff = dMassU;
            dMass -= dMassU;
        }

        // TODO: optimize away dataMassVec by a simple transient double within the loop
        DataArray1D< double > dataMassVec( nrows );  //vector of mass redistribution
        if( dMassDiff > 0.0 )
        {
            for( size_t i = 0; i < nrows; i++ )
            {
                dataMassVec[i] = ( dataUpperBound[i] - dataCorrection[i] ) / dMassCorrectU;
                dataCorrection[i] += dMassDiff * dataMassVec[i];
            }
        }
        else
        {
            for( size_t i = 0; i < nrows; i++ )
            {
                dataMassVec[i] = ( dataCorrection[i] - dataLowerBound[i] ) / dMassCorrectL;
                dataCorrection[i] += dMassDiff * dataMassVec[i];
            }
        }

        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
    }

    return;
}

std::pair< double, double > moab::TempestOnlineMap::ApplyBoundsLimiting( std::vector< double >& dataInDouble,
                                                                         std::vector< double >& dataOutDouble,
                                                                         CAASType caasType,
                                                                         int caasIteration,
                                                                         double mismatch )
{
    // Currently only implemented for FV to FV remapping
    // We should generalize this to other types of remapping
    assert( !dataGLLNodesSrcCov.IsAttached() && !dataGLLNodesDest.IsAttached() );

    std::pair< double, double > massDefect( 0.0, 0.0 );

    // Check if the source and target data are of the same size
    const size_t nTargetCount                    = dataOutDouble.size();
    const DataArray1D< double >& m_dOverlapAreas = this->m_remapper->m_overlap->vecFaceArea;

    // Apply the offline map to the data
    double dMassDiff = 0.0;
    std::vector< double > x( nTargetCount );
    std::vector< double > dataLowerBound( nTargetCount );
    std::vector< double > dataUpperBound( nTargetCount );
    std::vector< double > massVector( nTargetCount );
    std::vector< std::unordered_set< int > > vecSourceOvTarget( nTargetCount );

#undef USE_ComputeAdjacencyRelations
    constexpr bool useMOABAdjacencies = true;
#ifdef USE_ComputeAdjacencyRelations
    // Compute the adjacent faces to the source face
    // However, calling MOAB to do this does not work correctly as we need ixS to be the index
    // Cannot just iterate over all entities in the source covering mesh
    if( caasType == CAAS_QLT || caasType == CAAS_LOCAL_ADJACENT )
    {
        if( useMOABAdjacencies )
        {
            moab::ErrorCode rval =
                ComputeAdjacencyRelations( vecSourceOvTarget, caasIteration, m_remapper->m_covering_source_entities,
                                           useMOABAdjacencies );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
        }
        else
        {
            moab::ErrorCode rval =
                ComputeAdjacencyRelations( vecSourceOvTarget, caasIteration, m_remapper->m_covering_source_entities,
                                           useMOABAdjacencies, m_meshInputCov );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
        }
    }
#else
    moab::MeshTopoUtil mtu( m_interface );
#endif

    // Initialize the bounds on the given source and target data
    double dSourceMin = dataInDouble[0];
    double dSourceMax = dataInDouble[0];
    double dTargetMin = dataOutDouble[0];
    double dTargetMax = dataOutDouble[0];
    for( size_t i = 0; i < m_meshOverlap->faces.size(); i++ )
    {
        const int ixS = m_meshOverlap->vecSourceFaceIx[i];
        const int ixT = m_meshOverlap->vecTargetFaceIx[i];

        if( ixT < 0 ) continue;  // skip ghost target faces

        assert( m_dOverlapAreas[i] > 0.0 );
        assert( ixS >= 0 );
        assert( ixT >= 0 );

#ifndef USE_ComputeAdjacencyRelations
        // Compute the adjacent faces to the target face
        vecSourceOvTarget[ixT].insert( ixS );  // map target face to source face
        if( ( caasType == CAAS_QLT || caasType == CAAS_LOCAL_ADJACENT ) )
        {
            if( useMOABAdjacencies )
            {
                moab::Range ents;
                ents.insert( m_remapper->m_covering_source_entities[ixS] );
                moab::Range adjEnts;
                moab::ErrorCode rval = mtu.get_bridge_adjacencies( ents, 0, 2, adjEnts, caasIteration );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
                for( moab::Range::iterator it = adjEnts.begin(); it != adjEnts.end(); ++it )
                {
                    int adjIndex = m_remapper->m_covering_source_entities.index( *it );
                    if( adjIndex >= 0 ) vecSourceOvTarget[ixT].insert( adjIndex );
                }
            }
            else
            {
                // Compute the adjacent faces to the target face
                AdjacentFaceVector vecAdjFaces;
                GetAdjacentFaceVectorByEdge( *m_meshInputCov, ixS,
                                             ( caasIteration ) * ( m_input_order + 1 ) * ( m_input_order + 1 ),
                                             vecAdjFaces );

                //Compute min/max over neighboring faces
                for( size_t iadj = 0; iadj < vecAdjFaces.size(); iadj++ )
                    vecSourceOvTarget[ixT].insert( vecAdjFaces[iadj].first );  // map target face to source face
            }
        }
#endif

        // Update the min and max values of the source data
        dSourceMax = fmax( dSourceMax, dataInDouble[ixS] );
        dSourceMin = fmin( dSourceMin, dataInDouble[ixS] );

        // Update the min and max values of the target data
        dTargetMin = fmin( dTargetMin, dataOutDouble[ixT] );
        dTargetMax = fmax( dTargetMax, dataOutDouble[ixT] );

        const double locMassDiff = ( dataInDouble[ixS] * m_dOverlapAreas[i] ) -  // source mass
                                   ( dataOutDouble[ixT] * m_dOverlapAreas[i] );  // target mass

        // Update the mass difference between source and target faces
        // linked to the overlap mesh element
        dMassDiff += locMassDiff;  // target mass
        massVector[ixT] += locMassDiff;
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localMinMaxDefects( 5, 0.0 ), globalMinMaxDefects( 5, 0.0 );
    localMinMaxDefects[0] = dSourceMin;
    localMinMaxDefects[1] = dTargetMin;
    localMinMaxDefects[2] = dSourceMax;
    localMinMaxDefects[3] = dTargetMax;
    localMinMaxDefects[4] = dMassDiff;

    if( caasType == CAAS_GLOBAL )
    {
        MPI_Allreduce( localMinMaxDefects.data(), globalMinMaxDefects.data(), 2, MPI_DOUBLE, MPI_MIN, m_pcomm->comm() );
        MPI_Allreduce( localMinMaxDefects.data() + 2, globalMinMaxDefects.data() + 2, 2, MPI_DOUBLE, MPI_MAX,
                       m_pcomm->comm() );
        dSourceMin = globalMinMaxDefects[0];
        dSourceMax = globalMinMaxDefects[2];
        dTargetMin = globalMinMaxDefects[1];
        dTargetMax = globalMinMaxDefects[3];
    }
    if( caasIteration == 1 )
        MPI_Allreduce( localMinMaxDefects.data() + 4, globalMinMaxDefects.data() + 4, 1, MPI_DOUBLE, MPI_SUM,
                       m_pcomm->comm() );
    else
        globalMinMaxDefects[4] = mismatch;

    dMassDiff = localMinMaxDefects[4];
    // massDefect.first = localMinMaxDefects[4];
    massDefect.first = globalMinMaxDefects[4];
#else

    // massDefect.first = fabs( dMassDiff / ( dSourceMax - dSourceMin ) );
    massDefect.first = dMassDiff;
#endif

    // Early exit if the values are monotone already.
    // if( ( dTargetMax <= dSourceMax && dTargetMin <= dSourceMin ) || fabs( massDefect.first ) < 1e-16 )
    if( fabs( massDefect.first ) > 1e-20 )
    {
        if( caasType == CAAS_GLOBAL )
        {
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = dSourceMin - dataOutDouble[i];
                dataUpperBound[i] = dSourceMax - dataOutDouble[i];
            }
        }     // if( caasType == CAAS_GLOBAL )
        else  // caasType == CAAS_LOCAL
        {
            // Compute the local min and max values of the target data
            std::vector< double > vecLocalUpperBound( nTargetCount );
            std::vector< double > vecLocalLowerBound( nTargetCount );
            // Loop over the target faces and compute the min and max values
            // of the source data linked to the target faces
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                assert( vecSourceOvTarget[i].size() );

                double dMinI = 1E10;   // dataInDouble[vecSourceOvTarget[i][0]];
                double dMaxI = -1E10;  // dataInDouble[vecSourceOvTarget[i][0]];

                // Compute max over intersecting source faces
                for( const auto& srcElem : vecSourceOvTarget[i] )
                {
                    dMinI = fmin( dMinI, dataInDouble[srcElem] );  // min over intersecting source faces
                    dMaxI = fmax( dMaxI, dataInDouble[srcElem] );  // max over intersecting source faces
                }

                // Update the min and max values of the target data
                vecLocalLowerBound[i] = dMinI;
                vecLocalUpperBound[i] = dMaxI;
            }

            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = vecLocalLowerBound[i] - dataOutDouble[i];
                dataUpperBound[i] = vecLocalUpperBound[i] - dataOutDouble[i];
            }
        }  // caasType == CAAS_LOCAL

        // Invoke CAAS or QLT application on the map
        if( fabs( dMassDiff ) > 1e-20 )
        {
            if( caasType == CAAS_QLT )
                dMassDiff = QLTLimiter( caasIteration, dataOutDouble, dataLowerBound, dataUpperBound, massVector );
            else
                CAASLimiter( dataOutDouble, dataLowerBound, dataUpperBound, dMassDiff );
        }

        // Announce output mass
        double dMassDiffPost = 0.0;
        for( size_t i = 0; i < m_meshOverlap->faces.size(); i++ )
        {
            const int ixS = m_meshOverlap->vecSourceFaceIx[i];
            const int ixT = m_meshOverlap->vecTargetFaceIx[i];

            if( ixT < 0 ) continue;  // skip ghost target faces

            // Update the mass difference between source and target faces
            // linked to the overlap mesh element
            dMassDiffPost += ( dataInDouble[ixS] * m_dOverlapAreas[i] ) -  // source mass
                             ( dataOutDouble[ixT] * m_dOverlapAreas[i] );  // target mass
        }
        // massDefect.second = fabs( dMassDiffPost / ( dSourceMax - dSourceMin ) );
        massDefect.second = dMassDiffPost;
    }

    // Ideally should perform an AllReduce here to get the global mass difference across all processors
    // But if we satisfy the constraint on every task, essentially, the global mass difference should be zero!
    return massDefect;
}

///////////////////////////////////////////////////////////////////////////////

// ** Kahan Summation Algorithm for improved numerical accuracy **
struct KahanSum
{
    double sum        = 0.0;
    double correction = 0.0;

    void add( double value )
    {
        double y   = value - correction;  // Correct the input
        double t   = sum + y;             // Perform the sum
        correction = ( t - sum ) - y;     // Update correction
        sum        = t;                   // Store the new sum
    }

    double result() const
    {
        return sum;
    }
};

// Pairwise summation helper function
inline double pairwiseSum( const std::set< double >& sorted )
{
    if( sorted.empty() ) return 0.0;
    if( sorted.size() == 1 ) return *sorted.begin();

    // Accumulate pairwise to minimize rounding error
    double sum = 0.0;
    for( double val : sorted )
        sum += val;
    return sum;
}

// Pairwise summation helper function
inline double pairwiseKahanSum( const std::set< double >& sorted )
{
    if( sorted.empty() ) return 0.0;
    if( sorted.size() == 1 ) return *sorted.begin();

    // Accumulate pairwise to minimize rounding error
    // Apply pairwise summation with Kahan correction
    KahanSum kahan;
    for( double val : sorted )
        kahan.add( val );
    return kahan.result();
}

// Sparse matrix-vector multiplication using pairwise summation
inline void deterministicSparseMatVecMul( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                          const typename moab::TempestOnlineMap::WeightColVector& x,
                                          typename moab::TempestOnlineMap::WeightRowVector& result )
{
    constexpr bool useKahanSum    = false;
    constexpr bool usePairwiseSum = false;

    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        std::set< double > accumulators;
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            // accumulators contains the sorted values of the product: A(row, col) * x(col)
            accumulators.insert( it.value() * x( it.col() ) );
        }
        if( usePairwiseSum ) result( row ) = pairwiseSum( accumulators );
        if( useKahanSum ) result( row ) = pairwiseKahanSum( accumulators );

        if( !usePairwiseSum && !useKahanSum )
        {
            double sum = 0.0;
            for( double val : accumulators )
                sum += val;
            result( row ) = sum;
        }
    }
}

//
// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatVecMulKahan( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                               const typename moab::TempestOnlineMap::WeightColVector& x,
                                               typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        KahanSum kahan;
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            double product = it.value() * x( it.col() );  // Compute product
            kahan.add( product );
        }

        result( row ) = kahan.result();
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatVecMulClean( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                               const typename moab::TempestOnlineMap::WeightColVector& x,
                                               typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            const double product = it.value() * x( it.col() );  // Compute product
            result( it.row() ) += product;
        }
    }
}

inline void deterministicSparseMatVecMulNative( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                const typename moab::TempestOnlineMap::WeightColVector& x,
                                                typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result = A * x;  // Perform the matrix-vector multiplication using Eigen3
}

// Deterministic sparse matrix-vector multiplication with A^T * x using pairwise summation
inline void deterministicSparseMatTransposeVecMul( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                   const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                   typename moab::TempestOnlineMap::WeightColVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Temporary storage for pairwise summation
    std::vector< std::set< double > > accumulators( A.cols() );

    // Iterate over A row-wise, but accumulate into result as if computing A^T * x
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            accumulators[it.col()].insert( it.value() * x( row ) );
        }
    }

    // Compute final sum using pairwise summation for each entry
    for( int col = 0; col < A.cols(); ++col )
    {
        // result( col ) = pairwiseSum( accumulators[col] );
        result( col ) = pairwiseKahanSum( accumulators[col] );
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatTransposeVecMulClean( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                        const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                        typename moab::TempestOnlineMap::WeightColVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate over A row-wise, but accumulate into result as if computing A^T * x
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            const double product = it.value() * x( row );  // Compute product
            result( it.col() ) += product;                 // Accumulate contributions to the corresponding row in A^T
        }
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatTransposeVecMulNative( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                         const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                         typename moab::TempestOnlineMap::WeightColVector& result )
{
    result = A.adjoint() * x;  // Perform the adjoint.matrix-vector multiplication using Eigen3
}
///////////////////////////////////////////////////////////////////////////////
// #define VERBOSE
moab::ErrorCode moab::TempestOnlineMap::ApplyWeights( std::vector< double >& srcVals,
                                                      std::vector< double >& tgtVals,
                                                      bool transpose )
{
    // Reset the source and target data first
    m_rowVector.setZero();
    m_colVector.setZero();

#ifdef VERBOSE
    std::stringstream sstr;
    static int callId = 0;
    callId++;
    sstr << "projection_id_" << callId << "_s_" << size << "_rk_" << rank << ".txt";
    std::ofstream output_file( sstr.str() );
#endif
    // Perform the actual projection of weights: application of weight matrix onto the source
    // solution vector

    if( transpose )
    {
        // Permute the source data first
        for( unsigned i = 0; i < srcVals.size(); ++i )
        {
            if( row_dtoc_dofmap[i] >= 0 )
                m_rowVector( row_dtoc_dofmap[i] ) = srcVals[i];  // permute and set the row (source) vector properly
        }
        // deterministicSparseMatTransposeVecMulClean( m_weightMatrix, m_rowVector, m_colVector );
        //deterministicSparseMatTransposeVecMul( m_weightMatrix, m_rowVector, m_colVector );
        // deterministicSparseMatTransposeVecMulNative( m_weightMatrix, m_rowVector, m_colVector );
        m_colVector = m_weightMatrix.adjoint() * m_rowVector;

        // Permute the resulting target data back
        for( unsigned i = 0; i < tgtVals.size(); ++i )
        {
            if( col_dtoc_dofmap[i] >= 0 )
                tgtVals[i] = m_colVector( col_dtoc_dofmap[i] );  // permute and set the row (source) vector properly
        }
    }
    else
    {
#ifdef VERBOSE
        output_file << "ColVector: " << m_colVector.size() << ", SrcVals: " << srcVals.size()
                    << ", Sizes: " << m_nTotDofs_SrcCov << ", " << col_dtoc_dofmap.size() << "\n";
#endif
        for( size_t i = 0; i < srcVals.size(); ++i )
        {
            if( col_dtoc_dofmap[i] >= 0 )
                m_colVector( col_dtoc_dofmap[i] ) = srcVals[i];  // permute and set the row (source) vector properly
#ifdef VERBOSE
            output_file << i << " " << col_gdofmap[col_dtoc_dofmap[i]] + 1 << "  " << srcVals[i] << "\n";
#endif
        }
        // deterministicSparseMatVecMulClean( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMul( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMulNative( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMulKahan( m_weightMatrix, m_colVector, m_rowVector );
        m_rowVector = m_weightMatrix * m_colVector;

        // Permute the resulting target data back
#ifdef VERBOSE
        output_file << "RowVector: " << m_rowVector.size() << ", TgtVals:" << tgtVals.size()
                    << ", Sizes: " << m_nTotDofs_Dest << ", " << row_gdofmap.size() << "\n";
#endif
        for( size_t i = 0; i < tgtVals.size(); ++i )
        {
            if( row_dtoc_dofmap[i] >= 0 )
            {
                tgtVals[i] = m_rowVector( row_dtoc_dofmap[i] );  // permute and set the row (source) vector properly
#ifdef VERBOSE
                output_file << i << " " << row_gdofmap[row_dtoc_dofmap[i]] + 1 << "  " << tgtVals[i] << "\n";
#endif
            }
        }
    }

    // if( caasType != CAAS_NONE )
    // {
    //     constexpr int nmax_caas_iterations = 5;
    //     double mismatch                    = 1.0;
    //     int caasIteration                  = 0;
    //     while( mismatch > 1e-15 &&
    //            caasIteration++ < nmax_caas_iterations )  // iterate until convergence or a maximum of 5 iterations
    //     {
    //         std::pair< double, double > mDefect = this->ApplyCAASLimiting( srcVals, tgtVals, caasType );
    //         if( m_remapper->verbose )
    //             printf( "Rank %d: -- Iteration: %d, Net original mass defect: %3.4e, mass defect post-CAAS: %3.4e\n",
    //                     m_remapper->rank, caasIteration, mDefect.first, mDefect.second );
    //         mismatch = mDefect.second;
    //     }
    // }

#ifdef VERBOSE
    output_file.flush();  // required here
    output_file.close();
#endif

    // All done with matvec application
    return moab::MB_SUCCESS;
}
// #undef VERBOSE
#endif

///////////////////////////////////////////////////////////////////////////////

extern void ForceConsistencyConservation3( const DataArray1D< double >& vecSourceArea,
                                           const DataArray1D< double >& vecTargetArea,
                                           DataArray2D< double >& dCoeff,
                                           bool fMonotone,
                                           bool fSparseConstraints = false );

///////////////////////////////////////////////////////////////////////////////

extern void ForceIntArrayConsistencyConservation( const DataArray1D< double >& vecSourceArea,
                                                  const DataArray1D< double >& vecTargetArea,
                                                  DataArray2D< double >& dCoeff,
                                                  bool fMonotone );

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapSE4_Tempest_MOAB( const DataArray3D< int >& dataGLLNodes,
                                                          const DataArray3D< double >& dataGLLJacobian,
                                                          int nMonotoneType,
                                                          bool fContinuousIn,
                                                          bool fNoConservation )
{
    // Order of the polynomial interpolant
    int nP = dataGLLNodes.GetRows();

    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( TriQuadRuleOrder );

    int TriQuadraturePoints = triquadrule.GetPoints();

    const DataArray2D< double >& TriQuadratureG = triquadrule.GetG();

    const DataArray1D< double >& TriQuadratureW = triquadrule.GetW();

    // Sample coefficients
    DataArray2D< double > dSampleCoeff( nP, nP );

    // GLL Quadrature nodes on quadrilateral elements
    DataArray1D< double > dG;
    DataArray1D< double > dW;
    GaussLobattoQuadrature::GetPoints( nP, 0.0, 1.0, dG, dW );

    // Announcements
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapSE4_Tempest_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Volume Projection\n" );
        dbgprint.printf( 0, "Triangular quadrature rule order %i\n", TriQuadRuleOrder );
        dbgprint.printf( 0, "Order of the FE polynomial interpolant: %i\n", nP );
    }

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // NodeVector from m_meshOverlap
    const NodeVector& nodesOverlap = m_meshOverlap->nodes;
    const NodeVector& nodesFirst   = m_meshInputCov->nodes;

    // Vector of source areas
    DataArray1D< double > vecSourceArea( nP * nP );

    DataArray1D< double > vecTargetArea;
    DataArray2D< double > dCoeff;

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "remapdata_" << rank << ".txt";
    std::ofstream output_file( sstr.str() );
#endif

    // Current Overlap Face
    int ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    // generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

    // Loop over all input Faces
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        if( faceFirst.edges.size() != 4 )
        {
            _EXCEPTIONT( "Only quadrilateral elements allowed for SE remapping" );
        }
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so
        // on wrt the input mesh data Then the overlap_end and overlap_begin will be correct.
        // However, the relation with MOAB and Tempest will go out of the roof

        // Determine how many overlap Faces and triangles are present
        int nOverlapFaces    = 0;
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // if( m_meshOverlap->vecTargetFaceIx[ixOverlapTemp] < 0 ) continue;  // skip ghost target faces
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 ) break;

            nOverlapFaces++;
        }

        // No overlaps
        if( nOverlapFaces == 0 ) continue;

        // Allocate remap coefficients array for meshFirst Face
        DataArray3D< double > dRemapCoeff( nP, nP, nOverlapFaces );

        // Find the local remap coefficients
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            const Face& faceOverlap = m_meshOverlap->faces[ixOverlap + j];
            if( m_meshOverlap->vecFaceArea[ixOverlap + j] < 1.e-16 )  // machine precision
            {
                Announce( "Very small overlap at index %i area polygon: (%1.10e )", ixOverlap + j,
                          m_meshOverlap->vecFaceArea[ixOverlap + j] );
                int n = faceOverlap.edges.size();
                Announce( "Number nodes: %d", n );
                for( int k = 0; k < n; k++ )
                {
                    Node nd = nodesOverlap[faceOverlap[k]];
                    Announce( "Node %d  %d  : %1.10e  %1.10e %1.10e ", k, faceOverlap[k], nd.x, nd.y, nd.z );
                }
                continue;
            }

            // #ifdef VERBOSE
            // if ( is_root )
            //     Announce ( "\tLocal ID: %i/%i = %i, areas = %2.8e", j + ixOverlap, nOverlapFaces,
            //     m_remapper->lid_to_gid_covsrc[m_meshOverlap->vecSourceFaceIx[ixOverlap + j]],
            //     m_meshOverlap->vecFaceArea[ixOverlap + j] );
            // #endif

            int nbEdges           = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center;  // not used if nbEdges == 3
            if( nbEdges > 3 )
            {  // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for( int k = 0; k < nbEdges; k++ )
                {
                    const Node& node = nodesOverlap[faceOverlap[k]];
                    center           = center + node;
                }
                center = center / nbEdges;
                center = center.Normalized();  // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriangleArea;

            // Loop over all sub-triangles of this Overlap Face
            for( int k = 0; k < nOverlapTriangles; k++ )
            {
                if( nbEdges == 3 )  // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0         = nodesOverlap[faceOverlap[0]];
                    node1         = nodesOverlap[faceOverlap[1]];
                    node2         = nodesOverlap[faceOverlap[2]];
                    dTriangleArea = CalculateFaceArea( faceOverlap, nodesOverlap );
                }
                else  // decompose polygon in triangles around the center
                {
                    node0         = center;
                    node1         = nodesOverlap[faceOverlap[k]];
                    int k1        = ( k + 1 ) % nbEdges;
                    node2         = nodesOverlap[faceOverlap[k1]];
                    nodes[0]      = center;
                    nodes[1]      = node1;
                    nodes[2]      = node2;
                    dTriangleArea = CalculateFaceArea( faceTri, nodes );
                }
                // Coordinates of quadrature Node
                for( int l = 0; l < TriQuadraturePoints; l++ )
                {
                    Node nodeQuadrature;
                    nodeQuadrature.x = TriQuadratureG[l][0] * node0.x + TriQuadratureG[l][1] * node1.x +
                                       TriQuadratureG[l][2] * node2.x;

                    nodeQuadrature.y = TriQuadratureG[l][0] * node0.y + TriQuadratureG[l][1] * node1.y +
                                       TriQuadratureG[l][2] * node2.y;

                    nodeQuadrature.z = TriQuadratureG[l][0] * node0.z + TriQuadratureG[l][1] * node1.z +
                                       TriQuadratureG[l][2] * node2.z;

                    nodeQuadrature = nodeQuadrature.Normalized();

                    // Find components of quadrature point in basis
                    // of the first Face
                    double dAlpha;
                    double dBeta;

                    ApplyInverseMap( faceFirst, nodesFirst, nodeQuadrature, dAlpha, dBeta );

                    // Check inverse map value
                    if( ( dAlpha < -1.0e-13 ) || ( dAlpha > 1.0 + 1.0e-13 ) || ( dBeta < -1.0e-13 ) ||
                        ( dBeta > 1.0 + 1.0e-13 ) )
                    {
                        _EXCEPTION4( "Inverse Map for element %d and subtriangle %d out of range "
                                     "(%1.5e %1.5e)",
                                     j, l, dAlpha, dBeta );
                    }

                    // Sample the finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nP, dAlpha, dBeta, dSampleCoeff );

                    // Add sample coefficients to the map if m_meshOverlap->vecFaceArea[ixOverlap + j] > 0
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dRemapCoeff[p][q][j] += TriQuadratureW[l] * dTriangleArea * dSampleCoeff[p][q] /
                                                    m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
            }
        }

#ifdef VERBOSE
        output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    output_file << dRemapCoeff[p][q][j] << " ";
                }
            }
        }
        output_file << std::endl;
#endif

        // Force consistency and conservation
        if( !fNoConservation )
        {
            double dTargetArea = 0.0;
            for( int j = 0; j < nOverlapFaces; j++ )
            {
                dTargetArea += m_meshOverlap->vecFaceArea[ixOverlap + j];
            }

            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    vecSourceArea[p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                }
            }

            const double areaTolerance = 1e-10;
            // Source elements are completely covered by target volumes
            if( fabs( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea ) <= areaTolerance )
            {
                vecTargetArea.Allocate( nOverlapFaces );
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }

                dCoeff.Allocate( nOverlapFaces, nP * nP );

                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }

                // Target volumes only partially cover source elements
            }
            else if( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea > areaTolerance )
            {
                double dExtraneousArea = m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea;

                vecTargetArea.Allocate( nOverlapFaces + 1 );
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }
                vecTargetArea[nOverlapFaces] = dExtraneousArea;

#ifdef VERBOSE
                Announce( "Partial volume: %i (%1.10e / %1.10e)", ixFirst, dTargetArea,
                          m_meshInputCov->vecFaceArea[ixFirst] );
#endif
                if( dTargetArea > m_meshInputCov->vecFaceArea[ixFirst] )
                {
                    _EXCEPTIONT( "Partial element area exceeds total element area" );
                }

                dCoeff.Allocate( nOverlapFaces + 1, nP * nP );

                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                    }
                }
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[nOverlapFaces][p * nP + q] -=
                                dRemapCoeff[p][q][j] * m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] /= dExtraneousArea;
                    }
                }

                // Source elements only partially cover target volumes
            }
            else
            {
                Announce( "Coverage area: %1.10e, and target element area: %1.10e)", ixFirst,
                          m_meshInputCov->vecFaceArea[ixFirst], dTargetArea );
                _EXCEPTIONT( "Target grid must be a subset of source grid" );
            }

            ForceConsistencyConservation3( vecSourceArea, vecTargetArea, dCoeff, ( nMonotoneType > 0 )
                                           /*, m_remapper->lid_to_gid_covsrc[ixFirst]*/ );

            for( int j = 0; j < nOverlapFaces; j++ )
            {
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dRemapCoeff[p][q][j] = dCoeff[j][p * nP + q];
                    }
                }
            }
        }

#ifdef VERBOSE
        // output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        // for ( int j = 0; j < nOverlapFaces; j++ )
        // {
        //     for ( int p = 0; p < nP; p++ )
        //     {
        //         for ( int q = 0; q < nP; q++ )
        //         {
        //             output_file << dRemapCoeff[p][q][j] << " ";
        //         }
        //     }
        // }
        // output_file << std::endl;
#endif

        // Put these remap coefficients into the SparseMatrix map
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];

            // signal to not participate, because it is a ghost target
            if( ixSecondFace < 0 ) continue;  // do not do anything

            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    if( fContinuousIn )
                    {
                        int ixFirstNode = dataGLLNodes[p][q][ixFirst] - 1;

                        smatMap( ixSecondFace, ixFirstNode ) += dRemapCoeff[p][q][j] *
                                                                m_meshOverlap->vecFaceArea[ixOverlap + j] /
                                                                m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                    else
                    {
                        int ixFirstNode = ixFirst * nP * nP + p * nP + q;

                        smatMap( ixSecondFace, ixFirstNode ) += dRemapCoeff[p][q][j] *
                                                                m_meshOverlap->vecFaceArea[ixOverlap + j] /
                                                                m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                }
            }
        }
        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }
#ifdef VERBOSE
    output_file.flush();  // required here
    output_file.close();
#endif

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_MOAB( const DataArray3D< int >& dataGLLNodesIn,
                                                        const DataArray3D< double >& dataGLLJacobianIn,
                                                        const DataArray3D< int >& dataGLLNodesOut,
                                                        const DataArray3D< double >& dataGLLJacobianOut,
                                                        const DataArray1D< double >& dataNodalAreaOut,
                                                        int nPin,
                                                        int nPout,
                                                        int nMonotoneType,
                                                        bool fContinuousIn,
                                                        bool fContinuousOut,
                                                        bool fNoConservation )
{
    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( 8 );

    const DataArray2D< double >& dG = triquadrule.GetG();
    const DataArray1D< double >& dW = triquadrule.GetW();

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D< double > dSampleCoeffIn( nPin, nPin );
    DataArray2D< double > dSampleCoeffOut( nPout, nPout );

    // Announcemnets
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapGLLtoGLL2_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Element Projection\n" );
        dbgprint.printf( 0, "Order of the input FE polynomial interpolant: %i\n", nPin );
        dbgprint.printf( 0, "Order of the output FE polynomial interpolant: %i\n", nPout );
    }

    // Build the integration array for each element on m_meshOverlap
    DataArray3D< double > dGlobalIntArray( nPin * nPin, m_meshOverlap->faces.size(), nPout * nPout );

    // Number of overlap Faces per source Face
    DataArray1D< int > nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        // Determine how many overlap Faces and triangles are present
        int nOverlapFaces    = 0;
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 )
            {
                break;
            }

            nOverlapFaces++;
        }

        nAllOverlapFaces[ixFirst] = nOverlapFaces;

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Geometric area of each output node
    DataArray2D< double > dGeometricOutputArea( m_meshOutput->faces.size(), nPout * nPout );

    // Area of each overlap element in the output basis
    DataArray2D< double > dOverlapOutputArea( m_meshOverlap->faces.size(), nPout * nPout );

    // Loop through all faces on m_meshInputCov
    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    if( is_root ) dbgprint.printf( 0, "Building conservative distribution maps\n" );

    // generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Quantities from the First Mesh
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector& nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        if( !nOverlapFaces ) continue;

        // // Calculate total element Jacobian
        // double dTotalJacobian = 0.0;
        // for (int s = 0; s < nPin; s++) {
        //     for (int t = 0; t < nPin; t++) {
        //         dTotalJacobian += dataGLLJacobianIn[s][t][ixFirst];
        //     }
        // }

        // Loop through all Overlap Faces
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // Quantities from the overlap Mesh
            const Face& faceOverlap = m_meshOverlap->faces[ixOverlap + i];

            const NodeVector& nodesOverlap = m_meshOverlap->nodes;

            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            // signal to not participate, because it is a ghost target
            if( ixSecond < 0 ) continue;  // do not do anything

            const NodeVector& nodesSecond = m_meshOutput->nodes;

            const Face& faceSecond = m_meshOutput->faces[ixSecond];

            int nbEdges           = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center;  // not used if nbEdges == 3
            if( nbEdges > 3 )
            {  // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for( int k = 0; k < nbEdges; k++ )
                {
                    const Node& node = nodesOverlap[faceOverlap[k]];
                    center           = center + node;
                }
                center = center / nbEdges;
                center = center.Normalized();  // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriArea;

            // Loop over all sub-triangles of this Overlap Face
            for( int j = 0; j < nOverlapTriangles; j++ )
            {
                if( nbEdges == 3 )  // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0    = nodesOverlap[faceOverlap[0]];
                    node1    = nodesOverlap[faceOverlap[1]];
                    node2    = nodesOverlap[faceOverlap[2]];
                    dTriArea = CalculateFaceArea( faceOverlap, nodesOverlap );
                }
                else  // decompose polygon in triangles around the center
                {
                    node0    = center;
                    node1    = nodesOverlap[faceOverlap[j]];
                    int j1   = ( j + 1 ) % nbEdges;
                    node2    = nodesOverlap[faceOverlap[j1]];
                    nodes[0] = center;
                    nodes[1] = node1;
                    nodes[2] = node2;
                    dTriArea = CalculateFaceArea( faceTri, nodes );
                }

                for( int k = 0; k < triquadrule.GetPoints(); k++ )
                {
                    // Get the nodal location of this point
                    double dX[3];

                    dX[0] = dG( k, 0 ) * node0.x + dG( k, 1 ) * node1.x + dG( k, 2 ) * node2.x;
                    dX[1] = dG( k, 0 ) * node0.y + dG( k, 1 ) * node1.y + dG( k, 2 ) * node2.y;
                    dX[2] = dG( k, 0 ) * node0.z + dG( k, 1 ) * node1.z + dG( k, 2 ) * node2.z;

                    double dMag = sqrt( dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2] );

                    dX[0] /= dMag;
                    dX[1] /= dMag;
                    dX[2] /= dMag;

                    Node nodeQuadrature( dX[0], dX[1], dX[2] );

                    // Find the components of this quadrature point in the basis
                    // of the first Face.
                    double dAlphaIn;
                    double dBetaIn;

                    ApplyInverseMap( faceFirst, nodesFirst, nodeQuadrature, dAlphaIn, dBetaIn );

                    // Find the components of this quadrature point in the basis
                    // of the second Face.
                    double dAlphaOut;
                    double dBetaOut;

                    ApplyInverseMap( faceSecond, nodesSecond, nodeQuadrature, dAlphaOut, dBetaOut );

                    /*
                                        // Check inverse map value
                                        if ((dAlphaIn < 0.0) || (dAlphaIn > 1.0) ||
                                            (dBetaIn  < 0.0) || (dBetaIn  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaIn, dBetaIn);
                                        }

                                        // Check inverse map value
                                        if ((dAlphaOut < 0.0) || (dAlphaOut > 1.0) ||
                                            (dBetaOut  < 0.0) || (dBetaOut  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaOut, dBetaOut);
                                        }
                    */
                    // Sample the First finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPin, dAlphaIn, dBetaIn, dSampleCoeffIn );

                    // Sample the Second finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPout, dAlphaOut, dBetaOut, dSampleCoeffOut );

                    // Overlap output area
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            double dNodeArea = dSampleCoeffOut[s][t] * dW[k] * dTriArea;

                            dOverlapOutputArea[ixOverlap + i][s * nPout + t] += dNodeArea;

                            dGeometricOutputArea[ixSecond][s * nPout + t] += dNodeArea;
                        }
                    }

                    // Compute overlap integral
                    int ixp = 0;
                    for( int p = 0; p < nPin; p++ )
                    {
                        for( int q = 0; q < nPin; q++ )
                        {
                            int ixs = 0;
                            for( int s = 0; s < nPout; s++ )
                            {
                                for( int t = 0; t < nPout; t++ )
                                {
                                    // Sample the Second finite element at this point
                                    dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
                                        dSampleCoeffOut[s][t] * dSampleCoeffIn[p][q] * dW[k] * dTriArea;

                                    ixs++;
                                }
                            }

                            ixp++;
                        }
                    }
                }
            }
        }

        // Coefficients
        DataArray2D< double > dCoeff( nOverlapFaces * nPout * nPout, nPin * nPin );

        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            dCoeff[i * nPout * nPout + ixs][ixp] = dGlobalIntArray[ixp][ixOverlap + i][ixs] /
                                                                   dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Source areas
        DataArray1D< double > vecSourceArea( nPin * nPin );

        for( int p = 0; p < nPin; p++ )
        {
            for( int q = 0; q < nPin; q++ )
            {
                vecSourceArea[p * nPin + q] = dataGLLJacobianIn[p][q][ixFirst];
            }
        }

        // Target areas
        DataArray1D< double > vecTargetArea( nOverlapFaces * nPout * nPout );

        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];
            int ixs = 0;
            for( int s = 0; s < nPout; s++ )
            {
                for( int t = 0; t < nPout; t++ )
                {
                    vecTargetArea[i * nPout * nPout + ixs] = dOverlapOutputArea[ixOverlap + i][nPout * s + t];

                    ixs++;
                }
            }
        }

        // Force consistency and conservation
        if( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation( vecSourceArea, vecTargetArea, dCoeff, ( nMonotoneType != 0 ) );
        }

        // Update global coefficients
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            dGlobalIntArray[ixp][ixOverlap + i][ixs] =
                                dCoeff[i * nPout * nPout + ixs][ixp] * dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

#ifdef VVERBOSE
        // Check column sums (conservation)
        for( int i = 0; i < nPin * nPin; i++ )
        {
            double dColSum = 0.0;
            for( int j = 0; j < nOverlapFaces * nPout * nPout; j++ )
            {
                dColSum += dCoeff[j][i] * vecTargetArea[j];
            }
            printf( "Col %i: %1.15e\n", i, dColSum / vecSourceArea[i] );
        }

        // Check row sums (consistency)
        for( int j = 0; j < nOverlapFaces * nPout * nPout; j++ )
        {
            double dRowSum = 0.0;
            for( int i = 0; i < nPin * nPin; i++ )
            {
                dRowSum += dCoeff[j][i];
            }
            printf( "Row %i: %1.15e\n", j, dRowSum );
        }
#endif

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Build redistribution map within target element
    if( is_root ) dbgprint.printf( 0, "Building redistribution maps on target mesh\n" );
    DataArray1D< double > dRedistSourceArea( nPout * nPout );
    DataArray1D< double > dRedistTargetArea( nPout * nPout );
    std::vector< DataArray2D< double > > dRedistributionMaps;
    dRedistributionMaps.resize( m_meshOutput->faces.size() );

    for( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        dRedistributionMaps[ixSecond].Allocate( nPout * nPout, nPout * nPout );

        for( int i = 0; i < nPout * nPout; i++ )
        {
            dRedistributionMaps[ixSecond][i][i] = 1.0;
        }

        for( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistSourceArea[s] = dGeometricOutputArea[ixSecond][s];
        }

        for( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistTargetArea[s] = dataGLLJacobianOut[s / nPout][s % nPout][ixSecond];
        }

        if( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation( dRedistSourceArea, dRedistTargetArea, dRedistributionMaps[ixSecond],
                                                  ( nMonotoneType != 0 ) );

            for( int s = 0; s < nPout * nPout; s++ )
            {
                for( int t = 0; t < nPout * nPout; t++ )
                {
                    dRedistributionMaps[ixSecond][s][t] *= dRedistTargetArea[s] / dRedistSourceArea[t];
                }
            }
        }
    }

    // Construct the total geometric area
    DataArray1D< double > dTotalGeometricArea( dataNodalAreaOut.GetRows() );
    for( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        for( int s = 0; s < nPout; s++ )
        {
            for( int t = 0; t < nPout; t++ )
            {
                dTotalGeometricArea[dataGLLNodesOut[s][t][ixSecond] - 1] +=
                    dGeometricOutputArea[ixSecond][s * nPout + t];
            }
        }
    }

    // Compose the integration operator with the output map
    ixOverlap = 0;

    if( is_root ) dbgprint.printf( 0, "Assembling map\n" );

    // Map from source DOFs to target DOFs with redistribution applied
    DataArray2D< double > dRedistributedOp( nPin * nPin, nPout * nPout );

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        if( !nOverlapFaces ) continue;

        // Put composed array into map
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];

            // signal to not participate, because it is a ghost target
            if( ixSecondFace < 0 ) continue;  // do not do anything

            dRedistributedOp.Zero();
            for( int p = 0; p < nPin * nPin; p++ )
            {
                for( int s = 0; s < nPout * nPout; s++ )
                {
                    for( int t = 0; t < nPout * nPout; t++ )
                    {
                        dRedistributedOp[p][s] +=
                            dRedistributionMaps[ixSecondFace][s][t] * dGlobalIntArray[p][ixOverlap + j][t];
                    }
                }
            }

            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixFirstNode;
                    if( fContinuousIn )
                    {
                        ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                    }
                    else
                    {
                        ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
                    }

                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            int ixSecondNode;
                            if( fContinuousOut )
                            {
                                ixSecondNode = dataGLLNodesOut[s][t][ixSecondFace] - 1;

                                if( !fNoConservation )
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dataNodalAreaOut[ixSecondNode];
                                }
                                else
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dTotalGeometricArea[ixSecondNode];
                                }
                            }
                            else
                            {
                                ixSecondNode = ixSecondFace * nPout * nPout + s * nPout + t;

                                if( !fNoConservation )
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dataGLLJacobianOut[s][t][ixSecondFace];
                                }
                                else
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dGeometricOutputArea[ixSecondFace][s * nPout + t];
                                }
                            }

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_Pointwise_MOAB( const DataArray3D< int >& dataGLLNodesIn,
                                                                  const DataArray3D< double >& /*dataGLLJacobianIn*/,
                                                                  const DataArray3D< int >& dataGLLNodesOut,
                                                                  const DataArray3D< double >& /*dataGLLJacobianOut*/,
                                                                  const DataArray1D< double >& dataNodalAreaOut,
                                                                  int nPin,
                                                                  int nPout,
                                                                  int nMonotoneType,
                                                                  bool fContinuousIn,
                                                                  bool fContinuousOut )
{
    // Gauss-Lobatto quadrature within Faces
    DataArray1D< double > dGL;
    DataArray1D< double > dWL;

    GaussLobattoQuadrature::GetPoints( nPout, 0.0, 1.0, dGL, dWL );

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D< double > dSampleCoeffIn( nPin, nPin );

    // Announcemnets
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapGLLtoGLL2_Pointwise_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Element (Pointwise) Projection\n" );
        dbgprint.printf( 0, "Order of the input FE polynomial interpolant: %i\n", nPin );
        dbgprint.printf( 0, "Order of the output FE polynomial interpolant: %i\n", nPout );
    }

    // Number of overlap Faces per source Face
    DataArray1D< int > nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];

            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 ) break;

            nAllOverlapFaces[ixFirst]++;
        }

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Number of times this point was found
    DataArray1D< bool > fSecondNodeFound( dataNodalAreaOut.GetRows() );

    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    // Loop through all faces on m_meshInputCov
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Quantities from the First Mesh
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector& nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        // Loop through all Overlap Faces
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            // signal to not participate, because it is a ghost target
            if( ixSecond < 0 ) continue;  // do not do anything

            const NodeVector& nodesSecond = m_meshOutput->nodes;
            const Face& faceSecond        = m_meshOutput->faces[ixSecond];

            // Loop through all nodes on the second face
            for( int s = 0; s < nPout; s++ )
            {
                for( int t = 0; t < nPout; t++ )
                {
                    size_t ixSecondNode;
                    if( fContinuousOut )
                    {
                        ixSecondNode = dataGLLNodesOut[s][t][ixSecond] - 1;
                    }
                    else
                    {
                        ixSecondNode = ixSecond * nPout * nPout + s * nPout + t;
                    }

                    if( ixSecondNode >= fSecondNodeFound.GetRows() ) _EXCEPTIONT( "Logic error" );

                    // Check if this node has been found already
                    if( fSecondNodeFound[ixSecondNode] ) continue;

                    // Check this node
                    Node node;
                    Node dDx1G;
                    Node dDx2G;

                    ApplyLocalMap( faceSecond, nodesSecond, dGL[t], dGL[s], node, dDx1G, dDx2G );

                    // Find the components of this quadrature point in the basis
                    // of the first Face.
                    double dAlphaIn;
                    double dBetaIn;

                    ApplyInverseMap( faceFirst, nodesFirst, node, dAlphaIn, dBetaIn );

                    // Check if this node is within the first Face
                    if( ( dAlphaIn < -1.0e-10 ) || ( dAlphaIn > 1.0 + 1.0e-10 ) || ( dBetaIn < -1.0e-10 ) ||
                        ( dBetaIn > 1.0 + 1.0e-10 ) )
                        continue;

                    // Node is within the overlap region, mark as found
                    fSecondNodeFound[ixSecondNode] = true;

                    // Sample the First finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPin, dAlphaIn, dBetaIn, dSampleCoeffIn );

                    // Add to map
                    for( int p = 0; p < nPin; p++ )
                    {
                        for( int q = 0; q < nPin; q++ )
                        {
                            int ixFirstNode;
                            if( fContinuousIn )
                            {
                                ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                            }
                            else
                            {
                                ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
                            }

                            smatMap( ixSecondNode, ixFirstNode ) += dSampleCoeffIn[p][q];
                        }
                    }
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Check for missing samples
    for( size_t i = 0; i < fSecondNodeFound.GetRows(); i++ )
    {
        if( !fSecondNodeFound[i] )
        {
            _EXCEPTION1( "Can't sample point %i", i );
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////
// ----------------------
// KD-tree point cloud
// ----------------------
struct MOABCentroidCloud
{
    constexpr static int DIM = 3;
    std::vector< std::array< double, DIM > > points;
    std::vector< moab::EntityHandle > elements;

    inline void init( size_t length )
    {
        points.reserve( length );
        elements.reserve( length );
    }

    inline size_t kdtree_get_point_count() const
    {
        return points.size();
    }

    inline double kdtree_get_pt( const size_t idx, const size_t dim ) const
    {
        return points[idx][dim];
    }

    template < class BBOX >
    bool kdtree_get_bbox( BBOX& ) const
    {
        return false;
    }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, MOABCentroidCloud >,
                                                    MOABCentroidCloud,  // DatasetAdaptor
                                                    3,                  // 3D
                                                    size_t              // IndexType
                                                    >;

// ----------------------
// Radius search wrapper
// ----------------------
std::vector< size_t > radius_search_kdtree( const KDTree& tree,
                                            const MOABCentroidCloud& cloud,
                                            const std::array< double, 3 >& query_pt,
                                            double radius )
{
    double radius_sq = radius * radius;
    // double radius_sq = radius;
    std::vector< nanoflann::ResultItem< size_t, double > > matches;
    nanoflann::SearchParameters params;
    params.sorted = true;

    const double query_pt_sq = query_pt[0] * query_pt[0] + query_pt[1] * query_pt[1] + query_pt[2] * query_pt[2];
    assert( query_pt_sq > 1.0 - 1e-10 && query_pt_sq < 1.0 + 1.0e-10 ); // Check if the point is on the unit sphere

    tree.radiusSearch( query_pt.data(), radius_sq, matches, params );

    std::vector< size_t > found_elements;
    if( !matches.empty() )
    {
        // Return all matches within the radius
        for( const auto& match : matches )
            found_elements.emplace_back( cloud.elements[match.first] );
    }
    else
    {
        // Fallback to nearest neighbor
        size_t nearest_index;
        double nearest_dist_sq;
        tree.knnSearch( query_pt.data(), 1, &nearest_index, &nearest_dist_sq );
        found_elements.emplace_back( cloud.elements[nearest_index] );
    }

    return found_elements;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::LinearRemapFVtoGLL_Averaged( const DataArray3D< int >& dataGLLNodes,
                                                                     const DataArray3D< double >& dataGLLJacobian,
                                                                     const DataArray1D< double >& /*dataGLLNodalArea*/,
                                                                     int nOrder,
                                                                     bool fContinuous )
{
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if( m_meshInputCov->revnodearray.size() == 0 )
    {
        m_meshInputCov->ConstructReverseNodeArray();
    }
    if( m_meshInputCov->edgemap.size() == 0 )
    {
        m_meshInputCov->ConstructEdgeMap( false );
    }

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // Order of the finite element method
    int nP = dataGLLNodes.GetRows();

    // GLL nodes
    DataArray1D< double > dG;
    DataArray1D< double > dW;

    GaussLobattoQuadrature::GetPoints( nP, 0.0, 1.0, dG, dW );

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( TriQuadRuleOrder );

    // Number of elements needed
#ifdef RECTANGULAR_TRUNCATION
    int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION
    int nCoefficients = nOrder * ( nOrder + 1 ) / 2;
#endif

    // Announcements
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapFVtoSE_Averaged]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Volume to Spectral Element Projection\n" );
        dbgprint.printf( 0, "Triangular quadrature rule order %i\n", TriQuadRuleOrder );
        dbgprint.printf( 0, "Number of coefficients: %i\n", nCoefficients );
    }

    // Loop through all faces on meshInput
// #ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshOutput->faces.size() / 10 ) + 1;
    // const unsigned outputFrequency = 1;
// #endif

    // kd-tree for nearest neighbor search
    // kdtree* kdSource = kd_create( 3 );

    Range& source_vertices = m_remapper->m_covering_source_vertices;

    MOABCentroidCloud cloud;
    // {
    //     const moab::Range& elems = m_remapper->m_covering_source_entities;

    //     cloud.points( elems.size() * 3 );
    //     // Get coordinates of the elements
    //     MB_CHK_ERR( m_interface->get_coords( elems, cloud.points.data() ) );

    //     // Loop through all elements and add to the tree
    //     cloud.elements.reserve( elems.size() );
    //     for( auto elem : elems )
    //     {
    //         Node nodeRef = GetFaceCentroid( faceSecond, m_meshOutput->nodes );

    //         // Get the centroid of the element
    //         cloud.elements.emplace_back( elem );
    //     }
    // }
    {
        // Initialize the kd-tree
        cloud.init( source_vertices.size() );

        std::vector<double> srccoords(source_vertices.size() * 3);
        MB_CHK_ERR( m_interface->get_coords( source_vertices, srccoords.data() ) );

        // Loop through all elements and add to the tree
        for( size_t ielem = 0; ielem < source_vertices.size(); ielem++ )
        {
            const size_t offset = ielem * 3;
            const double query_pt_sq =
                std::sqrt( srccoords[offset] * srccoords[offset] + srccoords[offset+1] * srccoords[offset+1] + srccoords[offset+2] * srccoords[offset+2] );

            // Rescale the coordinates to the unit sphere and add to the point cloud
            cloud.points.emplace_back( std::array< double, 3 >( { srccoords[offset]/query_pt_sq, srccoords[offset+1]/query_pt_sq, srccoords[offset+2]/query_pt_sq } ) );

            // Get the vertex index
            cloud.elements.emplace_back( ielem );
        }
    }

    // if( is_root ) dbgprint.printf( 0, "Building Kd-tree now..." );
    KDTree tree( 3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams( 10 ) );
    tree.buildIndex();
    // if( is_root ) dbgprint.printf( 0, "Finished building Kd-tree index..." );
    for( size_t ixOutput = 0; ixOutput < m_meshOutput->faces.size(); ixOutput++ )
    {
        // Output every 1000 elements
// #ifdef VERBOSE
        if( ixOutput % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixOutput, m_meshOutput->faces.size() );
        }
// #endif
        // This Face
        const Face& faceSecond = m_meshOutput->faces[ixOutput];

        // Area of the First Face
        // double dSecondArea = m_meshOutput->vecFaceArea[ixOutput];
        // Node nodecenter = GetFaceCentroid( faceSecond, m_meshOutput->nodes );

        for( int p = 0; p < nP; p++ )
        {
            for( int q = 0; q < nP; q++ )
            {
                int ixOutputGlobal;
                if( fContinuous )
                    ixOutputGlobal = dataGLLNodes[p][q][ixOutput] - 1;
                else
                    ixOutputGlobal = ixOutput * nP * nP + p * nP + q;

                // Coordinate of the GLL point
                // Node nodeRef = m_meshOutput->nodes[ixOutputGlobal];

                Node nodeRef;
                ApplyLocalMap( faceSecond, m_meshOutput->nodes, dG[p], dG[q], nodeRef );

                const std::array< double, 3 > query = { nodeRef.x, nodeRef.y, nodeRef.z };

                // The radius for the search is the sqrt of the Jacobian
                const double radius   = std::sqrt( dataGLLJacobian[p][q][ixOutput] );

                std::cout << "Searching elements within radius " << radius
                          << " for query point: " << query[0] << ", " << query[1] << ", " << query[2] << "\n";
                // Now let us search for the nearest elements within search radius
                auto results          = radius_search_kdtree( tree, cloud, query, radius );

                // Find how many elements were found
                size_t nResults = results.size();

                // Newton-Cotes equal-weight method
                // const double dWeight  = dataGLLNodalArea[ixOutput] / nResults;
                const double dWeight = 1.0 / nResults;

                // Find nearest source mesh face and add its contribution
                // to the inverse distance.
                // kdres* kdresSource = kd_nearest_range( kdSource, query, radius );

                // nResults = kd_res_size( kdresSource );
                // if ( nResults == 0 )
                // {
                //     kd_res_free( kdresSource );
                //     // Find nearest source mesh face and add its contribution
                //     // to the inverse distance.
                //     kdresSource = kd_nearest3( kdSource, query[0], query[1], query[2] );
                //     nResults    = kd_res_size( kdresSource );
                // }
                // double pos[3];
                // while( !kd_res_end( kdresSource ) )
                // {
                //     /* get the data and position of the current result item */
                //     int* pFace = (int*)kd_res_item( kdresSource, pos );

                //     /* compute the distance of the current result from the pt */

                //     std::cout << "\t[ " << ixOutputGlobal << "] Found association of point "
                //               << query[0] << ", " << query[1] << ", " << query[2] << " to " << *pFace << "\n";
                //     smatMap( ixOutputGlobal, *pFace ) += dWeight;

                //     /* go to the next entry */
                //     kd_res_next( kdresSource );
                // }

                if (nResults == 0)
                {
                    _EXCEPTION4( "No elements found within radius %f for query point: %f, %f, %f",
                                 radius, query[0], query[1], query[2] );
                }

                std::cout << "\tFound " << nResults << " elements within radius " << radius
                          << " for query point: " << query[0] << ", " << query[1] << ", " << query[2] << "\n";
                for( auto ixFirstElement : results )
                {
                    // std::cout << "\t[ " << ixOutputGlobal << "] Found association of point " << query[0] << ", "
                    //           << query[1] << ", " << query[2] << " to " << ixFirstElement << "\n";
                    if ( ixFirstElement < 0 || ixFirstElement >= source_vertices.size() )
                    {
                        _EXCEPTION3( "Logic error: source element has to be between 0 and %d, but received %d for row %d\n",
                                     source_vertices.size(), ixFirstElement, ixOutputGlobal );
                    }
                    smatMap( ixOutputGlobal, ixFirstElement ) += dWeight;
                }

                // kd_res_free( kdresSource );
            }
        }
    }

    // kd_free( kdSource );
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
