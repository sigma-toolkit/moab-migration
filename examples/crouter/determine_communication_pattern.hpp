#ifndef DETERMINE_COMMUNICATION_PATTERN_HPP
#define DETERMINE_COMMUNICATION_PATTERN_HPP

#include "moab/ParallelComm.hpp"
#include "moab/TupleList.hpp"
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <set>

struct CommunicationPattern {
    std::map<int, std::vector<int>> send_map; // rank -> [global_rows]
    std::map<int, std::vector<int>> recv_map; // rank -> [global_rows]
};

// Function to determine the communication pattern for the initial reduction phase
moab::ErrorCode determine_communication_pattern(
    moab::ParallelComm* pcomm,
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& local_matrix_data,
    CommunicationPattern& comm_pattern,
    std::map<int, std::set<int>>& all_row_sharers);

// Function to perform the two-stage distributed reduction and synchronization
moab::ErrorCode perform_distributed_reduction(
    moab::ParallelComm* pcomm,
    const CommunicationPattern& comm_pattern,
    const std::map<int, std::set<int>>& all_row_sharers,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& local_matrix_data);

// Function to determine communication pattern from MOAB tag data
moab::ErrorCode determine_communication_pattern_from_tag(
    moab::ParallelComm* pcomm,
    const char* tag_name,
    moab::Tag dof_tag,
    CommunicationPattern& comm_pattern,
    std::map<int, std::set<int>>& all_dof_sharers);

// Function to clear log files at the start of a run
void clear_log_files(int rank, int size);

#endif // DETERMINE_COMMUNICATION_PATTERN_HPP
