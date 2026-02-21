#cython: language_level=3

"""Implements parallel communication functionality for MOAB."""

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.set cimport set as cpp_set

from pymoab cimport moab
from pymoab.moab cimport EntityHandle, Tag, ErrorCode, Range, Interface
from mpi4py.MPI cimport Comm as MPIComm
from mpi4py.libmpi cimport MPI_Comm, MPI_Op
from .core cimport Core

# C++ ParallelComm class declarations (directly inlined here because
# Cython's include directive doesn't export declarations to companion .pyx,
# and the Range type conflicts between C++ cppclass and Python cdef class).
# This file is only compiled when MPI is enabled (CMake excludes parallelcomm.pyx otherwise).
cdef extern from "moab/ParallelComm.hpp" namespace "moab":

    cdef cppclass cParallelComm "moab::ParallelComm":
        # Constructors
        cParallelComm(Interface* impl, MPI_Comm cm) except +
        cParallelComm(Interface* impl, MPI_Comm cm, int* id) except +

        # Rank, size, communicator
        unsigned int rank() except +
        unsigned int size() except +
        MPI_Comm comm() except +

        # Shared entity resolution
        ErrorCode resolve_shared_ents(EntityHandle this_set, int resolve_dim,
                                      int shared_dim, const Tag* id_tag) except +
        ErrorCode resolve_shared_sets(EntityHandle this_set, Tag id_tag) except +

        # Shared entity queries
        ErrorCode get_shared_entities(int other_proc, Range& shared_ents,
                                      int dim, bool iface, bool owned_filter) except +
        # Ghost cell exchange
        ErrorCode exchange_ghost_cells(int ghost_dim, int bridge_dim, int num_layers,
                                       int addl_ents, bool store_remote_handles,
                                       bool wait_all, EntityHandle* file_set) except +

        # Parallel status
        ErrorCode filter_pstatus(Range& entities, unsigned char pstatus,
                                  unsigned char operation, int target_proc,
                                  Range* returned_ents) except +
        ErrorCode get_pstatus(EntityHandle entity, unsigned char& pstatus_val) except +
        ErrorCode get_pstatus_entities(int dim, unsigned char pstatus_val,
                                       Range& pstatus_ents) except +

        # Ownership
        ErrorCode get_owner(EntityHandle entity, int& owner) except +
        ErrorCode get_owner_handle(EntityHandle entity, int& owner,
                                    EntityHandle& handle) except +

        # Sharing data
        ErrorCode get_sharing_data(EntityHandle entity, int* ps, EntityHandle* hs,
                                    unsigned char& pstat, unsigned int& num_ps) except +
        ErrorCode get_interface_procs(cpp_set[unsigned int]& procs, bool get_buffs) except +
        ErrorCode get_comm_procs(cpp_set[unsigned int]& procs) except +

        # Global IDs
        ErrorCode assign_global_ids(EntityHandle this_set, int dimension, int start_id,
                                     bool largest_dim_only, bool parallel,
                                     bool owned_only) except +
        ErrorCode check_global_ids(EntityHandle this_set, int dimension, int start_id,
                                    bool largest_dim_only, bool parallel,
                                    bool owned_only) except +

        # Tag exchange and reduction
        ErrorCode exchange_tags(vector[Tag]& src_tags, vector[Tag]& dst_tags,
                                Range& entities) except +
        ErrorCode reduce_tags(vector[Tag]& src_tags, vector[Tag]& dst_tags,
                               MPI_Op mpi_op, Range& entities) except +

        # Partitioning
        ErrorCode get_part_entities(Range& ents, int dim) except +
        EntityHandle get_partitioning() except +
        ErrorCode set_partitioning(EntityHandle h) except +
        ErrorCode get_global_part_count(int& count) except +
        ErrorCode get_part_owner(int part_id, int& owner) except +
        ErrorCode get_part_id(EntityHandle part, int& id) except +
        ErrorCode create_part(EntityHandle& part) except +
        ErrorCode destroy_part(EntityHandle part) except +

        # Broadcast, scatter, gather
        ErrorCode broadcast_entities(int from_proc, Range& entities,
                                      bool adjacencies, bool tags) except +
        ErrorCode scatter_entities(int from_proc, vector[Range]& entities,
                                    bool adjacencies, bool tags) except +
        ErrorCode gather_data(Range& gather_ents, Tag& tag_handle, Tag id_tag,
                               EntityHandle gather_set, int root_proc_rank) except +

        # Entity deletion
        ErrorCode delete_entities(Range& to_delete) except +

        # Interface entities
        ErrorCode get_iface_entities(int other_proc, int dim,
                                      Range& iface_ents) except +

        # Debugging
        void set_debug_verbosity(int verb) except +
        int get_debug_verbosity() except +

        # Shared handle verification
        ErrorCode check_all_shared_handles(bool print_em) except +

# Python class declaration
cdef class ParallelComm:
    """Parallel communication class for MOAB."""
    cdef cParallelComm* inst
    cdef moab.Core* core
    cdef MPIComm _comm
    cdef bint _mpi_basic
    cdef bint _mpi_io
    cdef bint _mpi_full
