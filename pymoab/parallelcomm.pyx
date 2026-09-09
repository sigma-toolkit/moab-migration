# cython: language_level=3
# cython: embedsignature=True

from cython.operator cimport dereference as deref

cimport numpy as np
import numpy as np

from pymoab cimport moab
from pymoab cimport eh
from .rng cimport Range
from .core cimport Core
from .types import check_error, np_tag_type, validate_type, _convert_array, _eh_array
from . import types

from pymoab.rng import intersect, subtract, unite

from libcpp.vector cimport vector
from libcpp.string cimport string as std_string
from libcpp.set cimport set as cpp_set
from libc.stdlib cimport malloc, free
from libcpp cimport bool

# Import EntityHandle and Tag from moab
from pymoab.moab cimport EntityHandle, Tag, ErrorCode, Interface

# Import MPI types from mpi4py
from mpi4py import MPI
from mpi4py.MPI cimport Comm as MPIComm
from mpi4py.libmpi cimport MPI_Comm, MPI_Op, MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD
from mpi4py.libmpi cimport MPI_LAND, MPI_LOR, MPI_BAND, MPI_BOR, MPI_BXOR

# Include PSTATUS macros
include "moab_defs.pxd"

cdef void* null = NULL


cdef class ParallelComm(object):
    """Parallel communication class for MOAB.

    This class provides parallel I/O and communication capabilities for MOAB meshes.
    It handles loading and writing mesh files in parallel, resolving shared entities
    between processes, and managing ghost entities.
    """

    def __cinit__(self, Core core not None, MPIComm comm not None):
        """Constructor

        Parameters
        ----------
        core : Core
            The MOAB Core instance to use
        comm : MPI.Comm
            The MPI communicator to use for parallel operations
        """
        # Get the C MPI communicator handle from mpi4py
        cdef MPI_Comm c_comm = comm.ob_mpi

        self.core = core.inst
        self._core_owner = core
        self._comm = comm
        self._mpi_basic = True
        self._mpi_io = True
        self._mpi_full = True

        # Create MOAB ParallelComm with C communicator
        self.inst = new cParallelComm(<moab.Interface*>(self.core), c_comm)

    def __dealloc__(self):
        if self.inst != null:
            del self.inst
            self.inst = NULL

    @property
    def has_basic_mpi(self):
        """Check if basic MPI support is available."""
        return self._mpi_basic

    @property
    def has_mpi_io(self):
        """Check if MPI I/O support is available."""
        return self._mpi_io

    @property
    def has_full_mpi(self):
        """Check if full MPI support is available."""
        return self._mpi_full

    def load_file(self, str file_name, str read_opts="", EntityHandle file_set=0):
        """Load a mesh file in parallel.

        Parameters
        ----------
        file_name : str
            Name of the file to load
        read_opts : str, optional
            Options for reading the file (e.g., "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")
        file_set : EntityHandle, optional
            Handle to the file set to load into

        Raises
        ------
        RuntimeError
            If file loading fails
        """
        cdef moab.ErrorCode err
        cdef bytes b_file_name = file_name.encode('utf-8')
        cdef bytes b_read_opts
        cdef eh.EntityHandle* ptr = NULL
        cdef const char* c_file_name = b_file_name
        cdef const char* c_read_opts = NULL

        if file_set != 0:
            fset = file_set
            ptr = &fset
        if read_opts:
            b_read_opts = read_opts.encode('utf-8')
            c_read_opts = b_read_opts

        err = self.core.load_file(c_file_name, ptr, c_read_opts)
        check_error(err)

    def write_file(self, str file_name, str write_opts="", EntityHandle file_set=0,
                  list tag_list=None):
        """Write a mesh file in parallel.

        Parameters
        ----------
        file_name : str
            Name of the file to write
        write_opts : str, optional
            Options for writing the file (e.g., "PARALLEL=WRITE_PART")
        file_set : EntityHandle, optional
            Handle to the file set to write
        tag_list : list of str, optional
            List of tag names to write

        Raises
        ------
        RuntimeError
            If file writing fails
        """
        cdef moab.ErrorCode err
        cdef bytes b_file_name = file_name.encode('utf-8')
        cdef bytes b_write_opts
        cdef const char* c_file_name = b_file_name
        cdef const char* c_write_opts = NULL
        cdef const char* c_format_type = NULL
        cdef vector[moab.Tag] tag_ptrs
        cdef eh.EntityHandle* ptr = NULL
        cdef eh.EntityHandle fset
        cdef int num_sets = 1
        cdef int num_tags = 0
        cdef moab.Tag tag_ptr = NULL

        if file_set != 0:
            fset = file_set
            ptr = &fset
        else:
            num_sets = 0

        if write_opts:
            b_write_opts = write_opts.encode('utf-8')
            c_write_opts = b_write_opts

        if tag_list is not None:
            for tag_name in tag_list:
                err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_ptr)
                check_error(err)
                if err == moab.MB_SUCCESS and tag_ptr != NULL:
                    tag_ptrs.push_back(tag_ptr)
            num_tags = tag_ptrs.size()

        err = self.core.write_file(c_file_name, c_format_type, c_write_opts, ptr, num_sets,
                                   tag_ptrs.data() if num_tags > 0 else NULL,
                                   num_tags)
        check_error(err)

    def get_shared_entities(self, int other_proc=-1, int dim=-1, bint iface=False,
                            bint owned_filter=False):
        """Get entities shared with another process.

        Parameters
        ----------
        other_proc : int, optional
            Process rank to get shared entities with (-1 for all processes)
        dim : int, optional
            Dimension of entities to get (-1 for all dimensions)
        iface : bool, optional
            If True, only return interface entities
        owned_filter : bool, optional
            If True, only return owned entities

        Returns
        -------
        Range
            Range containing shared entities
        """
        cdef moab.ErrorCode err
        cdef Range shared_ents = Range()

        err = self.inst.get_shared_entities(other_proc, deref(shared_ents.inst),
                                            dim, iface, owned_filter)
        check_error(err)
        return shared_ents

    def resolve_shared_ents(self, EntityHandle this_set=0, int resolve_dim=3,
                            int shared_dim=-1):
        """Resolve shared entities between processors.

        Parameters
        ----------
        this_set : EntityHandle, optional
            Set containing entities to resolve
        resolve_dim : int, optional
            Dimension of entities to resolve (default 3)
        shared_dim : int, optional
            Dimension of shared entities (-1 for resolve_dim - 1)

        Raises
        ------
        RuntimeError
            If entity resolution fails
        """
        cdef moab.ErrorCode err
        err = self.inst.resolve_shared_ents(this_set, resolve_dim, shared_dim, NULL)
        check_error(err)

    def exchange_ghost_cells(self, int ghost_dim, int bridge_dim, int num_layers=1,
                             int addl_ents=0, bint store_remote_handles=True,
                             bint wait_all=True, EntityHandle file_set=0):
        """Exchange ghost cells between processors.

        Parameters
        ----------
        ghost_dim : int
            Dimension of ghost entities
        bridge_dim : int
            Dimension of bridge entities
        num_layers : int, optional
            Number of ghost layers
        addl_ents : int, optional
            Additional entities to include
        store_remote_handles : bool, optional
            Whether to store remote handles
        wait_all : bool, optional
            Whether to wait for all communications
        file_set : EntityHandle, optional
            File set containing entities

        Raises
        ------
        RuntimeError
            If ghost cell exchange fails
        """
        cdef moab.ErrorCode err
        cdef EntityHandle* fset_ptr = NULL
        cdef EntityHandle fset
        if file_set != 0:
            fset = file_set
            fset_ptr = &fset
        err = self.inst.exchange_ghost_cells(ghost_dim, bridge_dim, num_layers,
                                             addl_ents, store_remote_handles,
                                             wait_all, fset_ptr)
        check_error(err)

    @property
    def rank(self):
        """Get the MPI rank of this processor."""
        return self.inst.rank()

    @property
    def size(self):
        """Get the total number of MPI processes."""
        return self.inst.size()

    @property
    def comm(self):
        """Get the mpi4py communicator object."""
        return self._comm

    def get_pstatus(self, EntityHandle entity):
        """Get parallel status of an entity.

        Parameters
        ----------
        entity : EntityHandle
            The entity to query

        Returns
        -------
        int
            Parallel status value (PSTATUS constants)
        """
        cdef moab.ErrorCode err
        cdef unsigned char pstatus_val
        err = self.inst.get_pstatus(entity, pstatus_val)
        check_error(err)
        return pstatus_val

    def get_pstatus_entities(self, int dim, int pstatus_val):
        """Get entities with specific parallel status.

        Parameters
        ----------
        dim : int
            Dimension of entities to query
        pstatus_val : int
            Parallel status value to filter by

        Returns
        -------
        Range
            Entities matching the parallel status
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_pstatus_entities(dim, <unsigned char>pstatus_val, deref(ents.inst))
        check_error(err)
        return ents

    def filter_pstatus(self, Range entities, int pstatus_val, int op, int to_proc=-1):
        """Filter entities by parallel status in-place.

        Parameters
        ----------
        entities : Range
            Entities to filter (modified in place)
        pstatus_val : int
            Parallel status value to filter by (PSTATUS_* constants)
        op : int
            Filter operation (PSTATUS_AND, PSTATUS_OR, PSTATUS_NOT)
        to_proc : int, optional
            Target processor (-1 for all)
        """
        cdef moab.ErrorCode err
        err = self.inst.filter_pstatus(deref(entities.inst), <unsigned char>pstatus_val,
                                        <unsigned char>op, to_proc, NULL)
        check_error(err)

    def get_owner(self, EntityHandle entity):
        """Get the owning processor rank for an entity.

        Parameters
        ----------
        entity : EntityHandle
            The entity to query

        Returns
        -------
        int
            The owning processor rank
        """
        cdef moab.ErrorCode err
        cdef int owner
        err = self.inst.get_owner(entity, owner)
        check_error(err)
        return owner

    def get_owner_handle(self, EntityHandle entity):
        """Get owning processor rank and remote handle for an entity.

        Parameters
        ----------
        entity : EntityHandle
            The entity to query

        Returns
        -------
        tuple
            (owner_rank, remote_handle)
        """
        cdef moab.ErrorCode err
        cdef int owner
        cdef EntityHandle remote_handle
        err = self.inst.get_owner_handle(entity, owner, remote_handle)
        check_error(err)
        return owner, remote_handle

    def get_sharing_data(self, EntityHandle entity):
        """Get all sharing processors and handles for an entity.

        Parameters
        ----------
        entity : EntityHandle
            The entity to query

        Returns
        -------
        dict
            Dictionary with 'procs', 'handles', 'pstatus', 'num_procs'
        """
        cdef moab.ErrorCode err
        cdef int procs[64]
        cdef EntityHandle handles[64]
        cdef unsigned char pstat
        cdef unsigned int num_ps
        err = self.inst.get_sharing_data(entity, procs, handles, pstat, num_ps)
        check_error(err)
        return {
            'procs': [procs[i] for i in range(num_ps)],
            'handles': [handles[i] for i in range(num_ps)],
            'pstatus': pstat,
            'num_procs': num_ps
        }

    def get_interface_procs(self):
        """Get all interface processor ranks.

        Returns
        -------
        list
            List of interface processor ranks
        """
        cdef moab.ErrorCode err
        cdef cpp_set[unsigned int] procs
        err = self.inst.get_interface_procs(procs, False)
        check_error(err)
        return list(procs)

    def get_comm_procs(self):
        """Get all processors in the communicator.

        Returns
        -------
        list
            List of all processor ranks
        """
        cdef moab.ErrorCode err
        cdef cpp_set[unsigned int] procs
        err = self.inst.get_comm_procs(procs)
        check_error(err)
        return list(procs)

    def get_owned_entities(self, int dim=-1, EntityHandle file_set=0):
        """Get entities owned by this process.

        Parameters
        ----------
        dim : int, optional
            Dimension of entities to get (-1 for all dimensions)
        file_set : EntityHandle, optional
            Handle to the file set to get entities from

        Returns
        -------
        Range
            Range containing owned entities
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()

        if dim == -1:
            err = self.core.get_entities_by_handle(file_set, deref(ents.inst), False)
        else:
            err = self.core.get_entities_by_dimension(file_set, dim, deref(ents.inst), False)
        check_error(err)

        # Filter to get only owned entities (remove NOT_OWNED)
        err = self.inst.filter_pstatus(deref(ents.inst), PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, NULL)
        check_error(err)

        return ents

    def get_ghost_entities(self, int dim=-1, EntityHandle file_set=0):
        """Get ghost entities on this process.

        Parameters
        ----------
        dim : int, optional
            Dimension of entities to get (-1 for all dimensions)
        file_set : EntityHandle, optional
            Handle to the file set to get entities from

        Returns
        -------
        Range
            Range containing ghost entities
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        cdef Range ownedents = Range()

        if dim == -1:
            err = self.core.get_entities_by_handle(file_set, deref(ents.inst), False)
        else:
            err = self.core.get_entities_by_dimension(file_set, dim, deref(ents.inst), False)
        check_error(err)

        ownedents.inst.merge(deref(ents.inst))
        err = self.inst.filter_pstatus(deref(ownedents.inst), PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, NULL)
        check_error(err)

        return subtract(ents, ownedents)

    def assign_global_ids(self, EntityHandle this_set=0, int dimension=3,
                          int start_id=1, bint largest_dim_only=True,
                          bint parallel=True, bint owned_only=False):
        """Assign global IDs to mesh entities.

        Parameters
        ----------
        this_set : EntityHandle, optional
            Set containing entities
        dimension : int, optional
            Dimension of entities to assign IDs
        start_id : int, optional
            Starting ID value
        largest_dim_only : bool, optional
            Only assign to highest dimension entities
        parallel : bool, optional
            Do in parallel
        owned_only : bool, optional
            Only assign to owned entities
        """
        cdef moab.ErrorCode err
        err = self.inst.assign_global_ids(this_set, dimension, start_id,
                                          largest_dim_only, parallel, owned_only)
        check_error(err)

    def check_global_ids(self, EntityHandle this_set=0, int dimension=3,
                         int start_id=1, bint largest_dim_only=True,
                         bint parallel=True, bint owned_only=False):
        """Check and create global IDs if missing.

        Parameters
        ----------
        this_set : EntityHandle, optional
            Set containing entities
        dimension : int, optional
            Dimension of entities to check
        start_id : int, optional
            Starting ID value
        largest_dim_only : bool, optional
            Only check highest dimension entities
        parallel : bool, optional
            Do in parallel
        owned_only : bool, optional
            Only check owned entities
        """
        cdef moab.ErrorCode err
        err = self.inst.check_global_ids(this_set, dimension, start_id,
                                         largest_dim_only, parallel, owned_only)
        check_error(err)

    def exchange_tags(self, list src_tag_names, list dst_tag_names, Range entities):
        """Exchange tag values between processes for shared entities.

        Parameters
        ----------
        src_tag_names : list of str
            Source tag names
        dst_tag_names : list of str
            Destination tag names
        entities : Range
            Entities to exchange tags for
        """
        cdef moab.ErrorCode err
        cdef vector[moab.Tag] src_tags
        cdef vector[moab.Tag] dst_tags
        cdef moab.Tag tag_handle

        for tag_name in src_tag_names:
            err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_handle)
            check_error(err)
            src_tags.push_back(tag_handle)

        for tag_name in dst_tag_names:
            err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_handle)
            check_error(err)
            dst_tags.push_back(tag_handle)

        err = self.inst.exchange_tags(src_tags, dst_tags, deref(entities.inst))
        check_error(err)

    def reduce_tags(self, list src_tag_names, list dst_tag_names,
                    str mpi_op_name, Range entities):
        """Reduce tags using MPI operations across shared entities.

        Parameters
        ----------
        src_tag_names : list of str
            Source tag names
        dst_tag_names : list of str
            Destination tag names
        mpi_op_name : str
            MPI operation name ('SUM', 'MAX', 'MIN', 'PROD', 'LAND', 'LOR', etc.)
        entities : Range
            Entities to reduce tags for

        Raises
        ------
        ValueError
            If mpi_op_name is not recognized
        """
        cdef moab.ErrorCode err
        cdef vector[moab.Tag] src_tags
        cdef vector[moab.Tag] dst_tags
        cdef moab.Tag tag_handle
        cdef MPI_Op mpi_op

        op_name = mpi_op_name.upper()
        if op_name == 'SUM':
            mpi_op = MPI_SUM
        elif op_name == 'MAX':
            mpi_op = MPI_MAX
        elif op_name == 'MIN':
            mpi_op = MPI_MIN
        elif op_name == 'PROD':
            mpi_op = MPI_PROD
        elif op_name == 'LAND':
            mpi_op = MPI_LAND
        elif op_name == 'LOR':
            mpi_op = MPI_LOR
        elif op_name == 'BAND':
            mpi_op = MPI_BAND
        elif op_name == 'BOR':
            mpi_op = MPI_BOR
        elif op_name == 'BXOR':
            mpi_op = MPI_BXOR
        else:
            raise ValueError(f"Unknown MPI operation: {mpi_op_name}. "
                             f"Valid operations: SUM, MAX, MIN, PROD, LAND, LOR, BAND, BOR, BXOR")

        for tag_name in src_tag_names:
            err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_handle)
            check_error(err)
            src_tags.push_back(tag_handle)

        for tag_name in dst_tag_names:
            err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_handle)
            check_error(err)
            dst_tags.push_back(tag_handle)

        err = self.inst.reduce_tags(src_tags, dst_tags, mpi_op, deref(entities.inst))
        check_error(err)

    def get_partitioning(self):
        """Get the partitioning entity set.

        Returns
        -------
        EntityHandle
            The partitioning set handle
        """
        return self.inst.get_partitioning()

    def set_partitioning(self, EntityHandle h):
        """Set the partitioning entity set.

        Parameters
        ----------
        h : EntityHandle
            The partitioning set handle
        """
        cdef moab.ErrorCode err
        err = self.inst.set_partitioning(h)
        check_error(err)

    def get_part_entities(self, int dim=-1):
        """Get entities in the local partition.

        Parameters
        ----------
        dim : int, optional
            Dimension of entities (-1 for all)

        Returns
        -------
        Range
            Entities in the partition
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_part_entities(deref(ents.inst), dim)
        check_error(err)
        return ents

    def get_global_part_count(self):
        """Get total number of partitions.

        Returns
        -------
        int
            Total number of partitions
        """
        cdef moab.ErrorCode err
        cdef int count
        err = self.inst.get_global_part_count(count)
        check_error(err)
        return count

    def get_part_owner(self, int part_id):
        """Get the owning processor for a partition.

        Parameters
        ----------
        part_id : int
            Partition ID

        Returns
        -------
        int
            Owning processor rank
        """
        cdef moab.ErrorCode err
        cdef int owner
        err = self.inst.get_part_owner(part_id, owner)
        check_error(err)
        return owner

    def get_part_id(self, EntityHandle part):
        """Get the ID for a partition handle.

        Parameters
        ----------
        part : EntityHandle
            Partition handle

        Returns
        -------
        int
            Partition ID
        """
        cdef moab.ErrorCode err
        cdef int part_id
        err = self.inst.get_part_id(part, part_id)
        check_error(err)
        return part_id

    def create_part(self):
        """Create a new partition.

        Returns
        -------
        EntityHandle
            Handle to the new partition
        """
        cdef moab.ErrorCode err
        cdef EntityHandle part
        err = self.inst.create_part(part)
        check_error(err)
        return part

    def destroy_part(self, EntityHandle part):
        """Destroy a partition.

        Parameters
        ----------
        part : EntityHandle
            Handle to the partition to destroy
        """
        cdef moab.ErrorCode err
        err = self.inst.destroy_part(part)
        check_error(err)

    def broadcast_entities(self, int from_proc, Range entities,
                           bint adjacencies=False, bint tags=True):
        """Broadcast entities from one processor to all others.

        Parameters
        ----------
        from_proc : int
            Source processor rank
        entities : Range
            Entities to broadcast (modified in place)
        adjacencies : bool, optional
            Include adjacencies
        tags : bool, optional
            Include tags
        """
        cdef moab.ErrorCode err
        err = self.inst.broadcast_entities(from_proc, deref(entities.inst),
                                           adjacencies, tags)
        check_error(err)

    def gather_data(self, Range gather_ents, str tag_name,
                    str id_tag_name=None, int root_proc_rank=0):
        """Gather tag data to root processor.

        Parameters
        ----------
        gather_ents : Range
            Entities whose tag data to gather
        tag_name : str
            Name of the tag to gather
        id_tag_name : str, optional
            Name of the ID tag
        root_proc_rank : int, optional
            Root processor rank
        """
        cdef moab.ErrorCode err
        cdef moab.Tag tag_handle
        cdef moab.Tag id_tag = NULL

        err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_handle)
        check_error(err)

        if id_tag_name is not None:
            err = self.core.tag_get_handle(id_tag_name.encode('utf-8'), id_tag)
            check_error(err)

        err = self.inst.gather_data(deref(gather_ents.inst), tag_handle, id_tag, 0, root_proc_rank)
        check_error(err)

    def delete_entities(self, Range to_delete):
        """Delete entities across all processors.

        Parameters
        ----------
        to_delete : Range
            Entities to delete
        """
        cdef moab.ErrorCode err
        err = self.inst.delete_entities(deref(to_delete.inst))
        check_error(err)

    def check_all_shared_handles(self, bint print_em=False):
        """Verify consistency of shared entity handles across processors.

        Parameters
        ----------
        print_em : bool, optional
            If True, print shared handle information
        """
        cdef moab.ErrorCode err
        err = self.inst.check_all_shared_handles(print_em)
        check_error(err)

    def set_debug_verbosity(self, int verbosity):
        """Set debug verbosity level.

        Parameters
        ----------
        verbosity : int
            Verbosity level (0 = silent, higher = more output)
        """
        self.inst.set_debug_verbosity(verbosity)

    def get_debug_verbosity(self):
        """Get debug verbosity level.

        Returns
        -------
        int
            Current verbosity level
        """
        return self.inst.get_debug_verbosity()

    def scatter_entities(self, int from_proc, list entities=None,
                         bint adjacencies=False, bint tags=True):
        """Scatter entities from one processor to all others.

        On the sending rank (from_proc), entities[i] is the Range
        to send to rank i.  On receiving ranks the list is populated
        with the received entities.

        Parameters
        ----------
        from_proc : int
            Source processor rank
        entities : list of Range, optional
            List of Range objects, one per processor. Required on from_proc.
        adjacencies : bool, optional
            Include adjacencies
        tags : bool, optional
            Include tags

        Returns
        -------
        list of Range
            One Range per processor with received entities
        """
        cdef moab.ErrorCode err
        cdef vector[moab.Range] c_entities
        cdef int i
        cdef Range r

        c_entities.resize(self.inst.size())
        if entities is not None:
            for i in range(min(len(entities), <int>self.inst.size())):
                r = <Range>entities[i]
                c_entities[i] = deref(r.inst)

        err = self.inst.scatter_entities(from_proc, c_entities,
                                          adjacencies, tags)
        check_error(err)

        result = []
        for i in range(<int>c_entities.size()):
            out = Range()
            out.inst.merge(c_entities[i])
            result.append(out)
        return result

    def resolve_shared_sets(self, EntityHandle this_set=0):
        """Resolve shared sets between processors.

        Matches entity sets across processes using global IDs and
        populates sharing data for sets.

        Parameters
        ----------
        this_set : EntityHandle, optional
            Set directly containing candidate sets (e.g. file set).
            Use 0 for root set.
        """
        cdef moab.ErrorCode err
        err = self.inst.resolve_shared_sets(this_set, NULL)
        check_error(err)

    def get_iface_entities(self, int other_proc, int dim=-1):
        """Get entities on interfaces shared with another processor.

        Parameters
        ----------
        other_proc : int
            Rank of the other processor sharing the interface
        dim : int, optional
            Dimension of entities to return (-1 for all dimensions)

        Returns
        -------
        Range
            Entities on the shared interface
        """
        cdef moab.ErrorCode err
        cdef Range iface_ents = Range()
        err = self.inst.get_iface_entities(other_proc, dim, deref(iface_ents.inst))
        check_error(err)
        return iface_ents
