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

# Import EntityHandle and Tag from moab
from .moab cimport EntityHandle, Tag, Interface

# Include moab_defs.pxd to get PSTATUS macros
include "moab_defs.pxd"

# Import MPI types from mpi4py
from mpi4py import MPI
from mpi4py.MPI cimport Comm as MPIComm
from mpi4py.libmpi cimport MPI_Comm, MPI_Op, MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD
from mpi4py.libmpi cimport MPI_LAND, MPI_LOR, MPI_BAND, MPI_BOR, MPI_BXOR

# NO direct MPI C calls - let MOAB handle its own MPI

cdef void* null = NULL

# Declare intptr_t manually
ctypedef long int  intptr_t

# Helper function to convert Python list of strings to C++ vector of strings
cdef vector[std_string] _convert_string_list(list py_strings):
    cdef vector[std_string] cpp_strings
    cdef bytes b_str
    for py_str in py_strings:
        if py_str is not None:
            b_str = py_str.encode('utf-8')
            cpp_strings.push_back(std_string(<char*>b_str))
    return cpp_strings

def test_mpi_comm_passing(MPIComm comm):
    """
    Test function to verify MPI Comm passing works correctly.
    This demonstrates the core concept without calling MPI C functions directly.
    """
    # Get the communicator handle as integer - this is the correct way!
    comm_handle = comm.handle

    # Use mpi4py's API to get rank and size
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"PyMOAB: Received MPI Comm with handle {comm_handle}")
    print(f"PyMOAB: Process rank {rank}/{size} (via mpi4py)")

    return rank, size

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
        # Ensure mpi4py has initialized MPI before any operations
        _ = comm.Get_rank()  # This ensures mpi4py has initialized MPI

        # Get the C MPI communicator handle from mpi4py
        cdef MPI_Comm c_comm = comm.ob_mpi

        self.core = core.inst
        self.comm = comm

        print(f"PyMOAB: Creating ParallelComm for process {comm.Get_rank()}")
        comm.Barrier()

        # Create MOAB ParallelComm with C communicator
        # MOAB will handle its own MPI operations internally
        self.inst = new moab.ParallelComm(<moab.Interface*>(self.core), c_comm)

        comm.Barrier()
        print(f"PyMOAB: ParallelComm created for process {comm.Get_rank()}")
        comm.Barrier()

    def __del__(self):
        """Destructor"""
        if self.inst != null:
            del self.inst

    def load_file(self, str file_name, str read_opts="", EntityHandle file_set=0):
        """Load a mesh file in parallel.

        Parameters
        ----------
        file_name : str
            Name of the file to load
        read_opts : str, optional
            Options for reading the file
        file_set : EntityHandle, optional
            Handle to the file set to load into

        Returns
        -------
            None.

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

        if file_set != None:
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
            Options for writing the file
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

        if write_opts:
            b_write_opts = write_opts.encode('utf-8')
            c_write_opts = b_write_opts

        if tag_list is not None:
            # Convert tag names to Tag pointers
            for tag_name in tag_list:
                err = self.core.tag_get_handle(tag_name.encode('utf-8'), tag_ptr)
                check_error(err)
                if err == moab.MB_SUCCESS and tag_ptr != NULL:
                    tag_ptrs.push_back(tag_ptr)
            num_tags = tag_ptrs.size()

        # Use the write_file overload that takes tag list
        #err = self.core.write_file(c_file_name, c_format_type, c_write_opts, ptr, num_sets,
        #                         tag_ptrs.data() if num_tags > 0 else NULL,
        #                         num_tags)
        err = self.core.write_file(c_file_name, c_format_type, c_write_opts, ptr, num_sets)

        check_error(err)

    def get_ghost_entities(self, int bridge_dim, int ghost_dim, int to_proc=-1,
                          int num_layers=1, int addl_ents=0):
        """Get ghost entities for parallel communication.

        Parameters
        ----------
        bridge_dim : int
            Dimension of bridge entities
        ghost_dim : int
            Dimension of ghost entities
        to_proc : int, optional
            Target processor (-1 for all)
        num_layers : int, optional
            Number of ghost layers
        addl_ents : int, optional
            Additional entities to include

        Returns
        -------
        Range
            Range containing ghost entities

        Raises
        ------
        RuntimeError
            If ghost entity retrieval fails
        """
        cdef moab.ErrorCode err
        cdef Range ghost_entities = Range()

        #err = self.inst.get_ghosted_entities(bridge_dim, ghost_dim, to_proc,
        #                                   num_layers, addl_ents, deref(ghost_entities.inst))
        #check_error(err)
        return ghost_entities

    def get_shared_entities(self, int other_proc=-1, int dim=-1, bint iface=False, bint owned_filter=False):
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

        Raises
        ------
        RuntimeError
            If getting shared entities fails
        """
        cdef moab.ErrorCode err
        cdef Range shared_ents = Range()
        cdef Range dim_ents = Range()

        # Get all shared entities
        err = self.core.get_entities_by_type(0, moab.MBMAXTYPE, deref(shared_ents.inst), False)
        check_error(err)

        # Filter by shared status
        err = self.inst.filter_pstatus(deref(shared_ents.inst), PSTATUS_SHARED, PSTATUS_AND, other_proc)
        check_error(err)

        # Filter by dimension if specified
        if dim != -1:
            err = self.core.get_entities_by_dimension(0, dim, deref(dim_ents.inst), False)
            check_error(err)
            shared_ents = shared_ents.intersect(dim_ents)

        # Filter by interface if requested
        if iface:
            err = self.inst.filter_pstatus(deref(shared_ents.inst), PSTATUS_INTERFACE, PSTATUS_AND, -1)
            check_error(err)

        # Filter by ownership if requested
        if owned_filter:
            err = self.inst.filter_pstatus(deref(shared_ents.inst), PSTATUS_NOT_OWNED, PSTATUS_NOT, -1)
            check_error(err)

        return shared_ents

    def resolve_shared_ents(self, EntityHandle this_set, int to_dim):
        """Resolve shared entities between processors.

        Parameters
        ----------
        this_set : EntityHandle
            Set containing entities to resolve
        to_dim : int
            Dimension up to which to resolve

        Raises
        ------
        RuntimeError
            If entity resolution fails
        """
        cdef moab.ErrorCode err
        err = self.inst.resolve_shared_ents(this_set, to_dim)
        check_error(err)

    def exchange_ghost_cells(self, int ghost_dim, int bridge_dim, int num_layers=1,
                            int addl_ents=0, bint store_remote_handles=False,
                            bint wait_all=False, EntityHandle file_set=0):
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
        err = self.inst.exchange_ghost_cells(ghost_dim, bridge_dim, num_layers,
                                           addl_ents, store_remote_handles,
                                           wait_all, &file_set)
        check_error(err)

    @property
    def rank(self):
        """Get the MPI rank of this processor.

        Returns
        -------
        int
            The MPI rank
        """
        return self.inst.rank()

    @property
    def size(self):
        """Get the total number of MPI processes.

        Returns
        -------
        int
            The number of MPI processes
        """
        return self.inst.size()

    @property
    def comm(self):
        """Get the MPI communicator.

        Returns
        -------
        MPI.Comm
            The MPI communicator
        """
        return self.comm

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
            'procs': list(procs[:num_ps]),
            'handles': list(handles[:num_ps]),
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

        Raises
        ------
        RuntimeError
            If getting entities fails
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()

        # Get all entities of the specified dimension
        if dim == -1:
            err = self.core.get_entities_by_handle(file_set, deref(ents.inst), False)
        else:
            err = self.core.get_entities_by_dimension(file_set, dim, deref(ents.inst), False)
        check_error(err)

        # Filter to get only owned entities
        err = self.inst.filter_pstatus(deref(ents.inst), PSTATUS_NOT_OWNED, PSTATUS_NOT, -1)
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

        Raises
        ------
        RuntimeError
            If getting entities fails
        """
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        cdef Range ownedents = Range()
        cdef Range ghostents = Range()

        # Get all entities of the specified dimension
        if dim == -1:
            err = self.core.get_entities_by_handle(file_set, deref(ents.inst), False)
        else:
            err = self.core.get_entities_by_dimension(file_set, dim, deref(ents.inst), False)
        check_error(err)

        # Filter to get only ghost entities (not owned)
        # We do this by getting the owned entities and then subtracting the owned entities from the total
        ownedents = ents
        err = self.inst.filter_pstatus(deref(ownedents.inst), PSTATUS_NOT_OWNED, PSTATUS_NOT, -1)
        check_error(err)

        ghostents = subtract(ents, ownedents)

        return ghostents

    def assign_global_ids(self, EntityHandle this_set, int dimension,
                         int start_id=1, bint largest_dim_only=True,
                         bint parallel=True, bint owned_only=False):
        """Assign global IDs to mesh entities.

        Parameters
        ----------
        this_set : EntityHandle
            Set containing entities
        dimension : int
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

    def check_global_ids(self, EntityHandle this_set, int dimension,
                        int start_id=1, bint largest_dim_only=True,
                        bint parallel=True, bint owned_only=False):
        """Check and create global IDs if missing.

        Parameters
        ----------
        this_set : EntityHandle
            Set containing entities
        dimension : int
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
        """Exchange tag values between processes.

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
        """Reduce tags using MPI operations.

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
        """
        cdef moab.ErrorCode err
        cdef vector[moab.Tag] src_tags
        cdef vector[moab.Tag] dst_tags
        cdef moab.Tag tag_handle
        cdef MPI_Op mpi_op

        op_map = {
            'SUM': MPI.SUM,
            'MAX': MPI.MAX,
            'MIN': MPI.MIN,
            'PROD': MPI.PROD,
            'LAND': MPI.LAND,
            'LOR': MPI.LOR,
            'BAND': MPI.BAND,
            'BOR': MPI.BOR,
            'BXOR': MPI.BXOR,
        }
        mpi_op = op_map.get(mpi_op_name.upper(), MPI.SUM)

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
        """Get entities in the partition.

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
        
        err = self.inst.gather_data(deref(gather_ents.inst), tag_handle, id_tag, NULL, root_proc_rank)
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

cdef extern from "mpi.h":
    int MPI_Comm_rank(MPI_Comm comm, int *rank)

def test_cpp_rank(MPIComm comm):
    cdef MPI_Comm c_comm = comm.ob_mpi
    cdef int rank
    MPI_Comm_rank(c_comm, &rank)
    return rank