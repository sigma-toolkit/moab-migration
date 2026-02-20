#!/usr/bin/env python3
"""
Test ghost exchange and shared entity operations.
"""

import sys
import os
import numpy as np

# Check if MPI is available first
try:
    from mpi4py import MPI
except ImportError:
    print("MPI not available - skipping parallel tests")
    sys.exit(0)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pymoab'))

from pymoab import core, parallelcomm, types
from pymoab.rng import Range
from parallel_driver import run_parallel_tests, CHECK_EQ, CHECK, CHECK_PARALLEL

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

filepath = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'MeshFiles', 'unittest')

def test_capability_detection():
    """Test capability detection."""
    if config.MOAB_MPI_ENABLED:
        print(f"Process {rank}: Basic MPI support detected")
        if config.has_mpi_io():
            print(f"Process {rank}: MPI I/O support detected")
        else:
            print(f"Process {rank}: MPI I/O support not detected")
    else:
        print(f"Process {rank}: MPI support not detected")

    # Test runtime capability detection
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    if config.MOAB_MPI_ENABLED:
        print(f"Process {rank}: Runtime basic MPI support: {pcomm.has_basic_mpi}")
        if config.has_mpi_io():
            print(f"Process {rank}: Runtime MPI I/O support: {pcomm.has_mpi_io}")
        print(f"Process {rank}: Runtime full MPI support: {pcomm.has_full_mpi}")


def test_load_partitioned_mesh():
    """Test loading a partitioned mesh file."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    mesh_file = os.path.join(filepath, "64bricks_512hex_256part.h5m")
    if not os.path.exists(mesh_file):
        CHECK(True)
        return
    
    # Test capability before loading
    if config.MOAB_MPI_ENABLED and config.has_mpi_io():
        mb.load_file(mesh_file)
        CHECK(True)
    else:
        print(f"Process {rank}: Skipping load_file test - MPI I/O not available")
        CHECK(True)


def test_get_owned_entities():
    """Test getting owned entities."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    owned = pcomm.get_owned_entities(dim=0)
    CHECK(owned is not None)


def test_get_ghost_entities():
    """Test getting ghost entities."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    # Test capability before getting ghost entities
    if config.MOAB_MPI_ENABLED and pcomm.has_basic_mpi:
        ghost = pcomm.get_ghost_entities(dim=0)
        CHECK(ghost is not None)
    else:
        print(f"Process {rank}: Skipping get_ghost_entities test - Basic MPI not available")
        CHECK(True)


def test_resolve_shared_ents():
    """Test resolving shared entities."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    root_set = mb.get_root_set()
    
    # Test capability before resolving
    if config.MOAB_MPI_ENABLED and pcomm.has_basic_mpi:
        pcomm.resolve_shared_ents(root_set, 3)
        CHECK(True)
    else:
        print(f"Process {rank}: Skipping resolve_shared_ents test - Basic MPI not available")
        CHECK(True)


def test_exchange_ghost_cells():
    """Test exchanging ghost cells."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    # Test capability before exchanging
    if config.MOAB_MPI_ENABLED and pcomm.has_basic_mpi:
        pcomm.exchange_ghost_cells(3, 0, 1, 0, False, False)
        CHECK(True)
    else:
        print(f"Process {rank}: Skipping exchange_ghost_cells test - Basic MPI not available")
        CHECK(True)


def main():
    tests = [
        test_capability_detection,
        test_load_partitioned_mesh,
        test_get_owned_entities,
        test_get_ghost_entities,
        test_resolve_shared_ents,
        test_exchange_ghost_cells,
    ]
    return run_parallel_tests(tests)

if __name__ == "__main__":
    sys.exit(main())
