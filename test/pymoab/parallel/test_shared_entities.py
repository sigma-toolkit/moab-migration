#!/usr/bin/env python3
"""
Test ownership and sharing query methods.
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


def test_get_pstatus():
    """Test getting parallel status."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    v = verts[0]
    pstatus = pcomm.get_pstatus(v)
    CHECK(pstatus is not None)


def test_get_pstatus_entities():
    """Test getting entities by parallel status."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    ents = pcomm.get_pstatus_entities(0, 0)
    CHECK(ents is not None)


def test_get_owner():
    """Test getting entity owner."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    v = verts[0]
    owner = pcomm.get_owner(v)
    CHECK(owner >= 0)


def test_get_owner_handle():
    """Test getting owner handle."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    v = verts[0]
    owner, handle = pcomm.get_owner_handle(v)
    CHECK(owner >= 0)


def test_get_sharing_data():
    """Test getting sharing data for entity."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    
    v = verts[0]
    data = pcomm.get_sharing_data(v)
    CHECK(data is not None)
    CHECK('procs' in data)
    CHECK('handles' in data)
    CHECK('num_procs' in data)


def test_get_interface_procs():
    """Test getting interface processors."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    procs = pcomm.get_interface_procs()
    CHECK(procs is not None)


def main():
    tests = [
        test_capability_detection,
        test_get_pstatus,
        test_get_pstatus_entities,
        test_get_owner,
        test_get_owner_handle,
        test_get_sharing_data,
        test_get_interface_procs,
    ]
    return run_parallel_tests(tests)

if __name__ == "__main__":
    sys.exit(main())
