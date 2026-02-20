#!/usr/bin/env python3
"""
Simple test script to verify PyMOAB parallel interface works with the new MPI pattern.
Run with: mpiexec -n 2 python test_simple_parallel.py
"""

import os
import sys
from mpi4py import MPI
import numpy as np

# Add the pymoab directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pymoab'))

try:
    from pymoab import core, parallelcomm
    print("✓ Successfully imported PyMOAB modules")
except ImportError as e:
    print(f"✗ Failed to import PyMOAB modules: {e}")
    sys.exit(1)

def test_capability_detection():
    """Test MPI capability detection."""
    try:
        # Test Python-level capability detection
        from pymoab import config
        if config.MOAB_MPI_ENABLED:
            print(f"Python: Basic MPI support detected")
            if config.has_mpi_io():
                print(f"Python: MPI I/O support detected")
            else:
                print(f"Python: MPI I/O support not detected")
        else:
            print(f"Python: MPI support not detected")
        
        # Test runtime capability detection
        from pymoab import core, parallelcomm
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(MPI.COMM_WORLD)
        
        if config.MOAB_MPI_ENABLED:
            print(f"Python: Runtime basic MPI support: {pcomm.has_basic_mpi}")
            if config.has_mpi_io():
                print(f"Python: Runtime MPI I/O support: {pcomm.has_mpi_io}")
            print(f"Python: Runtime full MPI support: {pcomm.has_full_mpi}")
        
        return True
    except Exception as e:
        print(f"✗ Capability detection failed: {e}")
        return False

def test_parallel_comm_creation():
    """Test ParallelComm creation."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    try:
        # Create MOAB Core
        mb = core.Core()
        print(f"✓ Created MOAB Core for rank {rank}")

        # Create ParallelComm
        pcomm = parallelcomm.ParallelComm(mb, comm)
        print(f"✓ Created ParallelComm for rank {rank}")

        # Test basic properties
        pcomm_rank = pcomm.rank
        pcomm_size = pcomm.size
        print(f"✓ ParallelComm properties: rank {pcomm_rank}, size {pcomm_size}")

        # Test capability properties
        if config.MOAB_MPI_ENABLED:
            print(f"✓ Runtime basic MPI support: {pcomm.has_basic_mpi}")
            if config.has_mpi_io():
                print(f"✓ Runtime MPI I/O support: {pcomm.has_mpi_io}")
            print(f"✓ Runtime full MPI support: {pcomm.has_full_mpi}")

        return True
    except Exception as e:
        print(f"✗ ParallelComm creation failed: {e}")
        return False

def test_simple_mesh_creation():
    """Test simple mesh creation in parallel."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    try:
        # Create MOAB Core and ParallelComm
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, comm)
        print(f"✓ Created ParallelComm for rank {rank}")

        # Test capability properties
        if config.MOAB_MPI_ENABLED:
            print(f"✓ Runtime basic MPI support: {pcomm.has_basic_mpi}")
            if config.has_mpi_io():
                print(f"✓ Runtime MPI I/O support: {pcomm.has_mpi_io}")
            print(f"✓ Runtime full MPI support: {pcomm.has_full_mpi}")

        # Create simple mesh data
        coords = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
        ], dtype='float64')

        # Offset by rank
        coords[:, 0] += rank * 2.0

        # Create vertices
        verts = mb.create_vertices(coords.flatten())
        print(f"✓ Created {len(verts)} vertices for rank {rank}")

        # Create hex elements
        hexes = mb.create_elements(4, [verts])  # 4 = MBHEX
        print(f"✓ Created {len(hexes)} hex elements for rank {rank}")

        return True
    except Exception as e:
        print(f"✗ Mesh creation failed: {e}")
        return False

def main():
    """Main test function."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"\n=== PyMOAB Parallel Interface Test ===")
    print(f"Process {rank}/{size}")

    # Test 1: Capability detection
    print(f"\n--- Test 1: Capability Detection ---")
    test1_passed = test_capability_detection()

    # Test 2: ParallelComm creation
    print(f"\n--- Test 2: ParallelComm Creation ---")
    test2_passed = test_parallel_comm_creation()

    # Test 3: Simple mesh creation
    print(f"\n--- Test 3: Simple Mesh Creation ---")
    test3_passed = test_simple_mesh_creation()

    # Summary
    comm.Barrier()
    if rank == 0:
        print(f"\n=== Test Summary ===")
        print(f"Test 1 (Capability Detection): {'✓ PASSED' if test1_passed else '✗ FAILED'}")
        print(f"Test 2 (ParallelComm Creation): {'✓ PASSED' if test2_passed else '✗ FAILED'}")
        print(f"Test 3 (Mesh Creation): {'✓ PASSED' if test3_passed else '✗ FAILED'}")

        all_passed = test1_passed and test2_passed and test3_passed
        print(f"\nOverall Result: {'✓ ALL TESTS PASSED' if all_passed else '✗ SOME TESTS FAILED'}")

if __name__ == "__main__":
    main()

def test_simple_mesh_creation():
    """Test simple mesh creation in parallel."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    try:
        # Create MOAB Core and ParallelComm
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, comm)

        # Create simple mesh data
        coords = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
        ], dtype='float64')

        # Offset by rank
        coords[:, 0] += rank * 2.0

        # Create vertices
        verts = mb.create_vertices(coords.flatten())
        print(f"✓ Created {len(verts)} vertices for rank {rank}")

        # Create hex elements
        hexes = mb.create_elements(4, [verts])  # 4 = MBHEX
        print(f"✓ Created {len(hexes)} hex elements for rank {rank}")

        return True
    except Exception as e:
        print(f"✗ Mesh creation failed: {e}")
        return False

def main():
    """Main test function."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"\n=== PyMOAB Parallel Interface Test ===")
    print(f"Process {rank}/{size}")

    # Test 1: MPI Comm passing
    print(f"\n--- Test 1: MPI Comm Passing ---")
    test1_passed = test_mpi_comm_passing()

    # Test 2: ParallelComm creation
    print(f"\n--- Test 2: ParallelComm Creation ---")
    test2_passed = test_parallel_comm_creation()

    # Test 3: Simple mesh creation
    print(f"\n--- Test 3: Simple Mesh Creation ---")
    test3_passed = test_simple_mesh_creation()

    # Summary
    comm.Barrier()
    if rank == 0:
        print(f"\n=== Test Summary ===")
        print(f"Test 1 (MPI Comm Passing): {'✓ PASSED' if test1_passed else '✗ FAILED'}")
        print(f"Test 2 (ParallelComm Creation): {'✓ PASSED' if test2_passed else '✗ FAILED'}")
        print(f"Test 3 (Mesh Creation): {'✓ PASSED' if test3_passed else '✗ FAILED'}")

        all_passed = test1_passed and test2_passed and test3_passed
        print(f"\nOverall Result: {'✓ ALL TESTS PASSED' if all_passed else '✗ SOME TESTS FAILED'}")

if __name__ == "__main__":
    main()