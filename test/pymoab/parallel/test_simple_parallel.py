#!/usr/bin/env python3

import sys
import numpy as np
from mpi4py import MPI

try:
    from pymoab import config, core, parallelcomm
except ImportError as e:
    print(f"✗ Failed to import PyMOAB modules: {e}")
    sys.exit(1)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def test_capability_detection():
    try:
        if config.MOAB_MPI_ENABLED:
            print("Python: Basic MPI support detected")
            if config.has_mpi_io():
                print("Python: MPI I/O support detected")
            else:
                print("Python: MPI I/O support not detected")
        else:
            print("Python: MPI support not detected")

        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, comm)

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
    try:
        mb = core.Core()
        print(f"✓ Created MOAB Core for rank {rank}")

        pcomm = parallelcomm.ParallelComm(mb, comm)
        print(f"✓ Created ParallelComm for rank {rank}")

        print(f"✓ ParallelComm properties: rank {pcomm.rank}, size {pcomm.size}")

        if config.MOAB_MPI_ENABLED:
            print(f"✓ Runtime basic MPI support: {pcomm.has_basic_mpi}")
            if config.has_mpi_io():
                print(f"✓ Runtime MPI I/O support: {pcomm.has_mpi_io}")
            print(f"✓ Runtime full MPI support: {pcomm.has_full_mpi}")

        return pcomm.rank == rank and pcomm.size == size
    except Exception as e:
        print(f"✗ ParallelComm creation failed: {e}")
        return False


def test_simple_mesh_creation():
    try:
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, comm)
        print(f"✓ Created ParallelComm for rank {rank}")

        coords = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ], dtype="float64")
        coords[:, 0] += rank * 2.0

        verts = mb.create_vertices(coords.flatten())
        print(f"✓ Created {len(verts)} vertices for rank {rank}")

        hexes = mb.create_elements(4, [verts])
        print(f"✓ Created {len(hexes)} hex elements for rank {rank}")
        return len(verts) == 8 and len(hexes) == 1
    except Exception as e:
        print(f"✗ Mesh creation failed: {e}")
        return False


def main():
    print(f"\n=== PyMOAB Parallel Interface Test ===")
    print(f"Process {rank}/{size}")

    print("\n--- Test 1: Capability Detection ---")
    test1_passed = test_capability_detection()

    print("\n--- Test 2: ParallelComm Creation ---")
    test2_passed = test_parallel_comm_creation()

    print("\n--- Test 3: Simple Mesh Creation ---")
    test3_passed = test_simple_mesh_creation()

    comm.Barrier()
    if rank == 0:
        print("\n=== Test Summary ===")
        print(f"Test 1 (Capability Detection): {'✓ PASSED' if test1_passed else '✗ FAILED'}")
        print(f"Test 2 (ParallelComm Creation): {'✓ PASSED' if test2_passed else '✗ FAILED'}")
        print(f"Test 3 (Mesh Creation): {'✓ PASSED' if test3_passed else '✗ FAILED'}")
        print(f"\nOverall Result: {'✓ ALL TESTS PASSED' if test1_passed and test2_passed and test3_passed else '✗ SOME TESTS FAILED'}")

    return 0 if test1_passed and test2_passed and test3_passed else 1


if __name__ == "__main__":
    sys.exit(main())
