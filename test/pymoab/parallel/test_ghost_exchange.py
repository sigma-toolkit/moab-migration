#!/usr/bin/env python3
"""
Test parallel I/O, ghost exchange, and shared entity operations.

Exercises the full parallel workflow:
  1. Load a partitioned mesh with PARALLEL=READ_PART
  2. Query owned and ghost entities
  3. Resolve shared entities
  4. Exchange ghost cells
  5. Verify entity counts across ranks
"""

import os
import sys
import numpy as np

try:
    from mpi4py import MPI
except ImportError:
    print("MPI not available - skipping parallel tests")
    sys.exit(0)

from pymoab import config, core, parallelcomm, types
from pymoab.rng import Range
from parallel_driver import (run_parallel_tests, CHECK_EQ, CHECK,
                             CHECK_PARALLEL, get_all_values)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MESH_DIR = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'MeshFiles', 'unittest')
PARTITIONED_MESH = os.path.join(MESH_DIR, "64bricks_512hex_256part.h5m")
READ_OPTS = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS"


def _have_parallel_io():
    return config.MOAB_MPI_ENABLED and config.has_mpi_io()


def _have_test_mesh():
    return os.path.exists(PARTITIONED_MESH)


def test_parallel_load():
    """Load a partitioned mesh in parallel and verify each rank got entities."""
    if not _have_parallel_io():
        print(f"Process {rank}: Skipping - MPI I/O not available")
        CHECK(True)
        return
    if not _have_test_mesh():
        print(f"Process {rank}: Skipping - test mesh not found: {PARTITIONED_MESH}")
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    local_hexes = mb.get_entities_by_type(0, types.MBHEX)
    local_verts = mb.get_entities_by_type(0, types.MBVERTEX)

    print(f"Process {rank}: loaded {len(local_hexes)} hexes, {len(local_verts)} vertices")

    CHECK(len(local_hexes) > 0)
    CHECK(len(local_verts) > 0)

    all_hex_counts = get_all_values(len(local_hexes))
    total_hexes = sum(all_hex_counts)
    if rank == 0:
        print(f"Total hexes across {size} ranks: {total_hexes}")
    CHECK(total_hexes >= 512)


def test_owned_vs_ghost_after_load():
    """After parallel load, owned + ghost = total local entities."""
    if not (_have_parallel_io() and _have_test_mesh()):
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    owned = pcomm.get_owned_entities(dim=3)
    ghost = pcomm.get_ghost_entities(dim=3)
    all_3d = mb.get_entities_by_dimension(0, 3)

    print(f"Process {rank}: owned={len(owned)}, ghost={len(ghost)}, total={len(all_3d)}")

    CHECK(len(owned) > 0)
    CHECK_EQ(len(owned) + len(ghost), len(all_3d))


def test_ghost_exchange_after_load():
    """Exchange ghost cells and verify ghost count increases or stays consistent."""
    if not (_have_parallel_io() and _have_test_mesh()):
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    ghost_before = pcomm.get_ghost_entities(dim=3)
    count_before = len(ghost_before)

    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    ghost_after = pcomm.get_ghost_entities(dim=3)
    count_after = len(ghost_after)

    print(f"Process {rank}: ghosts before={count_before}, after={count_after}")
    CHECK(count_after >= count_before)


def test_shared_entities_after_load():
    """After parallel load, shared entity queries should work."""
    if not (_have_parallel_io() and _have_test_mesh()):
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    shared = pcomm.get_shared_entities(-1, dim=0)
    print(f"Process {rank}: {len(shared)} shared vertices")
    CHECK(shared is not None)

    if size > 1:
        comm_procs = pcomm.get_comm_procs()
        print(f"Process {rank}: communicating with procs {comm_procs}")
        CHECK(isinstance(comm_procs, list))


def test_assign_global_ids():
    """Assign global IDs to a locally-created mesh in parallel."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    coords = np.array([
        0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,
        0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
    ], dtype="float64")
    coords[0::3] += rank * 2.0

    verts = mb.create_vertices(coords)
    hexes = mb.create_element(types.MBHEX, verts)

    pcomm.assign_global_ids(dimension=3, start_id=1)

    global_id_tag = mb.tag_get_handle("GLOBAL_ID")
    gids = mb.tag_get_data(global_id_tag, mb.get_entities_by_type(0, types.MBHEX))
    print(f"Process {rank}: assigned global IDs: {gids}")
    CHECK(len(gids) > 0)


def test_get_owned_entities_local():
    """Locally created entities should be owned by this rank."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)

    owned = pcomm.get_owned_entities(dim=0)
    CHECK(owned is not None)
    CHECK_EQ(len(owned), 1)


def test_get_ghost_entities_local():
    """Locally created mesh should have no ghost entities."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    coords = np.array([0.0, 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)

    ghost = pcomm.get_ghost_entities(dim=0)
    CHECK(ghost is not None)
    CHECK_EQ(len(ghost), 0)


def test_parallel_write(tmp_path=None):
    """Write a mesh in parallel and verify the file was created."""
    if not _have_parallel_io():
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    coords = np.array([
        0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,
        0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
    ], dtype="float64")
    coords[0::3] += rank * 2.0

    verts = mb.create_vertices(coords)
    mb.create_element(types.MBHEX, verts)

    import tempfile
    outfile = os.path.join(tempfile.gettempdir(), "pymoab_par_write_test.h5m")
    pcomm.write_file(outfile, "PARALLEL=WRITE_PART")

    comm.Barrier()
    if rank == 0:
        CHECK(os.path.exists(outfile))
        fsize = os.path.getsize(outfile)
        print(f"Written parallel file: {outfile} ({fsize} bytes)")
        CHECK(fsize > 0)
        os.unlink(outfile)
    else:
        CHECK(True)


def main():
    tests = [
        test_parallel_load,
        test_owned_vs_ghost_after_load,
        test_ghost_exchange_after_load,
        test_shared_entities_after_load,
        test_assign_global_ids,
        test_get_owned_entities_local,
        test_get_ghost_entities_local,
        test_parallel_write,
    ]
    return run_parallel_tests(tests)


if __name__ == "__main__":
    sys.exit(main())
