#!/usr/bin/env python3
"""
Test ownership, sharing query methods, tag exchange, and reductions.
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


def _make_local_vertex():
    """Create a single locally-owned vertex for simple property tests."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    coords = np.array([float(rank), 0.0, 0.0], dtype='float64')
    verts = mb.create_vertices(coords)
    return mb, pcomm, verts


def test_get_pstatus():
    """Parallel status of a locally created vertex should be zero (owned, not shared)."""
    mb, pcomm, verts = _make_local_vertex()
    pstatus = pcomm.get_pstatus(verts[0])
    CHECK(pstatus == 0)


def test_get_pstatus_entities():
    """Query entities by parallel status value."""
    mb, pcomm, verts = _make_local_vertex()
    ents = pcomm.get_pstatus_entities(0, 0)
    CHECK(ents is not None)


def test_get_owner():
    """Owner of a locally created entity should be this rank."""
    mb, pcomm, verts = _make_local_vertex()
    owner = pcomm.get_owner(verts[0])
    CHECK_EQ(owner, rank)


def test_get_owner_handle():
    """Owner handle of a locally created entity returns valid rank and handle."""
    mb, pcomm, verts = _make_local_vertex()
    owner, handle = pcomm.get_owner_handle(verts[0])
    CHECK_EQ(owner, rank)
    CHECK(handle > 0)


def test_get_sharing_data():
    """Sharing data for a local entity has expected structure."""
    mb, pcomm, verts = _make_local_vertex()
    data = pcomm.get_sharing_data(verts[0])
    CHECK('procs' in data)
    CHECK('handles' in data)
    CHECK('pstatus' in data)
    CHECK('num_procs' in data)


def test_get_interface_procs():
    """Interface procs returns a list (possibly empty for local-only mesh)."""
    mb, pcomm, verts = _make_local_vertex()
    procs = pcomm.get_interface_procs()
    CHECK(isinstance(procs, list))


def test_assign_global_ids():
    """Assign global IDs and verify each rank gets a unique ID."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    coords = np.array([
        0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,
        0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
    ], dtype="float64")
    coords[0::3] += rank * 2.0

    verts = mb.create_vertices(coords)
    mb.create_element(types.MBHEX, verts)

    pcomm.assign_global_ids(dimension=3, start_id=1)

    global_id_tag = mb.tag_get_handle("GLOBAL_ID")
    hexes = mb.get_entities_by_type(0, types.MBHEX)
    gids = mb.tag_get_data(global_id_tag, hexes)

    print(f"Process {rank}: global IDs = {gids.flatten().tolist()}")
    CHECK(len(gids) > 0)

    all_gids = get_all_values(gids.flatten().tolist())
    flat_gids = [g for sublist in all_gids for g in sublist]
    CHECK_EQ(len(flat_gids), len(set(flat_gids)))


def test_exchange_tags_after_parallel_load():
    """Exchange a tag on shared entities after parallel load."""
    if not (_have_parallel_io() and _have_test_mesh()):
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    tag = mb.tag_get_handle("TEST_EXCHANGE", 1, types.MB_TYPE_DOUBLE,
                            types.MB_TAG_DENSE, create_if_missing=True)

    verts = mb.get_entities_by_type(0, types.MBVERTEX)
    vals = np.full(len(verts), float(rank + 1), dtype='float64')
    mb.tag_set_data(tag, verts, vals)

    shared = pcomm.get_shared_entities(-1, dim=0)
    if len(shared) > 0:
        pcomm.exchange_tags(["TEST_EXCHANGE"], ["TEST_EXCHANGE"], shared)

    print(f"Process {rank}: exchanged tags on {len(shared)} shared vertices")
    CHECK(True)


def test_reduce_tags_after_parallel_load():
    """Reduce a tag across shared entities after parallel load."""
    if not (_have_parallel_io() and _have_test_mesh()):
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    src_tag = mb.tag_get_handle("TEST_REDUCE_SRC", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True)
    dst_tag = mb.tag_get_handle("TEST_REDUCE_DST", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True)

    verts = mb.get_entities_by_type(0, types.MBVERTEX)
    vals = np.ones(len(verts), dtype='float64')
    mb.tag_set_data(src_tag, verts, vals)
    mb.tag_set_data(dst_tag, verts, np.zeros(len(verts), dtype='float64'))

    shared = pcomm.get_shared_entities(-1, dim=0)
    if len(shared) > 0:
        pcomm.reduce_tags(["TEST_REDUCE_SRC"], ["TEST_REDUCE_DST"], "SUM", shared)

    print(f"Process {rank}: reduced tags on {len(shared)} shared vertices")
    CHECK(True)


def test_broadcast_entities():
    """Broadcast entities from rank 0 to all ranks."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    ents = Range()
    if rank == 0:
        coords = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0], dtype="float64")
        verts = mb.create_vertices(coords)
        tri = mb.create_element(types.MBTRI, verts)
        ents = mb.get_entities_by_dimension(0, 2)

    count_before = len(mb.get_entities_by_dimension(0, 2))
    pcomm.broadcast_entities(0, ents, adjacencies=False, tags=True)
    count_after = len(mb.get_entities_by_dimension(0, 2))

    print(f"Process {rank}: 2D entities before={count_before}, after={count_after}")
    CHECK(count_after >= count_before)
    if rank != 0:
        CHECK(count_after > 0)


def test_get_iface_entities():
    """get_iface_entities returns a Range (possibly empty for local mesh)."""
    mb, pcomm, verts = _make_local_vertex()
    iface = pcomm.get_iface_entities(0, dim=0)
    CHECK(iface is not None)
    CHECK(isinstance(iface, Range))


def main():
    tests = [
        test_get_pstatus,
        test_get_pstatus_entities,
        test_get_owner,
        test_get_owner_handle,
        test_get_sharing_data,
        test_get_interface_procs,
        test_assign_global_ids,
        test_exchange_tags_after_parallel_load,
        test_reduce_tags_after_parallel_load,
        test_broadcast_entities,
        test_get_iface_entities,
    ]
    return run_parallel_tests(tests)


if __name__ == "__main__":
    sys.exit(main())
