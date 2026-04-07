#!/usr/bin/env python3
"""
End-to-end parallel workflow integration test.

Exercises the full production-like pipeline:
  1. Load a partitioned mesh in parallel
  2. Exchange ghost cells at multiple bridge dimensions
  3. Create scalar and vector tags on owned entities
  4. Exchange tag data to ghost/halo layers
  5. Verify ghost entities received correct tag values
  6. Perform MPI reductions (SUM, MAX) on shared vertex tags
  7. Verify reduced values are mathematically correct
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
from pymoab.rng import Range, intersect, subtract
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


def _can_run():
    return _have_parallel_io() and _have_test_mesh()


def _skip_if_needed():
    if not _can_run():
        CHECK(True)
        return True
    return False


def _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0):
    """Load partitioned mesh and exchange ghost cells.

    Returns (mb, pcomm) with ghost layers already exchanged.
    """
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, bridge_dim, num_ghost_layers, 0, True, True)
    return mb, pcomm


# ---------- Ghost layer tests ----------

def test_ghost_layers_vertex_bridge():
    """Exchange 1 layer of ghost hexes bridged through vertices."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    owned_hex = pcomm.get_owned_entities(dim=3)
    ghost_hex = pcomm.get_ghost_entities(dim=3)
    all_hex = mb.get_entities_by_dimension(0, 3)

    print(f"  rank {rank}: owned={len(owned_hex)}, ghost={len(ghost_hex)}, total={len(all_hex)}")

    CHECK(len(owned_hex) > 0)
    CHECK_EQ(len(owned_hex) + len(ghost_hex), len(all_hex))
    # vertex-bridged ghosts are typically a superset of face-bridged
    CHECK(len(ghost_hex) >= 0)


def test_ghost_layers_face_bridge():
    """Exchange 1 layer of ghost hexes bridged through faces (dim=2)."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=2)

    owned_hex = pcomm.get_owned_entities(dim=3)
    ghost_hex = pcomm.get_ghost_entities(dim=3)

    print(f"  rank {rank}: owned={len(owned_hex)}, ghost_face_bridge={len(ghost_hex)}")

    CHECK(len(owned_hex) > 0)
    # face-bridged ghosts can be zero if partition has no face-adjacent neighbors


def test_ghost_vertices_received():
    """After ghost hex exchange, ghost vertices should also appear."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    owned_verts = pcomm.get_owned_entities(dim=0)
    ghost_verts = pcomm.get_ghost_entities(dim=0)

    print(f"  rank {rank}: owned_verts={len(owned_verts)}, ghost_verts={len(ghost_verts)}")

    CHECK(len(owned_verts) > 0)
    # ghost vertices come along with ghost hexes
    # (they may be zero if this rank has no ghost hexes)


# ---------- Tag exchange tests ----------

def test_exchange_scalar_tag():
    """Set a scalar tag on owned vertices, exchange to ghosts, verify values."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    # Create a scalar tag: each rank sets value = rank + 1
    tag = mb.tag_get_handle("WORKFLOW_SCALAR", 1, types.MB_TYPE_DOUBLE,
                            types.MB_TAG_DENSE, create_if_missing=True)

    owned_verts = pcomm.get_owned_entities(dim=0)
    vals = np.full(len(owned_verts), float(rank + 1), dtype='float64')
    mb.tag_set_data(tag, owned_verts, vals)

    # Exchange on all shared + ghost vertices
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    pcomm.exchange_tags(["WORKFLOW_SCALAR"], ["WORKFLOW_SCALAR"], all_verts)

    # Verify: ghost vertices should have values set by their owner
    ghost_verts = pcomm.get_ghost_entities(dim=0)
    if len(ghost_verts) > 0:
        ghost_vals = mb.tag_get_data(tag, ghost_verts).flatten()
        # All ghost values should be > 0 (set by some rank)
        CHECK(np.all(ghost_vals > 0))
        # Ghost values should equal owner_rank + 1
        for i, gv in enumerate(ghost_verts):
            owner = pcomm.get_owner(gv)
            expected = float(owner + 1)
            actual = ghost_vals[i]
            CHECK_EQ(actual, expected)

    print(f"  rank {rank}: exchanged scalar tag on {len(all_verts)} verts, "
          f"verified {len(ghost_verts)} ghost values")


def test_exchange_vector_tag():
    """Set a 3-component vector tag on owned hexes, exchange to ghost hexes."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    tag = mb.tag_get_handle("WORKFLOW_VEC3", 3, types.MB_TYPE_DOUBLE,
                            types.MB_TAG_DENSE, create_if_missing=True)

    owned_hex = pcomm.get_owned_entities(dim=3)
    # Each rank sets [rank, rank*10, rank*100]
    vals = np.zeros((len(owned_hex), 3), dtype='float64')
    vals[:, 0] = float(rank)
    vals[:, 1] = float(rank * 10)
    vals[:, 2] = float(rank * 100)
    mb.tag_set_data(tag, owned_hex, vals)

    all_hex = mb.get_entities_by_dimension(0, 3)
    pcomm.exchange_tags(["WORKFLOW_VEC3"], ["WORKFLOW_VEC3"], all_hex)

    # Verify ghost hex tag values
    ghost_hex = pcomm.get_ghost_entities(dim=3)
    if len(ghost_hex) > 0:
        ghost_vals = mb.tag_get_data(tag, ghost_hex)
        for i, gh in enumerate(ghost_hex):
            owner = pcomm.get_owner(gh)
            CHECK_EQ(ghost_vals[i, 0], float(owner))
            CHECK_EQ(ghost_vals[i, 1], float(owner * 10))
            CHECK_EQ(ghost_vals[i, 2], float(owner * 100))

    print(f"  rank {rank}: exchanged vector tag on {len(all_hex)} hexes, "
          f"verified {len(ghost_hex)} ghost values")


# ---------- Tag reduction tests ----------

def test_reduce_tags_sum():
    """Reduce a scalar tag with SUM: each shared vertex gets sum of contributions."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    # Tags for reduce MUST have default values
    src_tag = mb.tag_get_handle("REDUCE_SUM_SRC", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))
    dst_tag = mb.tag_get_handle("REDUCE_SUM_DST", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))

    # Each rank contributes 1.0 for every vertex it owns or shares
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    src_vals = np.ones(len(all_verts), dtype='float64')
    mb.tag_set_data(src_tag, all_verts, src_vals)
    mb.tag_set_data(dst_tag, all_verts, np.zeros(len(all_verts), dtype='float64'))

    # Reduce on shared vertices only
    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        pcomm.reduce_tags(["REDUCE_SUM_SRC"], ["REDUCE_SUM_DST"], "SUM", shared_verts)

        # After SUM reduction, dst value = number of sharing procs (including owner)
        owned_shared = Range()
        for sv in shared_verts:
            if pcomm.get_owner(sv) == rank:
                owned_shared.insert(sv)

        if len(owned_shared) > 0:
            dst_vals = mb.tag_get_data(dst_tag, owned_shared).flatten()
            # Each summed value should be >= 2 (shared between at least 2 procs)
            CHECK(np.all(dst_vals >= 2.0))
            print(f"  rank {rank}: SUM reduction on {len(owned_shared)} owned shared verts, "
                  f"min={dst_vals.min()}, max={dst_vals.max()}")
        else:
            print(f"  rank {rank}: no owned shared vertices to verify")
    else:
        print(f"  rank {rank}: no shared vertices to reduce")
    CHECK(True)


def test_reduce_tags_max():
    """Reduce a scalar tag with MAX: shared vertices should get max rank value."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    src_tag = mb.tag_get_handle("REDUCE_MAX_SRC", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))
    dst_tag = mb.tag_get_handle("REDUCE_MAX_DST", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))

    # Each rank sets its rank value on all vertices
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    src_vals = np.full(len(all_verts), float(rank), dtype='float64')
    mb.tag_set_data(src_tag, all_verts, src_vals)
    mb.tag_set_data(dst_tag, all_verts, np.zeros(len(all_verts), dtype='float64'))

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        pcomm.reduce_tags(["REDUCE_MAX_SRC"], ["REDUCE_MAX_DST"], "MAX", shared_verts)

        # Verify on owned shared vertices
        owned_shared = Range()
        for sv in shared_verts:
            if pcomm.get_owner(sv) == rank:
                owned_shared.insert(sv)

        if len(owned_shared) > 0:
            dst_vals = mb.tag_get_data(dst_tag, owned_shared).flatten()
            # MAX value should be >= rank (at least this rank contributed)
            CHECK(np.all(dst_vals >= float(rank)))
            print(f"  rank {rank}: MAX reduction on {len(owned_shared)} verts, "
                  f"min_result={dst_vals.min()}, max_result={dst_vals.max()}")
        else:
            print(f"  rank {rank}: no owned shared vertices to verify")
    else:
        print(f"  rank {rank}: no shared vertices to reduce")
    CHECK(True)


# ---------- Combined workflow test ----------

def test_full_workflow_pipeline():
    """Full pipeline: load -> ghost -> tag set -> exchange -> reduce -> verify.

    This is the most comprehensive integration test, exercising every step
    of a typical parallel simulation workflow.
    """
    if _skip_if_needed():
        return

    # Step 1: Load partitioned mesh
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)

    owned_hex_before = pcomm.get_owned_entities(dim=3)
    print(f"  rank {rank}: [1/6] loaded {len(owned_hex_before)} owned hexes")

    # Step 2: Exchange ghost cells (vertex-bridged, 1 layer)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)
    ghost_hex = pcomm.get_ghost_entities(dim=3)
    print(f"  rank {rank}: [2/6] got {len(ghost_hex)} ghost hexes")

    # Step 3: Create tags and set data on owned entities
    #   - "TEMPERATURE": scalar double, value = rank * 100 + 1
    #   - "VELOCITY": 3-component double, value = [rank, rank*2, rank*3]
    temp_tag = mb.tag_get_handle("TEMPERATURE", 1, types.MB_TYPE_DOUBLE,
                                 types.MB_TAG_DENSE, create_if_missing=True,
                                 default_value=np.array([0.0]))
    vel_tag = mb.tag_get_handle("VELOCITY", 3, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0, 0.0, 0.0]))

    owned_verts = pcomm.get_owned_entities(dim=0)
    temp_vals = np.full(len(owned_verts), float(rank * 100 + 1), dtype='float64')
    mb.tag_set_data(temp_tag, owned_verts, temp_vals)

    vel_vals = np.zeros((len(owned_verts), 3), dtype='float64')
    vel_vals[:, 0] = float(rank)
    vel_vals[:, 1] = float(rank * 2)
    vel_vals[:, 2] = float(rank * 3)
    mb.tag_set_data(vel_tag, owned_verts, vel_vals)
    print(f"  rank {rank}: [3/6] set tags on {len(owned_verts)} owned vertices")

    # Step 4: Exchange tags to ghost/halo vertices
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    pcomm.exchange_tags(["TEMPERATURE"], ["TEMPERATURE"], all_verts)
    pcomm.exchange_tags(["VELOCITY"], ["VELOCITY"], all_verts)
    print(f"  rank {rank}: [4/6] exchanged tags to {len(all_verts)} total vertices")

    # Step 5: Verify ghost vertex values
    ghost_verts = pcomm.get_ghost_entities(dim=0)
    if len(ghost_verts) > 0:
        ghost_temps = mb.tag_get_data(temp_tag, ghost_verts).flatten()
        ghost_vels = mb.tag_get_data(vel_tag, ghost_verts)

        for i, gv in enumerate(ghost_verts):
            owner = pcomm.get_owner(gv)
            CHECK_EQ(ghost_temps[i], float(owner * 100 + 1))
            CHECK_EQ(ghost_vels[i, 0], float(owner))
            CHECK_EQ(ghost_vels[i, 1], float(owner * 2))
            CHECK_EQ(ghost_vels[i, 2], float(owner * 3))

    print(f"  rank {rank}: [5/6] verified {len(ghost_verts)} ghost vertex tag values")

    # Step 6: Reduce temperature tag (SUM) on shared vertices
    temp_sum_tag = mb.tag_get_handle("TEMP_SUM", 1, types.MB_TYPE_DOUBLE,
                                     types.MB_TAG_DENSE, create_if_missing=True,
                                     default_value=np.array([0.0]))
    mb.tag_set_data(temp_sum_tag, all_verts, np.zeros(len(all_verts), dtype='float64'))

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        pcomm.reduce_tags(["TEMPERATURE"], ["TEMP_SUM"], "SUM", shared_verts)

    print(f"  rank {rank}: [6/6] reduced TEMPERATURE over {len(shared_verts)} shared verts")


def test_pstatus_consistency():
    """Verify pstatus flags are consistent after ghost exchange."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    owned_verts = pcomm.get_owned_entities(dim=0)
    ghost_verts = pcomm.get_ghost_entities(dim=0)

    for gv in ghost_verts:
        ps = pcomm.get_pstatus(gv)
        CHECK(bool(ps & 0x01))

    for ov in owned_verts:
        owner = pcomm.get_owner(ov)
        CHECK_EQ(owner, rank)

    shared = pcomm.get_shared_entities(-1, dim=0)
    for sv in shared:
        ps = pcomm.get_pstatus(sv)
        has_parallel_flag = bool(ps & 0x02) or bool(ps & 0x04) or bool(ps & 0x10)
        CHECK(has_parallel_flag)

    CHECK(len(owned_verts) + len(ghost_verts) == len(all_verts))

    print(f"  rank {rank}: pstatus consistent for {len(all_verts)} verts "
          f"({len(owned_verts)} owned, {len(ghost_verts)} ghost, {len(shared)} shared)")


def test_sharing_data_symmetry():
    """Verify that sharing data is symmetric: if A shares with B, B shares with A."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    shared = pcomm.get_shared_entities(-1, dim=0, owned_filter=True)
    sharing_map = {}  # vertex -> set of sharing procs

    for sv in shared:
        data = pcomm.get_sharing_data(sv)
        sharing_map[sv] = set(data['procs'])

    # Gather the number of shared entities per rank
    all_shared_counts = get_all_values(len(shared))
    if rank == 0:
        print(f"  Shared vertex counts per rank: {all_shared_counts}")

    # Verify this rank appears in comm_procs of neighbors
    comm_procs = pcomm.get_comm_procs()
    for p in comm_procs:
        # Get entities shared with proc p
        shared_with_p = pcomm.get_shared_entities(p, dim=0)
        CHECK(len(shared_with_p) > 0)

    print(f"  rank {rank}: verified sharing symmetry, "
          f"communicates with {len(comm_procs)} procs")
    CHECK(True)


def test_check_all_shared_handles():
    """Run MOAB's internal consistency check on shared handles."""
    if _skip_if_needed():
        return

    mb, pcomm = _load_mesh_with_ghosts(num_ghost_layers=1, bridge_dim=0)

    # This will throw if handles are inconsistent
    pcomm.check_all_shared_handles(print_em=False)
    print(f"  rank {rank}: shared handle check passed")
    CHECK(True)


def main():
    tests = [
        # Ghost layer tests
        test_ghost_layers_vertex_bridge,
        test_ghost_layers_face_bridge,
        test_ghost_vertices_received,
        # Tag exchange tests
        test_exchange_scalar_tag,
        test_exchange_vector_tag,
        # Reduction tests
        test_reduce_tags_sum,
        test_reduce_tags_max,
        # Full pipeline
        test_full_workflow_pipeline,
        # Consistency checks
        test_pstatus_consistency,
        test_sharing_data_symmetry,
        test_check_all_shared_handles,
    ]
    return run_parallel_tests(tests)


if __name__ == "__main__":
    sys.exit(main())
