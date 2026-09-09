#!/usr/bin/env python3
"""
End-to-end parallel mesh creation and shared entity resolution test.

Exercises the in-memory parallel workflow:
  1. Each rank creates a portion of a hex mesh with overlapping boundary vertices
  2. resolve_shared_ents() merges duplicate vertices and discovers shared edges/faces
  3. Verify shared entity counts and ownership consistency
  4. Assign global IDs
  5. Write the mesh in parallel
  6. Read it back and verify consistency
"""

import os
import sys
import tempfile
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


def _have_parallel_io():
    return config.MOAB_MPI_ENABLED and config.has_mpi_io()


def _create_hex_strip(mb, nx, ny, nz, x_offset=0.0, gid_offset=0):
    """Create an nx x ny x nz hex mesh starting at x_offset.

    Assigns GLOBAL_ID to each vertex based on its (i,j,k) position
    in a global grid: gid = k*(total_nx+1)*(ny+1) + j*(total_nx+1) + (i + gid_offset) + 1.
    Shared boundary vertices on adjacent ranks get identical global IDs,
    which is required for resolve_shared_ents to match them.
    """
    num_verts = (nx + 1) * (ny + 1) * (nz + 1)
    coords = np.zeros(num_verts * 3, dtype='float64')

    total_nx = nx * size
    gids = np.zeros(num_verts, dtype='int32')

    idx = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                coords[idx * 3] = x_offset + float(i)
                coords[idx * 3 + 1] = float(j)
                coords[idx * 3 + 2] = float(k)
                gids[idx] = k * (total_nx + 1) * (ny + 1) + j * (total_nx + 1) + (gid_offset + i) + 1
                idx += 1

    verts = mb.create_vertices(coords)

    gid_tag = mb.tag_get_handle("GLOBAL_ID", 1, types.MB_TYPE_INTEGER,
                                types.MB_TAG_DENSE, create_if_missing=True)
    mb.tag_set_data(gid_tag, verts, gids)

    hexes = Range()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                v = [0] * 8
                v[0] = verts[k * (ny + 1) * (nx + 1) + j * (nx + 1) + i]
                v[1] = verts[k * (ny + 1) * (nx + 1) + j * (nx + 1) + i + 1]
                v[2] = verts[k * (ny + 1) * (nx + 1) + (j + 1) * (nx + 1) + i + 1]
                v[3] = verts[k * (ny + 1) * (nx + 1) + (j + 1) * (nx + 1) + i]
                v[4] = verts[(k + 1) * (ny + 1) * (nx + 1) + j * (nx + 1) + i]
                v[5] = verts[(k + 1) * (ny + 1) * (nx + 1) + j * (nx + 1) + i + 1]
                v[6] = verts[(k + 1) * (ny + 1) * (nx + 1) + (j + 1) * (nx + 1) + i + 1]
                v[7] = verts[(k + 1) * (ny + 1) * (nx + 1) + (j + 1) * (nx + 1) + i]
                h = mb.create_element(types.MBHEX, v)
                hexes.insert(h)

    return verts, hexes, coords


def _create_overlapping_mesh(mb):
    """Create a strip of hexes per rank with shared boundary vertices.

    Rank r creates nx=2 hexes in x, shifted by r*2. Adjacent ranks share
    the boundary plane at x = r*2. Vertex global IDs ensure matching
    vertices on adjacent ranks have identical IDs.

    Returns (verts, hexes, coords, file_set) where file_set is a meshset
    containing all entities, required by resolve_shared_ents when no
    partition sets exist.
    """
    nx, ny, nz = 2, 2, 2
    x_offset = float(rank * nx)
    gid_offset = rank * nx
    verts, hexes, coords = _create_hex_strip(mb, nx, ny, nz, x_offset, gid_offset)

    file_set = mb.create_meshset()
    mb.add_entities(file_set, verts)
    mb.add_entities(file_set, hexes)
    return verts, hexes, coords, file_set


def _create_partition_set(mb, pcomm):
    """Create a PARALLEL_PARTITION entity set and populate it with owned entities.

    This creates a proper partition set tagged with PARALLEL_PARTITION = rank,
    which is required for the output file to contain the standard MOAB partition
    metadata. Without this, parallel files cannot be read back with
    PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION.

    Returns the partition EntityHandle.
    """
    part_set = pcomm.create_part()
    owned_verts = pcomm.get_owned_entities(dim=0)
    owned_hexes = pcomm.get_owned_entities(dim=3)
    mb.add_entities(part_set, owned_verts)
    mb.add_entities(part_set, owned_hexes)
    return part_set


def test_create_mesh_per_rank():
    """Each rank creates its portion of the mesh."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)

    print(f"  rank {rank}: created {len(hexes)} hexes, {len(verts)} vertices")

    CHECK_EQ(len(hexes), 8)   # 2x2x2
    CHECK_EQ(len(verts), 27)  # 3x3x3


def test_resolve_shared_ents():
    """Resolve shared entities and verify shared vertices appear."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)

    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    print(f"  rank {rank}: {len(shared_verts)} shared vertices after resolve")

    if size > 1:
        CHECK(len(shared_verts) > 0)
        if size == 2:
            CHECK_EQ(len(shared_verts), 9)


def test_shared_vertex_coordinates_match():
    """Shared vertices on adjacent ranks have identical coordinates."""
    if size < 2:
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        shared_coords = mb.get_coords(shared_verts)
        shared_coords_3d = shared_coords.reshape(-1, 3)

        # Shared boundary vertices have x = rank*2 (left boundary)
        # or x = (rank+1)*2 (right boundary)
        x_vals = shared_coords_3d[:, 0]
        boundary_x_left = float(rank * 2)
        boundary_x_right = float((rank + 1) * 2)

        for x in x_vals:
            CHECK(x == boundary_x_left or x == boundary_x_right)

    print(f"  rank {rank}: verified {len(shared_verts)} shared vertex coordinates")


def test_shared_vertex_ownership():
    """Exactly one rank owns each shared vertex."""
    if size < 2:
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    for sv in shared_verts:
        owner = pcomm.get_owner(sv)
        CHECK(0 <= owner < size)

    owned_shared = Range()
    for sv in shared_verts:
        if pcomm.get_owner(sv) == rank:
            owned_shared.insert(sv)

    all_owned_counts = get_all_values(len(owned_shared))
    total_owned = sum(all_owned_counts)

    # With 2 ranks: 9 shared vertices, each owned by exactly one rank
    if size == 2:
        CHECK_EQ(total_owned, 9)

    print(f"  rank {rank}: owns {len(owned_shared)} of {len(shared_verts)} shared vertices")


def test_assign_global_ids_after_resolve():
    """Assign global IDs after resolve and verify uniqueness."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)
    pcomm.assign_global_ids(dimension=3, start_id=1)

    gid_tag = mb.tag_get_handle("GLOBAL_ID")
    all_hexes = mb.get_entities_by_type(0, types.MBHEX)
    gids = mb.tag_get_data(gid_tag, all_hexes).flatten().tolist()

    all_gids = get_all_values(gids)
    flat_gids = [g for sublist in all_gids for g in sublist]

    CHECK_EQ(len(flat_gids), len(set(flat_gids)))
    print(f"  rank {rank}: global hex IDs = {gids}")


def test_check_shared_handles_after_resolve():
    """MOAB's internal consistency check should pass after resolve."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    pcomm.check_all_shared_handles(print_em=False)
    print(f"  rank {rank}: shared handle check passed")
    CHECK(True)


def test_exchange_tags_after_resolve():
    """Tag exchange works on a resolved in-memory mesh."""
    if size < 2:
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    tag = mb.tag_get_handle("RESOLVE_TEST_TAG", 1, types.MB_TYPE_DOUBLE,
                            types.MB_TAG_DENSE, create_if_missing=True)

    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    vals = np.full(len(all_verts), float(rank + 1), dtype='float64')
    mb.tag_set_data(tag, all_verts, vals)

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        pcomm.exchange_tags(["RESOLVE_TEST_TAG"], ["RESOLVE_TEST_TAG"], shared_verts)

        for sv in shared_verts:
            val = mb.tag_get_data(tag, sv).flatten()[0]
            owner = pcomm.get_owner(sv)
            CHECK_EQ(val, float(owner + 1))

    print(f"  rank {rank}: tag exchange verified on {len(shared_verts)} shared vertices")


def test_parallel_write_after_resolve():
    """Write the resolved mesh in parallel with PARALLEL_PARTITION sets."""
    if not _have_parallel_io():
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)
    pcomm.assign_global_ids(dimension=3, start_id=1)

    part_set = _create_partition_set(mb, pcomm)

    outfile = os.path.join(tempfile.gettempdir(), "pymoab_resolve_write_test.h5m")
    pcomm.write_file(outfile, "PARALLEL=WRITE_PART")

    comm.Barrier()
    if rank == 0:
        CHECK(os.path.exists(outfile))
        fsize = os.path.getsize(outfile)
        print(f"  Written resolved mesh: {outfile} ({fsize} bytes)")
        CHECK(fsize > 0)

        # Verify the file contains PARALLEL_PARTITION tag
        mb2 = core.Core()
        mb2.load_file(outfile)
        pp_tag = mb2.tag_get_handle("PARALLEL_PARTITION")
        CHECK(pp_tag is not None)
        part_sets = mb2.get_entities_by_type_and_tag(
            0, types.MBENTITYSET, pp_tag, [None])
        CHECK_EQ(len(part_sets), size)
        print(f"  File contains {len(part_sets)} PARALLEL_PARTITION sets")
        os.unlink(outfile)
    else:
        CHECK(True)


def test_parallel_write_read_roundtrip():
    """Write mesh in parallel with partition sets, read back in parallel."""
    if not _have_parallel_io():
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)
    pcomm.assign_global_ids(dimension=3, start_id=1)

    original_hex_count = len(mb.get_entities_by_type(0, types.MBHEX))
    original_owned_vert_count = len(pcomm.get_owned_entities(dim=0))
    all_orig_hexes = sum(get_all_values(original_hex_count))
    all_orig_verts = sum(get_all_values(original_owned_vert_count))

    part_set = _create_partition_set(mb, pcomm)

    outfile = os.path.join(tempfile.gettempdir(), "pymoab_roundtrip_test.h5m")
    pcomm.write_file(outfile, "PARALLEL=WRITE_PART")
    comm.Barrier()

    # Re-read in parallel using PARALLEL_PARTITION sets
    mb2 = core.Core()
    pcomm2 = parallelcomm.ParallelComm(mb2, comm)
    pcomm2.load_file(outfile,
                     "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION"
                     ";PARALLEL_RESOLVE_SHARED_ENTS")

    reread_hexes = mb2.get_entities_by_type(0, types.MBHEX)
    reread_owned_hexes = pcomm2.get_owned_entities(dim=3)
    reread_owned_verts = pcomm2.get_owned_entities(dim=0)

    all_reread_hexes = sum(get_all_values(len(reread_owned_hexes)))
    all_reread_verts = sum(get_all_values(len(reread_owned_verts)))

    CHECK_EQ(all_reread_hexes, all_orig_hexes)
    CHECK_EQ(all_reread_verts, all_orig_verts)

    print(f"  rank {rank}: roundtrip owned hexes {len(reread_owned_hexes)}, "
          f"owned verts {len(reread_owned_verts)}")
    print(f"  rank {rank}: global totals: hexes {all_orig_hexes} -> {all_reread_hexes}, "
          f"owned verts {all_orig_verts} -> {all_reread_verts}")

    comm.Barrier()
    if rank == 0:
        os.unlink(outfile)


def test_reduce_after_resolve():
    """Reduce tags on a resolved in-memory mesh."""
    if size < 2:
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)

    src_tag = mb.tag_get_handle("RESOLVE_REDUCE_SRC", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))
    dst_tag = mb.tag_get_handle("RESOLVE_REDUCE_DST", 1, types.MB_TYPE_DOUBLE,
                                types.MB_TAG_DENSE, create_if_missing=True,
                                default_value=np.array([0.0]))

    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    mb.tag_set_data(src_tag, all_verts, np.ones(len(all_verts), dtype='float64'))
    mb.tag_set_data(dst_tag, all_verts, np.zeros(len(all_verts), dtype='float64'))

    shared_verts = pcomm.get_shared_entities(-1, dim=0)
    if len(shared_verts) > 0:
        pcomm.reduce_tags(["RESOLVE_REDUCE_SRC"], ["RESOLVE_REDUCE_DST"], "SUM", shared_verts)

        owned_shared = Range()
        for sv in shared_verts:
            if pcomm.get_owner(sv) == rank:
                owned_shared.insert(sv)

        if len(owned_shared) > 0:
            dst_vals = mb.tag_get_data(dst_tag, owned_shared).flatten()
            # With 2 ranks sharing a boundary plane: each shared vertex has SUM=2
            if size == 2:
                CHECK(np.all(dst_vals == 2.0))

    print(f"  rank {rank}: reduced tags on {len(shared_verts)} shared vertices")


def test_full_resolve_pipeline():
    """Full pipeline: create -> resolve -> assign IDs -> partition -> tag -> exchange -> write.

    Exercises the complete in-memory parallel mesh workflow including
    proper PARALLEL_PARTITION sets for standard MOAB parallel file output.
    """
    if not _have_parallel_io():
        CHECK(True)
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)

    verts, hexes, coords, file_set = _create_overlapping_mesh(mb)
    print(f"  rank {rank}: [1/6] created {len(hexes)} hexes, {len(verts)} verts")

    pcomm.resolve_shared_ents(this_set=file_set, resolve_dim=3, shared_dim=-1)
    shared = pcomm.get_shared_entities(-1, dim=0)
    print(f"  rank {rank}: [2/6] resolved {len(shared)} shared vertices")

    pcomm.assign_global_ids(dimension=3, start_id=1)
    pcomm.assign_global_ids(dimension=0, start_id=1)
    print(f"  rank {rank}: [3/6] assigned global IDs")

    part_set = _create_partition_set(mb, pcomm)
    print(f"  rank {rank}: [4/6] created PARALLEL_PARTITION set")

    tag = mb.tag_get_handle("PIPELINE_DATA", 1, types.MB_TYPE_DOUBLE,
                            types.MB_TAG_DENSE, create_if_missing=True,
                            default_value=np.array([0.0]))
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    vert_coords = mb.get_coords(all_verts).reshape(-1, 3)
    dist_vals = np.sqrt(np.sum(vert_coords ** 2, axis=1))
    mb.tag_set_data(tag, all_verts, dist_vals)

    if len(shared) > 0:
        pcomm.exchange_tags(["PIPELINE_DATA"], ["PIPELINE_DATA"], shared)
    print(f"  rank {rank}: [5/6] set and exchanged distance tag")

    outfile = os.path.join(tempfile.gettempdir(), "pymoab_full_pipeline.h5m")
    pcomm.write_file(outfile, "PARALLEL=WRITE_PART")
    comm.Barrier()
    if rank == 0:
        CHECK(os.path.exists(outfile))
        fsize = os.path.getsize(outfile)
        print(f"  rank {rank}: [6/6] wrote {fsize} bytes to {outfile}")
        os.unlink(outfile)
    else:
        print(f"  rank {rank}: [6/6] write complete")
    CHECK(True)


def main():
    tests = [
        test_create_mesh_per_rank,
        test_resolve_shared_ents,
        test_shared_vertex_coordinates_match,
        test_shared_vertex_ownership,
        test_assign_global_ids_after_resolve,
        test_check_shared_handles_after_resolve,
        test_exchange_tags_after_resolve,
        test_parallel_write_after_resolve,
        test_parallel_write_read_roundtrip,
        test_reduce_after_resolve,
        test_full_resolve_pipeline,
    ]
    return run_parallel_tests(tests)


if __name__ == "__main__":
    sys.exit(main())
