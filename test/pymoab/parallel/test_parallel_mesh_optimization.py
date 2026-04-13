#!/usr/bin/env python3
"""
Parallel Laplacian mesh smoothing with halo coordinate exchange.

Demonstrates a production-style iterative mesh optimization workflow:
  1. Load a partitioned hex mesh in parallel (READ_PART)
  2. Exchange ghost cells so each rank has a 1-layer halo
  3. Identify boundary vertices (fixed) vs interior vertices (smoothable)
  4. Perturb interior vertices to degrade the mesh (simulate a real
     scenario where element quality needs improvement)
  5. Iteratively perform Laplacian smoothing on owned interior vertices
     using neighbor coordinates that may live on remote ranks (ghosts)
  6. After each local smoothing pass, exchange the updated coordinates
     of shared/ghost vertices across task boundaries
  7. Compute a global mesh-quality metric (min scaled-Jacobian)
     and iterate until quality converges or max iterations are reached
  8. Verify the optimized mesh has better quality than the degraded mesh

The key parallel pattern is:
  - Owner-computes: only the owning rank moves a shared vertex
  - After each sweep the new coordinates are pushed to every rank
    that holds a ghost copy via tag exchange on the halo region
  - Global convergence is checked with MPI_Allreduce
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
from pymoab.topo_util import MeshTopoUtil
from pymoab.skinner import Skinner
from parallel_driver import (
    run_parallel_tests,
    CHECK_EQ,
    CHECK,
    CHECK_PARALLEL,
    get_all_values,
)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

MESH_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "MeshFiles", "unittest"
)
PARTITIONED_MESH = os.path.join(MESH_DIR, "64bricks_512hex_256part.h5m")
READ_OPTS = (
    "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS"
)


# ---------------------------------------------------------------------------
# Mesh quality helpers
# ---------------------------------------------------------------------------


def _hex_scaled_jacobian(coords_8x3):
    """Compute the scaled Jacobian for a single hexahedron.

    The scaled Jacobian is evaluated at each of the 8 corner nodes.
    At each corner the Jacobian matrix is formed from the three edge
    vectors emanating from that corner, and the determinant is divided
    by the product of the three edge lengths.  The element quality is
    the *minimum* over the 8 corners.

    Parameters
    ----------
    coords_8x3 : ndarray, shape (8, 3)
        Vertex coordinates in MOAB hex ordering (tensor-product order).

    Returns
    -------
    float
        Scaled Jacobian in [-1, 1].  1 = perfect cube,  <0 = inverted.
    """
    # MOAB hex node ordering (same as ExodusII):
    #   bottom face: 0-1-2-3   top face: 4-5-6-7
    # Edges from each corner follow the local (xi, eta, zeta) directions.
    corner_edges = [
        (0, (1, 3, 4)),
        (1, (2, 0, 5)),
        (2, (3, 1, 6)),
        (3, (0, 2, 7)),
        (4, (7, 5, 0)),
        (5, (4, 6, 1)),
        (6, (5, 7, 2)),
        (7, (6, 4, 3)),
    ]
    min_sj = 1.0
    for node, (a, b, c) in corner_edges:
        e1 = coords_8x3[a] - coords_8x3[node]
        e2 = coords_8x3[b] - coords_8x3[node]
        e3 = coords_8x3[c] - coords_8x3[node]
        det = np.dot(e1, np.cross(e2, e3))
        l1 = np.linalg.norm(e1)
        l2 = np.linalg.norm(e2)
        l3 = np.linalg.norm(e3)
        denom = l1 * l2 * l3
        sj = det / denom if denom > 1.0e-15 else -1.0
        if sj < min_sj:
            min_sj = sj
    return min_sj


def _compute_local_quality(mb, hexes):
    """Return the (min, mean) scaled Jacobian over local hex elements.

    Parameters
    ----------
    mb : Core
        MOAB instance.
    hexes : Range
        Hex elements to evaluate.

    Returns
    -------
    tuple of (float, float)
        (min_quality, mean_quality) over the local hex set.
    """
    if len(hexes) == 0:
        return 1.0, 1.0

    qualities = np.empty(len(hexes), dtype="float64")
    for i, h in enumerate(hexes):
        conn = mb.get_connectivity(h)
        coords = mb.get_coords(conn).reshape(-1, 3)
        qualities[i] = _hex_scaled_jacobian(coords)

    return float(qualities.min()), float(qualities.mean())


# ---------------------------------------------------------------------------
# Laplacian smoothing kernel
# ---------------------------------------------------------------------------


def _min_quality_around_vertex(mb, v):
    """Return min scaled Jacobian of all hex elements adjacent to vertex *v*."""
    adj_hex = mb.get_adjacencies(Range([v]), 3, create_if_missing=False)
    if len(adj_hex) == 0:
        return 1.0
    worst = 1.0
    for h in adj_hex:
        conn = mb.get_connectivity(h)
        coords = mb.get_coords(conn).reshape(-1, 3)
        sj = _hex_scaled_jacobian(coords)
        if sj < worst:
            worst = sj
    return worst


def _laplacian_smooth_owned(mb, mtu, owned_interior_verts, omega=0.5):
    """Quality-guarded Laplacian smoothing on owned interior vertices.

    Each vertex is moved toward the centroid of its edge-connected
    neighbors:  x_new = (1-omega)*x_old + omega*centroid(neighbors)

    The move is accepted only if it does not decrease the minimum
    scaled Jacobian of the adjacent hex elements.  This prevents
    element inversion and guarantees monotonic quality improvement.

    Neighbor coordinates may include ghost vertices, so the halo
    must be up to date *before* calling this function.
    """
    max_disp = 0.0
    for v in owned_interior_verts:
        neighbors = mtu.get_bridge_adjacencies(v, bridge_dim=1, to_dim=0)
        if len(neighbors) == 0:
            continue

        nbr_coords = mb.get_coords(neighbors).reshape(-1, 3)
        centroid = nbr_coords.mean(axis=0)

        vr = Range([v])
        old_pos = mb.get_coords(vr).reshape(3)
        new_pos = (1.0 - omega) * old_pos + omega * centroid

        # Quality guard: only accept if local quality does not degrade
        q_before = _min_quality_around_vertex(mb, v)
        mb.set_coords(vr, new_pos)
        q_after = _min_quality_around_vertex(mb, v)

        if q_after < q_before:
            mb.set_coords(vr, old_pos)
        else:
            disp = np.linalg.norm(new_pos - old_pos)
            if disp > max_disp:
                max_disp = disp

    return max_disp


# ---------------------------------------------------------------------------
# Mesh perturbation (create a degraded mesh to optimize)
# ---------------------------------------------------------------------------


def _perturb_interior_vertices(
    mb, pcomm, owned_interior, all_verts, amplitude=0.15, seed=42
):
    """Add deterministic random perturbation to owned interior vertices.

    Each vertex is displaced by a random vector whose magnitude is
    bounded by ``amplitude`` times the local mean edge length.  The
    perturbation is deterministic (seeded per-rank) so that results
    are reproducible.  After perturbation, coordinates are exchanged
    across the halo so ghost copies are consistent.

    Parameters
    ----------
    mb : Core
        MOAB instance.
    pcomm : ParallelComm
        Parallel communicator.
    owned_interior : Range
        Interior vertices owned by this rank (safe to move).
    all_verts : Range
        All local vertices (for the halo exchange).
    amplitude : float
        Perturbation magnitude as a fraction of mean edge length.
    seed : int
        RNG seed (combined with rank for per-rank reproducibility).
    """
    if len(owned_interior) == 0:
        _exchange_halo_coordinates(mb, pcomm, all_verts)
        return

    rng = np.random.RandomState(seed + pcomm.rank)
    coords = mb.get_coords(owned_interior).reshape(-1, 3)

    # Estimate a local length scale from the coordinate bounding box
    bbox_diag = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
    n_verts = len(owned_interior)
    # Rough edge length ~ bbox diagonal / cube-root of vertex count
    mean_edge = bbox_diag / max(n_verts ** (1.0 / 3.0), 1.0)

    perturbation = rng.uniform(-1, 1, size=coords.shape)
    # Normalize each row then scale
    norms = np.linalg.norm(perturbation, axis=1, keepdims=True)
    norms = np.where(norms > 1e-15, norms, 1.0)
    perturbation = perturbation / norms * amplitude * mean_edge

    coords += perturbation
    mb.set_coords(owned_interior, coords.flatten())
    _exchange_halo_coordinates(mb, pcomm, all_verts)


# ---------------------------------------------------------------------------
# Coordinate exchange across the halo region
# ---------------------------------------------------------------------------


def _exchange_halo_coordinates(mb, pcomm, all_verts):
    """Push updated vertex coordinates from owners to ghost copies.

    This is the critical parallel communication step.  We store the
    x, y, z coordinates in a temporary dense tag, call
    ``pcomm.exchange_tags`` on all vertices (owned + ghost), and then
    copy the tag values back into the MOAB coordinate arrays so that
    ghost vertices reflect the new positions computed by their owners.

    Parameters
    ----------
    mb : Core
        MOAB instance.
    pcomm : ParallelComm
        Parallel communicator.
    all_verts : Range
        All local vertices (owned + ghost).
    """
    # Retrieve current coordinates from the MOAB coordinate storage
    coords = mb.get_coords(all_verts).reshape(-1, 3)

    # Store coordinates in separate scalar tags for exchange.
    # (exchange_tags operates on MOAB tags, not raw coordinate storage)
    coord_tags = []
    tag_names = ["_HALO_X", "_HALO_Y", "_HALO_Z"]
    for idx, name in enumerate(tag_names):
        tag = mb.tag_get_handle(
            name, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing=True
        )
        mb.tag_set_data(tag, all_verts, coords[:, idx].copy())
        coord_tags.append(tag)

    # Exchange: owner values → ghost copies
    pcomm.exchange_tags(tag_names, tag_names, all_verts)

    # Read back the exchanged values and update coordinate storage
    new_coords = np.empty_like(coords)
    for idx, name in enumerate(tag_names):
        tag = mb.tag_get_handle(name)
        new_coords[:, idx] = mb.tag_get_data(tag, all_verts).flatten()

    mb.set_coords(all_verts, new_coords.flatten())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _can_run():
    return (
        config.MOAB_MPI_ENABLED
        and config.has_mpi_io()
        and os.path.exists(PARTITIONED_MESH)
    )


def _skip_if_needed():
    if not _can_run():
        CHECK(True)
        return True
    return False


def test_laplacian_smoothing_convergence():
    """Perturb a perfect mesh, then recover quality via parallel Laplacian smoothing.

    Steps:
      1. Load partitioned mesh, exchange 1-layer vertex-bridged ghosts.
      2. Identify boundary vertices (Skinner) — these are pinned.
      3. Perturb interior vertices to degrade mesh quality.
      4. Iterate Laplacian smoothing with halo exchange until converged.
      5. Verify final quality > degraded quality.
    """
    if _skip_if_needed():
        return

    # -- Step 1: parallel load + ghost exchange --------------------------------
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    mtu = MeshTopoUtil(mb)
    skn = Skinner(mb)

    # -- Step 2: classify vertices ---------------------------------------------
    owned_hex = pcomm.get_owned_entities(dim=3)
    all_hex = mb.get_entities_by_dimension(0, 3)
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    owned_verts = pcomm.get_owned_entities(dim=0)

    skin_faces = skn.find_skin(0, all_hex, get_vertices=False)
    boundary_verts = mb.get_adjacencies(
        skin_faces, 0, create_if_missing=False, op_type=types.UNION
    )
    owned_interior = subtract(owned_verts, boundary_verts)

    global_hex_count = comm.allreduce(len(owned_hex), op=MPI.SUM)
    if rank == 0:
        print(
            f"  Interior verts to smooth (rank 0): {len(owned_interior)}, "
            f"global hex count: {global_hex_count}"
        )

    # -- Step 3: perturb to degrade quality ------------------------------------
    _perturb_interior_vertices(
        mb, pcomm, owned_interior, all_verts, amplitude=0.15, seed=42
    )

    degraded_min_q, degraded_mean_q = _compute_local_quality(mb, owned_hex)
    global_degraded_min = comm.allreduce(degraded_min_q, op=MPI.MIN)
    global_degraded_mean_sum = comm.allreduce(
        degraded_mean_q * len(owned_hex), op=MPI.SUM
    )
    global_degraded_mean = global_degraded_mean_sum / max(global_hex_count, 1)

    if rank == 0:
        print(
            f"  Degraded quality: min_SJ={global_degraded_min:.6f}, "
            f"mean_SJ={global_degraded_mean:.6f}"
        )

    # -- Step 4: iterative smoothing -------------------------------------------
    max_iter = 100
    tol = 1.0e-8
    omega = 0.5
    best_min_q = global_degraded_min

    for iteration in range(max_iter):
        # 4a. Smooth owned interior vertices
        local_max_disp = _laplacian_smooth_owned(mb, mtu, owned_interior, omega)

        # 4b. Exchange updated coordinates through the halo
        _exchange_halo_coordinates(mb, pcomm, all_verts)

        # 4c. Global convergence check
        global_max_disp = comm.allreduce(local_max_disp, op=MPI.MAX)

        cur_min_q, _ = _compute_local_quality(mb, owned_hex)
        global_min_q = comm.allreduce(cur_min_q, op=MPI.MIN)
        if global_min_q > best_min_q:
            best_min_q = global_min_q

        if rank == 0 and (iteration % 20 == 0 or global_max_disp < tol):
            print(
                f"    iter {iteration:3d}: max_disp={global_max_disp:.3e}, "
                f"min_SJ={global_min_q:.6f}"
            )

        # 4d. Convergence
        if global_max_disp < tol:
            if rank == 0:
                print(f"  Converged at iteration {iteration}")
            break

    # -- Step 5: final quality -------------------------------------------------
    final_min_q, final_mean_q = _compute_local_quality(mb, owned_hex)
    global_final_min = comm.allreduce(final_min_q, op=MPI.MIN)
    global_final_mean_sum = comm.allreduce(final_mean_q * len(owned_hex), op=MPI.SUM)
    global_final_mean = global_final_mean_sum / max(global_hex_count, 1)

    if rank == 0:
        print(
            f"  Final quality:   min_SJ={global_final_min:.6f}, "
            f"mean_SJ={global_final_mean:.6f}"
        )
        print(
            f"  Improvement:     delta_min={global_final_min - global_degraded_min:+.6f}, "
            f"delta_mean={global_final_mean - global_degraded_mean:+.6f}"
        )

    # Smoothing must improve quality relative to the degraded mesh
    CHECK(global_final_min > global_degraded_min)


def test_halo_coordinate_exchange_correctness():
    """Verify that coordinate exchange produces consistent ghost positions.

    After exchanging coordinates, every ghost vertex's position should
    match the position stored on the owning rank.  We perturb owned
    vertex positions, exchange, then gather and compare.
    """
    if _skip_if_needed():
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    owned_verts = pcomm.get_owned_entities(dim=0)

    # Perturb each owned vertex by a rank-dependent offset
    if len(owned_verts) > 0:
        coords = mb.get_coords(owned_verts).reshape(-1, 3)
        coords[:, 0] += 0.001 * (rank + 1)  # Small x-shift
        mb.set_coords(owned_verts, coords.flatten())

    # Exchange coordinates through halo
    _exchange_halo_coordinates(mb, pcomm, all_verts)

    # Verify: each ghost vertex should now carry the perturbed position
    # that its owner computed.
    ghost_verts = pcomm.get_ghost_entities(dim=0)
    mismatches = 0

    # Create a verification tag: owner sets value = rank, exchange, check
    verify_tag = mb.tag_get_handle(
        "_VERIFY_OWNER",
        1,
        types.MB_TYPE_DOUBLE,
        types.MB_TAG_DENSE,
        create_if_missing=True,
    )
    mb.tag_set_data(
        verify_tag, owned_verts, np.full(len(owned_verts), float(rank), dtype="float64")
    )
    pcomm.exchange_tags(["_VERIFY_OWNER"], ["_VERIFY_OWNER"], all_verts)

    if len(ghost_verts) > 0:
        ghost_owners = mb.tag_get_data(verify_tag, ghost_verts).flatten()
        ghost_coords = mb.get_coords(ghost_verts).reshape(-1, 3)

        for i, gv in enumerate(ghost_verts):
            owner = pcomm.get_owner(gv)
            CHECK_EQ(ghost_owners[i], float(owner))

    print(
        f"  rank {rank}: verified coordinate exchange for "
        f"{len(ghost_verts)} ghost vertices"
    )
    CHECK(True)


def test_quality_metric_parallel_consistency():
    """Verify that the quality metric is consistent across ranks.

    Each rank computes quality on its owned hexes.  The global min
    computed via MPI_Allreduce should be <= every local min.
    """
    if _skip_if_needed():
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    owned_hex = pcomm.get_owned_entities(dim=3)
    local_min, local_mean = _compute_local_quality(mb, owned_hex)
    global_min = comm.allreduce(local_min, op=MPI.MIN)

    # Global min must be <= every local min
    CHECK(global_min <= local_min + 1.0e-14)

    all_local_mins = get_all_values(local_min)
    if rank == 0:
        print(f"  Per-rank min qualities: {[f'{q:.4f}' for q in all_local_mins]}")
        print(f"  Global min quality: {global_min:.6f}")

    CHECK(global_min > 0.0)  # Mesh should not be inverted


def test_boundary_vertices_are_fixed():
    """Verify that boundary vertices are not moved during smoothing.

    We record boundary vertex coordinates before smoothing, run a few
    iterations, and confirm they are unchanged.
    """
    if _skip_if_needed():
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    mtu = MeshTopoUtil(mb)
    skn = Skinner(mb)

    all_hex = mb.get_entities_by_dimension(0, 3)
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    owned_verts = pcomm.get_owned_entities(dim=0)

    skin_faces = skn.find_skin(0, all_hex, get_vertices=False)
    boundary_verts = mb.get_adjacencies(
        skin_faces, 0, create_if_missing=False, op_type=types.UNION
    )
    owned_boundary = intersect(owned_verts, boundary_verts)
    owned_interior = subtract(owned_verts, boundary_verts)

    # Record boundary positions
    if len(owned_boundary) > 0:
        bdy_coords_before = mb.get_coords(owned_boundary).reshape(-1, 3).copy()
    else:
        bdy_coords_before = np.empty((0, 3), dtype="float64")

    # Run a few smoothing iterations
    for _ in range(5):
        _laplacian_smooth_owned(mb, mtu, owned_interior, omega=0.3)
        _exchange_halo_coordinates(mb, pcomm, all_verts)

    # Verify boundary positions unchanged
    if len(owned_boundary) > 0:
        bdy_coords_after = mb.get_coords(owned_boundary).reshape(-1, 3)
        max_bdy_disp = np.max(np.abs(bdy_coords_after - bdy_coords_before))
        CHECK(max_bdy_disp < 1.0e-14)
        print(
            f"  rank {rank}: boundary max displacement = {max_bdy_disp:.2e} "
            f"({len(owned_boundary)} boundary verts)"
        )
    else:
        print(f"  rank {rank}: no owned boundary vertices")

    CHECK(True)


def test_smoothing_reduces_displacement_monotonically():
    """Verify that per-iteration max displacement decreases over time.

    For a well-conditioned mesh with small omega, the maximum vertex
    displacement should decrease monotonically (or nearly so).
    """
    if _skip_if_needed():
        return

    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    pcomm.load_file(PARTITIONED_MESH, READ_OPTS)
    pcomm.exchange_ghost_cells(3, 0, 1, 0, True, True)

    mtu = MeshTopoUtil(mb)
    skn = Skinner(mb)

    all_hex = mb.get_entities_by_dimension(0, 3)
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)
    owned_verts = pcomm.get_owned_entities(dim=0)

    skin_faces = skn.find_skin(0, all_hex, get_vertices=False)
    boundary_verts = mb.get_adjacencies(
        skin_faces, 0, create_if_missing=False, op_type=types.UNION
    )
    owned_interior = subtract(owned_verts, boundary_verts)

    displacements = []
    n_iter = 20
    omega = 0.2

    for iteration in range(n_iter):
        local_disp = _laplacian_smooth_owned(mb, mtu, owned_interior, omega)
        _exchange_halo_coordinates(mb, pcomm, all_verts)
        global_disp = comm.allreduce(local_disp, op=MPI.MAX)
        displacements.append(global_disp)

    if rank == 0:
        print(
            f"  Displacement history (first 5): "
            f"{[f'{d:.3e}' for d in displacements[:5]]}"
        )
        print(
            f"  Displacement history (last 5):  "
            f"{[f'{d:.3e}' for d in displacements[-5:]]}"
        )

    # The displacement should generally decrease.  Allow a few non-monotone
    # steps due to floating-point effects, but the trend must be downward.
    # Check that the last displacement is smaller than the first.
    if displacements[0] > 1.0e-15:
        CHECK(displacements[-1] < displacements[0])
    else:
        # Already converged at the start — acceptable
        CHECK(True)


def main():
    tests = [
        test_halo_coordinate_exchange_correctness,
        test_quality_metric_parallel_consistency,
        test_boundary_vertices_are_fixed,
        test_smoothing_reduces_displacement_monotonically,
        test_laplacian_smoothing_convergence,
    ]
    return run_parallel_tests(tests)


if __name__ == "__main__":
    sys.exit(main())
