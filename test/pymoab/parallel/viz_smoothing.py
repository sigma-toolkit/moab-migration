#!/usr/bin/env python3
"""
Generate an animated GIF of Laplacian mesh smoothing on a z-midplane slice.

This is a serial visualization companion to test_parallel_mesh_optimization.py.
It loads the mesh on a single rank, perturbs interior vertices, runs
quality-guarded Laplacian smoothing, and captures a frame every few
iterations showing the quad wireframe colored by scaled-Jacobian quality.

Usage:
    python3 viz_smoothing.py                     # defaults
    python3 viz_smoothing.py --frames 30 -o smoothing.gif

Requires: matplotlib (with pillow for GIF writing).
"""

import argparse
import os
import sys
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from pymoab import core, types
from pymoab.rng import Range, subtract
from pymoab.topo_util import MeshTopoUtil
from pymoab.skinner import Skinner


MESH_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "MeshFiles", "unittest"
)
PARTITIONED_MESH = os.path.join(MESH_DIR, "64bricks_512hex_256part.h5m")


# ---------------------------------------------------------------------------
# Reuse quality / smoothing kernels from the test (serial versions)
# ---------------------------------------------------------------------------


def _hex_scaled_jacobian(c):
    corners = [
        (0, (1, 3, 4)),
        (1, (2, 0, 5)),
        (2, (3, 1, 6)),
        (3, (0, 2, 7)),
        (4, (7, 5, 0)),
        (5, (4, 6, 1)),
        (6, (5, 7, 2)),
        (7, (6, 4, 3)),
    ]
    worst = 1.0
    for n, (a, b, d) in corners:
        e1, e2, e3 = c[a] - c[n], c[b] - c[n], c[d] - c[n]
        det = np.dot(e1, np.cross(e2, e3))
        denom = np.linalg.norm(e1) * np.linalg.norm(e2) * np.linalg.norm(e3)
        sj = det / denom if denom > 1e-15 else -1.0
        if sj < worst:
            worst = sj
    return worst


def _quality_of_hex(mb, h):
    conn = mb.get_connectivity(h)
    return _hex_scaled_jacobian(mb.get_coords(conn).reshape(-1, 3))


def _min_quality_around_vertex(mb, v):
    adj = mb.get_adjacencies(Range([v]), 3, create_if_missing=False)
    if len(adj) == 0:
        return 1.0
    return min(_quality_of_hex(mb, h) for h in adj)


def _smooth_once(mb, mtu, interior, omega):
    max_disp = 0.0
    for v in interior:
        nbrs = mtu.get_bridge_adjacencies(v, bridge_dim=1, to_dim=0)
        if len(nbrs) == 0:
            continue
        centroid = mb.get_coords(nbrs).reshape(-1, 3).mean(axis=0)
        vr = Range([v])
        old = mb.get_coords(vr).reshape(3)
        new = (1.0 - omega) * old + omega * centroid
        q_before = _min_quality_around_vertex(mb, v)
        mb.set_coords(vr, new)
        if _min_quality_around_vertex(mb, v) < q_before:
            mb.set_coords(vr, old)
        else:
            d = np.linalg.norm(new - old)
            if d > max_disp:
                max_disp = d
    return max_disp


# ---------------------------------------------------------------------------
# Slice extraction: find hexes whose centroid is near z=z_slice, draw the
# quad cross-section as the 4 vertices closest to z_slice.
# ---------------------------------------------------------------------------


def _extract_z_slice(mb, hexes, z_val, tol=0.3):
    """Return list of (quad_xy, quality) for hexes crossing z=z_val."""
    quads = []
    for h in hexes:
        conn = mb.get_connectivity(h)
        coords = mb.get_coords(conn).reshape(-1, 3)
        zc = coords[:, 2]
        if zc.min() > z_val + tol or zc.max() < z_val - tol:
            continue
        # 4 vertices nearest z_val
        dist = np.abs(zc - z_val)
        idx = np.argsort(dist)[:4]
        xy = coords[idx, :2]
        # Sort by angle around centroid for proper polygon ordering
        cx, cy = xy.mean(axis=0)
        angles = np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx)
        order = np.argsort(angles)
        quads.append((xy[order], _hex_scaled_jacobian(coords)))
    return quads


def _render_frame(quads, ax, vmin, vmax, title):
    ax.clear()
    cmap = plt.cm.RdYlGn
    norm = Normalize(vmin=vmin, vmax=vmax)
    patches = []
    colors = []
    for xy, q in quads:
        patches.append(Polygon(xy, closed=True))
        colors.append(cmap(norm(q)))
    coll = mc.PatchCollection(
        patches, edgecolors="k", linewidths=0.4, facecolors=colors
    )
    ax.add_collection(coll)
    all_xy = np.vstack([q[0] for q in quads])
    margin = 0.2
    ax.set_xlim(all_xy[:, 0].min() - margin, all_xy[:, 0].max() + margin)
    ax.set_ylim(all_xy[:, 1].min() - margin, all_xy[:, 1].max() + margin)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        default="smoothing.gif",
        help="Output GIF path (default: smoothing.gif)",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=40,
        help="Total frames in the animation (default: 40)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=80,
        help="Max smoothing iterations (default: 80)",
    )
    parser.add_argument(
        "--omega", type=float, default=0.5, help="Relaxation factor (default: 0.5)"
    )
    parser.add_argument(
        "--amplitude",
        type=float,
        default=0.15,
        help="Perturbation amplitude (default: 0.15)",
    )
    parser.add_argument(
        "--z-slice",
        type=float,
        default=0.0,
        help="Z coordinate for the cross-section (default: 0.0)",
    )
    args = parser.parse_args()

    if not os.path.exists(PARTITIONED_MESH):
        print(f"Mesh not found: {PARTITIONED_MESH}")
        sys.exit(1)

    # Load mesh serially
    mb = core.Core()
    mb.load_file(PARTITIONED_MESH)
    mtu = MeshTopoUtil(mb)
    skn = Skinner(mb)

    hexes = mb.get_entities_by_type(0, types.MBHEX)
    all_verts = mb.get_entities_by_type(0, types.MBVERTEX)

    skin_faces = skn.find_skin(0, hexes, get_vertices=False)
    boundary_verts = mb.get_adjacencies(
        skin_faces, 0, create_if_missing=False, op_type=types.UNION
    )
    interior = subtract(all_verts, boundary_verts)

    print(
        f"Mesh: {len(hexes)} hexes, {len(all_verts)} verts, "
        f"{len(interior)} interior verts"
    )

    # Perturb
    rng = np.random.RandomState(42)
    if len(interior) > 0:
        coords = mb.get_coords(interior).reshape(-1, 3)
        bbox_diag = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
        mean_edge = bbox_diag / max(len(interior) ** (1.0 / 3.0), 1.0)
        pert = rng.uniform(-1, 1, size=coords.shape)
        norms = np.linalg.norm(pert, axis=1, keepdims=True)
        norms = np.where(norms > 1e-15, norms, 1.0)
        coords += pert / norms * args.amplitude * mean_edge
        mb.set_coords(interior, coords.flatten())

    # Decide which iterations to capture as frames
    iters_per_frame = max(1, args.max_iter // args.frames)
    capture_iters = set(range(0, args.max_iter, iters_per_frame))
    capture_iters.add(0)

    # Quality color range
    vmin, vmax = 0.3, 1.0

    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=plt.cm.RdYlGn)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8, label="Scaled Jacobian")

    frames = []
    print(
        f"Running {args.max_iter} smoothing iterations, "
        f"capturing ~{len(capture_iters)} frames..."
    )

    # Frame 0: degraded mesh
    quads = _extract_z_slice(mb, hexes, args.z_slice)
    qvals = [q for _, q in quads]
    _render_frame(
        quads,
        ax,
        vmin,
        vmax,
        f"iter 0  |  min SJ = {min(qvals):.3f}  mean = {np.mean(qvals):.3f}",
    )
    fig.canvas.draw()
    frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))
    frames.append(frame.copy())

    for it in range(1, args.max_iter + 1):
        disp = _smooth_once(mb, mtu, interior, args.omega)
        if it in capture_iters or disp < 1e-8:
            quads = _extract_z_slice(mb, hexes, args.z_slice)
            qvals = [q for _, q in quads]
            _render_frame(
                quads,
                ax,
                vmin,
                vmax,
                f"iter {it}  |  min SJ = {min(qvals):.3f}  mean = {np.mean(qvals):.3f}",
            )
            fig.canvas.draw()
            frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
            frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))
            frames.append(frame.copy())
            print(
                f"  frame {len(frames):3d}  iter {it:3d}  "
                f"disp={disp:.3e}  min_SJ={min(qvals):.4f}"
            )
        if disp < 1e-8:
            print(f"Converged at iteration {it}")
            break

    plt.close(fig)

    # Hold last frame a bit longer
    for _ in range(8):
        frames.append(frames[-1])

    # Write GIF via pillow (installed as matplotlib dependency)
    from PIL import Image

    pil_frames = [Image.fromarray(f[:, :, :3]) for f in frames]
    pil_frames[0].save(
        args.output,
        save_all=True,
        append_images=pil_frames[1:],
        duration=150,
        loop=0,
    )
    print(
        f"\nWrote {len(frames)} frames to {args.output} "
        f"({os.path.getsize(args.output) / 1024:.0f} KB)"
    )


if __name__ == "__main__":
    main()
