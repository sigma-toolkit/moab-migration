# MOAB + Eigen3 parallel mesh optimization examples

Date: 2026-09-15

## Purpose

Replace `examples/advanced/smoothing/OptimizeMeshMesquite.cpp`, which drives
Mesquite through the iMesh bindings, with examples that depend on nothing but
MOAB and Eigen3. The iTAPS bindings are not in this source tree and no build
system here detects Mesquite, so that example has never been buildable.

Two examples come out of this:

1. `MeshOptimization` — shape optimization of quad/tri surface meshes on a
   plane and hex/tet volume meshes in a box, in parallel.
2. `SphericalCVTOptimization` — generate a centroidal Voronoi tessellation on
   an analytic sphere, then optimize its Delaunay triangulation. This is the
   MPAS/E3SM mesh family, so it doubles as an earth-system use case.

Both share a templated metric header and a parallel optimizer core.

## Non-goals

- No Mesquite, no iMesh, no MsqIGeom, no CGM.
- No general CAD geometry. Domains are analytic: plane, axis-aligned box,
  sphere.
- No Hessian. First-order AutoDiff plus L-BFGS. Nesting `AutoDiffScalar` would
  give exact Hessians, but a Newton solve needs a distributed linear solver,
  which is a larger example than this should be.

## Components

### `MeshOptimizationMetrics.hpp`

Header-only, templated on the scalar type, no MOAB dependency beyond
`EntityType`. It must compile under `double`, under
`Eigen::AutoDiffScalar`, and under the finite-difference checker.

Per element, build the corner Jacobian `T` against the ideal-element weight
`W`, then:

- `inverse_mean_ratio = |T|_F^2 / (n * reg_det(T)^(2/n))`
- `condition_number   = |T|_F * |T^-1|_F / n`

with `n` the topological dimension.

**Escobar regularization.** `reg_det(d) = (d + sqrt(d^2 + 4*delta^2)) / 2`.
This keeps the metric finite and differentiable through element inversion, so
a single L-BFGS pass can untangle and then improve shape. At `delta = 0` it
reduces exactly to `d` for `d > 0`, which is what makes validation against
Verdict possible on a valid mesh.

**Sampling.** Quad and hex are sampled at every corner (4 and 8) and combined
by RMS. Not by max: L-BFGS needs a differentiable objective and `max` is not.

**Known identity.** For a 2x2 Jacobian the adjugate gives
`|T^-1|_F = |T|_F / det T`, so `cond = |T|_F^2 / (2 det T)`, which is the
inverse mean ratio. The two metrics coincide for tri and quad and differ only
for tet and hex. This is documented in the header rather than hidden.

**Surface elements in 3D.** Use the metric tensor `M = T^T T`, with
`det := sqrt(det M)` and `|T|_F^2 := tr M`. No local frame construction, so
AutoDiff propagates cleanly.

### `MeshOptimizerCore.hpp`

Parallel objective, gradient, constraints and L-BFGS.

**Objective.** `F(x) = (1/N) * sum over elements of metric_e(x)^p`, `p`
default 2. The mean rather than the sum, so the value does not depend on rank
count or mesh size.

**Gradient assembly.** Each element is evaluated exactly once, by its owner:

1. Zero a 3-component `GRADIENT` tag on all local vertices.
2. Loop owned elements. `AutoDiffScalar` gives the element gradient. Scatter-add
   into incident vertices, ghosts included.
3. `pcomm->reduce_tags(grad, MPI_SUM, shared_verts)` — owners hold the full sum.
4. `pcomm->exchange_tags(grad, shared_verts)` — sharers pick up the owner's
   value. Done unconditionally rather than relying on `reduce_tags` to
   back-propagate.
5. `MPI_Allreduce` the scalar objective and the element count.

L-BFGS state lives only on owned free vertices. After each step,
`exchange_tags(coords)` republishes positions before the next evaluation.

**Constraints.** Each vertex carries a symmetric 3x3 projector `P`. The search
direction and the gradient are both premultiplied by `P`, so a vertex moves
only within its allowed subspace.

`P = I - sum_i v_i v_i^T` over an orthonormal basis of the constrained
directions. Ranks:

| constrained directions | meaning | DOF |
|---|---|---|
| 0 | interior of a volume | 3 |
| 1 | on a plane or box face | 2 |
| 2 | on a box edge | 1 |
| 3 | box corner | 0 |

The projector is built analytically from the domain, **not** from skinning.
Parallel skinning reports partition-interface facets as boundary, and the
workarounds are fragile. Since the domain is a box by assumption, a vertex's
constraints follow from which of the bounding planes it lies on, within a
tolerance scaled to the mesh. That is purely local and therefore identical on
every rank with no communication.

A planar surface mesh gets the plane normal as a constraint for *every*
vertex, so the same machinery gives interior vertices 2 DOF and boundary
vertices 1.

For the sphere, `P = I - n n^T` with `n` radial, recomputed each iteration,
plus an explicit position re-projection onto the sphere after each step.

**L-BFGS.** Two-loop recursion, 7 history pairs, backtracking Armijo line
search. Every inner product is an `MPI_Allreduce` over owned DOFs. Armijo
alone does not guarantee the curvature condition, so pairs with `y^T s <= 0`
are skipped rather than silently corrupting the history.

### `MeshOptimization.cpp`

**Self-generating test problem.** `ScdInterface::construct_box` with
`ScdParData` (`SQIJ` for the plane, `SQIJK` for the box) and
`resolve_shared_ents` builds a distributed, sharing-resolved quad or hex mesh
with no input file and no HDF5. Interior vertices are then perturbed with a
seeded RNG and optimized back.

This matters for verification: the structured lattice **is** the exact
optimum of both metrics on a plane or box, so the example has a known answer.
The test is "quality returns to 1.0 and vertices return to the lattice", not
merely "it ran".

`-f <file>` still reads a tri/tet mesh. There the constraint planes come from
the global bounding box, which assumes a box-shaped domain; the help text says
so.

**Reporting.** `VerdictWrapper` before/after on `MB_CONDITION`, `MB_SHAPE`,
`MB_SCALED_JACOBIAN`, reduced across ranks. Independent of the templated
metrics, so it cross-checks them.

**Verification modes.**

- `--verify-metric`: templated condition number at `delta = 0` against
  `VerdictWrapper::quality_measure`, element by element. Catches a wrong metric.
- `--verify-gradient`: AutoDiff gradient against central finite differences.
  Catches a wrong derivative.

### `SphericalCVTOptimization.cpp`

1. Generate N generators on the unit sphere (Fibonacci spiral, or seeded
   random).
2. Spherical Delaunay triangulation via 3D convex hull. Every point on a
   sphere is a hull vertex, so an incremental insertion with visible-face
   removal suffices; MOAB has no hull utility, so this example carries its own.
3. Lloyd iteration to a CVT: each generator moves to the normalized centroid
   of its spherical Voronoi cell, which is the dual of the Delaunay
   triangulation. Retriangulate each sweep.
4. Emit both the Delaunay triangulation and the dual Voronoi polygonal mesh.
   The latter is the MPAS mesh form.
5. Run the shared optimizer on the triangulation with the spherical projector.

The CVT construction is serial (the hull is global); the optimization phase
runs in parallel after the mesh is distributed.

## Testing

| case | mesh | ranks | check |
|---|---|---|---|
| plane-quad | generated, perturbed | 1, 4 | quality -> 1.0, verts -> lattice |
| box-hex | generated, perturbed | 1, 8 | same |
| gradient | generated | 1 | AutoDiff vs central FD |
| metric | generated | 1 | delta=0 vs VerdictWrapper |
| determinism | generated | 1 vs 4 | identical objective history |
| sphere-cvt | generated | 1 | CVT energy decreases, quality improves |

The determinism row is the important one. A rank-dependent objective would
mean the gradient assembly is double-counting or dropping contributions, which
is the defect this design is most exposed to.

## Build wiring

- Both examples gate on MPI and Eigen3. Neither `moab.make` template currently
  defines an Eigen flag, so `MOAB_EIGEN3_ENABLED` is added to both —
  `@enableeigen@` for autotools, `@MOAB_MAKE_EIGEN3_ENABLED@` for CMake.
- `examples/CMakeLists.txt` does not reference `advanced/smoothing` at all, so
  the installed CMake examples tree builds none of it. Add the directory.
- `OptimizeMeshMesquite.cpp` is deleted from the tree, from `Makefile.am` and
  from the `MOAB_EXAMPLES` list in `CMakeLists.txt`.
