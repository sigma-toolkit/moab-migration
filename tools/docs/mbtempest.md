# mbtempest

A unified command-line tool for **climate-style mesh generation** and **conservative remapping weight computation**, built on top of MOAB and the TempestRemap library. Used heavily in the E3SM coupler workflow and as the offline counterpart to the iMOAB remap routines.

## What it does

`mbtempest` operates in a few distinct modes selected via `--type`:

- **Generate a sphere mesh from scratch** — Cubed-Sphere (CS), Regular Lat-Lon (RLL), or Icosahedral (ICO).
- **Compute the overlap (intersection) mesh** of two sphere meshes and optionally save it.
- **Compute conservative remapping weights** between a source and target mesh, with optional high-order, monotonic, and CAAS (Clip-And-Assured-Sum) variants.

Output is written in formats appropriate to the mode: `.h5m` / `.exo` / `.nc` for meshes, `.nc` (SCRIP) for remap-weight files.

## Quick-start examples

```bash
# Cubed-Sphere mesh, resolution 25
mbtempest --type 0 --res 25 --file cs25.h5m

# RLL mesh, 90 × 180 (lon × lat)
mbtempest --type 1 --res 90 --file rll90x180.h5m

# Icosahedral dual mesh
mbtempest --type 2 --res 25 --dual --file ico25_dual.h5m

# Overlap mesh of two existing sphere meshes
mbtempest --type 5 --load source.h5m --load target.h5m --intx intx.h5m

# FV → FV conservative weight map (defaults)
mbtempest --type 5 --load source.h5m --load target.h5m --file weights_fv_fv.nc

# SE → FV weights, explicitly: cgll order-4 source → fv order-1 target
mbtempest --type 5 --load source_se.h5m --load target_fv.h5m \
          --order 4 --method cgll --global_id GLOBAL_DOFS \
          --order 1 --method fv   --global_id GLOBAL_ID    \
          --file weights_se4_to_fv1.nc
```

Run `mbtempest --help` for the full option list; `mbtempest --manual` prints longer-form documentation with more examples.

## Mode reference (`-t`/`--type`)

| Value | Mode | Meaning |
|---|---|---|
| 0 | `CS` | Generate a Cubed-Sphere mesh. `--res N` controls the per-face resolution (N × N quads per panel). |
| 1 | `RLL` | Generate a Regular Lat-Lon mesh. `--res N` is the half-resolution; the mesh is `2N × N` (lon × lat). |
| 2 | `ICO` | Generate an Icosahedral mesh. Add `--dual` for the hex-pentagon dual. |
| 3 | `OVERLAP_FILES` | Compute overlap from two on-disk meshes via TempestRemap's file path. |
| 4 | `OVERLAP_MEMORY` | Compute overlap, holding the source/target in memory (avoid file I/O round-trip). |
| 5 | `OVERLAP_MOAB` | Compute overlap using MOAB-native intersection; recommended default for parallel runs and the path used by E3SM. |

For modes 3–5, you must supply both meshes with `--load src.h5m --load tgt.h5m`.

## Options

### Mesh generation

| Flag | Type | Default | Description |
|---|---|---|---|
| `-t`, `--type <0..5>` | int | 0 (CS) | Mesh type selector — see table above. |
| `-r`, `--res <N>` | int | 5 | Mesh resolution (mode-dependent meaning). |
| `-d`, `--dual` | flag | off | Output the dual mesh (ICO only). |
| `-f`, `--file <path>` | string | — | Output filename for the generated mesh or weights. |

### Inputs for overlap / remap workflows

| Flag | Type | Description |
|---|---|---|
| `-l`, `--load <path>` | string (repeatable) | Input mesh filenames. For overlap/weights modes, give exactly two: first is source, second is target. |
| `-i`, `--intx <path>` | string | Output filename for the intersection mesh (useful for caching). |
| `--baseline <path>` | string | Output baseline file (regression-test artifact). |

### Discretization

| Flag | Type | Default | Description |
|---|---|---|---|
| `-m`, `--method <name>` | string (repeatable) | `fv` | Discretization method for source then target. One of `fv`, `cgll`, `dgll`, `pcloud`. Give the flag twice (one per side) for asymmetric pairs. |
| `-o`, `--order <N>` | int (repeatable) | `1` | Discretization order for source then target. FV is order 1; SE/spectral methods take any `1..N`. Give twice for asymmetric pairs. |
| `-g`, `--global_id <tag>` | string (repeatable) | `GLOBAL_ID` | Tag name carrying global DoF IDs for source then target. Spectral element meshes typically use `GLOBAL_DOFS` instead of `GLOBAL_ID`. |
| `--fvmethod <name>` | string | `none` | FV-FV sub-method. One of `invdist`, `delaunay`, `bilin`, `intbilin`, `intbilingb`, `none`. |
| `--nobubble` | flag | off | Do not use bubble functions on interior of spectral element nodes. |

### Intersection / coverage tuning

| Flag | Type | Default | Description |
|---|---|---|---|
| `-a`, `--advfront` | flag | off | Use advancing-front intersection instead of the Kd-tree algorithm. |
| `--gnomonic` | flag | off | Project to gnomonic plane when computing the coverage mesh. |
| `--enforce_convexity` | flag | off | Check convexity of input meshes before intersection. |
| `--rrmgrids` | flag | off | At least one input is a regionally-refined mesh (RRM); enables an accelerated intersection path. |
| `--ghost <N>` | int | auto | Number of ghost layers in the coverage mesh. Auto-selects `0` for FV order 1, `p+1` for FV order `p > 1`. Override only if you know you need more. |
| `--boxeps <eps>` | double | `1e-7` | Tolerance used for bounding-box checks during the intersection. |

### Weights & filters

| Flag | Type | Default | Description |
|---|---|---|---|
| `-w`, `--weights` | flag | off | Compute and emit the remap-weights file alongside the overlap. |
| `--noconserve` | flag | off | Skip the post-projection conservation enforcement. |
| `--volumetric` | flag | off | Use the volumetric-projection variant when computing weights. |
| `--monotonicity <0\|1\|2\|3>` | int | 0 | Monotonicity constraint: 0 = none, 1 = local, 2 = global, 3 = strict CAAS. |
| `--limiter <N>` | int | 0 | Apply a non-linear post-filter (CAAS-style) after the linear map application. |
| `--sparseconstraints` | flag | off | Use a sparse solver for constraints — recommended when the source/target has high-valence cells (typical of high-res RLL). |

### Validation / verification

| Flag | Description |
|---|---|
| `--checkmap` | After computing the weights, verify conservation and consistency on the resulting map. |
| `--verify` | Project analytical reference functions through the map and report error metrics. |
| `--var <tagname>` | Tag name of the variable to use in the verification study. Built-in metrics are computed for the standard analytical fields; user-defined variables may not have error metrics. |

### Performance / I/O / misc

| Flag | Description |
|---|---|
| `--skip_intersection` | Skip the mesh-intersection step (e.g. when you already have a cached `--intx` file). |
| `--skip_output` | Skip all I/O (for performance studies). |
| `-v`, `--verbose` | Verbose diagnostic output during intersection and map computation. |
| `--manual` | Print extended documentation and examples and exit. |
| `--version` | Print version info and exit. |

## Discretization-method cheat sheet

| Source | Target | Typical CLI |
|---|---|---|
| FV (order 1) | FV (order 1) | `--method fv --order 1` *(default)* |
| SE / CGLL (order p) | FV (order 1) | `--method cgll --order p --global_id GLOBAL_DOFS` source, `--method fv --order 1` target |
| FV (order 1) | SE / CGLL (order p) | symmetric of the above |
| DGLL (order p) | DGLL (order q) | `--method dgll --order p` source, `--method dgll --order q` target |
| Point cloud | FV | `--method pcloud` source, `--method fv --order 1` target |

When pairing two different discretizations you must give `--method`, `--order`, and `--global_id` **twice**: the first occurrence applies to the source, the second to the target.

## Recommended monotonicity / filter combinations

| Goal | Recipe |
|---|---|
| Plain conservative high-order (no monotonicity) | `--order 4 --method cgll --monotonicity 0` |
| Bounded high-order with CAAS clipping at apply time | `--order 4 --monotonicity 1 --limiter 2` |
| MCT-style nonlinear "dual map" CAAS, BfB across rank counts | Compute a low-order monotone map and a high-order map separately, then apply with `--limiter 2`. The iMOAB driver applies these as a `dualMapPID` — see `test/parallel/imoab_dualmap_caas.cpp`. |

## Parallel execution

Just prefix with an MPI launcher:

```bash
mpiexec -n 8  mbtempest --type 5 \
    --load source.h5m --load target.h5m \
    --weights --file weights.nc
```

The intersection, weight computation, and output write all parallelize via MPI. For very large grids combine with:

```bash
mbtempest --type 5 --load ... --load ... \
    --weights --file w.nc \
    --rrmgrids --sparseconstraints --ghost 3 \
    --intx /tmp/intx_cache.h5m
```

`--intx` caches the intersection so subsequent runs (e.g. trying different weight options) can use `--skip_intersection`.

## Output formats

- Mesh generation: `.h5m` (MOAB native, recommended), `.exo` (Exodus), `.nc` (NetCDF SCRIP grid).
- Overlap meshes: same options as above.
- Weight maps: `.nc` (SCRIP). Reload with [`compareMaps`](compareMaps.md) for diffing, or via `iMOAB_LoadMapFile` in the E3SM coupler.

## See also

- [`mbIntxCheck`](mbIntxCheck.md) — verify a generated intersection mesh
- [`compareMaps`](compareMaps.md) — diff two weight maps in SCRIP format
- [`mbconvert`](mbconvert.md) — for converting outputs to other formats (e.g. emitting a generated CS mesh as SCRIP via `-o WRITE_FORMAT=SCRIP`)
- TempestRemap upstream: <https://github.com/ClimateGlobalChange/tempestremap>
