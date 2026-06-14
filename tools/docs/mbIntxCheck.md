# mbIntxCheck

Verify a precomputed mesh intersection against its source and target meshes. After a parallel run of `mbtempest` (or any other intersection generator) produces an intersection mesh, this tool checks that the per-element areas balance correctly between source / target / intersection.

For each source element, sum the areas of all intersection polygons that claim it as their source parent — that sum should equal the source element's own area, to within a configurable tolerance. Same check on the target side. Mismatches are written as element-level tags to per-side verification files.

## Synopsis

```
mpiexec -n N  mbintx_check  -s <source.h5m>  -t <target.h5m>  -i <intx.h5m> \
                            -v <source_verif.h5m> -w <target_verif.h5m> \
                            [options]
```

> Binary name on disk is `mbintx_check`.

## Options

| Flag | Default | Description |
|---|---|---|
| `-s <file>`, `--source <file>` | — | Source mesh (parent on one side of the intersection). |
| `-t <file>`, `--target <file>` | — | Target mesh (parent on the other side). |
| `-i <file>`, `--intersection <file>` | — | Precomputed intersection mesh. Each element carries `srcParent` / `tgtParent` tags pointing back into the source/target. |
| `-v <file>`, `--verif_source <file>` | — | Output: source mesh with per-cell `area_error` tag written. |
| `-w <file>`, `--verif_target <file>` | — | Output: target mesh with per-cell `area_error` tag written. |
| `-m <eps>`, `--threshold_source <eps>` | small | Pass/fail threshold for the per-cell source-area error. |
| `-q <eps>`, `--threshold_target <eps>` | small | Pass/fail threshold for the per-cell target-area error. |
| `-p {0\|1}`, `--sphere <flag>` | 1 | Mesh sits on a sphere (`1`) or in a plane (`0`); affects how areas are computed. |
| `-n {0\|1}`, `--old_convention <flag>` | 0 | Use the legacy `parent` tag names (pre-current convention). Set this if the intersection file was produced by an older mbtempest. |

## What it does

1. Reads the three meshes in parallel.
2. For each intersection element `i`:
   - Looks up its source parent `s = srcParent[i]` and target parent `t = tgtParent[i]`.
   - Accumulates `area(i)` into the running totals `src_sum[s]` and `tgt_sum[t]`.
3. After the gather, computes per-cell errors:
   ```
   src_error[s] = area(s) − src_sum[s]
   tgt_error[t] = area(t) − tgt_sum[t]
   ```
4. Writes the errors as element-level tags on the corresponding source / target mesh and saves them to the verification files.
5. Prints summary statistics (max error, number of cells above threshold) per side.

If the intersection is correct, both errors are zero to floating-point precision. Non-zero errors indicate missing or duplicated intersection polygons.

## Example

After a parallel mbtempest run that produced `intx.h5m`:

```bash
mpiexec -n 4  mbintx_check \
  -s source.h5m  -t target.h5m  -i intx.h5m \
  -v source_verif.h5m  -w target_verif.h5m \
  -m 1e-10  -q 1e-10
```

Open `source_verif.h5m` and `target_verif.h5m` in ParaView/VisIt and color by the `area_error` tag to spot regions where the intersection is incomplete.

## Sphere vs plane

The `-p` flag selects the geometric model:

- `-p 1` (default): areas are computed as spherical-polygon areas on the unit sphere. Use for climate / atmosphere / ocean grids.
- `-p 0`: areas are computed as planar-polygon areas. Use for 2-D regional or Cartesian meshes.

## See also

- [`mbtempest`](mbtempest.md) — generates the intersection meshes that `mbIntxCheck` verifies
- [`compareMaps`](compareMaps.md) — diff two remap weight files (different goal: matrix equivalence vs. area consistency)
