# compareMaps {#tools_compareMaps}

Compare two remap-weight map files in NetCDF (SCRIP) format. Loads each file's sparse weight matrix into an Eigen sparse representation, then reports differences in `row`, `col`, `S` entries and matrix norms. Useful for validating that a refactored remap workflow still produces the same matrix, or for diffing two implementations.

> The other map-file fields — `xc`, `yc`, `area`, `frac` — are best diffed with `ncdiff` from [NCO](http://nco.sourceforge.net/). `compareMaps` focuses on the sparse-matrix entries.

## Build requirements

Requires both **NetCDF** and **Eigen3**. Build via CMake or autotools; the tool is gated by `MOAB_HAVE_EIGEN3` (set by `-DENABLE_EIGEN3=ON` in CMake).

## Synopsis

```
mbcmpmaps -i <map1.nc> -j <map2.nc> [options]
```

> Binary name on disk is `mbcmpmaps`.

## Options

| Flag | Default | Description |
|---|---|---|
| `-i <map1.nc>`, `--firstMap <map1.nc>` | — | First input map file (SCRIP weights). |
| `-j <map2.nc>`, `--secondMap <map2.nc>` | — | Second input map file. |
| `-p <n>`, `--print_differences <n>` | 0 | Print up to `<n>` per-entry differences to stdout. Useful for spot-checking which weights diverged. |
| `-t <eps>`, `--tolerance <eps>` | — | Pass/fail threshold. The tool exits with code `1` if any matrix-norm difference exceeds `<eps>`; `0` otherwise. If unset, the tool only reports diagnostics and always exits 0. |

## What it reports

For each map file, the tool extracts the dimensions `n_a`, `n_b`, `n_s` and the `row[]`, `col[]`, `S[]` variables; it constructs an `Eigen::SparseMatrix<double>` of size `n_b × n_a` from them. Then between the two matrices it reports:

- **Structural agreement**: do they have the same shape and the same non-zero pattern?
- **Norms of the difference**: max absolute entry difference, Frobenius norm of the difference, etc.
- **Per-entry samples**: with `-p N`, prints up to `N` entries that differ, ordered by magnitude of the difference.

## Examples

```bash
# Diagnostic-only run
mbcmpmaps -i map1.nc -j map2.nc

# Show the 50 largest per-entry differences
mbcmpmaps -i map1.nc -j map2.nc -p 50

# CI gate: fail if any matrix-norm difference exceeds 1e-12
mbcmpmaps -i baseline.nc -j candidate.nc -t 1e-12
echo $?      # 0 if within tolerance, 1 if exceeded
```

## See also

- `ncdiff` (from NCO) — for diffing `xc / yc / area / frac` outside the sparse matrix
- [`mbtempest`](mbtempest.md) — generates map files in this same format
