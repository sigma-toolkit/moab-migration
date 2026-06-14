# mbconvert

Convert mesh files between formats supported by MOAB, optionally extracting subsets of the input.

## Synopsis

```
mbconvert [-h|-l]

mbconvert [-f format] [-a sat_file|-A] [-t] [-g]
          [-o option] [-O option] [-I dim] [-p|-P]
          [-v vol_ids] [-s surf_ids] [-c curve_ids] [-V vert_ids]
          [-m block_ids] [-d nodeset_ids] [-n sideset_ids] [-D part_ids]
          [-1] [-2] [-3]
          input_file [input_file2 ...] output_file
```

## Description

`mbconvert` uses the MOAB library to translate between mesh file formats or to extract subsets of a mesh file. The output file type is determined from the file extension when `-f` is not given. Use `-l` to list supported file formats.

`id_list` arguments accept individual IDs or ranges separated by commas (e.g. `1,5,10-15,20`). Lists may not contain spaces.

Any combination of subset options can be specified (repeating a flag with a different id-list is allowed); the output contains the union. With no subset options the entire input mesh is written.

Writer- and reader-specific option strings for the `-o` and `-O` flags are documented in [`../README.IO`](../README.IO).

## Options

### General

| Flag | Description |
|---|---|
| `-h` | Print help to stdout and exit. |
| `-l` | List supported file formats and exit. |
| `-f <format>` | Specify the output file format explicitly. |
| `-t` | Print read/write timing data. |
| `-g` | Enable verbose / debug output. |
| `-I <dim>` | Generate internal / adjacent entities of dimension `<dim>` (1 or 2). |
| `--` | Treat all subsequent arguments as file names (lets file names begin with `-`). |

### Format-specific

| Flag | Description |
|---|---|
| `-a <sat_file>` | ACIS SAT file dumped by the `.cub` reader. Equivalent to `-O SAT_FILE=<sat_file>`. |
| `-A` | Tell the `.cub` reader to **not** dump a SAT file (deprecated; current default). |

### I/O options

| Flag | Description |
|---|---|
| `-o <option>` | Pass a writer option. Repeat for multiple options. |
| `-O <option>` | Pass a reader option. Repeat for multiple options. |

> **Note:** lowercase `-o` is for **writes**; uppercase `-O` is for **reads**. (The older `mbconvert.man` page had these descriptions swapped; this doc and the in-binary help have it right.)

The most commonly used writer options:

- `WRITE_FORMAT={SCRIP|ESMF|DOMAIN}` — Force NetCDF output into a specific grid layout, regardless of the source mesh's `__MESH_TYPE` tag. See [Cross-format NetCDF conversion](#cross-format-netcdf-conversion) below.
- `PARALLEL=WRITE_PART` — Set automatically by `-M`; rarely needed by hand.

### Parallel (requires MPI build)

| Flag | Description |
|---|---|
| `-P` | Append `.<rank>` to the output file name (per-rank files). |
| `-p` | Replace every `%` in input and output file names with the world-comm rank. |
| `-M[0\|1\|2]` | Read/write in parallel; optional value selects `resolve_shared_ents` (1) or `exchange_ghosts` (2) on top. |
| `-z <file>` | Read Metis partition info for an MPAS grid file and emit an `.h5m` partition. |

### TempestRemap integration (when MOAB is built with TempestRemap)

| Flag | Description |
|---|---|
| `-B` | Use the TempestRemap Exodus reader, then convert to MOAB. |
| `-b` | Convert a MOAB mesh and write via the TempestRemap Exodus writer. |
| `-S` | Scale a climate mesh to the unit sphere. |
| `-i <tag>` | Name of the global-DoF tag to use with `mbtempest`. |
| `-r <order>` | Order of the field DoF (FV=1, SE=[1..N]). |

### Subset selection by geometric topology / sets

Each of the following flags is followed by an ID list (commas, hyphens for ranges; no spaces):

| Flag | Selects |
|---|---|
| `-v` | Geometric volumes |
| `-s` | Geometric surfaces |
| `-c` | Geometric curves |
| `-V` | Geometric vertices |
| `-m` | Material sets (blocks) |
| `-d` | Dirichlet sets (nodesets) |
| `-n` | Neumann sets (sidesets) |
| `-D` | Parallel partitioning sets (`PARALLEL_PARTITION`) |

### Subset selection by element dimension

| Flag | Description |
|---|---|
| `-1` | Edges |
| `-2` | Triangles, quads, polygons |
| `-3` | Tets, hexes, prisms, etc. (Vertices are always exported.) |

## Examples

### Format conversion (extension-driven)

```bash
mbconvert mesh.exo  mesh.h5m       # Exodus → MOAB HDF5
mbconvert mesh.h5m  mesh.vtk       # MOAB HDF5 → VTK
mbconvert mesh.nc   mesh.h5m       # NetCDF (MPAS/HOMME/Euler/...) → MOAB HDF5
```

### Cross-format NetCDF conversion

The new `WRITE_FORMAT` option lets you re-emit the mesh in a different NetCDF grid layout. Source can be any reader-supported format — NetCDF (MPAS/HOMME/CAM/SCRIP/ESMF/Domain), HDF5 (`.h5m`), Exodus (`.exo`), VTK, etc.

```bash
mbconvert -o WRITE_FORMAT=SCRIP   in_mpas.nc  out_scrip.nc
mbconvert -o WRITE_FORMAT=ESMF    in_mpas.nc  out_esmf.nc
mbconvert -o WRITE_FORMAT=DOMAIN  in_mpas.nc  out_domain.nc
```

`WRITE_FORMAT` overrides the source mesh's stored `__MESH_TYPE` tag, so the same input can be emitted in any of the three grid layouts above.

### Extracting subsets

```bash
# Just the material sets numbered 1, 3, 5
mbconvert -m 1,3,5  full.h5m  subset.h5m

# Surfaces 10 through 20
mbconvert -s 10-20  geom.h5m  surfaces.h5m

# Only the 3-D elements
mbconvert -3  mixed.h5m  volumes.h5m
```

### Parallel

```bash
mpiexec -n 4 mbconvert -M2 -O PARALLEL=READ_PART  big.h5m  big_out.h5m
```

### Convert + auto-generate adjacent entities

```bash
mbconvert -I 2  tets.exo  tets_with_faces.h5m    # add face adjacencies (dim 2)
```

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Usage error (bad flag or missing argument) |
| 2 | Read error (file open / parse failed) |
| 3 | Write error (writer rejected the file or set) |
| 4 | Other error (e.g. internal MOAB failure) |
| 5 | Entity / set not found for a subset filter |

## See also

- [`README.IO`](../README.IO) — full list of writer- and reader-specific options
- [`mbtempest`](mbtempest.md) — for climate mesh generation and remap-weight computation
- The legacy [`mbconvert.man`](../mbconvert.man) — `man` page version of this documentation
