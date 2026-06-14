# MOAB Tools

Command-line and library tools built on the MOAB library. Some are general-purpose mesh utilities (conversion, partitioning, inspection); others target specific application domains like climate remapping or hex meshing.

## Tool index

| Tool | Purpose | Reference |
|---|---|---|
| [mbconvert](mbconvert.md) | Convert meshes between formats; extract subsets | full options & examples |
| [mbtempest](mbtempest.md) | Generate climate meshes (CS/RLL/ICO) and compute conservative remap weights | full options & examples |
| [mbIntxCheck](mbIntxCheck.md) | Verify a precomputed mesh intersection against its source/target | full options |
| [compareMaps](compareMaps.md) | Diff two sparse remap-weight `.nc` map files | full options |
| [compareFiles](compareFiles.md) | Diff per-tag values between two `.h5m` files of the same mesh | full options |
| [measure](measure.md) | Library: length/area/volume of MOAB entities | API reference |
| [mbsize](mbsize.md) | List entity counts and per-set size statistics | full options |
| [mbskin](mbskin.md) | Generate the skin of a mesh; tag boundary vertices for Mesquite | usage |
| [mbsurfplot](mbsurfplot.md) | Plot a single geometric surface to GNUPlot / EPS / SVG | usage |
| [mbtagprop](mbtagprop.md) | Propagate tags from sets to their contained entities | full options |
| [hexmodops](hexmodops.md) | Build sample meshes for local hex modification operations | usage |
| [spheredecomp](spheredecomp.md) | All-hex mesh around vertex-located spheres in a substrate | overview |
| [tools (overview)](tools.md) | How tools are registered with the MOAB build system | build integration |

## Conventions across tools

- **Input format detection**: most tools determine reader/writer from the file extension. Override with `-f <format>` (mbconvert) where supported.
- **Option pass-through**: many tools accept `-o <opt>` (writer options) and `-O <opt>` (reader options) that flow through to MOAB's I/O layer. See [`README.IO`](../README.IO) for the full list of writer/reader option strings.
- **Parallel execution**: MPI-enabled tools detect `mpiexec` / `mpirun` automatically; some accept a `-M[0|1|2]` flag for parallel read/resolve/ghost layering.
- **Native format**: MOAB's native binary format is HDF5-based `.h5m`. Lossless round-tripping is only guaranteed through `.h5m`; other formats may not capture all of MOAB's data model.

## Build system integration

See [tools.md](tools.md) for how to add a new tool to MOAB's autotools / CMake build under `tools/`.
