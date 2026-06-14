# MOAB Tools Directory

This directory contains miscellaneous tools built on MOAB. Some, like `mbconvert` and `mbsize`, are useful examples of how to write MOAB applications. Others, like `mbzoltan`, are useful both as standalone tools and as building blocks to embed in other applications.

## Tools overview

| Tool | Purpose |
|---|---|
| `mbconvert` | Converts mesh between any two formats MOAB can read/write |
| `mbtempest` | Climate mesh generation and conservative remap weight computation |
| `mbIntxCheck` | Verify a precomputed mesh intersection |
| `compareMaps`, `compareFiles` | Diff weight maps / per-tag mesh files |
| `measure` | Library: length / area / volume of MOAB entities |
| `mbsize` | Reads mesh and prints size by geometric owner or block/boundary |
| `mbskin` | Computes and outputs the mesh skin (Mesquite test input) |
| `mbsurfplot` | Plots the mesh of a geometric surface projected to a plane |
| `mbtagprop` | Propagate tags from entity sets to their contained entities |
| `hexmodops` | Dual-based hex meshing research support |
| `spheredecomp` | Decompose a tet mesh into hexes that resolve vertex-located spheres |
| `mbpart` | Static partitioner — Zoltan / Metis / CHACO |
| `mcnpmit` | Interpolate `.unv` meshes on MCNP meshtal files |
| `qvdual` | Mesh visualization (Qt / VTK); hex-mesh dual support |
| `vtkMOABReader` | MOAB reader plugin for VTK |

## Adding a new tool

There are two ways a tool can be wired into the MOAB build system.

### Method 1 — driven by the main MOAB `configure`

1. Create a `Makefile.am` in the tool's directory. Use `tools/converter/Makefile.am` as a model.
2. Add a call to the `MB_OPTIONAL_TOOL` macro in the main MOAB `configure.in`, with three arguments:
   1. The name of the tool
   2. An empty string (two square brackets: `[]`)
   3. `[yes]` or `[no]` — whether the tool is built by default
3. Add the `Makefile` to the list passed to `AC_CONFIG_FILES` at the end of `configure.in` (e.g. `tools/<dir>/Makefile`).
4. Add a conditional to `tools/Makefile.am`:
   ```make
   if ENABLE_<name>
     SUBDIRS += <dir>
   endif
   ```
   where `<name>` matches `MB_OPTIONAL_TOOL`'s first argument and `<dir>` matches `AC_CONFIG_FILES`.

### Method 2 — tool owns its own `configure`

1. Create the tool's own makefiles, `configure`, etc.
2. Add a call to `MB_OPTIONAL_TOOL` in `configure.in` with:
   1. The name of the tool
   2. The path to the tool (e.g. `[tools/<dir>]`)
   3. `[yes]` or `[no]`
3. Add the same `tools/Makefile.am` conditional shown above.

## CMake

Add the tool's source to `tools/CMakeLists.txt` under the appropriate `add_executable` / target definition. Optional features should be gated with `MOAB_HAVE_*` flags emitted by `MOABConfig.h.in`.
