# measure (library)

A small utility module — **not** a standalone binary — providing geometric measurement (length, area, volume) for MOAB entity types. Used internally by other tools in this directory (`mbsize`, etc.) and available for embedding in client code.

## Source

- [`measure.hpp`](../measure.hpp) — public API
- [`measure.cpp`](../measure.cpp) — implementation

## API

```cpp
double edge_length(const double* start_vtx_coords,
                   const double* end_vtx_coords);

double measure(moab::EntityType type,
               int             num_vertices,
               const double*   vertex_coordinates);
```

### `edge_length(start, end)`

Returns the Euclidean distance between two 3-D points. The two arguments must each point to **three contiguous doubles** holding `{x, y, z}`.

### `measure(type, num_vertices, coords)`

Returns the geometric measure (length, area, or volume) of a single entity of the given `moab::EntityType`. `coords` is the flat XYZ coordinate array for the entity's vertices, in canonical MOAB ordering (3 × `num_vertices` doubles).

The function returns `0.0` for entity types it does not handle.

## Supported entity types and formulas

| Entity type | Measure | Formula |
|---|---|---|
| `MBEDGE` | length | `\|p1 − p0\|` |
| `MBTRI` | area | `½ \|(p1−p0) × (p2−p0)\|` |
| `MBQUAD` | area | Fan-triangulation around centroid (4 vertices) |
| `MBPOLYGON` | area | Fan-triangulation around centroid (`num_vertices` corners) |
| `MBTET` | volume | `⅙ ((p1−p0) × (p2−p0)) · (p3−p0)` |
| `MBPYRAMID` | volume | Two-tet decomposition through apex (vertex 4) |
| `MBPRISM` | volume | Three-tet decomposition |
| `MBHEX` | volume | Five-tet decomposition (canonical hex split) |
| anything else | — | returns `0.0` |

### Conventions

- All formulas assume **3-D Cartesian** coordinates. Spherical-polygon area on the unit sphere is not provided here — see `moab::IntxMesh::IntxUtils` for that.
- The polygon area formula uses fan triangulation around the centroid, which is signed-area-correct for convex polygons and a reasonable approximation for moderately non-convex ones. It does **not** handle self-intersecting polygons.
- Hex volume uses the canonical 5-tet decomposition; this matches Exodus / MOAB element ordering and is exact for trilinear hexes.

## Usage example

```cpp
#include "measure.hpp"
#include "moab/Core.hpp"

moab::Interface* mb = /* … */;
moab::Range cells;
mb->get_entities_by_dimension(0, 2, cells);

double total_area = 0.0;
for (moab::EntityHandle c : cells) {
    const moab::EntityHandle* conn = nullptr;
    int n = 0;
    mb->get_connectivity(c, conn, n);

    std::vector<double> coords(3 * n);
    mb->get_coords(conn, n, coords.data());

    total_area += measure(mb->type_from_handle(c), n, coords.data());
}
```

## License

Original code copyright Lawrence Livermore National Laboratory (LLNS contract B545069 with UW–Madison). BSD-2-Clause; see source header for the full notice.
