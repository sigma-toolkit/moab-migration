# spheredecomp

Construct an all-hex mesh for a set of spheres embedded in a substrate. Sphere locations and radii are specified as a set of vertices and tags on those vertices; each vertex implicitly defines one sphere.

## Algorithm

1. The input vertices are triangulated with a tet mesh.
2. Each tet is subdivided so the spherical surface around each vertex is resolved:
   - **Inside the spheres** — each "corner region" (the part of the tet that lies inside one of the four corner spheres) is split into 4 hexes using a standard tetrahedron primitive.
   - **Outside the spheres** — the remaining "interstices" region (the tet with the four corner spheres removed) is split using a midpoint subdivision.
3. The result is **28 hexes per input tet**: 16 hexes split across the 4 corner-sphere regions and 12 hexes filling the interstices.

This produces an all-hex mesh in which the spherical interfaces are explicitly captured by element faces.

## Use cases

- Particle / sphere-pack simulations where each sphere is a distinct material region.
- Multi-phase mesh setups where the spheres represent inclusions in a continuous substrate.
- Generating test inputs that exercise hex-meshing topology while still containing curved interfaces.
