# mbskin

Generate input files for testing **Mesquite** (MESh QUality Improvement Toolkit). Produces a mesh in which interior vertices are marked free and skin (boundary) vertices are marked fixed.

## Synopsis

```
mbskin  <input_file>  <output_file>
```

## What it does

1. Reads the input mesh.
2. Finds the set of all elements of the highest dimension present (volume elements if any exist; otherwise face elements).
3. Generates the **skin** of that set of elements.
4. Creates an integer tag named `fixed`.
5. Sets `fixed = 0` on all interior vertices and `fixed = 1` on all vertices in the skin.
6. Writes out the set of elements with the new tag.

## Example

```bash
mbskin  volumes.h5m  volumes_fixed.h5m
# Output has a `fixed` tag: 0 on interior verts, 1 on boundary verts.
```

## Notes

- The "fixed" naming is a Mesquite convention; vertices with `fixed=1` are held in place during quality improvement.
- Works on any mixed-dimension mesh by picking the highest dimension present.
