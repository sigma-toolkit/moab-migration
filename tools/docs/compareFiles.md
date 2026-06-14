# compareFiles

Compare two MOAB `.h5m` files that represent the **same mesh** and emit a third file that carries the per-tag differences. The inputs must share `GLOBAL_ID` values on their elements, but their `EntityHandle` values can differ arbitrarily (e.g. because the two files were partitioned differently or saved by different MPI rank counts).

Output is an `.h5m` file derived from one of the inputs, with an extra tag per compared tag whose value is the difference (`file1_value − file2_value`) per entity.

## Synopsis

```
mbcmpfiles -i <file1.h5m> -j <file2.h5m> -o <out.h5m> [options]
```

> Binary name on disk is `mbcmpfiles`.

## Options

| Flag | Default | Description |
|---|---|---|
| `-i <file1.h5m>`, `--input1 <file1.h5m>` | — | First input mesh. |
| `-j <file2.h5m>`, `--input2 <file2.h5m>` | — | Second input mesh. |
| `-o <out.h5m>`, `--outfile <out.h5m>` | — | Output file with the per-tag difference fields added. |
| `-n <tag>`, `--tagname <tag>` | *(all tags)* | Restrict the comparison to a single tag. If unset, every tag found on both files is compared. |
| `-d <dim>`, `--dimension <dim>` | — | Topological dimension of entities to compare (`0` vertices, `1` edges, `2` faces, `3` cells). |

## What it does

1. Loads both files independently.
2. Builds a map keyed on `GLOBAL_ID` so corresponding entities can be matched across the two files even when their `EntityHandle` values differ.
3. For each compared tag, computes `file1 − file2` per entity and stores the result as a new tag on the chosen entity range.
4. Writes the result to `<out.h5m>`. The output preserves the topology of file 1 plus the new difference tags.

## Examples

```bash
# Compare a single tag (e.g. a projected scalar field), only on cells
mbcmpfiles -i serial_run.h5m -j parallel_run.h5m  \
           -n AnalyticalSoln -d 2  -o diff.h5m

# Compare every shared tag, defaults to all entity dimensions
mbcmpfiles -i baseline.h5m -j candidate.h5m  -o all_diffs.h5m
```

After the run, open `diff.h5m` in your favorite viewer (ParaView via the MOAB plugin, VisIt) and look at the new difference tag — anywhere it's non-zero is where the two files disagree.

## Notes

- The two files must reference the same mesh topology — same GLOBAL_ID set on the elements being compared. Entity counts must match per-dimension.
- Tags that are present on only one of the two files are skipped (with a warning).
- Tag types must match between the two files (e.g. you can't compare a `DOUBLE` tag in file 1 against an `INTEGER` tag of the same name in file 2).

## See also

- [`compareMaps`](compareMaps.md) — for diffing remap-weight files in SCRIP format
