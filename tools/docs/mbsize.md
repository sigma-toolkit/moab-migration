# mbsize

List the entity types in a MOAB-readable mesh file, plus counts and size statistics. Useful for quick inspection: "how many cells does this mesh have, what types, by which set?"

## Synopsis

```
mbsize [-g] [-m] [-t|-l|-ll] <mesh_file> [<mesh_file2> ...]
```

## Modes (mutually exclusive)

| Flag | Output |
|---|---|
| *(none)* | Per-type entity measurement statistics (default) |
| `-t` | List tags and the count of entities of each type per tag |
| `-l` | Counts of entity types only |
| `-ll` | Verbose: list every entity |

## Grouping (mutually exclusive)

| Flag | Group by |
|---|---|
| *(none)* | One report per input file (default) |
| `-g` | Per geometric topology set |
| `-m` | Per material set and boundary-condition (Dirichlet/Neumann) set |

## Examples

```bash
mbsize  mesh.h5m                    # per-type stats
mbsize -l  mesh.h5m                 # raw counts only
mbsize -t  mesh.h5m                 # break down by tag
mbsize -m  mesh.h5m                 # one report per material/BC set
mbsize -g  mesh.h5m mesh2.h5m       # per geometric topology, across two files
```
