# mbsurfplot {#tools_mbsurfplot}

Plot the mesh of a single geometric surface projected to a plane. Output is written to stdout (redirect to a file).

## Synopsis

```
mbsurfplot [-g|-p|-s]  <Surface_ID>  <input_file>
```

## Output format (mutually exclusive)

| Flag | Format |
|---|---|
| `-g` | GNUPlot data file *(default)* |
| `-p` | Encapsulated PostScript |
| `-s` | SVG |

## Positional arguments

- `<Surface_ID>` — geometric-surface ID whose mesh is to be exported
- `<input_file>` — mesh file to read

## Examples

```bash
mbsurfplot -g  3  geom.h5m  > surface3.dat
gnuplot -e "plot 'surface3.dat' with lines" -p

mbsurfplot -p  3  geom.h5m  > surface3.eps
mbsurfplot -s  3  geom.h5m  > surface3.svg
```
