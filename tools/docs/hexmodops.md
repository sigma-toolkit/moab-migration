# hexmodops {#tools_hexmodops}

Build sample meshes that illustrate the local hex modification operations developed by Tautges. Useful for studying the operations themselves and for generating test inputs.

## Synopsis

```
hexmodops [-h5m] [-vtk] {-ap | -foc | -fs | -cp | -tcp | -thp}
```

## Output options

| Flag | Description |
|---|---|
| `-h5m` | Write the constructed mesh as a MOAB `.h5m` file. |
| `-vtk` | Write the constructed mesh as a VTK file. |

(One of `-h5m` or `-vtk` should be given. They can be combined.)

## Operation selectors (mutually exclusive)

| Flag | Operation |
|---|---|
| `-ap`  | Atomic pillow |
| `-foc` | Face open-collapse |
| `-fs`  | Face shrink |
| `-cp`  | Chord push |
| `-tcp` | Triple chord-push |
| `-thp` | Triple hex-push |

## Example

```bash
hexmodops -h5m -ap   # atomic-pillow sample → .h5m
hexmodops -vtk -tcp  # triple chord-push sample → .vtk
```
