# iMOAB parallel-test regression harness

A generic shell harness for the iMOAB parallel test binaries. Two regression
modes:

- **`bfb` mode** — strict byte-for-byte comparison of per-cell digests
  across rank counts. Requires the test binary to (a) write computed
  weight maps to disk and (b) emit `--digest_prefix` per-cell dumps.
  Currently: `dualmap_caas`, `spmv_bfb`.
- **`tol` mode** — tolerance-based comparison of `h5m` output artifacts
  via `mbcmpfiles` (L2-norm of every common double tag is checked
  against `TOLERANCE`, default `1e-9`). Required for tests whose
  online-computed weights are not BfB across rank counts (every
  iMOAB coupler test). Currently: `coupler` (more configs to follow).

## Layout

```
regression/
├── imoab_test_runner.sh    top-level dispatcher
├── lib_common.sh           shared helpers (logging, mbcmpfiles, mktemp)
├── lib_bfb_harness.sh      driver for TEST_MODE=bfb
├── lib_tol_harness.sh      driver for TEST_MODE=tol
├── configs/
│   ├── dualmap_caas.cfg
│   ├── spmv_bfb.cfg
│   └── coupler.cfg
└── README.md
```

## Quick start

```bash
cd /nfs/gce/projects/sigma/vijaysm/moab/test/parallel/regression

# List available configs
./imoab_test_runner.sh --list

# Run dual-map CAAS BfB sweep with config defaults (n=1,2,4 on stock meshes)
./imoab_test_runner.sh dualmap_caas

# Same test, user-supplied source/target meshes, custom rank sweep
./imoab_test_runner.sh dualmap_caas \
    --src /path/to/atm.h5m \
    --tgt /path/to/ocn.h5m \
    --ranks "1 2 4 8 16"

# Coupler tolerance regression (requires mbcmpfiles)
export MBCMPFILES=/nfs/gce/projects/sigma/vijaysm/moab/build/e3sm-gnu/tools/mbcmpfiles
./imoab_test_runner.sh coupler --tol 1e-10
```

## CLI options (all configs)

| Option | Effect |
|---|---|
| `--src <path>` | Override `DEFAULT_SRC` for this run |
| `--tgt <path>` | Override `DEFAULT_TGT` for this run |
| `--ranks "N1 N2 ..."` | Override `DEFAULT_RANKS` (must include `1`) |
| `--workdir <dir>` | Use `<dir>` for intermediate files (preserves them) |
| `--keep` | Keep mktemp workdir on exit (else auto-deleted) |
| `--mbcmp <path>` | Path to `mbcmpfiles` binary (sets `$MBCMPFILES`) |
| `--tol <value>` | Override `TOLERANCE` for `tol`-mode configs |
| `--extra <flag=value>` | Pass `EXTRA_VAL_1` for the config's `EXTRA_FLAG_1` |

## Adding a new config

Create `configs/<name>.cfg`. Required keys depend on `TEST_MODE`:

### `TEST_MODE=bfb`

```bash
# CONFIG_DESC: <one-line description shown by --list>
TEST_NAME=<binary-name>
TEST_MODE=bfb

EXE=$MOAB_BUILD/test/parallel/<binary>
SRC_FLAG="-t"; TGT_FLAG="-m"
DEFAULT_SRC="<absolute path>"
DEFAULT_TGT="<absolute path>"
ONLINE_FLAGS="--compute_online --write_maps ${TEST_NAME}_baseline"
LOAD_FLAGS="-l ${TEST_NAME}_baseline_lo.nc -h ${TEST_NAME}_baseline_hi.nc"
DIGEST_FLAG="--digest_prefix ${TEST_NAME}_digest"
DIGEST_PREFIX="${TEST_NAME}_digest"
DIGEST_KERNELS="lo hi dual"        # one ${TEST_NAME}_digest_<k>_<n>.txt per kernel
DEFAULT_RANKS="1 2 4"
H5M_ARTIFACTS=""                   # optional: also strict-cmp h5m outputs
EXTRA_ARGS=""                      # optional: extra flags for every run
```

> **Namespacing rule (mandatory).** Every persistent artifact a config writes
> — baselines, digests, renamed h5m references — **must** be prefixed with
> `${TEST_NAME}`. The harness namespaces h5m artifacts on rename
> (`<art> -> ${TEST_NAME}_<art>.n${n}`) automatically; configs are
> responsible for prefixing the map / digest filenames they pass on the
> command line. This guarantees that two configs running concurrently in a
> shared `WORKDIR` cannot overwrite each other's outputs.

### `TEST_MODE=tol`

```bash
# CONFIG_DESC: <one-liner>
TEST_NAME=<binary-name>
TEST_MODE=tol

EXE=$MOAB_BUILD/test/parallel/<binary>
SRC_FLAG="-t"; TGT_FLAG="-m"
DEFAULT_SRC="<absolute path>"
DEFAULT_TGT="<absolute path>"
H5M_ARTIFACTS="recvAtm.h5m recvOcn.h5m"   # files the test writes to CWD
DEFAULT_RANKS="1 2 4"
TOLERANCE=1e-9

# Optional: extra-mesh flag (e.g. -l for land)
EXTRA_FLAG_1="-l"
EXTRA_DEFAULT_1="<absolute path>"

# Optional: extra args every run gets
EXTRA_ARGS="--no_regression"
```

## Algorithm summaries

### `bfb` driver
1. Serial run computes maps online and writes them to disk.
2. Serial run reloads those maps and writes baseline digests.
3. Parallel runs at each `n>1` reload the same maps and write digests.
4. `cmp -s` every `(kernel, n>1)` digest against `<kernel>_1.txt`.
5. (optional) `mbcmpfiles` strict-compare each `H5M_ARTIFACTS` entry.

The serial reload (step 2) — not the online-compute (step 1) — is the
baseline because online computation populates per-cell areas via a
different code path than `iMOAB_LoadMapFile`, and we want the baseline
to use the same code path as the parallel runs.

### `tol` driver
1. Serial run produces every `H5M_ARTIFACTS` file; the harness renames
   each to `<file>.n1`.
2. Each parallel run produces the same files; the harness renames them
   to `<file>.n${n}`.
3. `mbcmpfiles -i <file>.n1 -j <file>.n${n}` checks that the L2 norm
   of every common double tag is `<= TOLERANCE`.

Online-computed weights in the coupler tests are not BfB across rank
counts by design (TempestRemap's polygon sweep depends on local mesh
ordering), so the harness only enforces tolerance, not byte equality.

## Environment variables

| Var | Default | Purpose |
|---|---|---|
| `MOAB_BUILD` | `/nfs/gce/projects/sigma/vijaysm/moab/build/e3sm-gnu` | Configs use this to locate `EXE` |
| `MOAB_MESH_DIR` | `/nfs/gce/projects/sigma/vijaysm/moab/MeshFiles` | Configs use this for `DEFAULT_SRC`/`DEFAULT_TGT` |
| `MPIRUN` | `mpirun` | Launcher to use |
| `MBCMPFILES` | (unset) | Path to `mbcmpfiles` (`tol` mode and BfB Phase E require this) |
| `TOLERANCE` | `1e-9` | Default L2-norm tolerance for `tol` mode |
| `WORKDIR` | (unset → mktemp) | Use this dir; implies `--keep` |
| `KEEP_WORKDIR` | `0` | Set to `1` (or pass `--keep`) to preserve mktemp dir |

## Migration from `imoab_dualmap_caas_bfb.sh`

The legacy script remains in place. The equivalent harness invocation is:

```bash
./regression/imoab_test_runner.sh dualmap_caas
```

Both produce the same digest files and the same pass/fail verdict.

## Roadmap (v2)

- Land + river + physgrid support (`coupler`, `apg2_ol_coupler`,
  `phatm_ocn_coupler`).
- Optional CTest integration: enumerate `configs/*.cfg` and register
  one `add_test` per config per rank-count.
- Per-config `tag-of-interest` filter for `mbcmpfiles` (currently it
  diffs every common double tag).
- `read_compute_map` BfB config (needs `--digest_prefix` added to
  the binary).
