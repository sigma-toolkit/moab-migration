# `install-moab.sh` — Handoff Document

> **Renamed in 2026-05.** Was `install-moab-e3sm.sh` and lived at
> `build-e3sm/install-moab-e3sm.sh`. The old path remains as a 3-line
> backwards-compat shim (`install/install-moab-e3sm.sh`) that exec's
> `install-moab.sh` with all args forwarded — existing E3SM workflows /
> docs that reference the old name keep working unchanged.
>
> **Push 1 (2026-05)** added a machine database and `--machine` /
> `--list-machines` / `--compiler` flags. Database entries are pure
> metadata (no hardcoded TPL paths).
>
> **Push 2 (this revision)** added `--profile={e3sm,standalone}`, a
> `--print` mode, `--e3sm-root=PATH` (or `$E3SM_ROOT`), and a Python helper
> (`install/scripts/e3sm_env.py`) that sources E3SM's
> `cime_config/machines/config_machines.xml` for the resolved
> (machine, compiler) pair. The e3sm profile auto-loads modules and env
> vars and even auto-fills `--with-blas=`/`--with-lapack=` from
> `$BLAS_ROOT`/`$LAPACK_ROOT` when E3SM exposes them. The legacy
> `suggest_configuration.sh` is now a 3-line shim onto
> `install-moab.sh --print --profile=standalone`.

End-to-end orchestrator that builds and installs MOAB plus its three required
TPLs (Eigen3, Zoltan, TempestRemap) on systems where E3SM is already
buildable. Designed for HPC clusters: portable across Linux/macOS, robust
against transient network and TPL build failures, resumable across runs,
and capable of self-correcting common configuration problems before they
cascade.

---

## Contents

1. [What the script does](#what-the-script-does)
2. [Prerequisites](#prerequisites)
3. [Quick start](#quick-start)
4. [Machine database](#machine-database)
5. [Directory layout](#directory-layout)
6. [All command-line options](#all-command-line-options)
6. [All environment variables](#all-environment-variables)
7. [Compiler resolution rules](#compiler-resolution-rules)
8. [MPI wrapper validation](#mpi-wrapper-validation)
9. [MOAB source management](#moab-source-management)
10. [TPL build pipeline](#tpl-build-pipeline)
11. [Build-system selection (autotools vs CMake)](#build-system-selection)
12. [Resume behavior and incremental rebuilds](#resume-behavior)
13. [Failure remediation](#failure-remediation)
14. [Cray PrgEnv detection](#cray-prgenv-detection)
15. [Dry-run mode](#dry-run-mode)
16. [Verification](#verification)
17. [Logging and live progress](#logging-and-live-progress)
18. [Exit codes](#exit-codes)
19. [Common workflows](#common-workflows)
20. [Troubleshooting](#troubleshooting)

---

## What the script does

A single invocation performs, in order:

1. Parse and validate arguments.
2. Resolve C/C++/Fortran compilers from `--cc/--cxx/--fc/--f77`, `MPI_ROOT`,
   PATH, or environment (in that priority).
3. Compile-and-link a tiny MPI program with each compiler to validate the
   wrapper (skippable).
4. Clone or update the MOAB source tree from `master` (or a chosen branch).
5. Run `autoreconf -fi` to bootstrap MOAB's autotools build (if missing).
6. Build the three TPLs into separate per-TPL prefixes:
   - `eigen3` (header-only): download tarball → copy headers
   - `zoltan` (autotools): download tarball → configure → make → install
   - `tempestremap` (autotools): shallow `git clone` of master → configure →
     make → install
   Each TPL has its own per-phase log file, retry budget, and remediation
   logic.
7. Configure MOAB (autotools or CMake) with `--with-eigen3=... --with-zoltan=...
   --with-tempestremap=...` (or the CMake equivalents) pointing at the
   pre-built TPL prefixes.
8. Build and install MOAB.
9. Verify the install by checking for known headers and libraries.

On a re-invocation, the script automatically skips work that's already
done and only rebuilds what changed (see [Resume behavior](#resume-behavior)).

---

## Prerequisites

### Required on the host

| Tool | Used for | How the script reacts if missing |
|---|---|---|
| `bash` (≥3.2) | the script itself | shebang failure |
| `git` | cloning MOAB & TempestRemap | hard error in `prepare_moab_source` |
| `curl` or `wget` | TPL tarball downloads | hard error during prefetch |
| `tar`, `make`, `autoreconf` | building TPLs and MOAB (autotools) | hard error at the relevant step |
| `cmake` (≥3.x) | only when `--build-system=cmake` | hard error in early validation |
| C / C++ / Fortran MPI wrappers | building everything | MPI validation step rejects, with the actual compiler error printed |

### Required environment

| Variable | Required? | What for |
|---|---|---|
| `HDF5_ROOT` | yes | parallel HDF5 install root (passed to MOAB and TempestRemap) |
| `NETCDF_C_PATH` | yes | NetCDF-C install root (passed to MOAB and TempestRemap) |
| `PNETCDF_PATH` | yes | Parallel-NetCDF install root (passed to MOAB) |
| `MPI_ROOT` | optional | MPI install root. If set, `$MPI_ROOT/bin/mpicc` etc. are preferred. If not set, the script falls back to PATH or, on Cray, `cc/CC/ftn`. |

The HDF5 / NetCDF / PNetCDF installs must already exist; this script does
not build them.

---

## Quick start

### Standard cluster (GCC + MPICH/OpenMPI + Spack-built TPLs)

```bash
# Load your environment first so wrappers are on PATH
module load gcc/12 mpich/4.1 hdf5/1.14 netcdf-c/4.9 parallel-netcdf/1.12

# Run the script — clones MOAB into $PWD/moab-src, installs to $PREFIX
mkdir -p /scratch/$USER/moab-build && cd /scratch/$USER/moab-build
/path/to/install-moab-e3sm.sh \
    --hdf5-root=$HDF5_ROOT \
    --netcdf-root=$NETCDF_C_PATH \
    --pnetcdf-root=$PNETCDF_PATH \
    --prefix=$HOME/install/MOAB
```

### Cray PrgEnv (Frontier, Perlmutter, Aurora, …)

```bash
module load PrgEnv-gnu cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf

# No --mpi-root needed; Cray detection picks cc/CC/ftn automatically
/path/to/install-moab-e3sm.sh \
    --hdf5-root=$HDF5_DIR \
    --netcdf-root=$NETCDF_DIR \
    --pnetcdf-root=$PNETCDF_DIR \
    --prefix=/lustre/orion/$USER/install/MOAB
```

### Re-running to update MOAB to latest master

```bash
# Same command as before — script does the right thing:
#   * git fetch + reset --hard origin/master
#   * skip TPLs (already installed)
#   * reconfigure MOAB if HEAD changed (fingerprint differs)
#   * incremental make + make install
/path/to/install-moab-e3sm.sh \
    --hdf5-root=$HDF5_ROOT \
    --netcdf-root=$NETCDF_C_PATH \
    --pnetcdf-root=$PNETCDF_PATH \
    --prefix=$HOME/install/MOAB
```

---

## Profiles (`--profile`)

The script supports two audiences via a profile knob:

| Profile | When to use | What changes |
|---|---|---|
| `e3sm` (default) | You have an E3SM checkout and want the same environment E3SM uses for the target machine. | Sources modules + env vars from `$E3SM_ROOT/cime_config/machines/config_machines.xml`. Requires `--e3sm-root=PATH` or `$E3SM_ROOT`. Auto-fills `--with-blas=`/`--with-lapack=` from `$BLAS_ROOT`/`$LAPACK_ROOT` if E3SM exposes them. |
| `standalone` | MOAB downstream user without E3SM, or you want full manual control over modules/env. | Trusts your loaded environment as-is. No CIME interaction. |

### `--profile=e3sm` walk-through

```bash
# 1. Tell the script where your E3SM checkout is (or set $E3SM_ROOT)
export E3SM_ROOT=$HOME/Code/E3SM

# 2. Load your minimum env (anything else needed beyond what config_machines.xml loads)
#    Typically nothing -- E3SM's XML lists its full module stack

# 3. Run with --profile=e3sm (the default; you can omit it)
./install/install-moab.sh \
    --machine=bebop \
    --hdf5-root=/lcrc/group/e3sm/soft/bebop/hdf5/1.12.3/gcc-13.2.0/openmpi-4.1.8 \
    --prefix=$PWD/installs
```

The script will:
1. Resolve `--machine=bebop` → `MACHINE_META_E3SM_NAME=bebop` (or whatever your entry maps to).
2. Run `python3 install/scripts/e3sm_env.py --e3sm-root=... --machine=bebop --compiler=gnu` and `source` its output (this happens under `--dry-run` too — the script runs in its own process, so it doesn't pollute your parent shell). This loads modules and exports env vars.
3. Run two adapter shims:
   - If `NETCDF_C_PATH` is unset but `NETCDF_PATH` is set (anlgce-ub22, crux), alias.
   - If `HDF5_ROOT` is unset but `CRAY_HDF5_PARALLEL_PREFIX` is set (Cray PrgEnv), alias.
4. If `BLAS_ROOT` and `LAPACK_ROOT` got set (bebop, improv pattern) and you didn't pass `--extra` / `--extra-tempestremap`, auto-fill them with the right `--with-blas=`/`--with-lapack=` specs (preferring static `lib*.a` when present).
5. **Show you the resolved env and prompt for y/N confirmation** before continuing. Bypass with `--yes` / `-y` / `ASSUME_YES=yes`. Rejects (anything other than `y` or `yes`) exit cleanly with a hint about re-running.
6. Continue with normal compiler resolution, MPI validation, TPL builds, MOAB build/install.

### The confirmation prompt

After step 4 you'll see something like:

```
[install-moab] E3SM env applied. Resolved roots:
[install-moab]   HDF5_ROOT     = /nfs/gce/projects/climate/.../hdf5/...
[install-moab]   NETCDF_C_PATH = /nfs/gce/projects/climate/.../netcdf/...
[install-moab]   PNETCDF_PATH  = /nfs/gce/projects/climate/.../pnetcdf/...
[install-moab]   BLAS_ROOT     = /lcrc/group/e3sm/soft/.../netlib-lapack/...
[install-moab]   LAPACK_ROOT   = /lcrc/group/e3sm/soft/.../netlib-lapack/...
[install-moab] Auto-set --extra-tempestremap from BLAS_ROOT/LAPACK_ROOT
[install-moab] Auto-set --extra (MOAB) from BLAS_ROOT/LAPACK_ROOT

[install-moab] Does this E3SM environment look correct? [y/N]
```

This is the safety net for "wrong machine name" or "wrong compiler family"
mistakes — you spot the bad path immediately rather than 5 minutes later
in a TPL configure log. Type `y` to continue, anything else to abort
cleanly.

Skip the prompt for CI / scripted runs:
```bash
./install-moab.sh --profile=e3sm --machine=bebop --yes ...
# or
ASSUME_YES=yes ./install-moab.sh --profile=e3sm --machine=bebop ...
```

If stdin isn't a tty (you piped the invocation, ran from a CI runner,
etc.) AND you didn't pass `--yes`, the script errors with a clear hint
rather than hanging on `read`.

### `--profile=standalone` walk-through

```bash
# Load your env yourself
module load gcc/13.2.0 openmpi hdf5 netcdf-c parallel-netcdf

./install/install-moab.sh \
    --machine=bebop --profile=standalone \
    --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH \
    --prefix=$PWD/installs
```

Identical behavior to Push 1 / pre-rename. The machine entry's
`standalone_hint` is printed informationally but no modules are loaded
automatically.

### `--print` mode

Quiet variant of `--dry-run`: emits **only** the resolved configure
(autotools) or `cmake` command on stdout, no banner, no resume state.
Useful when you want to inspect or hand-edit before running, or to feed
the configure line into a separate workflow.

```bash
./install/install-moab.sh --print --machine=bebop --profile=standalone
# /path/to/moab-src/configure \
#     --with-mpi \
#     CC=mpicc \
#     ...
```

The legacy `suggest_configuration.sh` at the repo root is now a 3-line
shim that delegates here:

```bash
./suggest_configuration.sh --machine=bebop      # equivalent to:
./install/install-moab.sh --print --profile=standalone --machine=bebop
```

---

## Machine database

The script ships with a registry of currently-supported HPC environments.
Each entry is **pure metadata**: a hostname/env-var match predicate, a default
compiler family, an E3SM machine name (used by Push 2 for CIME lookup), a
module-load hint string, and a last-validated date. **No TPL paths are
hardcoded** — the user is expected to load modules so that `HDF5_ROOT`,
`NETCDF_C_PATH`, and `PNETCDF_PATH` are set in the environment before
invoking the script. (Push 2 will automate this for `--profile=e3sm` by
sourcing `config_machines.xml`.)

### Listing entries

```bash
$ install-moab.sh --list-machines
Registered machines:

  NAME         DEFAULT SUPPORTED                      LAST-VALIDATED
  ----         ------- ---------                      --------------
 *bebop        gnu     gnu,intel                      2026-05-08
  improv       gnu     gnu,intel,aocc                 TBD
  crux         gnu     gnu,cray,nvhpc                 TBD
  gce          gnu     gnu,intel                      2026-05-08
  perlmutter   gnu     gnu,intel,nvidia,aocc          2026-05-08

  * = auto-detected on this host
```

`TBD` entries are stubs that have not yet been validated against a real build.
They are still selectable via `--machine=NAME`; please report failures so the
entry can be promoted (or fixed).

### Selecting an entry

| Invocation | Behavior |
|---|---|
| (no flag) | same as `--machine=auto` |
| `--machine=auto` | iterate registered `match()` predicates; use the first that returns 0 |
| `--machine=NAME` | force a specific entry; error if `NAME` isn't registered |
| `--compiler=NAME` | override the entry's `default_compiler`; warn if not in `supported_compilers` |

### What a match looks like

The orchestration banner reports the resolved machine before any TPL or MOAB
work begins:

```
[install-moab] Machine         : bebop (compiler family: gnu; last-validated: 2026-05-08)
[install-moab] Module hint     : module load gcc/13.2.0 openmpi/4.1.8 hdf5/1.12.3 netcdf-c parallel-netcdf
[install-moab] Machine notes   : Site-installed netlib BLAS/LAPACK at /lcrc/group/e3sm/soft/... (use --extra to point MOAB/TempestRemap at it; see below)
```

If no entry matches, the banner says so and the script continues in generic
mode (same behavior as before the database existed):

```
[install-moab] Machine         : (none detected; pass --machine=NAME or --list-machines to see options)
```

### Adding a new machine entry

In `install/install-moab.sh`, define two functions and append to `MACHINE_REGISTRY`.
Template:

```bash
machine_<name>_match() {
    [[ "${LMOD_SYSTEM_NAME:-}" == "<name>" ]] && return 0
    [[ "$(_hn)" == <hostname-prefix>* ]]
}
machine_<name>_meta() {
    MACHINE_META_E3SM_NAME="<name in config_machines.xml, or empty>"
    MACHINE_META_DEFAULT_COMPILER="gnu"           # or intel, cray, nvhpc, aocc, nvidia
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel"
    MACHINE_META_STANDALONE_HINT="module load gcc openmpi hdf5 netcdf-c parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="TBD"             # promote to YYYY-MM-DD after a successful build
    MACHINE_META_NOTES=""                          # optional, e.g. Cray PrgEnv quirks
}

# add to the registry list near the top of the database block:
MACHINE_REGISTRY="bebop improv crux gce perlmutter <name>"
```

Validate via `install-moab.sh --dry-run --machine=<name>` (no real build needed
for the resolution step). Promote `last_validated` to a date once a real build
on the target machine completes successfully.

### Why no paths

`suggest_configuration.sh` (the deprecated predecessor) had hardcoded paths for
six retired machines (vesta, mira, blogin, theta, cori, edison), all of which
silently misdirected users for years. The lesson: machine path data ages
faster than this script does. For Push 1, the only obligation is **detection
and informational hinting**. Push 2 delegates the path-resolution problem to
E3SM's `config_machines.xml`, which is maintained by the E3SM team and is
already the source of truth for E3SM builds.

---

## Directory layout

With `--build-dir=/scratch/me/build` and `--prefix=$HOME/install/MOAB`:

```
/scratch/me/build/                               BUILD_DIR
├── moab-src/                                    MOAB clone (script-managed)
├── tpl-archives/                                cached TPL tarballs (re-used across runs)
├── tpl-work/
│   ├── eigen3/src/                              extracted Eigen sources
│   ├── zoltan/{src,build}/                      Zoltan sources + out-of-tree build dir
│   └── tempestremap/{src,build}/                TempestRemap clone + out-of-tree build dir
├── tpl-logs/
│   ├── config_<tpl>.log                         per-phase per-TPL logs
│   ├── build_<tpl>.log
│   └── install_<tpl>.log
└── moab/                                        MOAB's own out-of-tree build dir
    ├── Makefile (autotools) or CMakeCache.txt (cmake)
    ├── config.status
    └── .install-moab.config_fingerprint         used for resume detection

$HOME/install/MOAB/                              PREFIX_PATH
├── tpls/                                        TPL_PREFIX (default: $PREFIX_PATH/tpls)
│   ├── eigen3/include/eigen3/Eigen/Dense
│   ├── zoltan/{include/zoltan.h, lib/libzoltan.*}
│   └── tempestremap/{include/TempestRemapAPI.h, lib/libTempestRemap.*}
├── include/moab/Core.hpp
└── lib/libMOAB.{a,so,dylib}
```

`BUILD_DIR` defaults are situational:

- **No `--src-dir` and no `--build-dir`** → `BUILD_DIR=$PWD/moab-build`,
  `MOAB_SRC_DIR=$BUILD_DIR/moab-src` (script-managed clone)
- **`--src-dir=PATH` set** → `BUILD_DIR=$MOAB_SRC_DIR/build-e3sm` (next to user's source)
- **`--build-dir=PATH` set** → `MOAB_SRC_DIR=$BUILD_DIR/moab-src` unless also overridden

---

## All command-line options

### Source and install paths

| Flag | Default | Description |
|---|---|---|
| `--prefix=PATH` | `$HOME/install/MOAB` | MOAB install prefix |
| `--tpl-prefix=PATH` | `$PREFIX_PATH/tpls` | Where TPLs install (each in its own subdir) |
| `--build-dir=PATH` | `$PWD/moab-build` (or `$MOAB_SRC_DIR/build-e3sm` if `--src-dir` set) | All build artifacts root |
| `--src-dir=PATH` | `$BUILD_DIR/moab-src` (script-managed) | Use this MOAB checkout. **No git operations performed.** |

### Source management

| Flag | Default | Description |
|---|---|---|
| `--moab-repo=URL` | `https://bitbucket.org/fathomteam/moab.git` | Git remote for MOAB |
| `--moab-branch=NAME` | `master` | Branch to track |
| `--no-source-update` | off | Skip `git fetch + reset --hard` on the script-managed clone (use whatever is on disk) |

### Build system & flags

| Flag | Default | Description |
|---|---|---|
| `--build-system=BS` | `autotools` | `autotools` or `cmake` |
| `--jobs=N` | nprocs | Parallel make jobs |
| `--shared` | off | Build shared libs |
| `--no-static` | off | Disable static (implies `--shared`) |
| `--no-debug` | off | Disable debug symbols |
| `--no-optimize` | off | Disable compiler optimization |
| `--extra=ARGS` | empty | Additional args appended verbatim to MOAB's configure or `cmake` invocation |
| `--extra-zoltan=ARGS` | empty | Additional args appended verbatim to **Zoltan's** `configure` invocation (e.g. `LIBS=-lopenblas`, `--with-gnumake`). Changing this between runs auto-triggers a Zoltan rebuild. |
| `--extra-tempestremap=ARGS` | empty | Additional args appended verbatim to **TempestRemap's** `configure` invocation. Use this for site-specific BLAS/LAPACK selection — e.g. `--extra-tempestremap="--with-blas=-L/opt/mkl/lib -lmkl_rt --with-lapack=-L/opt/mkl/lib"` or `--extra-tempestremap="LIBS=-lessl -lxlf90_r"`. Changing this between runs auto-triggers a TempestRemap rebuild. |

### TPL provider paths

| Flag | Description |
|---|---|
| `--mpi-root=PATH` | MPI install root. If set, `$MPI_ROOT/bin/mpi*` is preferred, and `--with-mpi=$MPI_ROOT` is passed to MOAB. |
| `--hdf5-root=PATH` | Required (or via `HDF5_ROOT`). |
| `--netcdf-root=PATH` | Required (or via `NETCDF_C_PATH`). |
| `--pnetcdf-root=PATH` | Required (or via `PNETCDF_PATH`). |

### Compiler overrides

| Flag | Description |
|---|---|
| `--cc=BIN` | Explicit C compiler (binding — overrides MPI auto-detection) |
| `--cxx=BIN` | Explicit C++ compiler |
| `--fc=BIN` | Explicit Fortran compiler |
| `--f77=BIN` | Explicit F77 compiler |

### Machine selection

| Flag | Default | Description |
|---|---|---|
| `--machine=NAME` | `auto` | Use the named entry from the machine database (`bebop`, `improv`, `crux`, `gce`, `perlmutter`, …). `auto` runs `detect_machine` against `LMOD_SYSTEM_NAME`/`NERSC_HOST`/hostname and silently falls through if no entry matches. |
| `--compiler=NAME` | entry's `default_compiler` | Compiler family on the chosen machine: `gnu`, `intel`, `cray`, `nvhpc`, `nvidia`, `aocc`. Warning (not error) if not in the entry's `supported_compilers`. With `--profile=e3sm`, this keys the `<modules compiler="X">` and `<environment_variables compiler="X">` filters in `config_machines.xml`. |
| `--list-machines` | — | Print the registry (with auto-detected entry marked) and exit. Requires no env vars; safe to run on a login node before any modules are loaded. |

### Profile + E3SM env

| Flag | Default | Description |
|---|---|---|
| `--profile=NAME` | `e3sm` | `e3sm` (load env from CIME) or `standalone` (trust user env). |
| `--e3sm-root=PATH` | `$E3SM_ROOT` env | Path to E3SM checkout; required for `--profile=e3sm`. Must contain `cime_config/machines/config_machines.xml`. |
| `--yes`, `-y` | — | Bypass the e3sm-profile env-confirmation prompt. Required when stdin is not a tty (CI). Equivalent to `ASSUME_YES=yes`. |
| `--print` | — | Emit only the resolved MOAB configure/cmake command (copy-pasteable). Implies `--dry-run`. |

### Resume / cleanup controls

| Flag | Default | Description |
|---|---|---|
| `--clean` | off | Wipe `BUILD_DIR` before doing anything (full restart) |
| `--clean-tpls` | off | Wipe `TPL_PREFIX` and force TPL rebuilds |
| `--no-reuse-tpls` | off | Always rebuild TPLs even if install marker exists |
| `--reconfigure` | off | Force MOAB re-configure even if build dir is already configured and fingerprint matches |

### Behavior toggles

| Flag | Default | Description |
|---|---|---|
| `--check` | off | Run `make check` (autotools) or `ctest` (CMake) after build |
| `--skip-mpi-validation` | off | Skip the MPI compile-and-link sanity test (only if you can't compile on the host, e.g. a login node with no compiler license) |
| `--no-tail` | off | Suppress live tail of TPL build logs |
| `--dry-run` | off | Print all planned commands and resume state, but do not execute |
| `-h`, `--help` | — | Show help and exit |

---

## All environment variables

Every flag has an env-var equivalent. Env vars are read once at startup; CLI
flags override them.

| Env var | Equivalent flag | Notes |
|---|---|---|
| `PREFIX_PATH` | `--prefix` | |
| `TPL_PREFIX` | `--tpl-prefix` | |
| `BUILD_DIR` | `--build-dir` | |
| `MOAB_SRC_DIR` | `--src-dir` | Setting this in env is treated as "user-managed" — no git ops. |
| `MOAB_REPO_URL` | `--moab-repo` | |
| `MOAB_BRANCH` | `--moab-branch` | |
| `NO_SOURCE_UPDATE` | `--no-source-update` | `yes`/`no` |
| `BUILD_SYSTEM` | `--build-system` | `autotools` / `cmake` |
| `JOBS` | `--jobs` | |
| `MPI_ROOT` | `--mpi-root` | Empty value (e.g. `--mpi-root=$MPI_ROOT` where the var is unset) is detected, warned about, and treated as unset. |
| `HDF5_ROOT` | `--hdf5-root` | Required. |
| `NETCDF_C_PATH` | `--netcdf-root` | Required. |
| `PNETCDF_PATH` | `--pnetcdf-root` | Required. |
| `CC`, `CXX`, `FC`, `F77` | `--cc/--cxx/--fc/--f77` | **Env-set values are HINTS only** — overridden by MPI wrappers from MPI_ROOT or PATH. Use `--cc=…` for a binding override. |
| `ENABLE_DEBUG` | `--no-debug` | `yes`/`no`, default `yes` |
| `ENABLE_OPTIMIZE` | `--no-optimize` | `yes`/`no`, default `yes` |
| `ENABLE_SHARED` | `--shared` | `yes`/`no`, default `no` |
| `ENABLE_STATIC` | `--no-static` | `yes`/`no`, default `yes` |
| `RUN_CHECK` | `--check` | |
| `CLEAN_BUILD` | `--clean` | |
| `CLEAN_TPLS` | `--clean-tpls` | |
| `REUSE_TPLS` | `--no-reuse-tpls` (negation) | default `yes` |
| `RECONFIGURE` | `--reconfigure` | |
| `SKIP_MPI_VALIDATION` | `--skip-mpi-validation` | |
| `TAIL_LOGS` | `--no-tail` (negation) | default `yes` |
| `DRY_RUN` | `--dry-run` | |
| `EXTRA_MOAB_ARGS` | `--extra` | |
| `EXTRA_ZOLTAN_ARGS` | `--extra-zoltan` | Verbatim extras for Zoltan configure. Stamped at `$TPL_PREFIX/zoltan/.install-moab.extra_args` for change detection. |
| `EXTRA_TEMPESTREMAP_ARGS` | `--extra-tempestremap` | Verbatim extras for TempestRemap configure. Stamped at `$TPL_PREFIX/tempestremap/.install-moab.extra_args`. |
| `DOWNLOAD_RETRIES` | — | Tarball download retry count, default 5 |
| `DOWNLOAD_RETRY_DELAY` | — | Seconds between download retries, default 10 |
| `TPL_MAX_ATTEMPTS` | — | Per-TPL configure/build/install retry budget, default 3 |
| `TPL_MIRROR` | — | Override default ANL TPL mirror |

---

## Compiler resolution rules

The script resolves each of `CC`, `CXX`, `FC`, `F77` independently using
the following priority:

1. **`--cc/--cxx/--fc/--f77` flag** — binding. If you pass an explicit flag,
   that compiler is used as-is, regardless of MPI auto-detection.
2. **`$MPI_ROOT/bin/<wrapper>`** — if `MPI_ROOT` is set and the canonical
   wrapper exists.
3. **`$MPI_ROOT/bin/<alternate>`** — alternates per language:
   - C: (none)
   - C++: `mpic++`, `mpiCC`, `CC`
   - F90/F77: `mpifort`, `mpiifort`, `ftn`
4. **PATH lookup of the canonical wrapper** (`mpicc/mpicxx/mpif90/mpif77`).
   This **overrides env-set non-MPI compilers**, with a log message
   explaining the override.
5. **PATH lookup of alternates**.
6. **Env-set `CC`/`CXX`/`FC`/`F77` value** — last resort, will likely fail
   the MPI compile test.
7. Hard error listing everything tried.

**Why env-set CC is demoted:** on Spack-managed systems, `module load gcc-12`
sets `CC=/path/to/gcc`, which is not an MPI wrapper. Treating env CC as a
binding override caused MPI auto-detection to be silently bypassed. Now,
env CC is a hint that the script may override when an MPI wrapper is
available; pass `--cc=…` if you really need a binding override.

### Cray exception

When the script detects a Cray PrgEnv environment **and** `MPI_ROOT` is not
set, the search order is reordered so that the canonical wrappers are
`cc/CC/ftn` and `mpicc/mpicxx/mpif90` are tried only as alternates. See
[Cray PrgEnv detection](#cray-prgenv-detection).

### F77 fallback

If no `mpif77`/`mpifort`/`mpiifort`/`ftn` is found, the script silently
reuses `FC` as `F77` (with a warning). Modern systems often don't ship a
distinct F77 compiler.

### Empty `MPI_ROOT`

If you invoke `--mpi-root=$MPI_ROOT` and your shell expands that to empty
(because `MPI_ROOT` isn't set in your env), the script emits:

```
WARNING: MPI_ROOT is set but empty (did you forget to load the MPI module
before invoking the script?). Treating as unset.
```

…and proceeds with PATH/Cray fallback resolution.

---

## MPI wrapper validation

After compiler resolution and **before any TPL or MOAB work**, each resolved
wrapper is exercised with a tiny MPI compile-and-link test:

- **C**: `mpi.h` + `MPI_Init` + `MPI_Comm_rank` + `MPI_Finalize`
- **C++**: same, with `<iostream>`
- **Fortran**: fixed-form, `include 'mpif.h'` + `MPI_Init/Finalize`

If any test fails, the script aborts with the actual compiler output
printed inline (not buried in a log file you have to find).

Why this matters: a misnamed wrapper, missing module, broken PrgEnv, or
unresolved `mpi.h` would otherwise only fail much later inside a TPL's
configure or MOAB's link step — wasting hours.

To bypass (e.g., login node with no compute license), pass
`--skip-mpi-validation`. The build will proceed and likely fail later in
the same way, but you'll have an explicit warning.

---

## MOAB source management

### Two modes

#### Script-managed (default)

If `--src-dir` is **not** set, the script owns the MOAB checkout at
`$BUILD_DIR/moab-src`:

- **First run**: `git clone --branch master --single-branch
  https://bitbucket.org/fathomteam/moab.git $BUILD_DIR/moab-src`
- **Subsequent runs**: `git fetch origin master && git checkout master &&
  git reset --hard origin/master` (forced fast-forward to upstream)
- After clone/update, if `configure` is missing, runs `autoreconf -fi`
  to generate it.

To opt out of the update (e.g., to pin a specific commit you've already
checked out manually): `--no-source-update`.

To use a custom remote or branch: `--moab-repo=URL --moab-branch=NAME`.

#### User-managed (when `--src-dir=PATH` is set)

The directory at `PATH` is used as-is. **Zero git operations are performed.**
This is the right mode for:

- Development checkouts on a feature branch with uncommitted changes
- Pre-fetched tarball extractions (non-git source trees)
- Pinned commits the user has already checked out

`autoreconf -fi` is still run if `configure` is missing.

### Refusal cases

If `$BUILD_DIR/moab-src` exists but isn't a `.git` directory and isn't
empty, the script refuses to proceed — won't nuke a directory it doesn't
own. Either remove it manually or use `--src-dir=PATH` to point elsewhere.

### Source HEAD in the configure fingerprint

The MOAB source HEAD's git hash is included in the configure fingerprint
(see [Resume behavior](#resume-behavior)), so a fresh `git pull` that
moves HEAD automatically forces a MOAB re-configure on the next run, even
if no flags changed.

---

## TPL build pipeline

Three TPLs are managed entirely by the orchestrator. MOAB's own
`--download-*` macros are **not** used.

### Eigen3 (3.4.0)

- **Source**: tarball from `$TPL_MIRROR/eigen-3_4_0.tar.gz`, with GitLab
  fallback at `https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz`
- **Build**: header-only — extract and copy `Eigen/`, `unsupported/` to
  `$TPL_PREFIX/eigen3/include/eigen3/`
- **Install marker**: `$TPL_PREFIX/eigen3/include/eigen3/Eigen/Dense`

### Zoltan (3.9.1)

- **Source**: tarball from `$TPL_MIRROR/zoltan-3_9_1.tar.gz`, with GitHub
  fallback at `https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz`
- **Build**: autotools out-of-tree:
  ```
  configure --prefix=$TPL_PREFIX/zoltan --libdir=$TPL_PREFIX/zoltan/lib \
            --with-pic=1 --enable-mpi \
            --enable-shared=$ENABLE_SHARED --enable-static=$ENABLE_STATIC \
            CC=$CC CXX=$CXX FC=$FC F77=$F77
  make -j$JOBS
  make install
  ```
- **Install marker**: `$TPL_PREFIX/zoltan/include/zoltan.h`

### TempestRemap

- **Source** (preferred): `git clone --depth 1 --branch master --single-branch
  https://github.com/E3SM-Project/tempestremap.git`
- **Source fallback** (only if `git` isn't on PATH): tarball download from
  `https://github.com/E3SM-Project/tempestremap/archive/refs/heads/master.tar.gz`
- **Build**: autotools out-of-tree:
  ```
  configure --prefix=$TPL_PREFIX/tempestremap \
            --libdir=$TPL_PREFIX/tempestremap/lib \
            --with-pic=1 --enable-shared=$ENABLE_SHARED --enable-static=$ENABLE_STATIC \
            --with-netcdf=$NETCDF_C_PATH --with-hdf5=$HDF5_ROOT \
            CC=$CC CXX=$CXX FC=$FC F77=$F77
  make -j$JOBS
  make install
  ```
- **Install marker**: `$TPL_PREFIX/tempestremap/include/TempestRemapAPI.h`

### Common TPL behaviors

- **Source preparation sentinel**: `$srcdir/.install-moab.source_ok` written
  only after a verified clone or extract. Partial extractions are detected
  and re-done.
- **Per-TPL retry loop**: each TPL gets up to `TPL_MAX_ATTEMPTS` (default 3)
  configure/build/install cycles. After each failure, the per-phase log is
  scanned for known patterns and a targeted remediation is applied (see
  [Failure remediation](#failure-remediation)).
- **Build dir wiped on retry** to avoid contamination from a previous
  failed attempt.
- **Reuse on subsequent runs**: by default (`REUSE_TPLS=yes`), if a TPL's
  install marker exists, the entire build/install phase is skipped. Force a
  rebuild with `--no-reuse-tpls`, or wipe with `--clean-tpls`.

---

## Build-system selection

### `--build-system=autotools` (default)

```
$MOAB_SRC_DIR/configure \
    --with-mpi=$MPI_ROOT                 # or bare --with-mpi if MPI_ROOT unset
    CC=$CC CXX=$CXX FC=$FC F77=$F77 \
    --with-hdf5=$HDF5_ROOT \
    --with-netcdf=$NETCDF_C_PATH \
    --with-pnetcdf=$PNETCDF_PATH \
    --with-eigen3=$TPL_PREFIX/eigen3/include/eigen3 \
    --with-zoltan=$TPL_PREFIX/zoltan \
    --with-tempestremap=$TPL_PREFIX/tempestremap \
    --enable-optimize --enable-debug \
    --disable-shared --enable-static \
    --prefix=$PREFIX_PATH

make -j$JOBS
[make -j$JOBS check]
make install
```

### `--build-system=cmake`

```
cmake -S $MOAB_SRC_DIR -B $BUILD_DIR/moab \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX=$PREFIX_PATH \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC \
    -DENABLE_MPI=ON -DMPI_HOME=$MPI_ROOT \
    -DENABLE_HDF5=ON -DHDF5_ROOT=$HDF5_ROOT \
    -DENABLE_NETCDF=ON -DNETCDF_DIR=$NETCDF_C_PATH \
    -DENABLE_PNETCDF=ON -DPNETCDF_DIR=$PNETCDF_PATH \
    -DENABLE_EIGEN3=ON \
        -DEIGEN3_DIR=$TPL_PREFIX/eigen3/include/eigen3 \
        -DEIGEN3_INCLUDE_DIR=$TPL_PREFIX/eigen3/include/eigen3 \
    -DENABLE_ZOLTAN=ON -DZOLTAN_DIR=$TPL_PREFIX/zoltan \
    -DENABLE_TEMPESTREMAP=ON -DTEMPESTREMAP_DIR=$TPL_PREFIX/tempestremap

cmake --build $BUILD_DIR/moab --parallel $JOBS
[ctest --output-on-failure -j$JOBS]
cmake --install $BUILD_DIR/moab
```

`CMAKE_BUILD_TYPE` is derived from `ENABLE_DEBUG`/`ENABLE_OPTIMIZE`:

| `ENABLE_DEBUG` | `ENABLE_OPTIMIZE` | `CMAKE_BUILD_TYPE` |
|:-:|:-:|:--|
| yes | yes | RelWithDebInfo |
| yes | no | Debug |
| no | yes | Release |

---

## Resume behavior

**Re-running the script with the same flags is a no-op for already-completed
work.** Four levels of granularity:

1. **TPL install markers** — if `$TPL_PREFIX/<tpl>/<known-header>` exists,
   the entire TPL build is skipped.
2. **TPL extraction sentinel** — `$srcdir/.install-moab.source_ok` written
   only after a verified extract; partial extractions get re-done on retry.
3. **MOAB configured?** — checked by presence of `Makefile` + `config.status`
   (autotools) or `CMakeCache.txt` (cmake).
4. **MOAB args fingerprint** — `sha256(<build-system>|src=<git-HEAD>|<args…>)`,
   stored at `$MOAB_BUILD_DIR/.install-moab.config_fingerprint`. Compared on
   each run; if it differs, MOAB is automatically re-configured.

The fingerprint includes the MOAB source git HEAD, so a `git pull` that
advances HEAD forces a re-configure even without any flag changes.

After re-configure is skipped, `make` (autotools) or `cmake --build` (CMake)
runs and **incrementally rebuilds only the changed translation units**.
Then `make install` / `cmake --install` runs.

### State banner

Before doing any work, the script prints a state report:

```
========== Resume state ==========

[install-moab]   eigen3:       installed at /home/me/install/MOAB/tpls/eigen3 -- WILL SKIP
[install-moab]   zoltan:       installed at /home/me/install/MOAB/tpls/zoltan -- WILL SKIP
[install-moab]   tempestremap: not installed -- WILL BUILD
[install-moab]   moab:         configured (fingerprint match) -- WILL SKIP CONFIGURE, resume make/install
```

### Forcing a re-do

| What you want | Flag |
|---|---|
| Force MOAB re-configure (TPLs unchanged) | `--reconfigure` |
| Force one TPL rebuild | Delete its install marker, e.g. `rm $TPL_PREFIX/zoltan/include/zoltan.h` |
| Force all TPL rebuilds | `--no-reuse-tpls` |
| Wipe TPL installs | `--clean-tpls` |
| Wipe everything in `BUILD_DIR` and start over | `--clean` |
| Skip the git pull on script-managed clone | `--no-source-update` |

---

## Failure remediation

When a TPL phase fails, the per-phase log (`$BUILD_DIR/tpl-logs/<phase>_<tpl>.log`)
is scanned for known failure fingerprints and the script applies a targeted
remediation, then retries.

| Pattern detected (regex on log) | Remediation | Common cause |
|---|---|---|
| `Could not resolve host` / `Connection (timed out\|refused)` / `Network is unreachable` | Delete cached archive, re-download | Transient network |
| `multiple definition of` / `undefined reference to Zoltan_*` / `relocation R_X86_64...` | Add `CFLAGS=-fcommon` and rebuild | GCC ≥10 default `-fno-common` breaks Zoltan 3.9.1 |
| `aclocal not found` / `cannot find install-sh` | Run `autoreconf -fi` in the TPL src dir | Fresh git clone missing autotools-generated files |
| `netcdf.h: No such file` / `Cannot find NetCDF` | Export `NETCDFROOT/NETCDF_DIR/NETCDF_PATH` for TempestRemap | TempestRemap doesn't see `--with-netcdf` because of bad pkg-config or env |
| C++14 errors (`requires C++14`, `range-based for loops are not allowed`) | Force `CXXFLAGS=-std=c++14` | Old default C++ standard on the system compiler |
| `cannot find -llapack` / undefined `dgemm_` / `dsyev_` | Append `-llapack -lblas` to `LDFLAGS` | TempestRemap needs LAPACK on systems without auto-detection |
| `Permission denied` / `Read-only file system` / `No space left` | Reported as unrecoverable; no retry | Filesystem problem — user must fix |

The retry budget is `TPL_MAX_ATTEMPTS` (default 3). After exhaustion, the
script aborts with the path to the per-phase log files.

If no known pattern matches a failure, the script aborts immediately
rather than retrying blindly.

---

## Cray PrgEnv detection

The script detects a Cray Programming Environment when **any** of these is
true:

- `PE_ENV` is set (set by `module load PrgEnv-{gnu,cray,intel,aocc,nvhpc}`)
- `CRAYPE_VERSION` is set (set by `craype` itself)
- `/opt/cray` exists

When detected **and** `MPI_ROOT` is not set, the compiler resolver swaps
the search order:

| | Default search order | Cray search order |
|---|---|---|
| C | `mpicc` → `cc` | `cc` → `mpicc` |
| C++ | `mpicxx` → `mpic++ mpiCC CC` | `CC` → `mpicxx mpic++ mpiCC` |
| Fortran | `mpif90` → `mpifort mpiifort ftn` | `ftn` → `mpif90 mpifort mpiifort` |
| F77 | `mpif77` → `mpifort mpiifort ftn` | `ftn` → `mpif77 mpifort mpiifort` |

The user is informed:
```
[install-moab] Cray PrgEnv detected (PE_ENV=GNU CRAYPE_VERSION=2.7.30); preferring cc/CC/ftn wrappers
```

If the user **does** pass `--mpi-root` on a Cray, that's interpreted as
"I'm using a non-craype MPI" and the standard order is used.

---

## Dry-run mode

`--dry-run` prints everything that would happen and exits without doing
work — useful for sanity-checking a new machine.

What `--dry-run` shows:

- The full orchestration banner (resolved compilers, paths, knobs)
- The MPI validation step is **skipped** (so dry-run works without an MPI
  install)
- The Resume state report
- Per-TPL recipes (source URL or git repo, install location)
- The exact MOAB invocation (`configure …` or `cmake …`) with all args

Example (truncated):

```
$ ./install-moab-e3sm.sh --dry-run --hdf5-root=/opt/hdf5 --netcdf-root=/opt/netcdf --pnetcdf-root=/opt/pnetcdf

========== MOAB orchestration ==========

[install-moab] MOAB source dir : /home/me/moab-build/moab-src (script-managed: https://bitbucket.org/fathomteam/moab.git @ master)
[install-moab] Build dir       : /home/me/moab-build
[install-moab] MOAB install    : /home/me/install/MOAB
[install-moab] TPL install     : /home/me/install/MOAB/tpls
[install-moab] Build system    : autotools
[install-moab] Compilers       : CC=/opt/mpich/bin/mpicc CXX=/opt/mpich/bin/mpicxx FC=/opt/mpich/bin/mpif90 F77=/opt/mpich/bin/mpif77

========== Resume state ==========

[install-moab]   eigen3:       not installed -- WILL BUILD
[install-moab]   zoltan:       not installed -- WILL BUILD
[install-moab]   tempestremap: not installed -- WILL BUILD
[install-moab]   moab:         not configured -- WILL CONFIGURE

[install-moab] Per-TPL recipes (would run, in order):
[install-moab]   eigen3:       https://web.cels.anl.gov/.../eigen-3_4_0.tar.gz -> .../tpls/eigen3
[install-moab]   zoltan:       https://web.cels.anl.gov/.../zoltan-3_9_1.tar.gz -> .../tpls/zoltan
[install-moab]   tempestremap: git --depth 1 https://github.com/E3SM-Project/tempestremap.git@master -> .../tpls/tempestremap

[install-moab] MOAB invocation (autotools):
  /home/me/moab-build/moab-src/configure \
      --with-mpi=/opt/mpich \
      CC=/opt/mpich/bin/mpicc \
      [...]
```

---

## Verification

After install, the script checks for known artifacts:

- `$PREFIX_PATH/include/moab/Core.hpp`
- `$PREFIX_PATH/lib/libMOAB.{a,so,dylib}` (any one suffices)
- `$TPL_PREFIX/eigen3/include/eigen3/Eigen/Dense`
- `$TPL_PREFIX/zoltan/include/zoltan.h`
- `$TPL_PREFIX/zoltan/lib/libzoltan.{a,so,dylib}`
- `$TPL_PREFIX/tempestremap/include/TempestRemapAPI.h`
- `$TPL_PREFIX/tempestremap/lib/libTempestRemap.{a,so,dylib}`

Each missing artifact is logged as a warning. The script does not abort on
verification failure (the install completed; just inform the user).

The final message lists environment variables to add to your downstream
shell:

```
[install-moab] MOAB installed at: /home/me/install/MOAB
[install-moab] TPLs installed under: /home/me/install/MOAB/tpls
[install-moab] Add the following to your environment for downstream use:
[install-moab]   export MOAB_ROOT=/home/me/install/MOAB
[install-moab]   export PATH=$MOAB_ROOT/bin:$PATH
[install-moab]   export PKG_CONFIG_PATH=$MOAB_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH
[install-moab]   export ZOLTAN_DIR=/home/me/install/MOAB/tpls/zoltan
[install-moab]   export TEMPESTREMAP_DIR=/home/me/install/MOAB/tpls/tempestremap
[install-moab]   export EIGEN3_DIR=/home/me/install/MOAB/tpls/eigen3/include/eigen3
```

---

## Logging and live progress

Each TPL phase writes its full output to a log file:

```
$BUILD_DIR/tpl-logs/config_<tpl>.log
$BUILD_DIR/tpl-logs/build_<tpl>.log
$BUILD_DIR/tpl-logs/install_<tpl>.log
```

By default (`TAIL_LOGS=yes`), while a phase is running the script tails the
log file in the background and prefixes each line with the log name:

```
  [eigen3-install] copying ...
  [zoltan-config]  checking for gcc... /opt/mpich/bin/mpicc
  [zoltan-build]   make[1]: Entering directory '...'
```

Disable with `--no-tail` if the live output is noisy.

MOAB's own configure output goes to `$MOAB_BUILD_DIR/config.log` (autotools)
or to stdout/stderr (CMake). On configure failure, the last 60 lines of
`config.log` are printed automatically.

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | usage / argument error (bad flag, missing src dir when `--src-dir` used, etc.) |
| 2 | environment / dependency missing (no `git`, no `cmake`, missing required env var, MPI validation failed, compiler not found) |
| 3 | configure / build / install / clone failure |

---

## Common workflows

### First-time install on a new machine

```bash
# Inspect what would happen first
./install-moab-e3sm.sh --dry-run \
    --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH

# Run for real
./install-moab-e3sm.sh \
    --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH \
    --prefix=$HOME/install/MOAB --jobs=16
```

### Update MOAB to latest master

Same command. The script:
1. `git fetch + reset --hard origin/master`
2. Skips TPLs (markers exist)
3. Reconfigures MOAB (HEAD changed → fingerprint differs)
4. `make` rebuilds only changed files
5. `make install`

### Pin MOAB to a specific commit

```bash
# Option 1: track a different branch
./install-moab-e3sm.sh --moab-branch=release-v5.6.0 ...

# Option 2: check out manually and use --src-dir
git clone https://bitbucket.org/fathomteam/moab.git /scratch/me/moab-src
cd /scratch/me/moab-src && git checkout 91b54bd8e
./install-moab-e3sm.sh --src-dir=/scratch/me/moab-src ...
```

### Develop on a feature branch

```bash
# Use --src-dir to keep your branch and uncommitted changes safe
./install-moab-e3sm.sh --src-dir=/path/to/my/checkout ...
# After editing source: just re-run -- make builds incrementally
./install-moab-e3sm.sh --src-dir=/path/to/my/checkout ...
```

### Rebuild MOAB only (TPLs already done)

```bash
# Just re-run; TPLs are skipped automatically
./install-moab-e3sm.sh ...

# Force MOAB re-configure (e.g., after editing configure.ac)
./install-moab-e3sm.sh --reconfigure ...
```

### Switch from autotools to CMake

```bash
./install-moab-e3sm.sh --build-system=cmake \
    --build-dir=$BUILD_DIR-cmake \      # use a separate build dir
    ...
```

(Both build dirs can coexist; TPLs are shared.)

### CI / scripted build

```bash
./install-moab-e3sm.sh \
    --hdf5-root=$HDF5_ROOT --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH \
    --prefix=$WORK/install/MOAB \
    --build-dir=$WORK/build/moab \
    --jobs=$(nproc) \
    --check \
    --no-tail \
    || { echo "MOAB build failed; see $WORK/build/moab/tpl-logs/ and $WORK/build/moab/moab/config.log"; exit 1; }
```

---

## Troubleshooting

### `--profile=e3sm` errors

**`--profile=e3sm requires --e3sm-root=PATH or $E3SM_ROOT`**
Set the env var or pass the flag. The path must contain
`cime_config/machines/config_machines.xml` — i.e., the root of an E3SM
checkout, not the CIME submodule.

**`machine 'X' has no E3SM config_machines.xml mapping`**
The machine entry in `install-moab.sh` has an empty
`MACHINE_META_E3SM_NAME`. Either fix the entry (look up `MACH=` in
`config_machines.xml`) or use `--profile=standalone`.

**`HDF5_ROOT is not set` after the e3sm env loaded**
config_machines.xml deliberately doesn't export `HDF5_ROOT` for some
machines (bebop, improv, pm-cpu, crux). The e3sm_env.py output prints a
`# NOTE` explaining this. Pass `--hdf5-root=PATH` explicitly. For Cray
PrgEnv, the script also auto-aliases `$CRAY_HDF5_PARALLEL_PREFIX` →
`HDF5_ROOT` after sourcing.

**`e3sm_env.py exited non-zero`**
Check `$BUILD_DIR/.../e3sm_env.*.sh` (the temp file path is logged). Most
common cause: machine name not in `config_machines.xml`. Verify with
`grep MACH= $E3SM_ROOT/cime_config/machines/config_machines.xml`.

**Tier 1 (CIME API) failed; falling back to direct XML parse**
Cosmetic warning. The direct XML parse handles all current schema features
correctly; verified byte-identical output for bebop. The only practical
difference is that CIME's API would catch some schema evolution
automatically, while the fallback may need an update if CIME ever
changes the `config_machines.xml` schema significantly.

**Modules fail to load when sourcing the e3sm env**
The emitted snippet runs `module load` for each module in
config_machines.xml. If one is unavailable on the target machine, check:
- That you're on the right machine (the entry in `config_machines.xml`
  may name modules that only exist on that specific cluster)
- That your shell has `module` available (the snippet sources Lmod's
  init script first; if that init script doesn't exist, the warning is
  printed but bash continues — the next `module load` will fail)
- That E3SM's master is the version you expect; module names get bumped
  over time

**Want to skip the e3sm env loading just this once?**
Pass `--profile=standalone` and load modules yourself. The machine entry's
`standalone_hint` printed in the orchestration banner shows what to load.

### "MPI_ROOT is set but empty" warning

You ran `--mpi-root=$MPI_ROOT` but `$MPI_ROOT` isn't set in your shell. The
script falls back to PATH-based MPI wrappers. If that's not what you want,
either:
- Load the module that sets `MPI_ROOT`
- Or pass an explicit path: `--mpi-root=/opt/mpich`

### MPI validation fails with "mpi.h: No such file or directory"

The compiler resolved is not an MPI wrapper. Check the `Compilers:` line in
the orchestration banner. Common causes:

1. **No MPI module loaded** — `module load mpich` (or your site's
   equivalent) and re-run.
2. **You explicitly passed `--cc=gcc`** — that's binding; remove it.
3. **`--mpi-root` was empty** — see above.

### "configure script missing"

For autotools mode, MOAB's `configure` is generated by `autoreconf -fi`.
The script runs that automatically — but `autoreconf` must be on PATH.
On Ubuntu: `apt install autoconf automake libtool m4`. On Cray:
`module load autoconf-archive` (or similar).

### Zoltan fails with "multiple definition of"

GCC ≥10 defaults to `-fno-common`, which breaks Zoltan 3.9.1. The script
detects this and adds `CFLAGS=-fcommon` automatically on the next attempt.
If the retry budget is exhausted before the fix lands, manually:

```bash
rm -rf $BUILD_DIR/tpl-work/zoltan
CFLAGS="-fcommon" ./install-moab-e3sm.sh ... --extra="CFLAGS=-fcommon"
```

### Site-specific BLAS / LAPACK (MKL, OpenBLAS, ESSL, libsci, AOCL, …)

Both TempestRemap and MOAB use **autoconf-archive's `AX_BLAS` / `AX_LAPACK`**
macros, which expose `--with-blas=<libs-string>` and `--with-lapack=<libs-string>`.
The `<libs-string>` is the LIBS spec used to link a tiny test program. The
script's automatic `lapack` remediation only adds `-llapack -lblas`, which
fails on systems where the host BLAS/LAPACK is MKL, ESSL, OpenBLAS at a
non-default path, Cray libsci, AMD AOCL, netlib-lapack at a custom prefix, etc.
Pass the right linker line explicitly with `--extra-tempestremap=` (TempestRemap)
and `--extra=` (MOAB).

**Quoting matters.** Each `--with-blas=...` / `--with-lapack=...` must arrive
as a single argv entry, even when the value contains spaces. Use shell-quoted
inner values inside outer single quotes:

```bash
# Static netlib-lapack at a custom prefix (e.g. Bebop/Improv at ANL)
BLAS_SPEC='/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/libblas.a -lgfortran'
LAPACK_SPEC='/lcrc/group/e3sm/soft/improv/netlib-lapack/3.12.0/gcc-12.3.0/liblapack.a -lm'
./install-moab-e3sm.sh \
    --extra-tempestremap="--with-blas=\"$BLAS_SPEC\" --with-lapack=\"$LAPACK_SPEC\"" \
    --extra="--with-blas=\"$BLAS_SPEC\" --with-lapack=\"$LAPACK_SPEC\"" \
    ...

# Intel MKL (sequential)
MKL_SPEC="-L$MKLROOT/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
./install-moab-e3sm.sh \
    --extra-tempestremap="--with-blas=\"$MKL_SPEC\" --with-lapack=\"$MKL_SPEC\"" \
    --extra="--with-blas=\"$MKL_SPEC\" --with-lapack=\"$MKL_SPEC\"" \
    ...

# OpenBLAS at a custom prefix
./install-moab-e3sm.sh \
    --extra-tempestremap='--with-blas="-L/opt/openblas/lib -lopenblas"' \
    --extra='--with-blas="-L/opt/openblas/lib -lopenblas"' \
    ...

# Cray libsci (PrgEnv: usually no extras — cc/CC/ftn auto-link libsci)
```

The script `eval`s the value of `--extra-*` so embedded shell quoting is
honored. Without inner quotes around a spaced value the words would be split
into broken half-tokens (e.g. `--with-blas=/path/lib.a` and a stray
`-lgfortran`), which is the failure mode you'll see if you forget them.

The script stamps the extras string at `$TPL_PREFIX/<tpl>/.install-moab.extra_args`
on a successful install. On a re-run, if the new `--extra-<tpl>` differs from
the stamp, that TPL is rebuilt automatically (no need to remember
`--no-reuse-tpls`). The Resume-state banner shows `WILL REBUILD` in that case.

### TempestRemap "Cannot find NetCDF"

The script auto-exports `NETCDFROOT/NETCDF_DIR/NETCDF_PATH` on detection.
If retries are exhausted:

```bash
NETCDFROOT=$NETCDF_C_PATH NETCDF_DIR=$NETCDF_C_PATH ./install-moab-e3sm.sh ...
```

### Re-running picks up no changes after a `git pull`

The fingerprint should detect that. If not, force it:

```bash
./install-moab-e3sm.sh --reconfigure ...
```

### "$BUILD_DIR/moab-src exists but is not a git clone"

You probably manually extracted MOAB into the script-managed location.
Either:
- Remove that directory: `rm -rf $BUILD_DIR/moab-src` and re-run (script
  will clone fresh)
- Or use `--src-dir=$BUILD_DIR/moab-src` to tell the script "this is
  user-managed, leave git alone"

### Want to test without committing to a long build

Use `--dry-run` to see exactly what would happen. To test just the source
prep without TPLs, you can manually edit the script and add `exit 0` after
the `prepare_moab_source` call.

### Logs

Everything is in `$BUILD_DIR/tpl-logs/` and `$MOAB_BUILD_DIR/config.log`
(or CMake stdout). Each TPL has three log files (`config_*`, `build_*`,
`install_*`).

---

## What this script does NOT do

- Does **not** build MPI, HDF5, NetCDF-C, or Parallel-NetCDF — those must
  be available already.
- Does **not** install any system packages (autotools, cmake, git, etc.) —
  those must be on PATH.
- Does **not** run `mpiexec` to test the binaries (login nodes typically
  can't run MPI). The MPI validation only compiles and links.
- Does **not** modify shell configuration files. The final summary prints
  the env vars you should add to your shell rc file.
- Does **not** push changes anywhere or modify the git remote of any
  existing checkout.

---

## File locations summary

```
/path/to/install-moab-e3sm.sh         the script
/path/to/INSTALL-MOAB-E3SM.md         this document

# All under $BUILD_DIR (default $PWD/moab-build):
$BUILD_DIR/moab-src/                  MOAB clone (script-managed)
$BUILD_DIR/tpl-archives/              cached tarballs
$BUILD_DIR/tpl-work/<tpl>/{src,build} per-TPL workspaces
$BUILD_DIR/tpl-logs/                  per-phase per-TPL logs
$BUILD_DIR/moab/                      MOAB out-of-tree build dir
    .install-moab.config_fingerprint  resume detection

# All under $PREFIX_PATH:
$PREFIX_PATH/include/                 MOAB headers
$PREFIX_PATH/lib/                     MOAB libraries
$PREFIX_PATH/tpls/<tpl>/              TPL installs
```
