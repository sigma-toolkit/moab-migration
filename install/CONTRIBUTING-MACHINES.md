# Adding a Machine Entry to `install-moab.sh`

This guide walks you through adding a new HPC machine to the `install-moab.sh`
machine database. Each entry is **pure metadata** (~15 lines). No TPL paths
get baked into the script — those come from E3SM's `config_machines.xml` (for
`--profile=e3sm`) or the user's loaded modules (`--profile=standalone`).

## When to add an entry

- You build MOAB regularly on a machine that isn't in `--list-machines`.
- A machine you use is in the registry but with `last_validated: TBD`, and
  you've confirmed a successful build.
- E3SM has added a new machine to `cime_config/machines/config_machines.xml`
  and you want install-moab.sh to support it.

You **don't** need an entry to use the script — `--profile=standalone` works
with any user-managed environment. Entries exist for convenience (auto-detection,
e3sm-profile shortcut, machine-specific notes).

## Anatomy of a machine entry

In `install/install-moab.sh`, each registered machine is two short bash
functions:

```bash
machine_<name>_match() {
    # Hostname/env-var predicate. Return 0 if this machine is the current host.
    [[ "${LMOD_SYSTEM_NAME:-}" == "<name>" ]] && return 0
    [[ "$(_hn)" == <hostname-prefix>* ]]
}

machine_<name>_meta() {
    MACHINE_META_E3SM_NAME="<name in config_machines.xml>"  # may be empty
    MACHINE_META_DEFAULT_COMPILER="gnu"           # gnu, intel, cray, nvhpc, nvidia, aocc
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel"  # comma-separated
    MACHINE_META_STANDALONE_HINT="module load gcc openmpi hdf5 netcdf-c parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="TBD"             # promote to YYYY-MM-DD after testing
    MACHINE_META_NOTES=""                         # optional: gotchas, GPU partition, etc.
}
```

Then add `<name>` to `MACHINE_REGISTRY` at the top of the database block.

### Field reference

| Field | Required? | Purpose |
|---|---|---|
| `MACHINE_META_E3SM_NAME` | Required for `--profile=e3sm` | Must equal the `MACH=` attribute in E3SM's `cime_config/machines/config_machines.xml`. Empty string disables the e3sm profile for this entry (force users to `--profile=standalone`). |
| `MACHINE_META_DEFAULT_COMPILER` | Required | Used when the user doesn't pass `--compiler=`. Must be one of the families CIME recognizes (`gnu`, `intel`, `cray`, `nvhpc`, `nvidia`, `aocc`). |
| `MACHINE_META_SUPPORTED_COMPILERS` | Required | Comma-separated list. The script warns (but doesn't error) if `--compiler=NAME` isn't in this list. Mirror the `<COMPILERS>` element in `config_machines.xml`. |
| `MACHINE_META_STANDALONE_HINT` | Recommended | Module-load command(s) printed informationally for `--profile=standalone` users. Not auto-executed. |
| `MACHINE_META_LAST_VALIDATED` | Required | `YYYY-MM-DD` if you've personally confirmed a clean build, else `TBD`. The `--list-machines` table surfaces this so users can spot stale entries. |
| `MACHINE_META_NOTES` | Optional | Free-form. Use for: HDF5_ROOT export gaps, Cray PrgEnv quirks, BLAS/LAPACK overrides, GPU partition vs CPU, etc. |

### Naming rules

- The registry name and function-name infix must be **bash identifier safe**:
  letters, digits, underscores. **No dashes.** (Function names with `-` cause
  parser ambiguities.)
- For machines that have dashes in their config_machines.xml MACH name (e.g.
  `pm-gpu`, `anlgce-ub22`), use a dashless registry name (`pmgpu`, `gce`)
  and put the original name in `MACHINE_META_E3SM_NAME`.

### Match-predicate guidance

Order, most specific to least:
1. `LMOD_SYSTEM_NAME` (Lmod sets this on most LCRC/NERSC systems)
2. `NERSC_HOST` (NERSC-specific)
3. `PE_ENV` / `CRAYPE_VERSION` / `/opt/cray` existence (Cray PrgEnv anywhere)
4. Hostname prefix match via `$(_hn)` helper

Avoid matches that could collide with another entry. If two machines share a
match (e.g., pm-cpu and pm-gpu both on perlmutter login nodes), make one
match unconditionally and require the other be selected explicitly via
`--machine=NAME` (return 1 from its match function). See `pmgpu` in
`install-moab.sh` for the pattern.

## Step-by-step

1. **Find the E3SM machine name.**
   ```bash
   grep MACH= $E3SM_ROOT/cime_config/machines/config_machines.xml | grep -i <your-machine>
   ```

2. **Sketch the match predicate.** On the target machine:
   ```bash
   echo "LMOD_SYSTEM_NAME=$LMOD_SYSTEM_NAME"
   echo "NERSC_HOST=$NERSC_HOST"
   echo "hostname=$(hostname)"
   ```
   Pick the most stable identifier (env vars first; hostname only if no env
   var is reliable).

3. **Edit `install/install-moab.sh`:**
   - Find the comment block `# Machine database` near line ~110.
   - Append `<name>` to `MACHINE_REGISTRY`.
   - Add `machine_<name>_match()` and `machine_<name>_meta()` near the
     other entries, alphabetically grouped.

4. **Validate without a build:**
   ```bash
   ./install/install-moab.sh --list-machines           # entry should appear
   ./install/install-moab.sh --dry-run --machine=<name> --profile=standalone
   ```
   The dry-run banner should report your entry. No env vars required.

5. **Validate the e3sm profile (if applicable):**
   ```bash
   ./install/install-moab.sh --dry-run --machine=<name> --profile=e3sm \
       --e3sm-root=$E3SM_ROOT --yes
   ```
   The script should source the env from config_machines.xml. Verify
   `Resolved roots` includes `HDF5_ROOT`, `NETCDF_C_PATH`, `PNETCDF_PATH`.
   If HDF5_ROOT is missing, you may need to pass `--hdf5-root=` (and add a
   note in `MACHINE_META_NOTES`).

6. **Real build:**
   ```bash
   ./install/install-moab.sh --machine=<name> --profile=e3sm \
       --e3sm-root=$E3SM_ROOT --prefix=$PWD/installs --yes
   ```

7. **Promote `last_validated`** from `TBD` to today's date (YYYY-MM-DD) in
   the same commit that fixes anything you discovered.

## Common gotchas

- **HDF5_ROOT not exported by config_machines.xml.** Many E3SM machines
  (bebop, improv, pm-cpu, crux) don't export `HDF5_ROOT` in their
  `<environment_variables>` block — they expect it to come from a separate
  HDF5 module. The script has two adapters for this:
  - Cray: `CRAY_HDF5_PARALLEL_PREFIX` → `HDF5_ROOT` (auto)
  - LCRC `.../soft/<machine>/netcdf-c/.../...` pattern: sniff sibling
    `.../soft/<machine>/hdf5/.../...` (auto)
  - Otherwise: pass `--hdf5-root=PATH` explicitly. Document in `MACHINE_META_NOTES`.

- **NETCDF_PATH vs NETCDF_C_PATH.** Some entries use `NETCDF_PATH`
  (anlgce-ub22, crux). The script auto-aliases. No action needed in your
  entry.

- **Compiler family mismatch.** If `--compiler=NAME` doesn't match the
  filters in `<modules compiler="X">` blocks, the e3sm_env.py helper will
  silently emit no module loads. Always test with the machine's actual
  default compiler first.

- **Module names changed in E3SM master.** E3SM rebases module versions
  occasionally. If your `last_validated` date is more than ~6 months old
  and the build now fails at `module load`, check for a recent change in
  `git log $E3SM_ROOT/cime_config/machines/config_machines.xml`.

## Removing an entry

When a machine is decommissioned:

1. Remove from `MACHINE_REGISTRY`.
2. Delete the two `machine_<name>_*` functions.
3. Search for any references in `install/INSTALL-MOAB.md` and update.
4. Note the removal in the commit message (so future maintainers can find
   the rationale via `git log`).

Don't leave `# DEPRECATED` stubs in the registry — the whole point of
removing the suggest_configuration.sh DB was that retired entries with
stale paths actively misdirect users.
