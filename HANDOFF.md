# iMOAB Rename And Build Handoff

## Scope

- Renamed and relocated all parallel iMOAB C++ and Fortran tests into `test/earthsystem/parallel` using the `iMOAB` spelling.
- Preserved the complete existing workflow and regression behavior in each moved test.
- Renamed the common support header to `iMOABCouplerUtils.hpp`; its header-defined functions are now `inline`, and immutable string arguments are const references.
- Added the self-contained C++17 example `examples/earthsystem/earth/ParallelRemapTemplate.cpp`.
- Updated CMake, Autotools, installed-example manifests, and regression executable paths.
- Raised the project C++ standard from C++14 to C++17 in both CMake and Autotools.
- Reworked CMake `distcheck` to package current Git-tracked working-tree files, configure an isolated build with the parent feature/toolchain settings, build, run all tests, install, build the installed examples as a downstream CMake project, uninstall, clean, verify the install prefix is file-free, verify the source snapshot is unchanged, and remove temporary state.
- Distribution archives use a lowercase `moab-<version>-<12-character-hash>[-dirty]` basename and contain a `GITHASH` file with the full HEAD hash, commit message, and tracked working-tree state. GPG signing is opt-in with `MOAB_DIST_SIGN=ON`.
- Fixed the CMake uninstall script for current CMake policy behavior and corrected its removal-result check.
- Added the missing `BoundaryDensity.hpp` to CMake and Autotools installed-example manifests.

## Verification

- Parallel CMake configuration completed in `build-cmake-modern`.
- CMake built all renamed C++ and Fortran iMOAB targets.
- CTest `^iMOAB`: 17/17 passed.
- Autotools regenerated successfully and built all renamed C++ and Fortran iMOAB targets.
- Standalone examples CMake build built `ParallelRemapTemplate`.
- `clang-tidy` passed for `ParallelRemapTemplate.cpp` with project checks promoted to errors; 330 third-party/system-header diagnostics were suppressed.
- macOS `leaks --atExit` on one-rank `iMOABPtest2`: 0 leaks.
- CMake `distcheck` completed the full tarball lifecycle with the MPI/HDF5/NetCDF/PNetCDF/Metis/ParMetis/Zoltan/TempestRemap configuration; all 117 tarball tests passed, every enabled installed C++ and Fortran example built against the installed package, uninstall left no files, source checksums were unchanged, and `_distcheck` was removed.

## Notes

- Regression config filenames and digest prefixes remain unchanged; executable paths use the renamed binaries.
- Existing unrelated untracked build artifacts in the worktree were not modified or removed.
- Distribution snapshots include current contents of files known to the Git index, including tracked unstaged edits. New source files must be added to Git before `distcheck`; this correctly exposed the currently untracked `ParallelRemapTemplate.cpp` during verification. The successful end-to-end run used a temporary alternate Git index to include that file without changing the real index.
