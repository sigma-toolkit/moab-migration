#!/bin/bash
# ==============================================================================
# generate_baselines.sh
# Generate golden reference files for mbtempest regression tests
#
# Usage:
#   cd <build_dir>
#   bash ../test/mbtempest_regression/generate_baselines.sh [build_dir]
#
# Or via CMake:
#   cmake -DMOAB_BUILD_MBTEMPEST_TESTS=ON -DMOAB_MBTEMPEST_GENERATE_BASELINES=ON ..
#   make -j$(nproc)
#   ctest -L mbtempest
#
# After generation, upload baselines from ${BASELINE_DIR} to the external
# server and set MBTEMPEST_BASELINE_URL in CMake.
# ==============================================================================

set -euo pipefail

# Paths — adjust for your build
BUILD_DIR="${1:-$(pwd)}"
MBTEMPEST="${BUILD_DIR}/bin/mbtempest"
MBCONVERT="${BUILD_DIR}/bin/mbconvert"
WORK_DIR="${BUILD_DIR}/test/mbtempest_regression/work"
BASELINE_DIR="${BUILD_DIR}/test/mbtempest_regression/baselines"
EXAMPLE_DIR="$(cd "${BUILD_DIR}/.." && pwd)/examples/advanced"

if [ ! -x "${MBTEMPEST}" ]; then
  echo "ERROR: mbtempest not found at ${MBTEMPEST}"
  echo "Usage: $0 [build_dir]"
  exit 1
fi

# MPI launcher
MPIRUN="${MPIRUN:-mpirun}"
if ! command -v "${MPIRUN}" &>/dev/null; then
  MPIRUN="mpiexec"
fi

mkdir -p "${WORK_DIR}" "${BASELINE_DIR}"
cd "${WORK_DIR}"

# Symlink pre-existing meshes
for f in ne30.h5m rll90.h5m ico30.h5m icod30.h5m; do
  [ -L "$f" ] || [ -f "$f" ] || ln -s "${EXAMPLE_DIR}/$f" "$f" 2>/dev/null || true
done

echo "============================================"
echo "Generating mbtempest regression baselines"
echo "Working directory: ${WORK_DIR}"
echo "Baseline directory: ${BASELINE_DIR}"
echo "============================================"

PASS=0
FAIL=0

# Helper: run mbtempest and copy output map as baseline
run_and_save() {
  local name="$1"
  local output="$2"
  shift 2
  echo ""
  echo "--- ${name} ---"
  if "${MPIRUN}" -np 1 "${MBTEMPEST}" "$@" > /dev/null 2>&1; then
    if [ -f "${output}" ]; then
      local gold_name="gold_${output#map_}"
      cp "${output}" "${BASELINE_DIR}/${gold_name}"
      echo "  PASS -> ${gold_name}"
      PASS=$((PASS + 1))
    else
      echo "  WARN: output ${output} not found"
      FAIL=$((FAIL + 1))
    fi
  else
    echo "  FAIL: mbtempest returned non-zero"
    FAIL=$((FAIL + 1))
  fi
}

# ==============================================================================
# Step 1: Generate SE meshes with DOF tags
# ==============================================================================
echo ""
echo "=== Generating SE meshes ==="

"${MPIRUN}" -np 1 "${MBTEMPEST}" -t 0 -r 30 -f ne30_se.g > /dev/null 2>&1
"${MPIRUN}" -np 1 "${MBTEMPEST}" -t 0 -r 25 -f ne25_se.g > /dev/null 2>&1
"${MBCONVERT}" -B ne30_se.g ne30_se_o4.h5m -i GLOBAL_DOFS -r 4 > /dev/null 2>&1
"${MBCONVERT}" -B ne30_se.g ne30_se_o2.h5m -i GLOBAL_DOFS -r 2 > /dev/null 2>&1
"${MBCONVERT}" -B ne25_se.g ne25_se_o4.h5m -i GLOBAL_DOFS -r 4 > /dev/null 2>&1
echo "  SE meshes ready"

# ==============================================================================
# Step 2: FV-FV remapping baselines
# ==============================================================================
echo ""
echo "=== FV-FV baselines ==="

run_and_save "ne30→rll90 FV-FV o1" map_ne30_rll90_fv_fv_o1.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_ne30_rll90_fv_fv_o1.nc

run_and_save "ne30→ico30 FV-FV o1" map_ne30_ico30_fv_fv_o1.nc \
  -t 5 -l ne30.h5m -l ico30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_ne30_ico30_fv_fv_o1.nc

run_and_save "ne30→icod30 FV-FV o1" map_ne30_icod30_fv_fv_o1.nc \
  -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_ne30_icod30_fv_fv_o1.nc

run_and_save "rll90→ne30 FV-FV o1" map_rll90_ne30_fv_fv_o1.nc \
  -t 5 -l rll90.h5m -l ne30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_rll90_ne30_fv_fv_o1.nc

run_and_save "ico30→ne30 FV-FV o1" map_ico30_ne30_fv_fv_o1.nc \
  -t 5 -l ico30.h5m -l ne30.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_ico30_ne30_fv_fv_o1.nc

# Higher order FV
run_and_save "ne30→icod30 FV-FV o2" map_ne30_icod30_fv_fv_o2.nc \
  -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 2 -o 2 --verify --var SH \
  -f map_ne30_icod30_fv_fv_o2.nc

run_and_save "ne30→icod30 FV-FV o3" map_ne30_icod30_fv_fv_o3.nc \
  -t 5 -l ne30.h5m -l icod30.h5m -w -m fv -m fv -o 3 -o 3 --verify --var SH \
  -f map_ne30_icod30_fv_fv_o3.nc

# ==============================================================================
# Step 3: CGLL-FV baselines
# ==============================================================================
echo ""
echo "=== CGLL-FV baselines ==="

run_and_save "CGLL→FV o2-o1" map_cgll_fv_o2_o1.nc \
  -t 5 -l ne30_se_o2.h5m -l rll90.h5m -w -m cgll -m fv -o 2 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV \
  -f map_cgll_fv_o2_o1.nc

run_and_save "CGLL→FV o4-o1" map_cgll_fv_o4_o1.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV \
  -f map_cgll_fv_o4_o1.nc

run_and_save "CGLL→FV o4-o1 mono1" map_cgll_fv_o4_o1_mono1.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 1 --verify --var SV \
  -f map_cgll_fv_o4_o1_mono1.nc

run_and_save "CGLL→FV o4-o1 mono2" map_cgll_fv_o4_o1_mono2.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 2 --verify --var SV \
  -f map_cgll_fv_o4_o1_mono2.nc

run_and_save "CGLL→FV o4-o1 mono3" map_cgll_fv_o4_o1_mono3.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --monotonicity 3 --verify --var SV \
  -f map_cgll_fv_o4_o1_mono3.nc

run_and_save "CGLL→FV o4-o1 nobubble" map_cgll_fv_o4_o1_nobubble.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --nobubble --verify --var SV \
  -f map_cgll_fv_o4_o1_nobubble.nc

# ==============================================================================
# Step 4: FV-CGLL and CGLL-CGLL baselines
# ==============================================================================
echo ""
echo "=== FV-CGLL and CGLL-CGLL baselines ==="

run_and_save "FV→CGLL o1-o4" map_fv_cgll_o1_o4.nc \
  -t 5 -l rll90.h5m -l ne30_se_o4.h5m -w -m fv -m cgll -o 1 -o 4 \
  -g GLOBAL_ID -g GLOBAL_DOFS --verify --var SV \
  -f map_fv_cgll_o1_o4.nc

run_and_save "CGLL→CGLL o4-o4" map_cgll_cgll_o4_o4.nc \
  -t 5 -l ne30_se_o4.h5m -l ne25_se_o4.h5m -w -m cgll -m cgll -o 4 -o 4 \
  -g GLOBAL_DOFS -g GLOBAL_DOFS --verify --var SV \
  -f map_cgll_cgll_o4_o4.nc

# ==============================================================================
# Step 4b: DGLL baselines (uses DGLOBAL_DOFS tag)
# ==============================================================================
echo ""
echo "=== DGLL baselines ==="

run_and_save "DGLL→FV o4-o1" map_dgll_fv_o4_o1.nc \
  -t 5 -l ne30_se_o4.h5m -l rll90.h5m -w -m dgll -m fv -o 4 -o 1 \
  -g DGLOBAL_DOFS -g GLOBAL_ID --verify --var SV \
  -f map_dgll_fv_o4_o1.nc

run_and_save "FV→DGLL o1-o4" map_fv_dgll_o1_o4.nc \
  -t 5 -l rll90.h5m -l ne30_se_o4.h5m -w -m fv -m dgll -o 1 -o 4 \
  -g GLOBAL_ID -g DGLOBAL_DOFS --verify --var SV \
  -f map_fv_dgll_o1_o4.nc

run_and_save "DGLL→DGLL o4-o4" map_dgll_dgll_o4_o4.nc \
  -t 5 -l ne30_se_o4.h5m -l ne25_se_o4.h5m -w -m dgll -m dgll -o 4 -o 4 \
  -g DGLOBAL_DOFS -g DGLOBAL_DOFS --verify --var SV \
  -f map_dgll_dgll_o4_o4.nc

run_and_save "DGLL→CGLL o4-o4" map_dgll_cgll_o4_o4.nc \
  -t 5 -l ne30_se_o4.h5m -l ne25_se_o4.h5m -w -m dgll -m cgll -o 4 -o 4 \
  -g DGLOBAL_DOFS -g GLOBAL_DOFS --verify --var SV \
  -f map_dgll_cgll_o4_o4.nc

run_and_save "CGLL→DGLL o4-o4" map_cgll_dgll_o4_o4.nc \
  -t 5 -l ne30_se_o4.h5m -l ne25_se_o4.h5m -w -m cgll -m dgll -o 4 -o 4 \
  -g GLOBAL_DOFS -g DGLOBAL_DOFS --verify --var SV \
  -f map_cgll_dgll_o4_o4.nc

# ==============================================================================
# Step 5: FV sub-method baselines
# ==============================================================================
echo ""
echo "=== FV sub-method baselines ==="

for method in bilin delaunay intbilin invdist; do
  run_and_save "ne30→rll90 ${method}" "map_ne30_rll90_fv_fv_${method}.nc" \
    -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod "${method}" \
    --verify --var SH -f "map_ne30_rll90_fv_fv_${method}.nc"
done

# ==============================================================================
# Step 6: Monotonicity (FV-FV) and CAAS baselines
# ==============================================================================
echo ""
echo "=== Monotonicity and CAAS baselines ==="

run_and_save "ne30→rll90 FV-FV mono1" map_ne30_rll90_fv_fv_mono1.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --monotonicity 1 \
  --verify --var SH -f map_ne30_rll90_fv_fv_mono1.nc

for caas in 1 2 3 4; do
  run_and_save "ne30→rll90 CAAS${caas}" "map_ne30_rll90_fv_fv_caas${caas}.nc" \
    -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --limiter "${caas}" \
    --verify --var SH -f "map_ne30_rll90_fv_fv_caas${caas}.nc"
done

# ==============================================================================
# Step 7: Advanced option baselines
# ==============================================================================
echo ""
echo "=== Advanced option baselines ==="

run_and_save "ne30→rll90 gnomonic" map_ne30_rll90_gnomonic.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --gnomonic \
  --verify --var SH -f map_ne30_rll90_gnomonic.nc

run_and_save "ne30→rll90 advfront" map_ne30_rll90_advfront.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --advfront \
  --verify --var SH -f map_ne30_rll90_advfront.nc

run_and_save "ne30→rll90 noconserve" map_ne30_rll90_noconserve.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --noconserve \
  --verify --var SH -f map_ne30_rll90_noconserve.nc

run_and_save "ne30→rll90 convexity" map_ne30_rll90_convexity.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --enforce_convexity \
  --verify --var SH -f map_ne30_rll90_convexity.nc

run_and_save "ne30→rll90 sparse" map_ne30_rll90_sparse.nc \
  -t 5 -l ne30.h5m -l rll90.h5m -w -m fv -m fv -o 1 -o 1 --sparseconstraints \
  --verify --var SH -f map_ne30_rll90_sparse.nc

run_and_save "FV→CGLL volumetric" map_fv_cgll_volumetric.nc \
  -t 5 -l rll90.h5m -l ne30_se_o4.h5m -w -m fv -m cgll -o 1 -o 4 \
  -g GLOBAL_ID -g GLOBAL_DOFS --volumetric --verify --var SV \
  -f map_fv_cgll_volumetric.nc

# ==============================================================================
# Step 8: Convergence baselines
# ==============================================================================
echo ""
echo "=== Convergence baselines ==="

run_and_save "ne30→icod30 CGLL-FV o4-o1 (convergence)" map_ne30_icod30_cgll_fv_o4_o1.nc \
  -t 5 -l ne30_se_o4.h5m -l icod30.h5m -w -m cgll -m fv -o 4 -o 1 \
  -g GLOBAL_DOFS -g GLOBAL_ID --verify --var SV \
  -f map_ne30_icod30_cgll_fv_o4_o1.nc

# ==============================================================================
# Done
# ==============================================================================

echo ""
echo "============================================"
echo "Baseline generation complete."
echo "  Passed: ${PASS}"
echo "  Failed: ${FAIL}"
echo ""
ls "${BASELINE_DIR}"/gold_*.nc 2>/dev/null | wc -l | xargs echo "  Files:"
echo ""
echo "Baselines saved to: ${BASELINE_DIR}"
echo "============================================"

[ "${FAIL}" -eq 0 ] || exit 1
