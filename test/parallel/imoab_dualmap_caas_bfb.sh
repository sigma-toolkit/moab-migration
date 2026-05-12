#!/bin/bash
# imoab_dualmap_caas_bfb.sh
#
# End-to-end BFB regression test for the dual-map CAAS path.
#
# Step 1: serially compute lo/hi weight maps online and persist them via
#         iMOAB_WriteMapFile (--compute_online --write_maps); also write
#         the serial baseline digests for the lo/hi/dual projections.
#
# Step 2: re-run at N=1, 2, 4 (configurable) loading those same map files
#         from disk and writing per-cell digests. Loading from disk gives
#         a byte-identical input matrix on every rank; coupled with the
#         partition-invariant source field on the on-disk ATM mesh, the
#         per-cell projected values must be byte-identical across rank
#         counts.
#
# Step 3: cmp -s every {kernel}_<n>.txt against {kernel}_1.txt. Any byte
#         difference is a regression. Exit code 0 iff all pairs match.
#
# Usage:
#   imoab_dualmap_caas_bfb.sh <path-to-imoab_dualmap_caas-binary> [N1 N2 ...]
# Default rank counts: 1 2 4
#
# Required environment: mpirun on PATH, write access to a working dir.

set -u

EXE=${1:-}
shift || true
RANKS=("$@")
if [ ${#RANKS[@]} -eq 0 ]; then
    RANKS=(1 2 4)
fi

if [ -z "$EXE" ] || [ ! -x "$EXE" ]; then
    echo "Usage: $0 <path-to-imoab_dualmap_caas> [rank-count ...]" >&2
    exit 2
fi

WORKDIR=$(mktemp -d -t dualmap_caas_bfb_XXXXXX)
trap "rm -rf $WORKDIR" EXIT
cd "$WORKDIR"

echo "=== imoab_dualmap_caas BFB regression ==="
echo " EXE     : $EXE"
echo " WORKDIR : $WORKDIR"
echo " RANKS   : ${RANKS[*]}"

# Step 1: serial compute online + write maps + serial baseline digest
echo ""
echo "--- Step 1: serial compute online, write maps + baseline digests ---"
if ! mpirun -n 1 "$EXE" --compute_online --write_maps baseline --digest_prefix digest >step1.log 2>&1 ; then
    echo "FAIL: serial compute step exited non-zero. See $WORKDIR/step1.log" >&2
    cat step1.log >&2
    exit 1
fi
for f in baseline_lo.nc baseline_hi.nc digest_lo_1.txt digest_hi_1.txt digest_dual_1.txt ; do
    if [ ! -s "$f" ]; then
        echo "FAIL: expected serial output file '$f' missing or empty" >&2
        exit 1
    fi
done
echo "  baseline maps:    baseline_{lo,hi}.nc"
echo "  serial digests:   digest_{lo,hi,dual}_1.txt"

# Step 2: parallel reload + parallel digests
echo ""
echo "--- Step 2: reload baseline maps at each rank count, write digests ---"
for n in "${RANKS[@]}" ; do
    if [ "$n" = "1" ]; then continue ; fi   # serial run already produced n=1
    if ! mpirun -n "$n" "$EXE" -l baseline_lo.nc -h baseline_hi.nc \
                                --digest_prefix digest >step2_n${n}.log 2>&1 ; then
        echo "FAIL: parallel run at n=$n exited non-zero. See $WORKDIR/step2_n${n}.log" >&2
        cat step2_n${n}.log >&2
        exit 1
    fi
    echo "  n=$n digests:     digest_{lo,hi,dual}_${n}.txt"
done

# Step 3: strict byte-for-byte diff against serial baseline
echo ""
echo "--- Step 3: strict cmp vs serial (n=1) baseline ---"
fail=0
for kernel in lo hi dual ; do
    for n in "${RANKS[@]}" ; do
        if [ "$n" = "1" ]; then continue ; fi
        a="digest_${kernel}_1.txt"
        b="digest_${kernel}_${n}.txt"
        if cmp -s "$a" "$b" ; then
            echo "  PASS  $b == $a"
        else
            echo "  FAIL  $b != $a"
            ndiff=$(diff "$a" "$b" 2>/dev/null | grep -c '^<')
            echo "        ($ndiff/$(wc -l < $a) lines differ; sample:)"
            diff "$a" "$b" | head -5 | sed 's/^/          /'
            fail=$((fail+1))
        fi
    done
done

echo ""
if [ $fail -eq 0 ]; then
    echo "OK: dual-map CAAS is BFB across rank counts ${RANKS[*]} for file-loaded maps."
    exit 0
else
    echo "FAIL: $fail digest pair(s) differ. See diffs above."
    exit 1
fi
