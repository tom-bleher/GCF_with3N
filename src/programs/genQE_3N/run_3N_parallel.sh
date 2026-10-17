#!/usr/bin/env bash
#
# run_3N_parallel.sh — generate 3N+FSI events across all CPU cores, then merge.
#
# This is purely an orchestration wrapper: it launches several independent
# `genQE_3N_FSI` processes (each generating a slice of the requested number of
# successful events) and merges their ROOT files with `hadd`. The generator
# binary and its physics are NOT modified, so the merged sample is
# STATISTICALLY EQUIVALENT to a single serial run of the same size (events are
# independent; see plan). It is NOT bit-for-bit identical to a serial run —
# each worker auto-seeds its own RNG (TRandom3(0)).
#
# Usage:
#   ./run_3N_parallel.sh [TOTAL] [OUTPUT] [EBEAM] [U] [-- <extra generator flags>]
#
#   TOTAL    total number of successful (weight>0) events   (default 5000000)
#   OUTPUT   merged output .root file                        (default events/3N_FSI_5M_12C.root)
#   EBEAM    beam energy in GeV                              (default 6.0)
#   U        interaction number (positional arg 2 of binary) (default 1)
#   extra    anything after `--` is passed verbatim to every worker, e.g.
#            -- -A 12 -Z 6 -f hN -K 0.220 -s 0.15
#            (with none, defaults reproduce the existing 1.5M-file settings)
#
# Environment overrides:
#   JOBS=N    number of parallel workers   (default: all cores)
#   FORCE=1   overwrite a non-empty OUTPUT  (same as passing -f)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ----------------------------------------------------------------------------
# Parse args (TOTAL OUTPUT EBEAM U  [-f]  [-- extra...])
# ----------------------------------------------------------------------------
TOTAL=5000000
OUTPUT="events/3N_FSI_5M_12C.root"
EBEAM="6.0"
U="1"
EXTRA=()
positional=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --) shift; EXTRA=("$@"); break ;;
        -f|--force) FORCE=1; shift ;;
        -h|--help)
            sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) positional+=("$1"); shift ;;
    esac
done
[[ ${#positional[@]} -ge 1 ]] && TOTAL="${positional[0]}"
[[ ${#positional[@]} -ge 2 ]] && OUTPUT="${positional[1]}"
[[ ${#positional[@]} -ge 3 ]] && EBEAM="${positional[2]}"
[[ ${#positional[@]} -ge 4 ]] && U="${positional[3]}"
FORCE="${FORCE:-0}"

# ----------------------------------------------------------------------------
# Pre-flight checks
# ----------------------------------------------------------------------------
BINARY="$SCRIPT_DIR/build/genQE_3N_FSI"
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: generator binary not found at $BINARY" >&2
    echo "       build it first:  cmake --build build" >&2
    exit 1
fi

HADD="$(command -v hadd || true)"
if [[ -z "$HADD" ]] && command -v root-config >/dev/null 2>&1; then
    HADD="$(root-config --bindir)/hadd"
fi
if [[ -z "$HADD" || ! -x "$HADD" ]]; then
    echo "ERROR: hadd not found (need ROOT in PATH)" >&2
    exit 1
fi

if ! [[ "$TOTAL" =~ ^[0-9]+$ ]] || [[ "$TOTAL" -le 0 ]]; then
    echo "ERROR: TOTAL must be a positive integer (got '$TOTAL')" >&2
    exit 1
fi

if [[ -s "$OUTPUT" && "$FORCE" != "1" ]]; then
    echo "ERROR: $OUTPUT already exists and is non-empty." >&2
    echo "       Pass -f (or FORCE=1) to overwrite." >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# GENIE environment (mirror ../genQE_FSI/run_SRC.sh, pointed from genQE_3N).
# Only set vars whose targets exist; otherwise leave the shell's profile env.
# ----------------------------------------------------------------------------
GENIE_DIR="$SCRIPT_DIR/../genQE_FSI/Generator-R-3_06_02"
if [[ -d "$GENIE_DIR" ]]; then
    export GENIE="$GENIE_DIR"
    OVERRIDE="$SCRIPT_DIR/../genQE_FSI/config/genie_override"
    if [[ -d "$OVERRIDE" ]]; then
        export GXMLPATH="$OVERRIDE:$GENIE_DIR/config/G18_10a:$GENIE_DIR/config"
    fi
    QUIET="$SCRIPT_DIR/../genQE_FSI/config/quiet_messenger.xml"
    [[ -f "$QUIET" ]] && export GMSGCONF="$QUIET"
    export LD_LIBRARY_PATH="$GENIE_DIR/lib:$HOME/local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
else
    echo "WARNING: $GENIE_DIR not found; relying on existing GENIE env from profile." >&2
fi

# ----------------------------------------------------------------------------
# Job plan: split TOTAL into JOBS chunks, remainder spread over first chunks.
# ----------------------------------------------------------------------------
JOBS="${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
[[ "$JOBS" -gt "$TOTAL" ]] && JOBS="$TOTAL"   # never more workers than events

base=$(( TOTAL / JOBS ))
rem=$(( TOTAL % JOBS ))

OUT_DIR="$(cd "$(dirname "$OUTPUT")" && pwd)"
OUT_BASE="$(basename "${OUTPUT%.root}")"
PARTS_DIR="$OUT_DIR/parts_${OUT_BASE}_$$"
mkdir -p "$PARTS_DIR"

echo "=================================================================="
echo " Parallel 3N generation"
echo "   total events : $TOTAL  (successful, weight>0)"
echo "   workers      : $JOBS  (each ~$base events; first $rem get +1)"
echo "   output       : $OUTPUT"
echo "   beam / u     : $EBEAM / $U"
echo "   extra flags  : ${EXTRA[*]:-(none -> defaults A=12 Z=6 hN kMin=0.25 sigCM=0.15 CM on FSI on)}"
echo "   parts dir    : $PARTS_DIR"
echo "   seeds        : auto (TRandom3(0)); launches staggered 1s to avoid collisions"
echo "=================================================================="

# ----------------------------------------------------------------------------
# Precompute the chunk plan (part-file + size per worker index).
# ----------------------------------------------------------------------------
declare -a PIDS PARTS SIZES
for (( i=0; i<JOBS; i++ )); do
    n=$base
    (( i < rem )) && n=$(( n + 1 ))
    PARTS+=("$(printf "%s/%s.part%02d.root" "$PARTS_DIR" "$OUT_BASE" "$i")")
    SIZES+=("$n")
done

# launch_chunk <idx> [attempt-label] — (re)start one worker on its chunk.
# Auto-seeded each time, so a retry draws a fresh RNG stream and dodges the
# rare per-event GENIE FSI segfault that may have killed the prior attempt.
launch_chunk() {
    local idx="$1" label="${2:-launched}"
    local part="${PARTS[$idx]}" log="${PARTS[$idx]%.root}.log" n="${SIZES[$idx]}"
    "$BINARY" "$EBEAM" "$U" "$part" "$n" -v ${EXTRA[@]+"${EXTRA[@]}"} >"$log" 2>&1 &
    PIDS[$idx]=$!
    echo "  $label worker $idx  pid=${PIDS[$idx]}  events=$n  -> $(basename "$part")"
    sleep 1   # stagger => distinct TRandom3(0) auto-seeds
}

# wait_and_collect <idx...> — wait on the listed workers; echo result; set
# the global `failed` array to the indices that exited non-zero.
wait_and_collect() {
    failed=()
    local idx rc
    for idx in "$@"; do
        if wait "${PIDS[$idx]}"; then
            echo "  worker $idx done ok"
        else
            rc=$?
            echo "  ERROR: worker $idx (pid ${PIDS[$idx]}) exited $rc — see ${PARTS[$idx]%.root}.log" >&2
            failed+=("$idx")
        fi
    done
}

# ----------------------------------------------------------------------------
# Launch all workers, then retry any that crash (up to RETRIES extra rounds).
# ----------------------------------------------------------------------------
MAX_RETRIES="${RETRIES:-2}"
declare -a failed

for (( i=0; i<JOBS; i++ )); do launch_chunk "$i"; done
echo "Waiting for $JOBS workers..."
wait_and_collect $(seq 0 $((JOBS-1)))

attempt=0
while [[ ${#failed[@]} -gt 0 && $attempt -lt $MAX_RETRIES ]]; do
    attempt=$(( attempt + 1 ))
    retry=("${failed[@]}")
    echo "Retry round $attempt/$MAX_RETRIES for ${#retry[@]} crashed worker(s): ${retry[*]}"
    for idx in "${retry[@]}"; do launch_chunk "$idx" "retry-launched"; done
    wait_and_collect "${retry[@]}"
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "ERROR: worker(s) ${failed[*]} still failing after $MAX_RETRIES retries;" \
         "NOT merging. Part files kept in $PARTS_DIR" >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# Merge and verify.
# ----------------------------------------------------------------------------
echo "Merging ${#PARTS[@]} files with hadd -> $OUTPUT"
"$HADD" -f "$OUTPUT" "${PARTS[@]}"

echo "------------------------------------------------------------------"
echo "Per-worker efficiency (filled / attempts):"
grep -h "^Done:" "$PARTS_DIR"/*.log 2>/dev/null || echo "  (no Done lines found)"
echo "------------------------------------------------------------------"

# Entry-count check (use project venv if present, else ROOT, else skip).
PYBIN=""
[[ -x "$SCRIPT_DIR/.venv/bin/python" ]] && PYBIN="$SCRIPT_DIR/.venv/bin/python"
if [[ -n "$PYBIN" ]]; then
    "$PYBIN" - "$OUTPUT" "$TOTAL" <<'PY'
import sys, uproot
path, total = sys.argv[1], int(sys.argv[2])
n = uproot.open(path)["genT"].num_entries
print(f"merged genT entries: {n}  (expected {total})")
sys.exit(0 if n == total else 3)
PY
    if [[ $? -ne 0 ]]; then
        echo "WARNING: merged entry count does not match expected $TOTAL" >&2
    fi
else
    echo "(skipping entry-count check: no .venv python found)"
fi

# Clean up part files only on full success.
rm -rf "$PARTS_DIR"
echo "Done. Output: $OUTPUT"
