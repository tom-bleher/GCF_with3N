#!/usr/bin/env bash
#
# run_convention_v2.sh — {absorbed events kept, thrown away} x {unique, all}
# on the existing 60M no-tail samples (30M SRC / 60M MF analysed).
#
# Common to all four: detector convention, CLAS fiducials + map at phi+180 deg,
# (e,e'p) with no recoil requirement.
#   absorption: thrown away = default (--absorption-veto, weight 0)
#               kept        = --no-absorption-veto
#   (e,e'pp) recoil:
#     unique  exactly one other accepted proton >= 0.35 GeV/c, else not (e,e'pp)
#             (the event stays (e,e'p)); softer extras allowed
#     all     every other accepted proton >= 0.35 GeV/c gives its own (e,e'pp)
#             entry with that proton as the recoil; (e,e'p) counted once
#
# Also: detector exactly2 (absorbed thrown away; (e,e'pp) only with exactly two
# accepted protons in total) and --convention wright (now with her phi rotation).
#
# Usage: run_convention_v2.sh [pid-to-wait-for]
# A cell whose run.log already ends in "Done." is skipped, so the script can be
# re-run safely. Cells run sequentially: each holds ~100 GB.
#
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
PY=/home/ilayws/Ilay-Generators/src/programs/genQE_3N/.venv/bin/python
EV=events/natalie
CAP_SRC=30000000
CAP_MF=60000000
WAIT_PID="${1:-}"
stamp() { date '+%F %T'; }

run() {   # name, extra args...
    local name=$1; shift
    local out="analysis/plots/natalie_paper/v2_$name"; mkdir -p "$out"
    if [[ -f "$out/run.log" ]] && tail -1 "$out/run.log" | grep -q "^Done\.$"; then
        echo "[$(stamp)] --- $name: already done, skipping"; return
    fi
    echo "[$(stamp)] --- $name -> $out   ($*)"
    "$PY" analysis/natalie_paper_plots.py \
        --pwia "$EV/pwia_SRC.root" --fsi "$EV/fsi_hA_SRC_rho1_notail.root" \
        --fsi-mf "$EV/fsi_hA_MF_localFG_notail.root" \
        --max-events $CAP_SRC --mf-max-events $CAP_MF --seed 12345 \
        --outdir "$out" "$@" > "$out/run.log" 2>&1
    echo "[$(stamp)]     rc=$?"
}

if [[ -n "$WAIT_PID" ]]; then
    echo "[$(stamp)] waiting for pid $WAIT_PID to finish"
    while kill -0 "$WAIT_PID" 2>/dev/null; do sleep 60; done
fi

echo "[$(stamp)] ======== {abs kept, abs thrown away} x {unique, all} ========"
free -g | sed -n '2p'
run detector_unique          --convention detector --epp-recoil unique
run detector_all             --convention detector --epp-recoil all
run detector_unique_keepabs  --convention detector --epp-recoil unique --no-absorption-veto
run detector_all_keepabs     --convention detector --epp-recoil all    --no-absorption-veto
run detector_exactly2       --convention detector --epp-recoil exactly2
run wright                  --convention wright

echo
echo "[$(stamp)] ===== SUMMARY (Full effective transparencies) ====="
printf "%-26s %-16s %10s %10s\n" "recoil rule" "absorbed events" "(e,e'p)" "(e,e'pp)"
for spec in "detector_unique:unique:thrown away" "detector_all:all:thrown away" \
            "detector_unique_keepabs:unique:kept" "detector_all_keepabs:all:kept" \
            "detector_exactly2:exactly2:thrown away" "wright:wright (hers):thrown away"; do
    IFS=: read -r name rule abs <<< "$spec"
    d="analysis/plots/natalie_paper/v2_$name/run.log"
    ep=$(grep  "(e,e'p) " "$d" 2>/dev/null | grep paper | tail -1 | awk '{print $3}')
    epp=$(grep "(e,e'pp)" "$d" 2>/dev/null | grep paper | tail -1 | awk '{print $3}')
    printf "%-26s %-16s %10s %10s\n" "$rule" "$abs" "${ep:--}" "${epp:--}"
done
printf "%-26s %-16s %10s %10s\n" "PAPER (Full)" "" "0.58" "0.73"
for name in detector_unique detector_all detector_unique_keepabs detector_all_keepabs detector_exactly2; do
    echo; echo "## $name"
    grep -E "\[detector\]|channels|possible recoils" "analysis/plots/natalie_paper/v2_$name/run.log"
done
echo "[$(stamp)] ======== DONE ========"
