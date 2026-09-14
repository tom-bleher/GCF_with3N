#!/usr/bin/env bash
#
# run_reqp1_test.sh — measure the effect of Wright's SECOND SRC branch test.
#
# Her generator decides "is this an SRC pair?" twice: first on the relative
# momentum k > kF(r), then AGAIN on |p1| = |k + 0.5*P_cm| > kF(r) after the
# c.m. motion has been added. Events failing the second test emit no recoil
# and fall into her mean-field class. Our `src_require_p1` switch reproduces
# that: a failing event gets weight 0 but still counts as a generation
# attempt, so the attempts-based luminosity normalisation stays exact.
#
# Everything else is pinned to the night-2 no-tail configuration so the only
# difference vs events/natalie/{pwia_SRC,fsi_hA_SRC_rho1_notail}.root is this
# one switch: local onset, kF 0.25, sigma_cm 0.150, 5.01 GeV, rho^1 vertex,
# hA, LocalFGM/Default (no SRC tail) cascade medium.
#
# 30M each, matching the 30M analysis cap already used for the night samples.
#
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PY=/home/ilayws/Ilay-Generators/src/programs/genQE_3N/.venv/bin/python
EV=events/natalie
N=30000000
export GENIE_OVERRIDE_DIR="$SCRIPT_DIR/config/genie_override"
stamp() { date '+%F %T'; }

echo "[$(stamp)] ===== src_require_p1 = 1 : generation ====="
#                                                            DOFSI MOD SIGCM WF EBEAM FG   KF  INDEP RHO ONSET REQP1
JOBS=150 FORCE=1 ./run_2N_parallel.sh $N "$EV/pwia_SRC_reqp1.root"      0 hA 0.150 2 5.01 global 0.25 0 1 local 1 || exit 1
JOBS=150 FORCE=1 ./run_2N_parallel.sh $N "$EV/fsi_hA_SRC_reqp1.root"    1 hA 0.150 2 5.01 global 0.25 0 1 local 1 || exit 1

echo "[$(stamp)] ===== analysis (both conventions) ====="
for conv in wright detector; do
    out="analysis/plots/natalie_paper/reqp1_$conv"; mkdir -p "$out"
    echo "[$(stamp)] --- $conv -> $out"
    "$PY" analysis/natalie_paper_plots.py \
        --pwia "$EV/pwia_SRC_reqp1.root" --fsi "$EV/fsi_hA_SRC_reqp1.root" \
        --fsi-mf "$EV/fsi_hA_MF_localFG_notail.root" --convention "$conv" \
        --max-events $N --mf-max-events 60000000 \
        --outdir "$out" --seed 12345 > "$out/run.log" 2>&1
    echo "[$(stamp)]     rc=$?"; grep -A3 "Channel  " "$out/run.log" | head -4
done

echo
echo "[$(stamp)] ===== SUMMARY (Full transparencies) ====="
printf "%-34s %10s %10s\n" "configuration" "(e,e'p)" "(e,e'pp)"
for d in analysis/plots/natalie_paper/{final_notail,reqp1}_*/run.log; do
    [[ -f "$d" ]] || continue
    name=$(basename "$(dirname "$d")")
    ep=$(grep  "(e,e'p)"  "$d" | tail -1 | awk '{print $3}')
    epp=$(grep "(e,e'pp)" "$d" | tail -1 | awk '{print $3}')
    printf "%-34s %10s %10s\n" "$name" "${ep:--}" "${epp:--}"
done
printf "%-34s %10s %10s\n" "PAPER (Full)" "0.58" "0.73"
echo "[$(stamp)] ============ DONE ============"
df -h /home | tail -1
