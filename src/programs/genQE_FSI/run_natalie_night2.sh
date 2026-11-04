#!/usr/bin/env bash
#
# run_natalie_night2.sh — analysis of the 60M/60M/150M samples, plus a
# matched no-tail-medium set, so both cascade media can be compared.
#
# WHY: the samples generated tonight use the LocalFGM/Correlated cascade
# medium (config/genie_override_srctail), which is what her
# ModelConfiguration.xml:37 selects. On those samples the (e,e'pp) effective
# transparency is 1.26. On the earlier 2M samples, which used the no-tail
# LocalFGM/Default medium, it is 0.74 -- matching her published 0.73 (Wright
# et al., doi:10.48550/arXiv.2104.05090). Her config and her published number
# disagree, so we need both at good statistics to decide.
#
# MEMORY: measured peak RSS is 39.1 GB for 45M events analysed (10M/10M/25M),
# i.e. ~0.85 GB per million. The analysis holds all three samples at once, so
# 270M would need ~220 GB against 178 GB available. Everything below is
# capped at 30M/30M/60M -> ~105 GB peak, a measured number rather than an
# estimate. The generated samples keep their full statistics on disk.
#
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PY=/home/ilayws/Ilay-Generators/src/programs/genQE_3N/.venv/bin/python
EV=events/natalie
CAP_SRC=30000000
CAP_MF=60000000
stamp() { date '+%F %T'; }

echo "[$(stamp)] ============ NIGHT 2: analysis + no-tail comparison ============"
echo "[$(stamp)] analysis caps: ${CAP_SRC} SRC / ${CAP_MF} MF  (peak RSS ~105 GB)"
free -g | sed -n '2p'

# ---------------------------------------------------------------- stage A
# Analyse the samples we already have (SRC-tail cascade medium).
echo
echo "[$(stamp)] ===== STAGE A: srctail medium, both conventions ====="
for conv in wright detector; do
    out="analysis/plots/natalie_paper/final_srctail_$conv"; mkdir -p "$out"
    echo "[$(stamp)] --- $conv -> $out"
    "$PY" analysis/natalie_paper_plots.py \
        --pwia "$EV/pwia_SRC.root" --fsi "$EV/fsi_hA_SRC_rho1.root" \
        --fsi-mf "$EV/fsi_hA_MF_localFG.root" --convention "$conv" \
        --max-events $CAP_SRC --mf-max-events $CAP_MF \
        --outdir "$out" --seed 12345 > "$out/run.log" 2>&1
    echo "[$(stamp)]     rc=$?"; grep -A3 "Channel  " "$out/run.log" | head -4
done

# ---------------------------------------------------------------- stage B
# Generate the matched no-tail-medium FSI and MF samples. PWIA needs no
# counterpart: with no cascade the medium is irrelevant, so pwia_SRC.root is
# reused for both.
echo
echo "[$(stamp)] ===== STAGE B: generate no-tail medium samples ====="
export GENIE_OVERRIDE_DIR="$SCRIPT_DIR/config/genie_override"
JOBS=150 FORCE=1 ./run_2N_parallel.sh 60000000 \
    "$EV/fsi_hA_SRC_rho1_notail.root" 1 hA 0.150 2 5.01 global 0.25 0 1 local 0
JOBS=150 FORCE=1 ./run_2N_parallel.sh 60000000 \
    "$EV/fsi_hA_MF_localFG_notail.root" 1 hA 0.150 4 5.01 local 0.25 0 1

# ---------------------------------------------------------------- stage C
echo
echo "[$(stamp)] ===== STAGE C: no-tail medium, both conventions ====="
for conv in wright detector; do
    out="analysis/plots/natalie_paper/final_notail_$conv"; mkdir -p "$out"
    echo "[$(stamp)] --- $conv -> $out"
    "$PY" analysis/natalie_paper_plots.py \
        --pwia "$EV/pwia_SRC.root" --fsi "$EV/fsi_hA_SRC_rho1_notail.root" \
        --fsi-mf "$EV/fsi_hA_MF_localFG_notail.root" --convention "$conv" \
        --max-events $CAP_SRC --mf-max-events $CAP_MF \
        --outdir "$out" --seed 12345 > "$out/run.log" 2>&1
    echo "[$(stamp)]     rc=$?"; grep -A3 "Channel  " "$out/run.log" | head -4
done

# ---------------------------------------------------------------- summary
echo
echo "[$(stamp)] ===== SUMMARY ====="
printf "%-34s %10s %10s\n" "configuration" "(e,e'p)" "(e,e'pp)"
for d in analysis/plots/natalie_paper/final_*/run.log; do
    [[ -f "$d" ]] || continue
    name=$(basename "$(dirname "$d")")
    ep=$(grep  "(e,e'p)"  "$d" | tail -1 | awk '{print $3}')
    epp=$(grep "(e,e'pp)" "$d" | tail -1 | awk '{print $3}')
    printf "%-34s %10s %10s\n" "$name" "${ep:--}" "${epp:--}"
done
printf "%-34s %10s %10s\n" "PAPER (Full)" "0.58" "0.73"
echo
echo "[$(stamp)] ============ ALL DONE ============"
df -h /home | tail -1
