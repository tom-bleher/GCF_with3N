#!/usr/bin/env bash
#
# run_natalie_overnight.sh — full production + all figures for the
# Wright et al. (arXiv:2104.05090, doi:10.48550/arXiv.2104.05090) reproduction,
# in one detached run.
#
# Launch with:
#   setsid nohup bash run_natalie_overnight.sh > natalie_overnight.log 2>&1 < /dev/null &
# It then survives logout / closing the laptop.
#
# Stage 1  generate 20M PWIA, 20M GCF+FSI, 60M mean-field
#          (run_natalie_samples.sh pins the paper conventions: hA, 5.01 GeV,
#           sigma_cm 0.150, rho^1 vertex, local-FG MF, SRC onset "local",
#           and the LocalFGM/Correlated cascade medium that matches hers)
# Stage 2  run the analysis TWICE:
#            --convention wright    reproduces her genie_gcf_loop.C
#            --convention detector  our own detector-style convention
# Stage 3  print a comparison summary of the numbers that matter
#
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PY=/home/ilayws/Ilay-Generators/src/programs/genQE_3N/.venv/bin/python
EV=events/natalie
stamp() { date '+%F %T'; }

N_SRC="${N_SRC:-20000000}"
N_MF="${N_MF:-60000000}"
export JOBS="${JOBS:-150}"
export SRC_ONSET="${SRC_ONSET:-local}"

echo "[$(stamp)] ================= NATALIE OVERNIGHT RUN ================="
echo "[$(stamp)] SRC samples: $N_SRC   MF sample: $N_MF   cores: $JOBS"
echo "[$(stamp)] SRC onset: $SRC_ONSET"
df -h /home | tail -1

# ---------------------------------------------------------------- stage 1
echo
echo "[$(stamp)] ===== STAGE 1: generation ====="
./run_natalie_samples.sh "$N_SRC" "$N_MF"
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "[$(stamp)] ERROR: generation failed (rc=$rc). Stopping before figures."
    exit 1
fi

for f in "$EV/pwia_SRC.root" "$EV/fsi_hA_SRC_rho1.root" "$EV/fsi_hA_MF_localFG.root"; do
    if [[ ! -s "$f" ]]; then
        echo "[$(stamp)] ERROR: missing or empty $f. Stopping before figures."
        exit 1
    fi
done
echo "[$(stamp)] generation complete:"
ls -lh "$EV"/*.root

# ---------------------------------------------------------------- stage 2
echo
echo "[$(stamp)] ===== STAGE 2: figures, both conventions ====="
for conv in wright detector; do
    out="analysis/plots/natalie_paper/final_$conv"
    mkdir -p "$out"
    echo "[$(stamp)] --- convention: $conv -> $out"
    "$PY" analysis/natalie_paper_plots.py \
        --pwia   "$EV/pwia_SRC.root" \
        --fsi    "$EV/fsi_hA_SRC_rho1.root" \
        --fsi-mf "$EV/fsi_hA_MF_localFG.root" \
        --convention "$conv" \
        --outdir "$out" --seed 12345 \
        > "$out/run.log" 2>&1
    if [[ $? -ne 0 ]]; then
        echo "[$(stamp)] ERROR: analysis failed for convention=$conv; see $out/run.log"
        tail -25 "$out/run.log"
    else
        echo "[$(stamp)] done. Transparencies:"
        grep -A3 "Channel  " "$out/run.log" | head -4
    fi
done

# ---------------------------------------------------------------- stage 3
echo
echo "[$(stamp)] ===== STAGE 3: summary ====="
for conv in wright detector; do
    log="analysis/plots/natalie_paper/final_$conv/run.log"
    [[ -f "$log" ]] || continue
    echo
    echo "########## convention = $conv"
    grep -E "Convention|wright\]|Q2 pre-cut|accepted nucleons|channels:" "$log" | head -12
    grep -A3 "Channel  " "$log" | head -4
    echo "  figures:"; ls "analysis/plots/natalie_paper/final_$conv"/*.png 2>/dev/null | xargs -n1 basename | sed 's/^/    /'
done

echo
echo "[$(stamp)] paper reference values:  (e,e'p) T+SCX 0.61 / Full 0.58"
echo "[$(stamp)]                          (e,e'pp) T+SCX 0.83 / Full 0.73"
echo "[$(stamp)] NOTE: her published (e,e'pp) T+SCX 0.83 is hard-coded in"
echo "[$(stamp)]       gst_plotter.C:1356 as (.44/.53), not simulated."
echo
df -h /home | tail -1
echo "[$(stamp)] ================= ALL DONE ================="
