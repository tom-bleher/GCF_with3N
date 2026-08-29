#!/usr/bin/env bash
# Detached paper-production pipeline: waits for the running 10M validation
# jobs, then generates the final 100M samples (sequentially, 150 cores each)
# and converts the 3N genT to the detection-based events tree.
# Launch with:  setsid nohup bash run_paper_production.sh > production_pipeline.log 2>&1 &
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEN3N="$ROOT_DIR/src/programs/genQE_3N"
GEN2N="$ROOT_DIR/src/programs/genQE_FSI"

# GENIE env (needed for the direct SRC_analysis_3N call; the run_*_parallel
# scripts set their own copies).
export GENIE="$GEN2N/Generator-R-3_06_02"
export GXMLPATH="$GEN2N/config/genie_override:$GENIE/config/G18_10a:$GENIE/config"
export GMSGCONF="$GEN2N/config/quiet_messenger.xml"
export LD_LIBRARY_PATH="$GENIE/lib:$HOME/local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

stamp() { date '+%F %T'; }

echo "[$(stamp)] pipeline start; waiting for val10M_* runs to finish..."
while pgrep -f "val10M_" > /dev/null; do sleep 60; done
echo "[$(stamp)] validation runs finished."

echo "[$(stamp)] === 3N+FSI(indep) 100M, hN, 12C, kF=0.22, sigCM=0.15, Q2>1 ==="
cd "$GEN3N"
JOBS=150 FORCE=1 ./run_3N_parallel.sh 100000000 \
    events/3N_FSIindep_hN_100M_12C_kF220_6GeV_Q2gt1.root 6.0 1 -- \
    -A 12 -Z 6 -f hN -K 0.22 -s 0.15 -I -q 1:1000000 \
    || { echo "[$(stamp)] ERROR: 3N production failed"; exit 1; }

echo "[$(stamp)] === convert 3N genT -> events (detection-based) ==="
.venv/bin/python plotting/make_events_from_genT.py \
    --input events/3N_FSIindep_hN_100M_12C_kF220_6GeV_Q2gt1.root \
    --output events/events_3N_100M.root --ebeam 6.0 --kF 0.22 \
    || { echo "[$(stamp)] ERROR: 3N conversion failed"; exit 1; }

echo "[$(stamp)] === 3N PWIA no-CM 20M (Fig 1/2/4) ==="
./build/SRC_analysis_3N 20000000 0 0 events/3N_PWIA_noCM_20M_events.root \
    || { echo "[$(stamp)] ERROR: PWIA no-CM run failed"; exit 1; }

echo "[$(stamp)] === 2N+FSI(indep) 100M, hN, 6 GeV, kF=0.22 ==="
cd "$GEN2N"
JOBS=150 FORCE=1 ./run_2N_parallel.sh 100000000 \
    events/misc/2N_hN_100M_6GeV_kF220_indep.root 1 hN 0.150 0 6.0 global 0.220 1 2 \
    || { echo "[$(stamp)] ERROR: 2N production failed"; exit 1; }

echo "[$(stamp)] === MF local-FG 100M, hN, 6 GeV, kF=0.22 (Fig 7) ==="
JOBS=150 FORCE=1 ./run_2N_parallel.sh 100000000 \
    events/misc/MF_localFG_hN_100M_6GeV_kF220.root 1 hN 0.150 4 6.0 local 0.220 0 1 \
    || { echo "[$(stamp)] ERROR: MF production failed"; exit 1; }

echo "[$(stamp)] ALL PRODUCTION DONE"
