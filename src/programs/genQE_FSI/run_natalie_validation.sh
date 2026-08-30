#!/usr/bin/env bash
#
# run_natalie_validation.sh — small (2M) samples for the staged validation of
# the Wright et al. recreation. Each stage changes exactly ONE thing so its
# effect on the effective transparencies and on Figs 5/8/10 is attributable:
#
#   fixed : k_rel > kF                (current behaviour)
#   local : k_rel > kF(r)             (Wright's onset)
#   MF    : local Fermi gas, unchanged by either switch
#
# PWIA and FSI share identical primary configuration within each stage so the
# attempts-based luminosity matching in natalie_paper_plots.py is exact.
#
# Usage:  ./run_natalie_validation.sh [TOTAL_PER_SAMPLE]
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TOTAL="${1:-2000000}"
OUTDIR="events/natalie_val"
mkdir -p "$OUTDIR"

SIGCM=0.150; EBEAM=5.01; KF=0.25; RHOPOW=1   # Wright et al. conventions
JOBS="${JOBS:-120}"
stamp() { date '+%F %T'; }

#      TOTAL  OUTPUT  DOFSI MODEL SIGCM WF EBEAM FG KF INDEP RHOPOW ONSET REQP1
run() {
  echo "[$(stamp)] === $2"
  JOBS="$JOBS" FORCE=1 ./run_2N_parallel.sh "$@" \
      || { echo "[$(stamp)] ERROR on $2"; return 1; }
}

for onset in fixed local; do
  run "$TOTAL" "$OUTDIR/pwia_$onset.root" 0 hA "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW" "$onset" 0
  run "$TOTAL" "$OUTDIR/fsi_$onset.root"  1 hA "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW" "$onset" 0
done
run "$TOTAL" "$OUTDIR/mf_localFG.root" 1 hA "$SIGCM" 4 "$EBEAM" local "$KF" 0 "$RHOPOW"

echo "[$(stamp)] VALIDATION SAMPLES DONE"
ls -lh "$OUTDIR"
