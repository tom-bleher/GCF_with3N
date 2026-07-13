#!/usr/bin/env bash
#
# run_natalie_samples.sh — produce the three matched samples for the
# Wright et al. (Natalie) paper comparison, with the paper's conventions
# PINNED:
#
#   * Ebeam = 5.01 GeV, 12C, sigma_CM = 0.150 GeV/c
#   * FSI model hA (paper page 3), shared record (standard eGENIE transport)
#   * SRC-pair FSI vertex ~ rho^1(r) r^2 (pair position follows the one-body
#     charge density, paper page 3) — NOT the rho^2 project default
#   * MF sample: single struck nucleon from a LOCAL Fermi gas (paper page 6)
#   * kF = 0.25 GeV/c (SRC p_rel > kF split; MF FG ceiling)
#
# PWIA and FSI SRC samples are generated with IDENTICAL primary
# configuration (only doFSI differs) so the attempts-based luminosity
# matching in natalie_paper_plots.py is exact.
#
# Usage:  ./run_natalie_samples.sh [TOTAL_PER_SAMPLE]
#   TOTAL_PER_SAMPLE   successful events per sample (default 20000000)
# Environment:
#   JOBS=N  parallel workers per sample (passed through to run_2N_parallel.sh)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TOTAL="${1:-20000000}"
OUTDIR="events/natalie"
mkdir -p "$OUTDIR"

SIGCM=0.150
EBEAM=5.01
KF=0.25
RHOPOW=1      # Wright et al. convention

echo "### natalie samples: $TOTAL events each -> $OUTDIR (hA, rho^1 vertex, local-FG MF)"

#                          TOTAL    OUTPUT                            DOFSI MODEL SIGCM  WF EBEAM  FG     KF  INDEP RHOPOW
FORCE=1 ./run_2N_parallel.sh "$TOTAL" "$OUTDIR/pwia_SRC.root"           0    hA  "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW"
FORCE=1 ./run_2N_parallel.sh "$TOTAL" "$OUTDIR/fsi_hA_SRC_rho1.root"    1    hA  "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW"
FORCE=1 ./run_2N_parallel.sh "$TOTAL" "$OUTDIR/fsi_hA_MF_localFG.root"  1    hA  "$SIGCM" 4 "$EBEAM" local  "$KF" 0 "$RHOPOW"

echo "### natalie samples complete:"
ls -lh "$OUTDIR"
