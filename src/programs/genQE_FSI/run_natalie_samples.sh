#!/usr/bin/env bash
#
# run_natalie_samples.sh — produce the matched samples for the
# Wright et al. (arXiv:2104.05090, doi:10.48550/arXiv.2104.05090) paper
# comparison, with the paper's conventions PINNED:
#
#   * Ebeam = 5.01 GeV, 12C, sigma_CM = 0.150 GeV/c
#   * FSI model hA (paper page 3), shared record (both nucleons in ONE
#     GENIE record, as in her eGENIE chain)
#   * SRC-pair FSI vertex ~ rho^1(r) r^2 (pair position follows the one-body
#     charge density, paper page 3) — NOT the rho^2 project default
#   * MF sample: single struck nucleon from a LOCAL Fermi gas (paper page 6)
#   * kF = 0.25 GeV/c (MF Fermi-gas ceiling; also the "fixed" SRC onset)
#   * Intranuke cascade targets from genie::LocalFGM/Correlated (20% SRC
#     tail) via config/genie_override_srctail — this is what her
#     config/G18_10a/ModelConfiguration.xml:37 selects. See that
#     directory's README.txt.
#
#   * SRC onset (see SRC_analysis_2N argv[12]): "local" by default, i.e.
#     k_rel > kF(r) at the pair vertex, matching LocalFGM::ProbDistro.
#     The mean-field sample is already capped at the local kF(r), so this
#     makes the two samples exactly complementary; a constant floor would
#     leave the nucleons between kF(r) and kF in neither sample (78% of
#     vertices have kF(r) < 0.25 GeV/c). Override with SRC_ONSET=fixed to
#     generate the systematic variant.
#
# PWIA and FSI SRC samples are generated with IDENTICAL primary
# configuration (only doFSI differs) so the attempts-based luminosity
# matching in natalie_paper_plots.py is exact.
#
# Usage:  ./run_natalie_samples.sh [TOTAL_SRC] [TOTAL_MF]
#   TOTAL_SRC  successful events per SRC sample  (default 20000000)
#   TOTAL_MF   successful events for the MF sample (default 60000000 — its
#              survival rate after the SRC cuts is far lower, and Figs 8/10
#              are visibly statistics-limited otherwise)
# Environment:
#   JOBS=N              parallel workers per sample
#   SRC_ONSET=local|fixed
#   GENIE_OVERRIDE_DIR  override the cascade-target config (see above)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TOTAL_SRC="${1:-20000000}"
TOTAL_MF="${2:-60000000}"
OUTDIR="events/natalie"
mkdir -p "$OUTDIR"

SIGCM=0.150
EBEAM=5.01
KF=0.25
RHOPOW=1                       # Wright et al. convention
SRC_ONSET="${SRC_ONSET:-local}"
REQP1=0                        # see SRC_analysis_2N argv[13]

# Cascade targets with the 20% SRC tail, matching her ModelConfiguration.
export GENIE_OVERRIDE_DIR="${GENIE_OVERRIDE_DIR:-$SCRIPT_DIR/config/genie_override_srctail}"

echo "### natalie samples -> $OUTDIR"
echo "###   SRC: $TOTAL_SRC events each   MF: $TOTAL_MF events"
echo "###   hA, rho^1 vertex, local-FG MF, src_onset=$SRC_ONSET"
echo "###   cascade override: $GENIE_OVERRIDE_DIR"

#                          TOTAL         OUTPUT                           DOFSI MODEL SIGCM  WF EBEAM  FG     KF  INDEP RHOPOW ONSET REQP1
FORCE=1 ./run_2N_parallel.sh "$TOTAL_SRC" "$OUTDIR/pwia_SRC.root"            0    hA  "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW" "$SRC_ONSET" "$REQP1"
FORCE=1 ./run_2N_parallel.sh "$TOTAL_SRC" "$OUTDIR/fsi_hA_SRC_rho1.root"     1    hA  "$SIGCM" 2 "$EBEAM" global "$KF" 0 "$RHOPOW" "$SRC_ONSET" "$REQP1"
FORCE=1 ./run_2N_parallel.sh "$TOTAL_MF"  "$OUTDIR/fsi_hA_MF_localFG.root"   1    hA  "$SIGCM" 4 "$EBEAM" local  "$KF" 0 "$RHOPOW"

echo "### natalie samples complete:"
ls -lh "$OUTDIR"
