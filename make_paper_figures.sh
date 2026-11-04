#!/usr/bin/env bash
# Regenerate the paper figures (all except Fig 3) from the production samples
# with the unified detection-based (e,e'pp) selection. Outputs go to
# Paper_Figures_eepp/ together with a README of the region fractions.
# Usage: bash make_paper_figures.sh   (after run_paper_production.sh finished)
# PY, OUT and the four sample paths (S3N_FSI, S3N_PWIA, S2N, SMF) can be
# overridden from the environment, e.g. to read a production made in another
# checkout (absolute paths):
#   S3N_FSI=/path/events_3N_100M.root ... bash make_paper_figures.sh
set -uo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
G3="$ROOT_DIR/src/programs/genQE_3N"
G2="$ROOT_DIR/src/programs/genQE_FSI"
PY="${PY:-$G3/.venv/bin/python}"
[[ -x "$PY" ]] || PY="$(command -v python3)"
OUT="${OUT:-$ROOT_DIR/Paper_Figures_eepp}"
mkdir -p "$OUT" && OUT="$(cd "$OUT" && pwd)"   # absolute: the script cd's below

S3N_FSI="${S3N_FSI:-$G3/events/events_3N_100M.root}"
S3N_PWIA="${S3N_PWIA:-$G3/events/3N_PWIA_noCM_20M_events.root}"
S2N="${S2N:-$G2/events/misc/2N_hN_100M_6GeV_kF220_indep.root}"
SMF="${SMF:-$G2/events/misc/MF_localFG_hN_100M_6GeV_kF220.root}"
for f in "$S3N_FSI" "$S3N_PWIA" "$S2N" "$SMF"; do
  [[ -s "$f" ]] || { echo "ERROR: missing sample $f"; exit 1; }
done
KF=0.22

cd "$G3"
echo "=== Fig 1/2/4 (PWIA no-CM) ==="
$PY plotting/plot_theta_heatmap_paper.py   --input "$S3N_PWIA" --output "$OUT/Fig1.pdf"  || exit 1
$PY plotting/plot_fig2_momenta_by_config.py --input "$S3N_PWIA" --output "$OUT/Fig2.pdf" || exit 1
$PY plotting/plot_fig4_xB_by_config.py      --input "$S3N_PWIA" --output "$OUT/Fig4.pdf" || exit 1

echo "=== Fig 5 (3N+FSI vs 2N+FSI, unified (e,e'pp)) ==="
$PY plotting/plot_fig5_theta_heatmap_pair.py --input-3N "$S3N_FSI" --input-2N "$S2N" \
    --kF $KF --output "$OUT/Fig5.pdf" | tee "$OUT/Fig5_summary.txt" || exit 1
R3N=$(grep '\[3N+FSI\]' "$OUT/Fig5_summary.txt" | sed 's/.*L\/R=//')
R2N=$(grep '\[2N+FSI\]' "$OUT/Fig5_summary.txt" | sed 's/.*L\/R=//')
echo "R_3N=$R3N R_2N=$R2N"

echo "=== Fig 6 ==="
$PY plotting/plot_fig6_R_xL.py --R-3N "$R3N" --R-2N "$R2N" --output "$OUT/Fig6.pdf" || exit 1

echo "=== Fig 7 (MF fake-3N, unified (e,e'pp)) ==="
cd "$G2"
$PY plotting/plot_mf_fake3n_heatmap.py --input "$SMF" --output "$OUT/Fig7.pdf" \
    --ebeam 6.0 --kf $KF --src-cuts --theta-e-max 45 --theta-e-min 0 | tee "$OUT/Fig7_summary.txt" || exit 1

echo "=== recoil-threshold scan ==="
cd "$G3"
$PY plotting/lr_recoil_scan.py "3N+FSI=$S3N_FSI" "2N+FSI=$S2N" "MF=$SMF" --kF $KF \
    | tee "$OUT/LR_recoil_scan.txt" || exit 1

cat > "$OUT/README.txt" <<TXT
PAPER FIGURES -- flat-list FSI, unified detection-based (e,e'pp) selection
========================================================================
12C, 6 GeV, GENIE hN FSI (independent transport, LocalFGM no-tail override),
kF=0.22 GeV/c. 3N: 100M (CM on, sigCM=0.15, Q2>1 at generation); Fig1/2/4 from
a 20M PWIA no-CM sample. 2N: 100M AV18, sigCM=0.150. MF: 100M local FG.

Selection (identical for 3N, 2N, MF):
  detected proton = final-state proton with |p| > kF
  lead   = detected proton closest in angle to q, with theta_pq < 8 deg and
           |p_lead| > 0.75 |q|
  recoil = the other detected proton (exactly 2 detected protons required;
           neutrons are not detected and not vetoed)
  third  = reconstructed, p3 = -(p_lead - q + p_recoil), |p3| > kF
  theta_e < 45 deg, Q2 > 1 GeV^2, xB < 1.2, 0.25 < pmiss < 0.9 GeV/c
  theta12 = angle(p_lead - q, p_recoil); theta23 = angle(p_recoil, p3)

Fig1  no-FSI 3N configuration map.   Fig2  ground-state momenta by config.
Fig3  NOT regenerated.               Fig4  xB by config.
Fig5  3N+FSI vs 2N+FSI:  R_3N = $R3N   R_2N = $R2N   (see Fig5_summary.txt)
Fig6  R(x_L) with those values.
Fig7  MF fake-3N (see Fig7_summary.txt / Fig7_region_fractions.txt)
LR_recoil_scan.txt: L/R vs an extra |p_recoil| threshold for all three samples.
TXT
echo "FIGURES DONE -> $OUT"
