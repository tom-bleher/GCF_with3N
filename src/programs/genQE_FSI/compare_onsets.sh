#!/usr/bin/env bash
# Run the paper analysis on each SRC-onset validation sample and print the
# numbers side by side. One variable changes between the two columns.
set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY=/home/ilayws/Ilay-Generators/src/programs/genQE_3N/.venv/bin/python
VAL=events/natalie_val
SEED="${SEED:-12345}"

for onset in fixed local; do
  out="analysis/plots/natalie_paper/val_$onset"
  mkdir -p "$out"
  echo "############### SRC onset = $onset"
  "$PY" analysis/natalie_paper_plots.py \
      --pwia "$VAL/pwia_$onset.root" --fsi "$VAL/fsi_$onset.root" \
      --fsi-mf "$VAL/mf_localFG.root" --outdir "$out" --seed "$SEED" \
      2>&1 | grep -vE "RuntimeWarning|np\.sqrt|return p0|^\s+mmiss" \
           | tee "analysis/plots/natalie_paper/val_$onset/run.log" \
           | grep -E "Transparenc|e,e'p\)|e,e'pp\)|Q2 pre-cut|channels|Saved"
  echo
done
