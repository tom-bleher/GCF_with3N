#!/usr/bin/env bash
# Waits for run_paper_production.sh to finish, then makes the figures.
cd "$(dirname "${BASH_SOURCE[0]}")"
until grep -q "ALL PRODUCTION DONE" production_pipeline.log 2>/dev/null; do
  if grep -q "ERROR" production_pipeline.log 2>/dev/null; then
    echo "[$(date '+%F %T')] production reported an error; not making figures"; exit 1
  fi
  sleep 120
done
echo "[$(date '+%F %T')] production done; making figures"
bash make_paper_figures.sh
