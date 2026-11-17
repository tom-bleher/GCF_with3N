#!/usr/bin/env bash
# Waits for run_paper_production.sh to finish, then makes the figures.
# Only a production_pipeline.log written after this script started counts:
# an older log (left by an earlier production, or the copy in git) may
# already say "ALL PRODUCTION DONE" or "ERROR" and must not trigger anything.
cd "$(dirname "${BASH_SOURCE[0]}")"
LOG=production_pipeline.log
START_MARKER="$(mktemp)"   # mtime = watcher start
trap 'rm -f "$START_MARKER"' EXIT
if grep -q "ALL PRODUCTION DONE\|ERROR" "$LOG" 2>/dev/null; then
  echo "[$(date '+%F %T')] ignoring $LOG from $(date -r "$LOG" '+%F %T') (older than this watcher);" \
       "waiting for a production run to rewrite it (if that log is the run you want," \
       "run make_paper_figures.sh directly)"
fi
fresh_log() { [[ "$LOG" -nt "$START_MARKER" ]]; }
until fresh_log && grep -q "ALL PRODUCTION DONE" "$LOG" 2>/dev/null; do
  if fresh_log && grep -q "ERROR" "$LOG" 2>/dev/null; then
    echo "[$(date '+%F %T')] production reported an error; not making figures"; exit 1
  fi
  sleep "${POLL_SECONDS:-120}"
done
echo "[$(date '+%F %T')] production done; making figures"
bash make_paper_figures.sh
