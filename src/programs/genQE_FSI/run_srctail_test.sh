#!/usr/bin/env bash
# Step 5: retest the LocalFGM/Correlated cascade medium under the corrected
# event definition. PWIA needs no counterpart (no cascade => medium irrelevant).
set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GENIE_OVERRIDE_DIR="$PWD/config/genie_override_srctail"
JOBS=120 FORCE=1 ./run_2N_parallel.sh 2000000 events/natalie_val/fsi_local_srctail.root 1 hA 0.150 2 5.01 global 0.25 0 1 local 0
JOBS=120 FORCE=1 ./run_2N_parallel.sh 2000000 events/natalie_val/mf_srctail.root      1 hA 0.150 4 5.01 local  0.25 0 1
echo "SRCTAIL TEST DONE"
