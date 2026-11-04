#!/usr/bin/env bash
# Python packages for the web interface and gcf quicklook / summary, in <repo>/.venv,
# which ./gcf then uses by itself. Not needed if your python3 already has them. Never
# fatal: the pipeline itself runs without them.
#
#   tools/gcf/setup_python.sh [--force]
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
VENV=$REPO/.venv
REQ=$REPO/tools/gcf/requirements.txt
IMPORTS='import uproot, awkward, numpy, matplotlib, scipy'

if ! command -v python3 >/dev/null 2>&1; then
    echo "No python3: ./gcf needs Python 3.11 or newer."
    exit 0
fi
if ! python3 -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
    echo "./gcf needs Python 3.11 or newer; python3 here is $(python3 -V 2>&1). The compiled generators work without it."
    exit 0
fi
if [[ "${1:-}" != "--force" ]]; then
    if [[ -x "$VENV/bin/python3" ]] && "$VENV/bin/python3" -c "$IMPORTS" 2>/dev/null; then
        echo "Python packages: ready in .venv"
        exit 0
    fi
    if python3 -c "$IMPORTS" 2>/dev/null; then
        echo "Python packages: python3 already has them"
        exit 0
    fi
fi

made=""
if python3 -m venv "$VENV" >/dev/null 2>&1; then
    made=pip
elif command -v uv >/dev/null 2>&1 && uv venv --quiet --allow-existing --python "$(command -v python3)" "$VENV" >/dev/null 2>&1; then
    made=uv
fi
if [[ -z "$made" ]]; then
    rm -rf "$VENV"
    echo "Could not make a virtual environment, so plots and run summaries won't work yet."
    echo "On Ubuntu or Debian: sudo apt install python3-venv, then run $0 again."
    exit 0
fi
if [[ $made == pip ]]; then
    "$VENV/bin/python3" -m pip install --quiet --upgrade pip >/dev/null 2>&1
    "$VENV/bin/python3" -m pip install --quiet -r "$REQ" >/dev/null 2>&1
else
    VIRTUAL_ENV=$VENV uv pip install --quiet -r "$REQ" >/dev/null 2>&1
fi
if "$VENV/bin/python3" -c "$IMPORTS" 2>/dev/null; then
    echo "Python packages: installed in .venv"
else
    rm -rf "$VENV"
    echo "Could not install the Python packages, so plots and run summaries won't work yet (no network?)."
    echo "Run $0 again, or: pip install -r tools/gcf/requirements.txt"
fi
exit 0
