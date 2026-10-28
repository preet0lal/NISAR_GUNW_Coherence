#!/usr/bin/env bash
set -euo pipefail
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ "$#" -ne 1 ]; then
    echo "Usage: run_gunw.sh <s3_url> (ONE positional text input)" >&2
    exit 2
fi
mkdir -p input output
echo "NISAR GUNW → NSIDC EASE2 M0p1 (~100 m); provisional only"
python "${BASEDIR}/nisar_gunw_dps_wrapper.py" --s3_url "$1"
