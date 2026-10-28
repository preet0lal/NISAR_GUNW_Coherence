#!/usr/bin/env bash
set -euo pipefail
basedir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# MAAP DPS File input: register ONE required input named gunw_file.
# MAAP stages it in input/ or inputs/; use that file, not its source URL.
mkdir -p output
h5file=""
for folder in input inputs /input /inputs; do
  if [ -d "$folder" ]; then
    h5file=$(find -L "$folder" -maxdepth 3 -type f \( -iname '*.h5' -o -iname '*.hdf5' \) -print -quit)
    [ -n "$h5file" ] && break
  fi
done
# Optional local direct use: bash run_gunw.sh /path/file.h5
if [ -z "$h5file" ] && [ "${1:-}" != "" ] && [ -f "$1" ]; then h5file="$1"; fi
if [ -z "$h5file" ]; then
  echo 'ERROR: No GUNW .h5 found in MAAP staged input directory.' >&2
  exit 2
fi
python "$basedir/nisar_gunw_coherence_ease2.py" --input "$h5file" --output-dir output --resampling nearest
