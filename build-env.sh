#!/usr/bin/env bash
set -euo pipefail

# MAAP's minimal maap_base image may not bundle the libraries needed by
# this worker. Install them explicitly during the Docker build.
echo "[GUNW BUILD] Python: $(python --version 2>&1)"
echo "[GUNW BUILD] Pip: $(python -m pip --version)"
echo "[GUNW BUILD] Installing dependencies..."
python -m pip install --no-cache-dir \
  'maap-py==5.1.0' \
  'numpy>=1.23,<3' \
  'h5py>=3.8' \
  'rasterio>=1.3' \
  'affine>=2.4' \
  's3fs>=2023.6'

echo "[GUNW BUILD] Checking imports..."
python - <<'PY'
from maap.maap import MAAP
import s3fs
import h5py
import rasterio
import numpy as np
from affine import Affine
print("[GUNW BUILD] All imports succeeded.")
print("[GUNW BUILD] rasterio", rasterio.__version__)
print("[GUNW BUILD] numpy", np.__version__)
PY
echo "[GUNW BUILD] COMPLETE"
