#!/usr/bin/env bash
set -euo pipefail
# maap-py is supplied by MAAP's maap_base image, matching the working GCOV job.
python -m pip install --no-cache-dir 'numpy>=1.23,<3' 'h5py>=3.8' 'rasterio>=1.3' 'affine>=2.4' 's3fs>=2023.6'
python -c 'from maap.maap import MAAP; import s3fs, rasterio, h5py; print("MAAP/NISAR dependencies OK")'
