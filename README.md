# NISAR GUNW coherence → 100 m EASE2, Australia (MAAP DPS)

Includes the **modified original GCOV MAAP submission notebook**, and a **new GUNW coherence DPS worker**. **The original `nisar-gcov-ease2` registration cannot be reused for GUNW.**

## Inputs and outputs

- CMR collection **`NISAR_L2_GUNW_PROVISIONAL_V1` only**. Beta products are never searched or submitted; granule names must begin `NISAR_L2_PR_GUNW_`.
- L-band GUNW `frequencyA/unwrappedInterferogram/HH or VV/coherenceMagnitude` (~80 m source). Picks the available co-polarization(s).
- Per-granule/per-interferometric-pair GeoTIFF, using native HDF5 pixel-center coordinates and native EPSG, projected to **EPSG:6933**, **100 m × 100 m** and snapped to multiples of 100 m (common global pixel-edge origin `0,0`).
- Default output resampling **nearest**, preserving individual native coherence magnitudes without spatial averaging. Optional `--resampling average` for area-weighted aggregation from 80 m to 100 m (not a speckle filter). The file's fill/nodata is preserved as NaN. No landcover mask, connected-component mask, coherence threshold, temporal composite, or smoothing filter.
- Retains all original GUNW filename tokens: two acquisition time periods, both frame fields, processing version, etc., plus `_frequencyA_HH_coherence_100m_EASE2.tif` (or VV).
- Output TIFF in DPS `output/` directory, which MAAP collects in individual job output folders.

## Register in MAAP DPS

1. Place `nisar_gunw_coherence_ease2.py`, `run_gunw.sh`, and `build-env.sh` in a **public repository** accessible to MAAP (MAAP registration requires a reachable repository).
2. Register a **new** algorithm (`nisar-gunw-coherence-ease2`, version `main`) using the MAAP Register Algorithm UI. Supply `bash build-env.sh` as the build command and `bash run_gunw.sh` as the run command.
3. Register **exactly one required File Input named `gunw_file`**, supplying an HTTPS Earthdata NISAR GUNW HDF5 URL per job. MAAP stages this file into the worker's `input/` (or `/inputs`) directory. The run script automatically discovers the HDF5.
4. Open `Australia_NISAR_GUNW_Coherence_EASE2_100m_MAAP_DPS.ipynb` in the MAAP environment. Run discovery, inspect `manifest` and confirm `GUNW algorithm registered: True`.
5. Set `SUBMIT_JOBS=True` with `MAX_NEW_SUBMISSIONS=1` and run one MAAP smoke-test job. Inspect its output for correct georeferencing and values.
6. Set `MAX_NEW_SUBMISSIONS=None` and rerun the controller to submit the rest, up to `MAX_IN_FLIGHT=500` jobs in flight. Ledger CSVs persist under the MAAP private bucket. Previously successful/failed job records are not automatically resubmitted.

**Note:** Registered MAAP algorithms may have different file input handling between deployment modes; the run script checks the common staging paths `input/`, `inputs/`, `/input`, and `/inputs`. The provided notebook assumes the submitted `gunw_file=` keyword matches the registered File Input's name.

## Local single-file run (before MAAP registration)

```bash
python -m pip install numpy h5py rasterio affine
python nisar_gunw_coherence_ease2.py --input /path/NISAR_L2_PR_GUNW_....h5 --output-dir output
```

For optional area-averaging instead of nearest:

```bash
python nisar_gunw_coherence_ease2.py --input /path/NISAR_L2_PR_GUNW_....h5 --output-dir output --resampling average
```

In the notebook, use `RUN_LOCAL_SMOKE_TEST=True` to download/process one GUNW file *locally* with `earthaccess` (requires the worker file next to the notebook and suitable Python packages installed).

## Processing logic

1. Read the coordinate arrays (`xCoordinates`, `yCoordinates`) and `projection` from the product's **unwrappedInterferogram** group (not GCOV and not the 20 m wrapped group).
2. Write the native coherence data in blocks to a temporary GeoTIFF, with HDF5 `_FillValue` as NoData. No additional quality masking is applied.
3. Compute an EPSG:6933 footprint and snap the resulting bounds to the 100 m common lattice.
4. Reproject the raster to EPSG:6933 with nearest resampling (or optional average) to 100 m, compressed float32 GeoTIFF.
5. Skip already completed files *only when the existing file's raster dimensions, transform and CRS match*.

No input images have been provided for live validation; a synthetic HDF5 test is included in the development verification, but the MAAP registration and downloads must be tested in your MAAP workspace.
