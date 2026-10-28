# Fixes for MAAP GUNW 401/AccessDenied and PROJ proj.db errors

**Latest correction:** If the worker reports `CRSError: Cannot find proj.db`, the
`nisar_gunw_proj.py` helper now locates Rasterio's compatible PROJ database,
configures the active GDAL/PROJ search path, and verifies EPSG:4326 and
EPSG:6933 **before downloading the GUNW HDF5**. The Docker build also tests
EPSG lookup, so missing PROJ data can be detected at registration time.
The notebook and positional `s3_url` input configuration are unchanged.
Rebuild/re-register in MAAP to pick up the updated GitHub commit.

Your previous GUNW algorithm registered `gunw_file` as a **File Input**. That made MAAP's *separate pre-job staging service* download NISAR files without the user's NISAR authorization. S3 gave AccessDenied, while HTTPS gave Earthdata 401.

The working GCOV repository `preet0lal/NISAR_GCOV_Processing` uses a different strategy: `s3_url` is a **Positional Input** (text). Its run script uses `MAAP().aws.earthdata_s3_credentials(...)`, `s3fs`, then opens/downloads the NISAR product inside the job. This package copies that proven access design for GUNW.

## Files to upload to https://github.com/preet0lal/NISAR_GUNW_Coherence

Overwrite the three existing files:

- `run_gunw.sh` — now expects the sole positional text argument `$1`, a direct S3 URL.
- `build-env.sh` — installs `s3fs` and checks `maap-py` imports in the MAAP base environment.
- `nisar_gunw_coherence_ease2.py` — reprojects only coherenceMagnitude; aligned to GCOV's **exact NSIDC M0p1 grid**.

Add the **new** file:

- `nisar_gunw_dps_wrapper.py` — requests temporary S3 credentials *inside the job*, downloads the HDF5, processes it, deletes the HDF5. No credentials saved or printed.

The old GUNW script used exact 100 m pixels and a `(0,0)` origin; your actual GCOV algorithm instead uses NSIDC M0p1, ~100.089502334956 m. This package corrects the mismatch.

## MUST re-register algorithm in MAAP (changing the notebook alone is NOT enough)

Use your same repository `https://github.com/preet0lal/NISAR_GUNW_Coherence.git`, branch `main`, algorithm name `nisar-gunw-coherence-ease2`, version `main`.

Run Command: `NISAR_GUNW_Coherence/run_gunw.sh`
Build Command: `NISAR_GUNW_Coherence/build-env.sh`

**Remove** the old File Input named `gunw_file`. There must be **zero File Inputs**.

Add **one required Positional Input**, Name `s3_url`, Description `NISAR L2 provisional GUNW S3 URL`, default empty.

Use the same `maap_base` container and compatible resource queue as GCOV (64 GB is a reasonable starting point).

A successful registration/build is NOT a successfully processed GUNW data job.

## Run the notebook

Upload and open `Australia_NISAR_GUNW_MAAP_DPS_Provisional_S3Auth.ipynb` in MAAP. Its batch identifier is new (the previous two failed job records are excluded), and defaults to `SUBMIT_JOBS=False`, `MAX_NEW_SUBMISSIONS=1`, `MAX_IN_FLIGHT=1`. Inspect manifest and registered algorithm; change only `SUBMIT_JOBS=True` for ONE job. Check output of the actual DPS job before increasing concurrency.

The worker will print that it is requesting credentials and downloading via S3, then save one TIFF per available HH/VV in `output/`. No quality filtering/averaging. Grid EPSG:6933 NSIDC M0p1, nominal 100 m, nearest resampling.

If the MAAP worker cannot acquire temporary credentials, investigate its MAAP authentication or deployment; do not embed temporary credentials, Earthdata passwords, or tokens in your GitHub repository or job arguments.

## Local syntax check

```bash
python -m py_compile nisar_gunw_dps_wrapper.py nisar_gunw_coherence_ease2.py
bash -n run_gunw.sh build-env.sh
```

## References

- MAAP GUNW/GCOV direct S3 with temporary credentials: https://docs.maap-project.org/en/develop/science/NISAR/NISAR_access.html
- MAAP DPS positional vs file input: https://docs.maap-project.org/en/latest/technical_tutorials/dps_tutorial/dps_tutorial_demo.html
- NISAR S3 credentials/region: https://nisar-docs.asf.alaska.edu/aws-s3-access/
