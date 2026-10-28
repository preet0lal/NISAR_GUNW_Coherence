#!/usr/bin/env python3
"""DPS runner: use temporary NISAR S3 credentials *inside* the MAAP worker.

MAAP input must be a POSITIONAL TEXT parameter, not a File Input.
Matches the user's working GCOV `nisar_dps_wrapper.py` access pattern.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
from urllib.parse import urlparse

ASF_CREDS = "https://nisar.asf.earthdatacloud.nasa.gov/s3credentials"
BUCKET = "sds-n-cumulus-prod-nisar-products"
COLLECTION = "NISAR_L2_GUNW_PROVISIONAL_V1"


def check_input_url(value: str) -> str:
    u = urlparse(value)
    if u.scheme != 's3' or u.netloc != BUCKET:
        raise ValueError('Expected S3 URL in the authorized NISAR production bucket')
    parts = u.path.lstrip('/').split('/')
    if len(parts) != 3 or parts[0] != COLLECTION:
        raise ValueError('Only provisional NISAR L2 GUNW product paths are permitted')
    granule, filename = parts[1], parts[2]
    if not granule.startswith('NISAR_L2_PR_GUNW_') or filename != granule + '.h5':
        raise ValueError('Input filename must match NISAR provisional GUNW product directory')
    if u.query or u.fragment:
        raise ValueError('S3 URL must not include query or fragment')
    return filename


def run(s3_url: str):
    filename = check_input_url(s3_url)
    try:
        import s3fs
        from maap.maap import MAAP
    except ImportError as exc:
        raise RuntimeError('MAAP base environment needs maap-py and s3fs installed') from exc
    # Validate Rasterio's PROJ database before downloading a multi-GB GUNW file.
    from nisar_gunw_proj import configure_proj_database
    configure_proj_database()
    from nisar_gunw_coherence_ease2 import convert

    input_dir = Path('input').resolve()
    output_dir = Path('output').resolve()
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    local_h5 = input_dir / filename

    print(f'DPS GUNW provisional input: {s3_url}', flush=True)
    print('Requesting temporary NISAR S3 credentials from MAAP...', flush=True)
    creds = MAAP().aws.earthdata_s3_credentials(ASF_CREDS)
    # Never print or persist temporary credentials.
    fs = s3fs.S3FileSystem(
        anon=False,
        key=creds['accessKeyId'],
        secret=creds['secretAccessKey'],
        token=creds['sessionToken'],
        client_kwargs={'region_name': 'us-west-2'},
    )
    try:
        if not local_h5.exists() or local_h5.stat().st_size == 0:
            print(f'Downloading to {local_h5}', flush=True)
            fs.get(s3_url, str(local_h5))
        if not local_h5.is_file() or local_h5.stat().st_size <= 0:
            raise RuntimeError('Downloaded HDF5 is missing or empty')
        print(f'Downloaded {local_h5.stat().st_size / (1024**3):.3f} GiB', flush=True)
        produced = convert(local_h5, output_dir, resampling='nearest')
        if not produced:
            raise RuntimeError('Conversion completed but produced no GeoTIFF')
        for p in produced:
            print(f'SAVED: {p}', flush=True)
    finally:
        # Always remove downloaded HDF5 so DPS collects only desired GeoTIFFs.
        if local_h5.exists():
            local_h5.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--s3_url', required=True)
    args = parser.parse_args()
    run(args.s3_url)


if __name__ == '__main__':
    main()
