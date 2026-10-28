#!/usr/bin/env python3
"""Convert NISAR L-band GUNW coherenceMagnitude (~80 m) to NSIDC EASE2 M0p1 (~100 m).

Per-GUNW interferometric pair/per-polarization output. No denoising, masking by
quality or land cover, thresholding, or temporal/spatial mosaicking. Only NODATA
handling and required reprojection/resampling are performed.

Usage:
    python nisar_gunw_coherence_ease2.py --input granule.h5 --output-dir output
    python nisar_gunw_coherence_ease2.py --input granule.h5 --output-dir output \
        --resampling nearest  # preserve individual sample magnitudes
"""
from __future__ import annotations
import argparse
import logging
import math
import tempfile
from pathlib import Path

import h5py
import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import Window

LOG = logging.getLogger('gunw_to_ease2')
OUT_EPSG = 6933
# EXACT grid definition used by the user's GCOV processing (NSIDC M0p1).
# "100 m" is nominal. Native posting is 100.089502334956 m.
PIXEL_M = 100.089502334956
EASE2_COLS = 347040
EASE2_ROWS = 146160
EASE2_LEFT = -0.5 * EASE2_COLS * PIXEL_M
EASE2_TOP = 0.5 * EASE2_ROWS * PIXEL_M
GROUP = 'unwrappedInterferogram'  # ~80 m input; no wrappedInterferogram (~20 m)
FREQUENCY = 'frequencyA'  # mission-generated GUNW is frequency A only
BLOCK_ROWS = 256


def _source_grid(hf: h5py.File):
    for radar in ('LSAR', 'SSAR'):
        path = f'/science/{radar}/GUNW/grids/{FREQUENCY}/{GROUP}'
        if path in hf:
            group = hf[path]
            xs = np.asarray(group['xCoordinates'][:], dtype='float64')
            ys = np.asarray(group['yCoordinates'][:], dtype='float64')
            # projection is the scalar EPSG code in NISAR L2 CF grid mapping
            epsg = int(group['projection'][()])
            if epsg <= 0:
                raise ValueError(f'Invalid input EPSG code: {epsg}')
            if xs.ndim != 1 or ys.ndim != 1 or len(xs) < 2 or len(ys) < 2:
                raise ValueError('Invalid xCoordinates/yCoordinates axes')
            dx = float(np.median(np.diff(xs)))
            dy = float(np.median(np.diff(ys)))
            if not dx or not dy:
                raise ValueError('Zero source pixel size')
            if not (np.allclose(np.diff(xs), dx, rtol=1e-5, atol=1e-4) and
                    np.allclose(np.diff(ys), dy, rtol=1e-5, atol=1e-4)):
                raise ValueError('Source axes are not regularly spaced')
            # NISAR coordinate vectors give pixel CENTERS, rasterio needs EDGES.
            transform = Affine(dx, 0, xs[0] - dx/2, 0, dy, ys[0] - dy/2)
            return group, CRS.from_epsg(epsg), transform, len(xs), len(ys)
    raise KeyError('GUNW frequencyA/unwrappedInterferogram group not found')


def _target_grid(src_crs, src_transform, width, height):
    """Crop the canonical NSIDC M0p1 global raster to the scene footprint.

    Uses integer global row/column indices; all frames align with GCOV output.
    """
    left, bottom, right, top = rasterio.transform.array_bounds(
        height, width, src_transform)
    xmin, ymin, xmax, ymax = transform_bounds(
        src_crs, CRS.from_epsg(OUT_EPSG), left, bottom, right, top,
        densify_pts=41)

    col0 = max(0, math.floor((xmin - EASE2_LEFT) / PIXEL_M))
    col1 = min(EASE2_COLS, math.ceil((xmax - EASE2_LEFT) / PIXEL_M))
    row0 = max(0, math.floor((EASE2_TOP - ymax) / PIXEL_M))
    row1 = min(EASE2_ROWS, math.ceil((EASE2_TOP - ymin) / PIXEL_M))
    ncols, nrows = col1-col0, row1-row0
    if ncols <= 0 or nrows <= 0:
        raise ValueError(f'Invalid output grid: {nrows} x {ncols}')
    origin_x = EASE2_LEFT + col0 * PIXEL_M
    origin_y = EASE2_TOP - row0 * PIXEL_M
    return Affine(PIXEL_M, 0, origin_x, 0, -PIXEL_M, origin_y), ncols, nrows


def _valid_fill_values(ds):
    values = []
    for key in ('_FillValue', 'missing_value'):
        if key in ds.attrs:
            arr = np.asarray(ds.attrs[key]).ravel()
            if arr.size:
                values.append(float(arr[0]))
    return values


def convert(input_path: str | Path, output_dir: str | Path, *, resampling='nearest',
            polarizations='all', overwrite=False):
    input_path = Path(input_path)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    method = {'average': Resampling.average, 'nearest': Resampling.nearest,
              'bilinear': Resampling.bilinear}[resampling]
    outputs = []
    with h5py.File(input_path, 'r') as hf:
        grp, src_crs, src_tr, width, height = _source_grid(hf)
        dst_tr, dst_width, dst_height = _target_grid(src_crs, src_tr, width, height)
        pols = [p for p in ('HH', 'VV') if p in grp and 'coherenceMagnitude' in grp[p]]
        if polarizations.upper() != 'ALL':
            pols = [p for p in pols if p == polarizations.upper()]
        if not pols:
            raise KeyError(f'No requested HH/VV coherenceMagnitude in {grp.name}; available: {list(grp.keys())}')
        LOG.info('Input %s: CRS=%s shape=%sx%s, output EPSG:%d %sx%s',
                 input_path.name, src_crs, width, height, OUT_EPSG, dst_width, dst_height)
        for pol in pols:
            ds = grp[pol]['coherenceMagnitude']
            if ds.ndim != 2 or ds.shape != (height, width):
                raise ValueError(f'{ds.name} shape {ds.shape} does not match coordinates {(height, width)}')
            # Full filename retained: both dates, both frame fields and product revision.
            out_file = output_dir / f'{input_path.stem}_{FREQUENCY}_{pol}_coherence_100m_EASE2.tif'
            if out_file.is_file() and not overwrite:
                with rasterio.open(out_file) as previous:
                    if (previous.count == 1 and previous.crs == CRS.from_epsg(OUT_EPSG)
                        and previous.width == dst_width and previous.height == dst_height
                        and previous.transform.almost_equals(dst_tr)):
                        LOG.info('SKIP existing valid output: %s', out_file)
                        outputs.append(out_file)
                        continue
                raise RuntimeError(f'Existing output has wrong grid; choose --overwrite: {out_file}')
            fills = _valid_fill_values(ds)
            # A temporary on-disk source raster enables GDAL to warp blockwise
            # without loading an entire frame into memory.
            with tempfile.TemporaryDirectory(prefix='gunw_') as tmp:
                source_tif = Path(tmp) / 'native_coherence.tif'
                with rasterio.open(
                    source_tif, 'w', driver='GTiff', width=width, height=height,
                    count=1, dtype='float32', crs=src_crs, transform=src_tr,
                    nodata=np.nan, tiled=True, blockxsize=256, blockysize=256,
                    compress='deflate', predictor=3, bigtiff='IF_SAFER') as src_writer:
                    for row in range(0, height, BLOCK_ROWS):
                        n = min(BLOCK_ROWS, height - row)
                        vals = np.asarray(ds[row:row+n, :], dtype='float32')
                        invalid = ~np.isfinite(vals)
                        for fill in fills:
                            if np.isfinite(fill): invalid |= (vals == fill)
                        if invalid.any():
                            vals[invalid] = np.nan
                        src_writer.write(vals, 1, window=Window(0, row, width, n))
                # Avoid partial files on errors; use .tif.part name with explicit GTiff driver.
                part = out_file.with_suffix('.tif.part')
                try:
                    with rasterio.open(source_tif) as src, rasterio.open(
                        part, 'w', driver='GTiff', width=dst_width, height=dst_height,
                        count=1, dtype='float32', crs=CRS.from_epsg(OUT_EPSG),
                        transform=dst_tr, nodata=np.nan, tiled=True,
                        blockxsize=256, blockysize=256, compress='deflate',
                        predictor=3, bigtiff='IF_SAFER') as dst:
                        reproject(
                            source=rasterio.band(src, 1),
                            destination=rasterio.band(dst, 1),
                            src_transform=src_tr, src_crs=src_crs,
                            src_nodata=np.nan, dst_transform=dst_tr,
                            dst_crs=CRS.from_epsg(OUT_EPSG), dst_nodata=np.nan,
                            resampling=method, init_dest_nodata=True,
                            num_threads=2, warp_mem_limit=256)
                        dst.update_tags(product='NISAR GUNW coherenceMagnitude',
                                        source_granule=input_path.name,
                                        source_dataset=ds.name, source_epsg=src_crs.to_epsg(),
                                        target_grid='EASE-Grid 2.0 NSIDC M0p1',
                                        resampling=resampling, unit='unitless',
                                        quality_filter='none', temporal_mosaic='none')
                    part.replace(out_file)
                finally:
                    if part.exists(): part.unlink()
            with rasterio.open(out_file) as check:
                if (check.crs != CRS.from_epsg(OUT_EPSG) or
                    not np.allclose(check.res, (PIXEL_M, PIXEL_M), atol=1e-8)):
                    raise AssertionError('Output grid verification failed')
            LOG.info('SAVED %s', out_file)
            outputs.append(out_file)
    return outputs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='Local GUNW HDF5 path (DPS stages file input)')
    p.add_argument('--output-dir', default='output')
    p.add_argument('--resampling', choices=['average','nearest','bilinear'], default='nearest')
    p.add_argument('--polarizations', choices=['all','HH','VV'], default='all')
    p.add_argument('--overwrite', action='store_true')
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    files = convert(args.input, args.output_dir, resampling=args.resampling,
                    polarizations=args.polarizations, overwrite=args.overwrite)
    print('Outputs:', *map(str, files), sep='\n  ')


if __name__ == '__main__':
    main()
