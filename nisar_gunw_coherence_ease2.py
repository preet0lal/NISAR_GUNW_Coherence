#!/usr/bin/env python3
"""Convert NISAR L-band GUNW coherenceMagnitude (~80 m) to NSIDC EASE2 M0p1 (~100 m).

Per-GUNW interferometric pair/per-polarization output. Mask invalid regions
using the native GUNW unwrappedInterferogram/mask and decoding the
reference/secondary inputDataExceptionMask flags packed into uint32 mask
values when present. Separate same-grid inputDataExceptionMask datasets
are also supported. A 2-pixel inward boundary buffer remains optional.
Valid water and low-coherence data are retained; no coherence threshold,
speckle filtering or mosaicking.

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
DEFAULT_EDGE_BUFFER_PIXELS = 2


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


def _gunw_mask_components(mask: np.ndarray):
    """Decode GUNW swath and input-data exception flags.

    GUNW mask layouts observed in products:
      * Legacy uint8: decimal WRS (water/ref-subswath/sec-subswath).
      * Packed uint32: low byte = 10*ref_subswath + sec_subswath;
        bits 8-15 = secondary inputDataExceptionMask;
        bits 16-23 = reference inputDataExceptionMask.

    Do not treat packed uint32 as a decimal WRS integer: doing so can
    incorrectly keep pixels flagged for sample slip/shift.
    """
    m = np.asarray(mask)
    if not np.issubdtype(m.dtype, np.integer):
        raise TypeError(f'GUNW mask must have integer dtype, got {m.dtype}')
    bits = m.astype(np.uint32, copy=False)
    low = bits & np.uint32(255)
    subswath_valid = ((low != 255) &
                      (((low // 10) % 10) != 0) &
                      ((low % 10) != 0))
    if m.dtype.itemsize > 1:
        # Decode both source RSLC exception channels. Nonzero means
        # reject, same criterion used for the standalone GCOV mask.
        secondary = (bits >> np.uint32(8)) & np.uint32(255)
        reference = (bits >> np.uint32(16)) & np.uint32(255)
        exception_valid = (secondary == 0) & (reference == 0)
    else:
        exception_valid = np.ones(m.shape, dtype=bool)
    return subswath_valid, exception_valid


def _gunw_mask_is_valid(mask: np.ndarray, use_exception_mask=True) -> np.ndarray:
    """Public native validity check; no low-coherence rejection."""
    swath, exception = _gunw_mask_components(mask)
    return swath & exception if use_exception_mask else swath


def _inward_buffer(valid: np.ndarray, radius: int) -> np.ndarray:
    """Square erosion of native validity; no scipy required."""
    if radius == 0:
        return valid
    rows, cols = valid.shape
    padded = np.pad(valid, ((radius, radius), (radius, radius)),
                    constant_values=False)
    buffered = np.ones((rows, cols), dtype=bool)
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            buffered &= padded[radius + dy:radius + dy + rows,
                               radius + dx:radius + dx + cols]
    return buffered


def _find_exception_dataset(grp: h5py.Group, pol: str, shape: tuple):
    """Only search same-grid parent/polarization groups; no cross-grid reuse."""
    for parent in (grp, grp[pol]):
        ds = parent.get('inputDataExceptionMask')
        if ds is None:
            continue
        if not isinstance(ds, h5py.Dataset):
            raise TypeError(f'{parent.name}/inputDataExceptionMask is not a dataset')
        if ds.shape != shape or not np.issubdtype(ds.dtype, np.integer):
            raise ValueError(
                f'{ds.name}: expected integer exception mask with shape {shape}, '
                f'got shape={ds.shape}, dtype={ds.dtype}')
        return ds
    return None


def _read_valid_mask_rows(mask_ds: h5py.Dataset, row: int, n: int,
                          height: int, buffer_pixels: int,
                          exception_ds=None, use_exception_mask=True):
    """Validity and statistics for a row block, with correct erosion halo."""
    row_start = max(0, row - buffer_pixels)
    row_end = min(height, row + n + buffer_pixels)
    raw = mask_ds[row_start:row_end, :]
    swath_ok, packed_exception_ok = _gunw_mask_components(raw)
    if use_exception_mask and exception_ds is not None:
        # Same as the GCOV code: exactly zero means no input anomalies.
        separate_exception_ok = (exception_ds[row_start:row_end, :] == 0)
    else:
        separate_exception_ok = np.ones(swath_ok.shape, dtype=bool)

    good = swath_ok.copy()
    if use_exception_mask:
        good &= packed_exception_ok & separate_exception_ok

    buffered = _inward_buffer(good, buffer_pixels)
    offset = row - row_start
    sl = slice(offset, offset+n)
    counts = {
        'invalid_swath': int(np.count_nonzero(~swath_ok[sl])),
        'invalid_packed_exception': (
            int(np.count_nonzero(~packed_exception_ok[sl]))
            if use_exception_mask else 0),
        'invalid_separate_exception': (
            int(np.count_nonzero(~separate_exception_ok[sl]))
            if use_exception_mask else 0),
        'removed_by_buffer': int(np.count_nonzero(good[sl] & ~buffered[sl])),
    }
    return buffered[sl], counts


def convert(input_path: str | Path, output_dir: str | Path, *, resampling='nearest',
            polarizations='all', overwrite=False,
            edge_buffer_pixels=DEFAULT_EDGE_BUFFER_PIXELS,
            use_exception_mask=True):
    if not isinstance(edge_buffer_pixels, int) or not (0 <= edge_buffer_pixels <= 20):
        raise ValueError('edge_buffer_pixels must be an integer between 0 and 20')
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
        if 'mask' not in grp:
            raise KeyError(
                f'GUNW WRS validity mask missing at {grp.name}/mask; '
                'refusing to export unmasked edge data')
        mask_ds = grp['mask']
        if not isinstance(mask_ds, h5py.Dataset) or mask_ds.shape != (height, width):
            raise ValueError(
                f'GUNW mask shape {getattr(mask_ds, "shape", None)} '
                f'does not match grid {(height, width)}')
        pols = [p for p in ('HH', 'VV') if p in grp and 'coherenceMagnitude' in grp[p]]
        print(f'[GUNW MASK] Source mask: {mask_ds.name}, dtype={mask_ds.dtype}, '
              f'packed_exception_bits={mask_ds.dtype.itemsize > 1}, '
              f'use_exception_mask={use_exception_mask}', flush=True)
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
                        tags = previous.tags()
                        if (tags.get('source_validity_mask') == 'GUNW_WRS_exception_v2' and
                            tags.get('mask_edge_buffer_pixels') == str(edge_buffer_pixels) and
                            tags.get('use_input_data_exception_mask') == str(use_exception_mask).lower()):
                            LOG.info('SKIP existing masked output: %s', out_file)
                            outputs.append(out_file)
                            continue
                        raise RuntimeError(
                            f'Output is from older/unmasked processing: {out_file}; '
                            'use --overwrite to regenerate with validity mask')
                raise RuntimeError(f'Existing output has wrong grid; choose --overwrite: {out_file}')
            fills = _valid_fill_values(ds)
            exception_ds = _find_exception_dataset(grp, pol, (height, width))
            if exception_ds is None and mask_ds.dtype.itemsize == 1 and use_exception_mask:
                print(f'[GUNW MASK] WARNING: {pol}: source uses uint8 WRS '
                      'and contains no separate inputDataExceptionMask; '
                      'input-data exception flags are unavailable. '
                      'Cannot guarantee sample-slip strips are removed.',
                      flush=True)
            elif exception_ds is not None:
                print(f'[GUNW MASK] {pol}: separate exception dataset = '
                      f'{exception_ds.name}', flush=True)
            rejected_mask = 0
            rejected_packed_exception = 0
            rejected_separate_exception = 0
            rejected_buffer = 0
            rejected_bad_coherence = 0
            kept = 0
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
                        good_swath, mask_counts = _read_valid_mask_rows(
                            mask_ds, row, n, height, edge_buffer_pixels,
                            exception_ds=exception_ds,
                            use_exception_mask=use_exception_mask)
                        good_value = np.isfinite(vals) & (vals >= 0) & (vals <= 1)
                        for fill in fills:
                            if np.isfinite(fill):
                                good_value &= (vals != fill)
                        rejected_mask += mask_counts['invalid_swath']
                        rejected_packed_exception += mask_counts['invalid_packed_exception']
                        rejected_separate_exception += mask_counts['invalid_separate_exception']
                        rejected_buffer += mask_counts['removed_by_buffer']
                        rejected_bad_coherence += int((~good_value).sum())
                        good = good_swath & good_value
                        kept += int(good.sum())
                        vals[~good] = np.nan
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
                                        quality_filter='native_GUNW_mask_plus_input_exceptions',
                                        source_validity_mask='GUNW_WRS_exception_v2',
                                        mask_dataset=mask_ds.name,
                                        mask_rule='subswath_R_and_S>0; packed_ref_sec_exception==0',
                                        use_input_data_exception_mask=str(use_exception_mask).lower(),
                                        packed_exception_mask_present=str(mask_ds.dtype.itemsize > 1).lower(),
                                        separate_exception_dataset=(exception_ds.name if exception_ds is not None else 'none'),
                                        mask_edge_buffer_pixels=str(edge_buffer_pixels),
                                        source_valid_pixels=str(kept),
                                        source_invalid_mask_pixels=str(rejected_mask),
                                        source_invalid_packed_exception_pixels=str(rejected_packed_exception),
                                        source_invalid_separate_exception_pixels=str(rejected_separate_exception),
                                        source_buffer_removed_pixels=str(rejected_buffer),
                                        source_invalid_coherence_pixels=str(rejected_bad_coherence),
                                        temporal_mosaic='none')
                    part.replace(out_file)
                finally:
                    if part.exists(): part.unlink()
            with rasterio.open(out_file) as check:
                if (check.crs != CRS.from_epsg(OUT_EPSG) or
                    not np.allclose(check.res, (PIXEL_M, PIXEL_M), atol=1e-8)):
                    raise AssertionError('Output grid verification failed')
            LOG.info(
                'GUNW %s: retained=%d invalid_swath=%d packed_exception=%d '
                'separate_exception=%d buffer_removed=%d invalid_coherence=%d '
                '(native buffer=%d px, water retained)',
                pol, kept, rejected_mask, rejected_packed_exception,
                rejected_separate_exception, rejected_buffer,
                rejected_bad_coherence, edge_buffer_pixels)
            LOG.info('SAVED %s', out_file)
            outputs.append(out_file)
    return outputs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='Local GUNW HDF5 path')
    p.add_argument('--output-dir', default='output')
    p.add_argument('--resampling', choices=['average','nearest','bilinear'], default='nearest')
    p.add_argument('--polarizations', choices=['all','HH','VV'], default='all')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--edge-buffer-pixels', type=int, default=DEFAULT_EDGE_BUFFER_PIXELS,
                   help='Inward trim of invalid-mask boundaries on native GUNW '
                        'grid; default 2; 0 keeps mask without buffer')
    p.add_argument('--use-exception-mask', choices=['yes','no'], default='yes',
                   help='Reject any flagged reference/secondary input data '
                        'exceptions in packed GUNW mask and separate mask (default yes)')
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    files = convert(args.input, args.output_dir, resampling=args.resampling,
                    polarizations=args.polarizations, overwrite=args.overwrite,
                    edge_buffer_pixels=args.edge_buffer_pixels,
                    use_exception_mask=(args.use_exception_mask == 'yes'))
    print('Outputs:', *map(str, files), sep='\n  ')


if __name__ == '__main__':
    main()
