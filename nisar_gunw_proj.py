"""Initialize the PROJ database used by Rasterio inside MAAP DPS images.

In mixed conda/pip images PROJ_DATA or PROJ_LIB can point at an absent or
incompatible proj.db. Prefer Rasterio's own bundled database and verify EPSG.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def configure_proj_database() -> Path:
    import rasterio
    from rasterio.crs import CRS

    candidates = [
        Path(rasterio.__file__).resolve().parent / 'proj_data',
        Path(sys.prefix) / 'share' / 'proj',
        Path('/opt/conda/share/proj'),
        Path('/usr/local/share/proj'),
        Path('/usr/share/proj'),
    ]
    for env_name in ('PROJ_DATA', 'PROJ_LIB'):
        value = os.environ.get(env_name)
        if value:
            candidates.append(Path(value))
    seen = set()
    errors = []
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        if not (candidate / 'proj.db').is_file():
            continue
        try:
            # Explicitly configure GDAL/PROJ after Rasterio has been imported.
            # Changing os.environ alone can be too late in mixed environments.
            if hasattr(rasterio._env, 'set_proj_data_search_path'):
                rasterio._env.set_proj_data_search_path(str(candidate))
            os.environ['PROJ_DATA'] = str(candidate)
            os.environ['PROJ_LIB'] = str(candidate)
            CRS.from_epsg(4326)
            CRS.from_epsg(6933)
            print(f'[GUNW PROJ] Using {candidate / "proj.db"}', flush=True)
            return candidate
        except Exception as exc:
            errors.append(f'{candidate}: {exc}')

    raise RuntimeError(
        'No compatible proj.db found for Rasterio/GDAL. '
        f'Rasterio: {rasterio.__file__}. '
        f'Checked: {list(map(str, seen))}. '
        f'Failures: {errors}'
    )
