#!/usr/bin/env python3
'''
Where a run's tiles have no data over land -- read from the finished tiles, not from footprints.

  gapReport /scratch/ianj/mosaics/cycle30ascending [--top 40]

For every tile in <work>/tiles/<product>/ (frequency A and B layers combined) this samples the
tile on a --cells x --cells grid, marks land (Natural Earth 50 m, as the coverage maps use) and
counts land cells without data. Writes <work>/quicklooks/gaps.txt (worst first: uncovered land
in km2 and as a share of the tile's land) and prints the top --top lines. Polar caps are sampled
in their own polar stereographic grid.

Also writes <work>/quicklooks/gaps.geojson: every no-data area inside each tile (land or water),
one feature per tile, in the tile's own CRS (property epsg) -- the input for
`globalGCOVTiles --fillGaps`, which adds granules only there.
'''
import argparse
import collections
import glob
import os
import sys

import numpy as np
from osgeo import gdal
from rasterio import features
from rasterio.transform import from_bounds
from shapely.geometry import box, mapping, shape
from shapely.ops import transform, unary_union

gdal.UseExceptions()
NODATA = -3000


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', help='run (work) directory')
    ap.add_argument('--product', default='gamma0')
    ap.add_argument('--cells', type=int, default=600, help='samples per tile side [600]')
    ap.add_argument('--top', type=int, default=40, help='lines to print [40]')
    args = ap.parse_args()
    import cartopy.feature as cf
    import pyproj
    land = unary_union(list(cf.LAND.with_scale('50m').geometries()))
    groups = collections.defaultdict(list)
    for tif in sorted(glob.glob(f'{args.work}/tiles/{args.product}/*.tif')):
        groups[os.path.basename(tif).split('.')[0]].append(tif)
    rows, gapFeats = [], []
    n = args.cells
    for name, tifs in groups.items():
        ds = gdal.Open(tifs[0])
        gt, nx, ny = ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize
        x0, y1 = gt[0], gt[3]
        x1, y0 = x0 + nx * gt[1], y1 + ny * gt[5]
        epsg = ds.GetSpatialRef().GetAuthorityCode(None)
        have = np.zeros((n, n), bool)
        for tif in tifs:
            t = gdal.Open(tif)          # keep the dataset alive while its band is read
            a = t.GetRasterBand(1).ReadAsArray(buf_xsize=n, buf_ysize=n)
            have |= a != NODATA
        # every no-data area of the tile, for globalGCOVTiles --fillGaps
        if not have.all():
            polys = [shape(g) for g, v in features.shapes((~have).astype(np.uint8), mask=~have,
                     transform=from_bounds(x0, y0, x1, y1, n, n)) if v]
            gapFeats.append({'type': 'Feature', 'geometry': mapping(unary_union(polys)),
                             'properties': {'tile': name, 'epsg': int(epsg),
                                            'gapFrac': round(float((~have).mean()), 4)}})
        if epsg == '4326':
            tileLand = land.intersection(box(x0, y0, x1, y1))
            lat = y0 + (np.arange(n)[::-1] + 0.5) * (y1 - y0) / n
            cellKm2 = (111.32 * (x1 - x0) / n) * (111.32 * (y1 - y0) / n) * np.cos(np.radians(lat))[:, None]
        else:
            fwd = pyproj.Transformer.from_crs(4326, int(epsg), always_xy=True).transform
            sign = -1 if y0 < 0 and int(epsg) == 3031 else 1
            capLand = land.intersection(box(-180, -90, 180, -60) if sign < 0 else box(-180, 60, 180, 90))
            tileLand = transform(fwd, capLand.segmentize(0.25)).buffer(0).intersection(box(x0, y0, x1, y1))
            cellKm2 = np.full((n, 1), (x1 - x0) / n * (y1 - y0) / n / 1e6)
        if tileLand.is_empty:
            continue
        isLand = features.rasterize([(tileLand, 1)], out_shape=(n, n), dtype='uint8',
                                    transform=from_bounds(x0, y0, x1, y1, n, n)).astype(bool)
        if not isLand.any():
            continue
        area = np.broadcast_to(cellKm2, (n, n))
        gapKm2 = float(area[isLand & ~have].sum())
        landKm2 = float(area[isLand].sum())
        rows.append((gapKm2, gapKm2 / landKm2, landKm2, name))
    rows.sort(reverse=True)
    os.makedirs(f'{args.work}/quicklooks', exist_ok=True)
    total = sum(r[0] for r in rows)
    landTotal = sum(r[2] for r in rows)
    lines = [f'{len(rows)} tiles with land; land without data {total / 1e6:.3f} M km2 '
             f'({total / max(landTotal, 1):.2%} of their land)',
             f'{"tile":<16}{"no-data land km2":>18}{"share of land":>15}']
    lines += [f'{name:<16}{g:>18,.0f}{f:>15.1%}' for g, f, _, name in rows]
    open(f'{args.work}/quicklooks/gaps.txt', 'w').write('\n'.join(lines) + '\n')
    import json
    with open(f'{args.work}/quicklooks/gaps.geojson', 'w') as fp:
        json.dump({'type': 'FeatureCollection', 'features': gapFeats}, fp)
    print(f'no-data areas of {len(gapFeats)} tiles: {args.work}/quicklooks/gaps.geojson')
    print('\n'.join(lines[:args.top + 2]))
    print(f'full list: {args.work}/quicklooks/gaps.txt')
    return 0


if __name__ == '__main__':
    sys.exit(main())
