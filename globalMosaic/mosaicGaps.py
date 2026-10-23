#!/usr/bin/env python3
'''
Find the holes in a tiled mosaic: places with no data that a granule DOES cover.

  mosaicGaps /scratch/ianj/mosaics/ease/global/tiles/gamma0 \
             --catalogue cycle30/catalogue.geojson \
             --fillCatalogues cycle30/catalogue0*.geojson --allowV --out gaps

The distinction that matters is between a hole and an absence. Open ocean away from the ice
is empty because NISAR never looked there, and reporting that as a gap buries the real ones.
So the expected area is the union of the granule FOOTPRINTS from the same catalogues the
tiler used, and a gap is a cell inside that union with no data in the mosaic -- something was
acquired and did not reach the product. Causes worth separating, which the report does by
attributing every gap to its tile:

  * the tile failed or has not run yet          -> the tile is absent or mostly empty
  * the tiler selected nothing there            -> tiles.geojson has no entry, or selected 0
  * geomosaic dropped it                        -> the tile exists and is covered around it

A footprint overstates: it is where the frame was acquired, not where it holds valid data, so
a partly unfocused frame shows up as a gap that nothing can fill. Treat the gap area as an
upper bound unless the tiler was given --validFootprints.

Outputs (--out PREFIX)
  PREFIX.txt       the summary printed to the terminal
  PREFIX.geojson   one polygon per gap, in lat/lon, with area and tile name
  PREFIX.png       the gaps in red over the mosaic, if --figure
'''
import argparse
import glob
import json
import os

import numpy as np
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import mapping, shape
from shapely.ops import transform, unary_union

from .easeTiles import isWagner7, wagner7Inverse
from .globalGCOVTiles import loadGranules

NODATA = -3000


def mosaicGrid(src, cellM):
    ''' The mosaic read at cellM spacing: (valid mask, transform, CRS). Averaged, not
    decimated, so a cell counts as covered if ANY of its full-resolution pixels had data. '''
    from osgeo import gdal
    gdal.UseExceptions()
    if os.path.isdir(src):
        tifs = sorted(glob.glob(f'{src}/*.tif'))
        if not tifs:
            raise SystemExit(f'mosaicGaps: no .tif in {src}')
        ds = gdal.BuildVRT('', tifs, VRTNodata=NODATA)
    else:
        ds = gdal.Open(src)
    # Translate, not Warp: Warp needs an inverse transform, which Wagner VII has not got
    out = gdal.Translate('', ds, format='MEM', xRes=cellM, yRes=cellM,
                         resampleAlg='average', noData=NODATA)
    a = out.ReadAsArray()
    if a.ndim == 3:
        a = a[0]
    gt = out.GetGeoTransform()
    crs = out.GetSpatialRef().ExportToWkt()
    del out, ds
    return (a != NODATA) & np.isfinite(a), from_origin(gt[0], gt[3], gt[1], -gt[5]), crs, a.shape


def expectedMask(cands, crs, shape2d, tr):
    ''' Where a granule footprint lands, on the mosaic's own grid. '''
    import pyproj
    fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
    geoms = []
    for g in cands:
        try:
            p = transform(fwd, shape(g['raw'])).buffer(0)
        except Exception:
            continue
        if not p.is_empty:
            geoms.append(p)
    if not geoms:
        raise SystemExit('mosaicGaps: no granule footprints projected successfully')
    return features.rasterize(((p, 1) for p in geoms), out_shape=shape2d, transform=tr,
                              dtype='uint8', all_touched=True).astype(bool)


def toLonLat(geom, crs):
    ''' Projected polygon -> lat/lon. Wagner VII has no PROJ inverse, so use ours. '''
    if isWagner7(crs):
        return transform(lambda x, y: wagner7Inverse(x, y), geom)
    import pyproj
    inv = pyproj.Transformer.from_crs(crs, 4326, always_xy=True).transform
    return transform(inv, geom)


def tileLookup(tilesGeojson):
    ''' name -> projected box, for attributing a gap to the tile it falls in. '''
    if not tilesGeojson or not os.path.exists(tilesGeojson):
        return []
    feats = json.load(open(tilesGeojson))['features']
    out = []
    for f in feats:
        p = f['properties']
        if all(k in p for k in ('xMin', 'yMin', 'xMax', 'yMax')):
            out.append((p['name'], p['xMin'], p['yMin'], p['xMax'], p['yMax'],
                        p.get('selected', 0)))
    return out


def whichTile(x, y, tiles):
    for name, x0, y0, x1, y1, sel in tiles:
        if x0 <= x <= x1 and y0 <= y <= y1:
            return name, sel
    return '-', 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('raster', help='mosaic VRT, GeoTIFF, or directory of tiles')
    ap.add_argument('--catalogue', required=True, help='the catalogue the tiler used')
    ap.add_argument('--fillCatalogues', nargs='+', default=[], help='and its fill catalogues')
    ap.add_argument('--direction', choices=['A', 'D', 'both'], default='A',
                    help='orbit direction, matching the run [A]')
    ap.add_argument('--allowV', action='store_true', help='V fallback, matching the run')
    ap.add_argument('--cellKm', type=float, default=10., help='gap raster cell, km [10]')
    ap.add_argument('--minAreaKm2', type=float, default=2500.,
                    help='ignore gaps smaller than this [2500, i.e. 50x50 km]')
    ap.add_argument('--tiles', default=None, help='tiles.geojson, to name the tile a gap is in')
    ap.add_argument('--top', type=int, default=25, help='gaps listed in the summary [25]')
    ap.add_argument('--out', default='gaps', help='output prefix [gaps]')
    ap.add_argument('--figure', action='store_true', help='also write PREFIX.png')
    args = ap.parse_args()

    cellM = args.cellKm * 1000.
    valid, tr, crs, shape2d = mosaicGrid(args.raster, cellM)
    cands = loadGranules(args.catalogue, allowV=args.allowV)
    for k, fc in enumerate(args.fillCatalogues, 1):
        cands += loadGranules(fc, allowV=args.allowV, fill=k)
    if args.direction != 'both':
        cands = [g for g in cands if g['direction'] == args.direction]
    expected = expectedMask(cands, crs, shape2d, tr)

    cell = (args.cellKm ** 2)
    gaps = expected & ~valid
    lines = [f'mosaic      {args.raster}',
             f'grid        {shape2d[1]} x {shape2d[0]} cells of {args.cellKm:g} km',
             f'granules    {len(cands)} ({args.direction}'
             + (', V fallback' if args.allowV else '') + ')',
             f'expected    {expected.sum() * cell:,.0f} km2 covered by a granule footprint',
             f'in product  {(expected & valid).sum() * cell:,.0f} km2 '
             f'({100. * (expected & valid).sum() / max(expected.sum(), 1):.1f}%)',
             f'GAPS        {gaps.sum() * cell:,.0f} km2 '
             f'({100. * gaps.sum() / max(expected.sum(), 1):.1f}%)',
             f'(data outside any footprint: {(valid & ~expected).sum() * cell:,.0f} km2)',
             '',
             'A footprint is the ACQUISITION polygon, not the valid-data polygon: a frame can',
             'be partly unfocused (mask 0), which geomosaic drops. Such areas appear here as',
             'gaps although nothing is wrong. globalGCOVTiles --validFootprints exists for',
             'exactly this, and passing those files to the tiler shrinks both the expected',
             'area and the gap list.']

    tiles = tileLookup(args.tiles)
    polys = []
    for geom, _ in features.shapes(gaps.astype('uint8'), mask=gaps, transform=tr):
        p = shape(geom)
        if p.area / 1e6 >= args.minAreaKm2:
            polys.append(p)
    polys.sort(key=lambda p: -p.area)
    lines.append(f'{len(polys)} gap(s) of at least {args.minAreaKm2:g} km2')
    lines.append('')
    lines.append(f'{"area km2":>12}  {"lon":>8} {"lat":>7}  {"tile":<10} selected')
    feats = []
    for p in polys:
        c = p.centroid
        name, sel = whichTile(c.x, c.y, tiles)
        ll = toLonLat(p, crs)
        cl = ll.centroid
        if len(feats) < args.top:
            lines.append(f'{p.area / 1e6:12,.0f}  {cl.x:8.2f} {cl.y:7.2f}  {name:<10} {sel}')
        feats.append(dict(type='Feature', geometry=mapping(ll),
                          properties=dict(areaKm2=round(p.area / 1e6, 1), tile=name,
                                          selected=sel, lon=round(cl.x, 3),
                                          lat=round(cl.y, 3))))
    if len(polys) > args.top:
        lines.append(f'... and {len(polys) - args.top} smaller')
    text = '\n'.join(lines)
    print(text)
    with open(f'{args.out}.txt', 'w') as fp:
        fp.write(text + '\n')
    with open(f'{args.out}.geojson', 'w') as fp:
        json.dump(dict(type='FeatureCollection', features=feats), fp)
    print(f'\nwrote {args.out}.txt and {args.out}.geojson ({len(feats)} gap polygons)')

    if args.figure:
        drawFigure(valid, gaps, tr, shape2d, crs, f'{args.out}.png')
        print(f'wrote {args.out}.png')


def drawFigure(valid, gaps, tr, shape2d, crs, out):
    ''' Grey where the mosaic has data, red where a gap is, on the projection's outline. '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .globeFigure import mapBoundary
    ny, nx = shape2d
    ext = [tr.c, tr.c + nx * tr.a, tr.f + ny * tr.e, tr.f]
    bx, by = mapBoundary(crs)
    fig, ax = plt.subplots(figsize=(18, 18 * np.ptp(by) / np.ptp(bx)))
    ax.fill(bx, by, facecolor='#1b4f72', edgecolor='none', zorder=1)
    rgba = np.zeros((ny, nx, 4))
    rgba[valid] = (0.85, 0.85, 0.85, 1.)
    rgba[gaps] = (0.9, 0.1, 0.1, 1.)
    ax.imshow(rgba, extent=ext, interpolation='nearest', zorder=2)
    ax.plot(bx, by, color='#0d2b3e', linewidth=0.8, zorder=3)
    ax.set_xlim(bx.min(), bx.max()), ax.set_ylim(by.min(), by.max())
    ax.set_aspect('equal'), ax.axis('off')
    ax.set_title('grey: in the product    red: covered by a granule but missing')
    fig.tight_layout(pad=0.2)
    fig.savefig(out, dpi=100)
    plt.close(fig)


if __name__ == '__main__':
    main()
