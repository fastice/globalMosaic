#!/usr/bin/env python3
'''
Lift the ocean part of a Google Earth superoverlay to sea level, so sea ice does not fade out.

  python3 -m globalMosaic.oceanAltitude TREE [TREE ...] [--processes 16] [--undo]

gdal2tiles' overlays are draped on the terrain, and over the ocean Google Earth's terrain is the
sea floor: the water surface covers them more and more as you zoom in. Tiles with data over the
ocean are changed to sit at --altitude m above sea level (altitudeMode absolute):
  - no land in the tile: its KML only gains <altitude>/<altitudeMode>
  - land and ocean: split in two overlays, <y>.land.png (draped, as before) and <y>.sea.png (at sea
    level); the original <y>.png is kept and the change is undone with --undo
Land = Natural Earth 10 m land + lakes + Antarctic ice shelves + the Caspian (as the tiered
product's land). Tiles already done are skipped, so reruns are cheap; works on any tree with
<z>/<x>/<y>.kml (the tiered product's base/, land/, detail/ each).
'''
import argparse
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from osgeo import gdal

gdal.UseExceptions()
MARK = '<altitudeMode>absolute</altitudeMode>'
LAND = None


def landParts():
    import cartopy.feature as cf
    from shapely.geometry import box
    geoms = list(cf.LAND.with_scale('10m').geometries()) + list(cf.LAKES.with_scale('10m').geometries()) + \
        list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '10m').geometries()) + \
        [box(46, 36, 55.5, 47.5)]                                         # the Caspian
    parts = []
    for g in geoms:
        parts.extend(getattr(g, 'geoms', [g]))
    return parts


def initWorker():
    global LAND
    from shapely.strtree import STRtree
    parts = landParts()
    LAND = (STRtree(parts), parts)


def overlay(href, box, altitude):
    alt = f'\n      <altitude>{altitude:g}</altitude>\n      {MARK}' if altitude is not None else ''
    return (f'    <GroundOverlay>\n      <drawOrder>{box["drawOrder"]}</drawOrder>\n      <Icon>\n'
            f'        <href>{href}</href>\n      </Icon>{alt}\n      <LatLonBox>\n'
            + ''.join(f'        <{k}>{box[k]}</{k}>\n' for k in ('north', 'south', 'east', 'west'))
            + '      </LatLonBox>\n    </GroundOverlay>\n')


def writePng(path, grey, alpha):
    m = gdal.GetDriverByName('MEM').Create('', grey.shape[1], grey.shape[0], 2, gdal.GDT_Byte)
    m.GetRasterBand(1).WriteArray(grey)
    m.GetRasterBand(2).WriteArray(alpha)
    m.GetRasterBand(2).SetColorInterpretation(gdal.GCI_AlphaBand)
    gdal.GetDriverByName('PNG').CreateCopy(path, m)


def doTile(job):
    ''' One tile: 'sea', 'split', or '' (unchanged). '''
    kml, altitude = job
    from shapely import prepared
    from shapely.geometry import box as sbox
    s = open(kml).read()
    if MARK in s:
        return ''
    g = re.search(r'    <GroundOverlay>.*?</GroundOverlay>\n', s, re.S)
    if g is None:
        return ''
    b = {k: re.search(f'<{k}>(.*?)</{k}>', g.group(0)).group(1) for k in ('drawOrder', 'href', 'north', 'south', 'east', 'west')}
    n, so, e, w = (float(b[k]) for k in ('north', 'south', 'east', 'west'))
    png = os.path.join(os.path.dirname(kml), b['href'])
    if not os.path.exists(png):
        return ''
    tree, parts = LAND
    tb = sbox(w, max(so, -90), e, min(n, 90))
    near = [parts[i] for i in tree.query(tb)]
    near = [p for p in near if p.intersects(tb)]
    if any(prepared.prep(p).contains(tb) for p in near):
        return ''                                                     # all land
    ds = gdal.Open(png)
    alpha = ds.GetRasterBand(ds.RasterCount).ReadAsArray()
    if not (alpha > 0).any():
        return ''
    if near:
        from rasterio import features
        from rasterio.transform import from_bounds
        ny, nx = alpha.shape
        land = features.rasterize([(p, 1) for p in near], out_shape=(ny, nx),
                                  transform=from_bounds(w, so, e, n, nx, ny), dtype='uint8').astype(bool)
        if not ((alpha > 0) & ~land).any():
            return ''                                                 # data only on land
    else:
        land = None
    if land is None or not ((alpha > 0) & land).any():
        new = overlay(b['href'], b, altitude)                         # data only over the ocean
        kind = 'sea'
    else:
        grey = ds.GetRasterBand(1).ReadAsArray()
        stem = os.path.splitext(png)[0]
        writePng(f'{stem}.land.png', grey, np.where(land, alpha, 0).astype('u1'))
        writePng(f'{stem}.sea.png', grey, np.where(land, 0, alpha).astype('u1'))
        base = os.path.splitext(b['href'])[0]
        new = overlay(f'{base}.land.png', b, None) + overlay(f'{base}.sea.png', b, altitude)
        kind = 'split'
    tmp = kml + '.tmp'
    with open(tmp, 'w') as fp:
        fp.write(s[:g.start()] + new + s[g.end():])
    os.replace(tmp, kml)
    return kind


def undoTile(kml):
    ''' Back to gdal2tiles' single draped overlay of <y>.png. '''
    s = open(kml).read()
    if MARK not in s:
        return ''
    s = re.sub(r'\n      <altitude>.*?</altitude>\n      ' + re.escape(MARK), '', s)
    m = re.search(r'    <GroundOverlay>.*?</GroundOverlay>\n    <GroundOverlay>.*?</GroundOverlay>\n', s, re.S)
    if m:                                                             # split: keep the first, point it at <y>.png
        first = re.search(r'    <GroundOverlay>.*?</GroundOverlay>\n', m.group(0), re.S).group(0)
        s = s[:m.start()] + first.replace('.land.png', '.png') + s[m.end():]
        stem = os.path.splitext(kml)[0]
        for p in (f'{stem}.land.png', f'{stem}.sea.png'):
            if os.path.exists(p):
                os.remove(p)
    with open(kml, 'w') as fp:
        fp.write(s)
    return 'undone'


def liftOcean(tree, altitude=10., processes=16, undo=False):
    kmls = sorted(glob.glob(f'{tree}/[0-9]*/*/*.kml'))
    print(f'ocean at sea level: {len(kmls)} tile KMLs under {tree}', flush=True)
    counts = {}
    if not undo:
        landParts()             # download the coastline once here: workers fetching it at once corrupt it
    with ProcessPoolExecutor(processes, initializer=None if undo else initWorker) as pool:
        results = pool.map(undoTile, kmls, chunksize=256) if undo else \
            pool.map(doTile, [(k, altitude) for k in kmls], chunksize=256)
        for i, r in enumerate(results, 1):
            counts[r] = counts.get(r, 0) + 1
            if i % 20000 == 0:
                print(f'  {i}/{len(kmls)}', flush=True)
    print(f'at sea level: {counts.get("sea", 0)}; split land/sea: {counts.get("split", 0)}; '
          f'undone: {counts.get("undone", 0)}; unchanged: {counts.get("", 0)}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tree', nargs='+', help='superoverlay directories (doc.kml + <z>/<x>/<y>.kml)')
    ap.add_argument('--altitude', type=float, default=10., help='ocean overlays, m above sea level [10]')
    ap.add_argument('--processes', type=int, default=16)
    ap.add_argument('--undo', action='store_true', help='restore the draped gdal2tiles overlays')
    args = ap.parse_args()
    for t in args.tree:
        liftOcean(t, args.altitude, args.processes, args.undo)
    return 0


if __name__ == '__main__':
    sys.exit(main())
