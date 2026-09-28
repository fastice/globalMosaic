#!/usr/bin/env python3
'''
No-data areas of a run read from its progress/monitor picture (pale-blue pixels = planned, no data),
as a gaps.geojson for `globalGCOVTiles --keepFrom ... --fillGaps` -- when the run's tiles are not
at hand (e.g. only the PNG was copied off the machine). gapReport on the tiles is more exact.

  gapsFromPicture progress.png cycle30/ascending --out gaps.geojson [--allWater]

Only land, Antarctic ice shelves, lakes and the Caspian are kept (open ocean along coasts is not
acquired and would add many granules for nothing) unless --allWater. Needs the picture layout of
runGeomosaicTiles.Progress.draw (1800 x 1500 px, dpi 100).
'''
import argparse
import json
import math
import sys

import numpy as np
import pyproj
from PIL import Image
from shapely.geometry import Point, box, mapping, shape
from shapely.ops import transform, unary_union
from shapely.strtree import STRtree

PLANNED = (170, 200, 235)
MID = (50, 1750, 89, 656)                               # frame x0, x1, y0, y1 of the 60S-60N panel
POLAR = (((57, 851), 3413, 1), ((949, 1743), 3031, -1))   # frame x0, x1 of the polar panels
POLAR_Y = (691, 1485)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('picture')
    ap.add_argument('tileRun', help='the tiling the run used (tiles.geojson)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--allWater', action='store_true', help='keep gaps over open ocean too')
    args = ap.parse_args()
    im = np.asarray(Image.open(args.picture).convert('RGB')).astype(int)
    if im.shape[:2] != (1500, 1800):
        sys.exit(f'{args.picture}: {im.shape[1]} x {im.shape[0]} px, expected 1800 x 1500')
    blue = np.abs(im - np.array(PLANNED)).sum(2) <= 6
    tiles = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    pts, polarPS = [], []
    x0, x1, y0, y1 = MID
    ys, xs = np.nonzero(blue[y0 + 1:y1, x0 + 1:x1])
    dx, dy = 360 / (x1 - x0), 120 / (y1 - y0)
    for a, b in zip(-180 + (xs + 1.5) * dx, 60 - (ys + 1.5) * dy):
        if abs(b) < 60:
            pts.append((a, b, dy / 2, dx / 2))
    for (px0, px1), epsg, sign in POLAR:
        fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True)
        inv = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True)
        r = float(np.hypot(*fwd.transform(0., sign * 60.)))
        py0, py1 = POLAR_Y
        ys, xs = np.nonzero(blue[py0 + 1:py1, px0 + 1:px1])
        res = 2 * r / (px1 - px0)
        X, Y = -r + (xs + 1.5) * res, r - (ys + 1.5) * res
        lo, la = inv.transform(X, Y)
        for a, b, x, y in zip(lo, la, X, Y):
            if sign * b >= 60:
                h = res / 2 / 111320
                pts.append((a, b, h, min(h / max(math.cos(math.radians(b)), 0.02), 30)))
                polarPS.append((epsg, x, y, res / 2))
    land = None
    if not args.allWater:
        import cartopy.feature as cf
        land = unary_union(list(cf.LAND.with_scale('50m').geometries()) +
                           list(cf.LAKES.with_scale('50m').geometries()) +
                           list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m').geometries()) +
                           [box(46, 36, 55.5, 47.5)])                     # the Caspian
    geoms = [shape(t['geometry']) for t in tiles]
    names = [t['properties']['name'] for t in tiles]
    tree = STRtree(geoms)
    per = {}
    for a, b, hy, hx in pts:
        p = Point(a, b)
        for i in tree.query(p):
            if tiles[i]['properties']['kind'] == 'latlon' and geoms[i].contains(p):
                per.setdefault(names[i], []).append(box(a - hx, b - hy, a + hx, b + hy))
    feats = []
    for n, v in per.items():
        g = unary_union(v).intersection(geoms[names.index(n)])
        if land is not None:
            g = g.intersection(land.buffer(0.05))
        if not g.is_empty:
            feats.append({'type': 'Feature', 'geometry': mapping(g), 'properties': {'tile': n, 'epsg': 4326}})
    for t in tiles:
        if t['properties']['kind'] != 'ps':
            continue
        h, e = t['properties']['halfWidthM'], t['properties']['epsg']
        g = unary_union([box(x - s, y - s, x + s, y + s) for ep, x, y, s in polarPS
                         if ep == e and abs(x) <= h and abs(y) <= h])
        if land is not None and not g.is_empty:
            fwd = pyproj.Transformer.from_crs(4326, e, always_xy=True).transform
            polarLand = land.intersection(box(-180, -90, 180, -70) if e == 3031 else box(-180, 70, 180, 90))
            g = g.intersection(transform(fwd, polarLand.segmentize(0.25)).buffer(5000))
        if not g.is_empty:
            feats.append({'type': 'Feature', 'geometry': mapping(g), 'properties': {'tile': t['properties']['name'], 'epsg': e}})
    json.dump({'type': 'FeatureCollection', 'features': feats}, open(args.out, 'w'))
    print(f'{len(pts)} no-data pixels -> gaps in {len(feats)} tiles: {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
