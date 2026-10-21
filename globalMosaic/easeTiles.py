#!/usr/bin/env python3
'''
Chunk NISAR GCOV granules into EASE-Grid 2.0 tiles -- the input lists for an equal-area tiled
mosaic, in the same form globalGCOVTiles writes for the lat/lon product.

  easeTiles --catalogue cycle30/catalogue.geojson --hemisphere south --out cycle30/easeSouth
  easeTiles --catalogue cycle30/catalogue.geojson --hemisphere north --tileKm 600 --out ease/north

Projection
  EPSG:6931 (north) / EPSG:6932 (south): Lambert azimuthal equal-area on WGS84, centred on the
  pole, zero false easting/northing. geomosaic accepts them -- grimpProjFromSRS classifies them
  as projected/non-conformal, which is usable for positions (a backscatter mosaic); only vx/vy
  decomposition needs a conformal grid. Each hemisphere's grid reaches the equator, so north and
  south together cover the globe and meet there.

  LAEA is NOT conformal: shape and the direction of local scale vary away from the pole, while
  area is preserved everywhere. That is the point -- use it for area-integrated quantities, not
  where shape matters.

Tiles
  A square grid of --tileKm tiles about the pole, out to --extentKm (default 9000 km, just past
  the equator at 9010 km). Tiles that miss the hemisphere disc, or that no granule touches, are
  dropped. --tileKm 600 is about 5.4 deg of latitude, close to the 6 x 6 deg tiles that work well
  in the lat/lon product.

Selection (per tile)
  As globalGCOVTiles, but rasterised in metres on the tile's own projected grid rather than in
  degrees: the tile padded by --marginKm is rasterised at --cellM, candidates are ranked in tiers
  by bandwidth (77 -> 40 -> 20 -> 5 MHz), and a greedy set cover adds the granule covering the
  most still-uncovered cells. No lat/lon distortion enters, which matters at high latitude.

Outputs (--out DIR)
  tiles.geojson         one feature per tile: the tile's lat/lon outline, with kind 'ease', the
                        epsg and the projected bounds xMin/yMin/xMax/yMax in metres
  tiles/<name>.csv      selected granules (same columns as the lat/lon tiler)
  tiles/<name>.*.yaml   geomosaic -gcov yaml(s), one per frequency+polarization

  Then mosaic it with runGeomosaicTiles, which reads kind 'ease' and builds the job grid at
  --psResM spacing on the tile's projected bounds.
'''
import argparse
import json
import math
import os

import pyproj
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import box, mapping, shape
from shapely.ops import transform

from .globalGCOVTiles import greedyCover, loadGranules, rasterTarget, writeTile

EPSG = {'north': 6931, 'south': 6932}
# LAEA radius of the equator: R * sqrt(2) on the authalic sphere, 9009.96 km on WGS84
EQUATOR_M = 9009964.76


def tileBoxes(extentM, tileM):
    ''' (name, box) for every tile of the square grid about the pole that meets the hemisphere. '''
    n = int(math.ceil(extentM / tileM))
    disc = EQUATOR_M
    for j in range(-n, n):
        for i in range(-n, n):
            x0, y0 = i * tileM, j * tileM
            b = box(x0, y0, x0 + tileM, y0 + tileM)
            # nearest corner inside the equator disc: a tile wholly beyond it has no hemisphere
            near = min(math.hypot(x, y) for x, y in b.exterior.coords)
            if near > disc:
                continue
            yield f'X{i:+03d}Y{j:+03d}', b


def selectForBox(epsg, tileBox, cands, marginKm, cellM, minNewFrac, forced=(), target=None):
    ''' Greedy cover of one projected tile, rasterised in metres on its own grid. '''
    fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
    pad = marginKm * 1000.
    x0, y0, x1, y1 = tileBox.bounds
    padded = box(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
    nx = max(int(round((x1 - x0 + 2 * pad) / cellM)), 1)
    ny = max(int(round((y1 - y0 + 2 * pad) / cellM)), 1)
    tr = from_origin(x0 - pad, y1 + pad, (x1 - x0 + 2 * pad) / nx, (y1 - y0 + 2 * pad) / ny)
    masks = []
    for g in cands:
        try:
            part = transform(fwd, shape(g['raw'])).buffer(0).intersection(padded)
        except Exception:
            continue
        if part.is_empty:
            continue
        m = features.rasterize([(part, 1)], out_shape=(ny, nx), transform=tr, dtype='uint8',
                               all_touched=True).astype(bool)
        if m.any():
            masks.append((g, m))
    if not masks:
        return [], 0., 0
    # greedyCover -> (picks, coverage, coveredCells, totalCells)
    picks, cover, _, _ = greedyCover(masks, (ny, nx), minNewFrac, forced,
                                     rasterTarget(target, (ny, nx), tr))
    return picks, cover, len(masks)


def outline(epsg, tileBox, nSeg=12):
    ''' The tile's lat/lon outline, densified so the LAEA curvature shows, clipped to its
    hemisphere (points beyond the equator have no meaning for this grid). '''
    inv = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True).transform
    try:
        ll = transform(inv, tileBox.segmentize(min(tileBox.bounds[2] - tileBox.bounds[0],
                                                   tileBox.bounds[3] - tileBox.bounds[1]) / nSeg))
        return ll.buffer(0)
    except Exception:
        return tileBox


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalogue', required=True, help='granule catalogue.geojson')
    ap.add_argument('--hemisphere', choices=['north', 'south'], required=True)
    ap.add_argument('--out', required=True, help='output directory')
    ap.add_argument('--tileKm', type=float, default=600.,
                    help='tile size, km (600 ~ 5.4 deg, like the 6x6 deg lat/lon tiles) [600]')
    ap.add_argument('--extentKm', type=float, default=9000.,
                    help='half-width of the grid about the pole, km; the equator is at 9010 [9000]')
    ap.add_argument('--marginKm', type=float, default=10., help='feather margin padding [10]')
    ap.add_argument('--cellM', type=float, default=2000., help='selection raster cell, m [2000]')
    ap.add_argument('--minNewFrac', type=float, default=0.001,
                    help='stop adding granules below this new-cell fraction [0.001]')
    ap.add_argument('--allowV', action='store_true', help='also use V-only granules')
    args = ap.parse_args()

    epsg = EPSG[args.hemisphere]
    tileM, extentM = args.tileKm * 1000., args.extentKm * 1000.
    os.makedirs(f'{args.out}/tiles', exist_ok=True)
    cands = loadGranules(args.catalogue, allowV=args.allowV)
    # keep only granules in this hemisphere (a footprint straddling the equator stays)
    sign = 1 if args.hemisphere == 'north' else -1
    cands = [g for g in cands
             if sign * shape(g['raw']).bounds[3] > 0 or sign * shape(g['raw']).bounds[1] > 0]
    print(f'{len(cands)} candidate granules in the {args.hemisphere}ern hemisphere', flush=True)

    feats, used, nTiles = [], set(), 0
    boxes = list(tileBoxes(extentM, tileM))
    print(f'{len(boxes)} tiles meet the hemisphere ({args.tileKm:g} km, out to {args.extentKm:g} km)',
          flush=True)
    for k, (name, b) in enumerate(boxes, 1):
        picks, cover, nAvail = selectForBox(epsg, b, cands, args.marginKm, args.cellM,
                                            args.minNewFrac)
        if not picks:
            continue
        nTiles += 1
        writeTile(args.out, name, picks, used)
        x0, y0, x1, y1 = b.bounds
        bw = {}
        for g, _ in picks:
            bw[str(g['bw'])] = bw.get(str(g['bw']), 0) + 1
        feats.append(dict(type='Feature', geometry=mapping(outline(epsg, b)),
                          properties=dict(name=name, kind='ease', epsg=epsg,
                                          xMin=x0, yMin=y0, xMax=x1, yMax=y1,
                                          lonMult=None, height=None, width=None,
                                          available=nAvail,
                                          selected=len(picks), coverage=round(cover, 4), bw=bw)))
        if k % 25 == 0 or k == len(boxes):
            print(f'  {k}/{len(boxes)} boxes, {nTiles} tiles with data', flush=True)
    with open(f'{args.out}/tiles.geojson', 'w') as fp:
        json.dump(dict(type='FeatureCollection', features=feats), fp)
    print(f'{nTiles} tiles, {len(used)} granules -> {args.out}/tiles.geojson')


if __name__ == '__main__':
    main()
