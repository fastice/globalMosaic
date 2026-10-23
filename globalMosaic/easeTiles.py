#!/usr/bin/env python3
'''
Chunk NISAR GCOV granules into EASE-Grid 2.0 tiles -- the input lists for an equal-area tiled
mosaic, in the same form globalGCOVTiles writes for the lat/lon product.

  easeTiles --catalogue cycle30/catalogue.geojson --hemisphere south --out cycle30/easeSouth
  easeTiles --catalogue cycle30/catalogue.geojson --hemisphere north --tileKm 600 --out ease/north

Projection
  --hemisphere global takes any equal-area CRS through --crs:
    ESRI:54076   Wagner VII, the lens-shaped whole-Earth map. PROJ has its forward and NO
                 inverse, so geomosaic needs the GP_WAGNER7 support in grimpProj.c (and the
                 outlines here use this module's own wagner7Inverse).
    EPSG:6933    EASE-Grid 2.0 Global, cylindrical. Has an EPSG code, so it reaches geomosaic
                 as -epsg and needs nothing special of it.
    +proj=eqearth  Equal Earth. No EPSG code: passed as -wkt, which needs grimpSRSString's
                 GP_GENERIC branch. This is the default for --hemisphere global.

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
  --direction A or D restricts to one orbit direction, as globalGCOVTiles' flag does. Worth
  using: every candidate is reprojected and rasterised for every box, so the tiling cost is
  proportional to the catalogue size (measured 2.1-2.7 s per box against 11,386 granules).


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
import csv
import functools
import glob
import json
import math
import os

import numpy as np
import pyproj
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import box, mapping, shape
from shapely.ops import transform, unary_union

from .globalGCOVTiles import greedyCover, loadGranules, rasterTarget, writeTile

EPSG = {'north': 6931, 'south': 6932}
# A global equal-area alternative with no EPSG code, passed to geomosaic as -wkt. Equal Earth is
# pseudocylindrical: straight parallels, curved meridians, area preserved. PROJ has its inverse,
# so geomosaic takes it through GP_GENERIC with no new C code -- unlike Wagner VII, which has no
# inverse at all. The extent below is the full map, 34488 x 16786 km.
GLOBAL_CRS = '+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs'
# Wagner VII (ESRI:54076) constants, the same ones mosaicSource/common/grimpProj.c uses.
# PROJ has the forward but NO INVERSE, so the inverse is written out here in closed form --
# (x/KX R, y/KY R) are equatorial Lambert azimuthal equal-area coordinates, which Snyder
# 24-16 inverts directly. PROJ treats the CRS as spherical with R = the semi-major axis.
W7 = dict(sin65=0.90630778703664996, kx=2.66723, ky=1.24104, R=6378137.0)
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


@functools.lru_cache(maxsize=8)
def isWagner7(crs):
    ''' Asked once per tile outline, so it is cached and does NOT go through to_proj4():
    that is lossy and pyproj warns about it every single call, which buried the tiling log
    in thousands of identical warnings. The projection method name says it directly. '''
    try:
        op = pyproj.CRS.from_user_input(crs).coordinate_operation
        return op is not None and 'wagner' in op.method_name.lower()
    except Exception:
        return False


def wagner7Inverse(x, y, lon0=0.):
    ''' Closed-form Wagner VII inverse, the twin of wagner7Inverse() in grimpProj.c. PROJ
    cannot do this direction at all, so pyproj returns inf and every outline would be lost. '''
    u = np.asarray(x, float) / (W7['R'] * W7['kx'])
    v = np.asarray(y, float) / (W7['R'] * W7['ky'])
    rho = np.hypot(u, v)
    c = 2. * np.arcsin(np.clip(rho / 2., -1., 1.))
    with np.errstate(invalid='ignore', divide='ignore'):
        st = np.where(rho > 0, v * np.sin(c) / np.where(rho > 0, rho, 1.) / W7['sin65'], 0.)
        lam = np.where(rho > 0, np.arctan2(u * np.sin(c), rho * np.cos(c)), 0.)
    lon = lon0 + 3. * np.degrees(lam)
    lat = np.degrees(np.arcsin(np.clip(st, -1., 1.)))
    # Outside the lens the formula gives |lon| > 180. WRAPPING that (lon % 360) is wrong and
    # was a real bug: it sends a point just past the edge to the far side of the world, which
    # is how data came to be painted beyond the map boundary. grimpProj.c rejects such pixels
    # outright; here they are clamped to the boundary meridian so an outline stays a closed
    # ring on the edge of the map.
    return np.clip(lon, lon0 - 180., lon0 + 180.), lat


def globalExtent(crs, step=1.):
    ''' Half-width and half-height of the whole-world map, from the projection itself rather
    than a table: the graticule is projected and the extremes taken. Hardcoding Equal Earth's
    numbers put an EPSG:6933 grid 371 rows past that projection's south edge. '''
    fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
    lon, lat = np.meshgrid(np.arange(-180., 180. + step, step),
                           np.arange(-90., 90. + step, step))
    x, y = fwd(lon.ravel(), lat.ravel())
    x, y = np.asarray(x), np.asarray(y)
    good = np.isfinite(x) & np.isfinite(y)
    if not good.any():
        raise SystemExit(f'easeTiles: cannot project the graticule into {crs}')
    return float(np.abs(x[good]).max()), float(np.abs(y[good]).max())


def globalBoxes(tileM, extent):
    ''' (name, box) over the whole map rectangle. Tiles outside the projection's lens are
    dropped later, when no granule footprint lands in them. '''
    nx = int(math.ceil(extent[0] / tileM))
    ny = int(math.ceil(extent[1] / tileM))
    for j in range(-ny, ny):
        for i in range(-nx, nx):
            yield f'X{i:+03d}Y{j:+03d}', box(i * tileM, j * tileM, (i + 1) * tileM, (j + 1) * tileM)


def projectCandidates(crs, cands):
    ''' Project every footprint ONCE, and keep its projected bounds for a cheap per-tile
    reject. selectForBox used to reproject all of them for every box: O(boxes x granules)
    shapely transforms, measured at 1.0-1.3 s per box against 5,631 granules, which the fill
    catalogues would take to ~42,000. '''
    fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
    out = []
    for g in cands:
        try:
            geom = transform(fwd, shape(g['raw'])).buffer(0)
        except Exception:
            continue
        if geom.is_empty:
            continue
        out.append(dict(g, projGeom=geom))
    bounds = np.array([g['projGeom'].bounds for g in out]) if out else np.zeros((0, 4))
    return out, bounds


def selectForBox(epsg, tileBox, cands, marginKm, cellM, minNewFrac, forced=(), target=None,
                 bounds=None):
    ''' Greedy cover of one projected tile, rasterised in metres on its own grid. '''
    pad = marginKm * 1000.
    x0, y0, x1, y1 = tileBox.bounds
    padded = box(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
    nx = max(int(round((x1 - x0 + 2 * pad) / cellM)), 1)
    ny = max(int(round((y1 - y0 + 2 * pad) / cellM)), 1)
    tr = from_origin(x0 - pad, y1 + pad, (x1 - x0 + 2 * pad) / nx, (y1 - y0 + 2 * pad) / ny)
    if bounds is not None and len(bounds):
        # bounding boxes that cannot touch this tile, dropped before any geometry work
        hit = ((bounds[:, 0] <= x1 + pad) & (bounds[:, 2] >= x0 - pad) &
               (bounds[:, 1] <= y1 + pad) & (bounds[:, 3] >= y0 - pad))
        cands = [g for g, h in zip(cands, hit) if h]
    fwd = None if cands and 'projGeom' in cands[0] else \
        pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
    masks = []
    for g in cands:
        try:
            geom = g['projGeom'] if fwd is None else transform(fwd, shape(g['raw'])).buffer(0)
            part = geom.intersection(padded)
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
    inv = (wagner7Inverse if isWagner7(epsg)
           else pyproj.Transformer.from_crs(epsg, 4326, always_xy=True).transform)
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
    ap.add_argument('--hemisphere', choices=['north', 'south', 'global'], required=True,
                    help="north/south: EASE-Grid 2.0 polar (EPSG 6931/6932). global: --crs over "
                         "the whole map, written as proj.wkt and passed to geomosaic with -wkt")
    ap.add_argument('--crs', default=GLOBAL_CRS,
                    help='--hemisphere global only: any CRS PROJ can invert [Equal Earth]')
    ap.add_argument('--out', required=True, help='output directory')
    ap.add_argument('--tileKm', type=float, default=600.,
                    help='tile size, km (600 ~ 5.4 deg, like the 6x6 deg lat/lon tiles) [600]')
    ap.add_argument('--extentKm', type=float, default=9000.,
                    help='half-width of the grid about the pole, km; the equator is at 9010 [9000]')
    ap.add_argument('--marginKm', type=float, default=10., help='feather margin padding [10]')
    ap.add_argument('--cellM', type=float, default=2000., help='selection raster cell, m [2000]')
    ap.add_argument('--minNewFrac', type=float, default=0.001,
                    help='stop adding granules below this new-cell fraction [0.001]')
    ap.add_argument('--keepFrom', default=None,
                    help='an earlier tiling directory: every tile keeps exactly the granules it '
                         'already had and only gains more, so re-mosaicking redoes just the '
                         'tiles that changed')
    ap.add_argument('--fillGaps', default=None,
                    help='a gaps GeoJSON in LAT/LON (mosaicGaps writes one): with --keepFrom, '
                         'tiles gain granules only over those areas. Lat/lon so that only the '
                         'forward transform is needed -- Wagner VII has no inverse in PROJ')
    ap.add_argument('--fillCatalogues', nargs='+', default=[],
                    help='catalogues of other cycles, used in the order given, only where the '
                         'main catalogue leaves a tile uncovered (as globalGCOVTiles)')
    ap.add_argument('--direction', choices=['A', 'D', 'both'], default='both',
                    help='orbit direction from the file name [both]')
    ap.add_argument('--allowV', action='store_true', help='also use V-only granules')
    args = ap.parse_args()

    isGlobal = args.hemisphere == 'global'
    crs = args.crs if isGlobal else f'EPSG:{EPSG[args.hemisphere]}'
    tileM, extentM = args.tileKm * 1000., args.extentKm * 1000.
    os.makedirs(f'{args.out}/tiles', exist_ok=True)
    wktName, epsg = None, 0
    if isGlobal:
        import pyproj
        c = pyproj.CRS.from_user_input(crs)
        code = c.to_authority()
        # Prefer the EPSG code when the CRS has one (EPSG:6933): geomosaic then takes -epsg, which
        # every build understands. Only a CRS with no code (Equal Earth) needs -wkt, and that
        # needs grimpSRSString's GP_GENERIC branch and gcovMosaic's outputSRS.
        if code is not None and code[0] == 'EPSG':
            epsg = int(code[1])
            print(f'{crs} -> EPSG:{epsg}, passed to geomosaic with -epsg', flush=True)
        else:
            wktName = 'proj.wkt'
            with open(f'{args.out}/{wktName}', 'w') as fp:
                fp.write(c.to_wkt() + '\n')
            print(f'{crs} has no EPSG code: wrote {args.out}/{wktName}, passed with -wkt\n'
                  f'  (needs a geomosaic with the GP_GENERIC SRS fixes)', flush=True)
    else:
        epsg = EPSG[args.hemisphere]
    cands = loadGranules(args.catalogue, allowV=args.allowV)
    for k, fc in enumerate(args.fillCatalogues, 1):
        cands += loadGranules(fc, allowV=args.allowV, fill=k)
    if args.fillCatalogues:
        print(f'{len(cands)} granules including {len(args.fillCatalogues)} fill catalogue(s)',
              flush=True)
    if args.direction != 'both':
        n = len(cands)
        cands = [g for g in cands if g['direction'] == args.direction]
        print(f'direction {args.direction}: {len(cands)} of {n} granules', flush=True)
    if args.validFootprints:
        # select on where a granule HAS data, not merely where it was acquired
        vf = {f['properties']['name']: f for path in args.validFootprints
              for f in json.load(open(path))['features']}
        keep, dropped, trimmed = [], 0, 0
        for g in cands:
            v = vf.get(g['name'])
            if v is not None:
                if v['geometry'] is None or v['properties'].get('validFrac', 1.) < args.minValid:
                    dropped += 1
                    continue
                g = dict(g, raw=v['geometry'])
                trimmed += 1
            keep.append(g)
        cands = keep
        print(f'valid footprints: {trimmed} granule(s) trimmed to their valid area, '
              f'{dropped} dropped as empty, {len(cands)} left', flush=True)
    if not isGlobal:
        # keep only granules in this hemisphere (a footprint straddling the equator stays)
        sign = 1 if args.hemisphere == 'north' else -1
        cands = [g for g in cands
                 if sign * shape(g['raw']).bounds[3] > 0 or sign * shape(g['raw']).bounds[1] > 0]
        print(f'{len(cands)} candidate granules in the {args.hemisphere}ern hemisphere', flush=True)
    else:
        print(f'{len(cands)} candidate granules, global', flush=True)

    cands, projBounds = projectCandidates(crs, cands)
    byName = {g['name']: g for g in cands}

    # --keepFrom: the granules each tile already had, forced into its new selection
    keep = {}
    if args.keepFrom:
        for path in glob.glob(f'{args.keepFrom}/tiles/*.csv'):
            name = os.path.basename(path)[:-4]
            keep[name] = [r['name'] for r in csv.DictReader(open(path))]
        print(f'--keepFrom: {len(keep)} existing tile(s)', flush=True)

    # --fillGaps: only these areas are open for new granules
    gapGeom = None
    if args.fillGaps:
        fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
        parts = []
        for f in json.load(open(args.fillGaps))['features']:
            try:
                g = transform(fwd, shape(f['geometry'])).buffer(0)
            except Exception:
                continue
            if not g.is_empty:
                parts.append(g)
        if not parts:
            raise SystemExit(f'easeTiles: no usable polygons in {args.fillGaps}')
        gapGeom = unary_union(parts)
        print(f'--fillGaps: {len(parts)} gap polygon(s), {gapGeom.area / 1e6:,.0f} km2',
              flush=True)

    feats, used, nTiles = [], set(), 0
    if isGlobal:
        ex = globalExtent(crs)
        print(f'{crs}: map half-extent {ex[0] / 1000.:.0f} x {ex[1] / 1000.:.0f} km', flush=True)
        boxes = list(globalBoxes(tileM, ex))
    else:
        boxes = list(tileBoxes(extentM, tileM))
    print(f'{len(boxes)} candidate tiles ({args.tileKm:g} km)' +
          ('' if isGlobal else f', out to {args.extentKm:g} km'), flush=True)
    for k, (name, b) in enumerate(boxes, 1):
        forced = [byName[n] for n in keep.get(name, []) if n in byName]
        picks, cover, nAvail = selectForBox(crs, b, cands, args.marginKm, args.cellM,
                                            args.minNewFrac, forced=forced, target=gapGeom,
                                            bounds=projBounds)
        if not picks:
            continue
        nTiles += 1
        writeTile(args.out, name, picks, used)
        x0, y0, x1, y1 = b.bounds
        bw = {}
        for g, _ in picks:
            bw[str(g['bw'])] = bw.get(str(g['bw']), 0) + 1
        feats.append(dict(type='Feature', geometry=mapping(outline(crs, b)),
                          properties=dict(name=name, kind='ease', epsg=epsg, wkt=wktName,
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
