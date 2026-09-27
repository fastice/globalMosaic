#!/usr/bin/env python3
'''
Chunk a list of NISAR GCOV granules into lat/lon tiles and pick, for each tile, a
non-redundant set of granules that covers it -- the input lists for a tiled global mosaic.

  globalGCOVTiles --cycle 030 --out cycle30/both            # search ASF, cache, tile the globe
  globalGCOVTiles --catalogue cycle30/catalogue.geojson --bbox -180 -60 180 77.5 --out test

Granules
  From an ASF search for one cycle (cached as catalogue.geojson), or from an existing
  catalogue. Only routine products (L2_PR) are used: urgent-response (L2_UR) products are always
  reprocessed as PR, so they are redundant and lower quality. Direction, mode and polarization
  come from the FILE NAME. Where a scene exists in more than one processing version, a P (production)
  CRID is used over an X (experimental) one; an X granule is kept only if no P version exists.
  Direction, mode and polarization (ASF's flightDirection is
  wrong for ~4% of granules). Only the latest product counter of each granule is kept.
  For an HH mosaic a granule is usable if frequency A carries H (SH/DH/QP) -- then A is used --
  or else frequency B does (5 MHz); V-only granules are dropped unless --allowV (then their V
  channel is used; yamls are split by frequency and polarization).

Tiles
  Lat/lon bands (--bands). Default: 6 x 6 deg to 60 deg; then 6 deg tall and 12 / 18 / 24 / 60 deg
  wide toward the poles, and one cap tile above 84 deg. All widths are multiples of 6 and divide
  360, so every tile edge lies on the 6 deg grid. Each band also carries a suggested longitude
  pixel-spacing multiplier (lonMult ~ 1/cos(band middle), integer) for the output grid.
  Tiles are made only where granules exist (no empty ocean tiles).

Selection (per tile)
  The tile, padded by --marginKm (the feather margin), is rasterised at --cellDeg. Candidates are
  ranked in tiers by the bandwidth of the channel used (77 -> 40 -> 20 -> 5 MHz). Within a tier a
  greedy set cover adds the granule covering the most still-uncovered cells, until nothing new is
  gained or a granule would add fewer than --minNewFrac of the tile; lower tiers only fill what is
  left. Mixed-mode (_M_) granules rank after non-mixed ones in the same tier.
  --exclude5MHzOcean drops 5 MHz granules whose footprint is < --landFrac land and lies within
  --oceanLat of the equator (mid-latitude ocean), keeping 5 MHz over land (e.g. the Sahara).

Outputs (--out DIR)
  catalogue.geojson     cached granule list (when searched)
  tiles.geojson         one feature per tile: bounds, band, lonMult, granules available/selected,
                        coverage of the padded tile
  tiles/<name>.csv      selected granules: name, url, frequency, pol, bandwidth, direction, newCells
  tiles/<name>.*.yaml   geomosaic -gcov yaml(s), one per frequency+polarization present (geomosaic
                        takes one of each per yaml): <name>.AHH.yaml, <name>.BHH.yaml, <name>.BVV.yaml
  summary.txt           counts
'''
import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union
from shapely.strtree import STRtree

NAME = re.compile(r'NISAR_L2_PR_GCOV_(\d{3})_(\d{3})_([AD])_(\d{3})_(\d{2})(\d{2})_'
                  r'([A-Z]{2})([A-Z]{2})_([A-Z])_.*_(\d{3})$')
HPOLS = ('SH', 'DH', 'QP')
VPOLS = ('SV', 'DV')
# lat0, lat1, tile height, tile width, lonMult
# poleward of +-84 deg a lat/lon tile degenerates (E-W pixel size -> 0 at the pole), so each cap
# is one polar stereographic tile instead (EPSG 3031 south, 3413 north)
DEFAULT_BANDS = [(-84, -78, 6, 60, 6), (-78, -72, 6, 24, 4),
                 (-72, -66, 6, 18, 3), (-66, -60, 6, 12, 2), (-60, 60, 6, 6, 1),
                 (60, 66, 6, 12, 2), (66, 72, 6, 18, 3), (72, 78, 6, 24, 4),
                 (78, 84, 6, 60, 6)]
CAP_LAT = 84.
CAPS = [('S84PS3031', -1, 3031), ('N84PS3413', 1, 3413)]


def parseName(stem, allowV=False):
    ''' Direction, bandwidths, pols, mixed flag from a GCOV file name (no extension).
    Picks the H channel (A preferred); with allowV, a granule with no H channel falls back to
    its V channel (A preferred) instead of being dropped. '''
    m = NAME.match(stem)
    if m is None:
        return None
    cyc, trk, d, frame, bwA, bwB, polA, polB, flag, ver = m.groups()
    if polA in HPOLS:
        freq, bw, pol = 'A', int(bwA), polA
    elif polB in HPOLS:
        freq, bw, pol = 'B', int(bwB), polB
    elif allowV and polA in VPOLS:
        freq, bw, pol = 'A', int(bwA), polA
    elif allowV and polB in VPOLS:
        freq, bw, pol = 'B', int(bwB), polB
    else:
        return None                     # no usable channel
    return dict(cycle=cyc, track=trk, direction=d, frame=frame, freq=freq, bw=bw, pol=pol,
                mixed=(flag == 'M'), version=ver)


def sceneKey(stem):
    ''' Identity of an acquisition regardless of processing version: everything up to the stop time
    (cycle, track, direction, frame, mode, pol, flag, start, stop). '''
    return '_'.join(stem.split('_')[:13])


def versionRank(stem):
    ''' Preference among copies of one scene: a P (production) CRID over an X (experimental) one --
    an X is kept only when no P exists -- then the higher CRID, then the higher product counter. '''
    parts = stem.split('_')
    crid = parts[13]
    return (crid[0] == 'P', int(crid[1:]) if crid[1:].isdigit() else 0, parts[-1])


def unwrap(geom):
    ''' Split a footprint that crosses the antimeridian into pieces within -180..180. '''
    xs = [x for poly in getattr(geom, 'geoms', [geom]) for x, _ in poly.exterior.coords]
    if max(xs) - min(xs) <= 180:
        return geom
    from shapely.affinity import translate
    from shapely.geometry import Polygon
    shifted = Polygon([(x + 360 if x < 0 else x, y) for x, y in geom.exterior.coords])
    west = shifted.intersection(box(0, -90, 180, 90))
    east = translate(shifted.intersection(box(180, -90, 540, 90)), xoff=-360)
    return unary_union([west, east])


def searchCatalogue(cycle, outFile):
    ''' All GCOVs of one cycle (by file name), latest product counter each, with footprints. '''
    import asf_search as asf
    asf.constants.INTERNAL.CMR_TIMEOUT = 600
    from datetime import datetime, timedelta
    # cycle 030 began 2026-09-06; cycles are 12 days and spill past the nominal end (030 has
    # acquisitions on 09-18), so search 3 days either side and keep the cycle by file name
    start = datetime(2026, 9, 6) + timedelta(days=12 * (int(cycle) - 30) - 3)
    res = asf.search(dataset=asf.DATASET.NISAR, processingLevel='GCOV',
                     start=start.isoformat(), end=(start + timedelta(days=18)).isoformat(),
                     maxResults=1000000)
    best = {}
    for r in res:
        stem = os.path.basename(r.properties.get('url', ''))[:-3]
        m = NAME.match(stem) if stem else None
        if m is None or m.group(1) != cycle:
            continue
        # one entry per scene: P over X processing version, then highest CRID and counter
        key, rank = sceneKey(stem), versionRank(stem)
        if key not in best or rank > best[key][3]:
            best[key] = (stem, r.properties['url'], r.geometry, rank)
    feats = [{'type': 'Feature', 'geometry': g, 'properties': {'name': s, 'url': u}}
             for s, u, g, _ in best.values()]
    with open(outFile, 'w') as fp:
        # some asf_search geometries hold shapely coordinate sequences: write them as lists
        json.dump({'type': 'FeatureCollection', 'features': feats}, fp,
                  default=lambda o: [list(c) if hasattr(c, '__iter__') else c for c in o])
    return outFile


def loadGranules(catalogue, allowV=False, fill=0):
    grans = []
    feats = json.load(open(catalogue))['features']
    # keep one copy per scene (P over X, then highest version) -- also for catalogues written
    # before scenes were keyed without the version
    best = {}
    for f in feats:
        name = f['properties']['name']
        k = sceneKey(name)
        if k not in best or versionRank(name) > versionRank(best[k]['properties']['name']):
            best[k] = f
    for f in best.values():
        name = f['properties']['name']
        p = parseName(name, allowV)
        if p is None:
            continue
        p.update(name=name, url=f['properties'].get('url', ''), raw=f['geometry'],
                 geom=unwrap(shape(f['geometry']).buffer(0)), fill=fill)
        grans.append(p)
    return grans


def tilesFor(bands, bbox):
    ''' (name, box, band) for every tile of the band table inside bbox. '''
    x0, y0, x1, y1 = bbox
    for lat0, lat1, h, w, mult in bands:
        for lat in np.arange(lat0, lat1, h):
            top = min(lat + h, lat1)
            if top <= y0 or lat >= y1:
                continue
            for lon in np.arange(-180, 180, w):
                if lon + w <= x0 or lon >= x1:
                    continue
                s, n = max(lat, y0), min(top, y1)
                name = f'{"N" if lat >= 0 else "S"}{abs(int(lat)):02d}' \
                       f'{"E" if lon >= 0 else "W"}{abs(int(lon)):03d}_{int(top - lat)}x{int(w)}'
                yield name, box(lon, s, lon + w, n), dict(lonMult=mult, height=h, width=w)


def selectForTile(tileBox, cands, marginKm, cellDeg, minNewFrac):
    ''' Greedy tiered set cover on a raster of the padded tile. Returns (picks, coverage). '''
    lat = tileBox.centroid.y
    dLat = marginKm / 111.32
    dLon = min(marginKm / (111.32 * max(math.cos(math.radians(lat)), 0.02)), 30)
    x0, y0, x1, y1 = tileBox.bounds
    padded = box(x0 - dLon, max(y0 - dLat, -90), x1 + dLon, min(y1 + dLat, 90))
    px0, py0, px1, py1 = padded.bounds
    nx, ny = max(int(round((px1 - px0) / cellDeg)), 1), max(int(round((py1 - py0) / cellDeg)), 1)
    tr = from_origin(px0, py1, (px1 - px0) / nx, (py1 - py0) / ny)
    masks = []
    for g in cands:
        part = g['geom'].intersection(padded)
        if part.is_empty:
            continue
        m = features.rasterize([(part, 1)], out_shape=(ny, nx), transform=tr, dtype='uint8',
                               all_touched=True).astype(bool)
        if m.any():
            masks.append((g, m))
    return greedyCover(masks, (ny, nx), minNewFrac)


def greedyCover(masks, shape2d, minNewFrac):
    ''' Tiered greedy set cover over rasterised footprints [(granule, bool mask)]. Tiers: the main
    catalogue (fill 0) through all bandwidths first, then each fill catalogue (--fillCatalogues)
    in turn, which therefore only covers what the main cycle leaves open. '''
    ny, nx = shape2d
    union = np.zeros((ny, nx), bool)
    for _, m in masks:
        union |= m
    covered = np.zeros((ny, nx), bool)
    picks = []
    minNew = max(int(minNewFrac * nx * ny), 1)
    for fill, bw in sorted({(g.get('fill', 0), -g['bw']) for g, _ in masks}):
        tier = sorted([(g, m) for g, m in masks if g.get('fill', 0) == fill and g['bw'] == -bw],
                      key=lambda gm: gm[0]['mixed'])
        while tier:
            gains = [(int((m & ~covered).sum()) - (0.5 if g['mixed'] else 0), i)
                     for i, (g, m) in enumerate(tier)]
            gain, i = max(gains)
            if gain < minNew:
                break
            g, m = tier.pop(i)
            picks.append((g, int((m & ~covered).sum())))
            covered |= m
    coverage = covered.sum() / (nx * ny)
    return picks, coverage, int(union.sum()), nx * ny


def writeTile(outDir, name, picks, used):
    ''' Selected granules of one tile: a CSV, and one geomosaic yaml per frequency/pol group. '''
    with open(f'{outDir}/tiles/{name}.csv', 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['name', 'url', 'freq', 'pol', 'bwMHz', 'direction', 'mixed', 'newCells'])
        for g, new in picks:
            w.writerow([g['name'], g['url'], g['freq'], g['pol'], g['bw'], g['direction'],
                        int(g['mixed']), new])
            used.add(g['name'])
    byFreq = defaultdict(list)
    for g, _ in picks:
        byFreq[(g['freq'], 'H' if g['pol'] in HPOLS else 'V')].append(g)
    for (freq, p), gs in byFreq.items():
        with open(f'{outDir}/tiles/{name}.{freq}{p}{p}.yaml', 'w') as fp:
            fp.write(f'# {name}: {len(gs)} granules, frequency {freq}\n'
                     f'polarization: {p}{p}\nfrequency: {freq}\nuseMask: true\nfiles:\n')
            for g in gs:
                fp.write(f'  - {g["name"]}.h5\n')


def selectForCap(sign, epsg, cands, marginKm, cellM, minNewFrac):
    ''' Greedy cover of the polar cap (|lat| > CAP_LAT) on a polar stereographic raster, from the
    footprints' own corners (no lat/lon distortion near the pole). '''
    import pyproj
    from shapely.ops import transform
    fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
    half = math.hypot(*fwd(0., sign * CAP_LAT)) + marginKm * 1000.
    n = max(int(round(2 * half / cellM)), 1)
    tr = from_origin(-half, half, 2 * half / n, 2 * half / n)
    masks = []
    for g in cands:
        part = transform(fwd, shape(g['raw'])).buffer(0).intersection(box(-half, -half, half, half))
        if part.is_empty:
            continue
        m = features.rasterize([(part, 1)], out_shape=(n, n), transform=tr, dtype='uint8',
                               all_touched=True).astype(bool)
        if m.any():
            masks.append((g, m))
    return greedyCover(masks, (n, n), minNewFrac), half - marginKm * 1000.


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--cycle', help='search ASF for this cycle (e.g. 030) and cache the catalogue')
    src.add_argument('--catalogue', help='existing catalogue.geojson')
    ap.add_argument('--out', required=True, help='output directory')
    ap.add_argument('--bbox', type=float, nargs=4, default=[-180, -90, 180, 77.5],
                    metavar=('W', 'S', 'E', 'N'), help='area [globe south of 77.5N]')
    ap.add_argument('--bands', help='JSON list of [lat0, lat1, height, width, lonMult] '
                    '[built-in 6-deg table]')
    ap.add_argument('--marginKm', type=float, default=10., help='feather margin padding [10]')
    ap.add_argument('--cellDeg', type=float, default=0.02, help='selection raster cell [0.02]')
    ap.add_argument('--minNewFrac', type=float, default=0.001,
                    help='skip a granule adding less than this fraction of the tile [0.001]')
    ap.add_argument('--exclude5MHzOcean', action='store_true',
                    help='drop 5 MHz granules over mid-latitude ocean (keep over land)')
    ap.add_argument('--landFrac', type=float, default=0.05,
                    help='with --exclude5MHzOcean: keep if at least this land fraction [0.05]')
    ap.add_argument('--oceanLat', type=float, default=60.,
                    help='with --exclude5MHzOcean: mid-latitude limit |lat| [60]')
    ap.add_argument('--direction', choices=['A', 'D', 'both'], default='both',
                    help='orbit direction by file name [both]')
    ap.add_argument('--allowV', action='store_true',
                    help='use the V channel of granules that have no H channel (e.g. 5 MHz SV-only '
                    'modes over land); the channel used is in each tile CSV')
    ap.add_argument('--noMixed', action='store_true', help='drop mixed-mode (_M_) granules')
    ap.add_argument('--fillCatalogues', nargs='+', default=[],
                    help='catalogues of other cycles (e.g. cycle30/catalogue029.geojson), used in the '
                    'order given, only where the main catalogue leaves a tile uncovered')
    ap.add_argument('--validFootprints', nargs='+', default=None,
                    help='scanMasks output: use each scanned granule\'s valid-data footprint (mask '
                    '1..254) instead of its acquisition footprint; drop it when too little is valid')
    ap.add_argument('--minValid', type=float, default=0.05,
                    help='with --validFootprints: drop granules with less valid data than this [0.05]')
    args = ap.parse_args()

    os.makedirs(f'{args.out}/tiles', exist_ok=True)
    catalogue = args.catalogue or searchCatalogue(args.cycle, f'{args.out}/catalogue.geojson')
    nCat = len(json.load(open(catalogue))['features'])
    grans = loadGranules(catalogue, args.allowV)
    for k, fc in enumerate(args.fillCatalogues, 1):
        grans += loadGranules(fc, args.allowV, fill=k)
    nAll = len(grans)
    dropped = Counter()
    if args.direction != 'both':
        dropped['direction'] = sum(g['direction'] != args.direction for g in grans)
        grans = [g for g in grans if g['direction'] == args.direction]
    if args.noMixed:
        dropped['mixed'] = sum(g['mixed'] for g in grans)
        grans = [g for g in grans if not g['mixed']]
    if args.validFootprints:
        # partially focused (mask 0) samples are dropped by geomosaic: select on what is valid
        vf = {f['properties']['name']: f for path in args.validFootprints
              for f in json.load(open(path))['features']}
        keep = []
        for g in grans:
            v = vf.get(g['name'])
            if v is not None:
                if v['geometry'] is None or v['properties']['validFrac'] < args.minValid:
                    dropped['no valid data (mask)'] += 1
                    continue
                g['geom'] = unwrap(shape(v['geometry']).buffer(0))
                g['raw'] = v['geometry']
            keep.append(g)
        grans = keep
    if args.exclude5MHzOcean:
        import cartopy.feature as cf
        land = unary_union(list(cf.LAND.with_scale('50m').geometries()))
        keep = []
        for g in grans:
            c = g['geom'].centroid
            if g['bw'] == 5 and abs(c.y) < args.oceanLat:
                frac = g['geom'].intersection(land).area / max(g['geom'].area, 1e-12)
                if frac < args.landFrac:
                    dropped['5MHz ocean'] += 1
                    continue
            keep.append(g)
        grans = keep
    tree = STRtree([g['geom'] for g in grans])
    bands = [tuple(b) for b in json.loads(args.bands)] if args.bands else DEFAULT_BANDS

    feats, stats = [], Counter()
    used = set()
    for name, tbox, band in tilesFor(bands, args.bbox):
        idx = tree.query(tbox, predicate='intersects')
        if len(idx) == 0:
            continue                    # no frames: no tile
        cands = [grans[i] for i in idx]
        picks, coverage, unionCells, nCells = selectForTile(tbox, cands, args.marginKm,
                                                            args.cellDeg, args.minNewFrac)
        if not picks:
            continue
        stats['tiles'] += 1
        stats['available'] += len(cands)
        stats['selected'] += len(picks)
        writeTile(args.out, name, picks, used)
        feats.append({'type': 'Feature', 'geometry': mapping(tbox), 'properties': dict(
            name=name, kind='latlon', epsg=4326, **band, available=len(cands), selected=len(picks),
            coverage=round(float(coverage), 4), unionFrac=round(unionCells / nCells, 4),
            bw=dict(Counter(str(g['bw']) for g, _ in picks)))})
    # polar caps as single polar stereographic tiles
    for name, sign, epsg in CAPS:
        if (sign < 0 and args.bbox[1] >= -CAP_LAT) or (sign > 0 and args.bbox[3] <= CAP_LAT):
            continue
        capBox = box(-180, -90, 180, -CAP_LAT) if sign < 0 else box(-180, CAP_LAT, 180, 90)
        idx = tree.query(capBox, predicate='intersects')
        if len(idx) == 0:
            continue
        cands = [grans[i] for i in idx]
        (picks, coverage, unionCells, nCells), halfM = selectForCap(
            sign, epsg, cands, args.marginKm, 2000., args.minNewFrac)
        if not picks:
            continue
        writeTile(args.out, name, picks, used)
        stats['tiles'] += 1
        stats['available'] += len(cands)
        stats['selected'] += len(picks)
        feats.append({'type': 'Feature', 'geometry': mapping(capBox), 'properties': dict(
            name=name, kind='ps', epsg=epsg, halfWidthM=round(halfM), lonMult=None, height=None,
            width=None, available=len(cands), selected=len(picks), coverage=round(float(coverage), 4),
            unionFrac=round(unionCells / nCells, 4),
            bw=dict(Counter(str(g['bw']) for g, _ in picks)))})
    with open(f'{args.out}/tiles.geojson', 'w') as fp:
        json.dump({'type': 'FeatureCollection', 'features': feats}, fp)
    lines = [f'catalogue {catalogue}: {nCat} granules, {nAll} usable '
             f'({"H, V fallback" if args.allowV else "H only"}); dropped {dict(dropped)}',
             f'tiles {stats["tiles"]}; granule-tile pairs available {stats["available"]}, '
             f'selected {stats["selected"]}; distinct granules used {len(used)} of {len(grans)}']
    open(f'{args.out}/summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    sys.exit(main())
