#!/usr/bin/env python3
'''
Area-weighted backscatter histograms of a runGeomosaicTiles run, globally and by region.

  backscatterHistograms /scratch/ianj/mosaics/cycle30ascending [--stride 8] [--processes 16]

Reads each tile's VRT (frequency B under A, as in the mosaic) at every --stride-th pixel, and
weights each sample by its ground area (lat/lon pixels shrink with cos(lat)), so the histograms are
of area, not pixel counts. Regions:
  global, land, ocean
  continents (Natural Earth admin-0 CONTINENT): Africa, Antarctica (+ ice shelves), Asia, Europe,
      North America (without Greenland), Oceania, South America
  Greenland
  sea ice: no ice-concentration mask is used; stand-ins are ocean south of --southIceLat (default
      55 S, near the September maximum) and ocean north of --northIceLat (default 60 N; the tiling
      stops at 77.5 N, so most of the Arctic pack is outside it). --iceMask FILE (polygons, lon/lat)
      replaces both with a real sea-ice extent.
Writes <work>/stats/: histograms.csv (0.1 dB bins, km2 per region), percentiles.csv (1-99 %),
histograms.png, and prints the percentile table.
'''
import argparse
import concurrent.futures
import json
import math
import os
import sys

import numpy as np
from osgeo import gdal
from rasterio import features
from rasterio.transform import Affine
from shapely.geometry import box, shape
from shapely.ops import transform, unary_union

gdal.UseExceptions()
NODATA = -3000
BINS = np.arange(-30.0, 10.001, 0.1)          # dB bin edges
PCTS = [1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99]
CONTINENTS = ['Africa', 'Antarctica', 'Asia', 'Europe', 'North America', 'Oceania', 'South America']


def regionGeometries(args):
    ''' {name: lon/lat geometry} of the land regions, plus land and the sea-ice zone. '''
    import cartopy.feature as cf
    import cartopy.io.shapereader as shpreader
    reader = shpreader.Reader(shpreader.natural_earth(resolution='50m', category='cultural',
                                                      name='admin_0_countries'))
    parts = {c: [] for c in CONTINENTS}
    greenland = []
    for rec in reader.records():
        a = rec.attributes
        if a.get('ADMIN') == 'Greenland':
            greenland.append(rec.geometry)
            continue
        c = a.get('CONTINENT')
        if c in parts:
            parts[c].append(rec.geometry)
    shelves = list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m').geometries())
    regions = {c: unary_union(g) for c, g in parts.items()}
    regions['Antarctica'] = unary_union([regions['Antarctica']] + shelves)
    regions['Greenland'] = unary_union(greenland)
    land = unary_union(list(cf.LAND.with_scale('50m').geometries()) +
                       list(cf.LAKES.with_scale('50m').geometries()) + shelves)
    if args.iceMask:
        feats = json.load(open(args.iceMask))['features']
        ice = unary_union([shape(f['geometry']) for f in feats]).difference(land)
        seaIce = {'sea ice (mask)': ice}
    else:
        seaIce = {f'ocean S of {args.southIceLat:g}S': box(-180, -90, 180, -args.southIceLat),
                  f'ocean N of {args.northIceLat:g}N': box(-180, args.northIceLat, 180, 90)}
    return regions, land, seaIce


def tileHist(job):
    ''' Area-weighted histograms of one tile: {region: counts per bin (km2)}. '''
    vrt, stride, geoms = job
    ds = gdal.Open(vrt)
    gt = ds.GetGeoTransform()
    nx, ny = ds.RasterXSize // stride, ds.RasterYSize // stride
    if nx < 1 or ny < 1:
        return {}
    a = ds.GetRasterBand(1).ReadAsArray(buf_xsize=nx, buf_ysize=ny, resample_alg=gdal.GRIORA_NearestNeighbour)
    valid = a != NODATA
    if not valid.any():
        return {}
    db = a / 100.
    gts = (gt[0], gt[1] * ds.RasterXSize / nx, 0., gt[3], 0., gt[5] * ds.RasterYSize / ny)
    tr = Affine.from_gdal(*gts)
    epsg = ds.GetSpatialRef().GetAuthorityCode(None)
    if epsg == '4326':
        lat = gts[3] + (np.arange(ny) + 0.5) * gts[5]
        area = np.repeat((abs(gts[1]) * 111.32 * np.cos(np.radians(lat)) * abs(gts[5]) * 111.32)[:, None], nx, 1)
    else:
        area = np.full((ny, nx), abs(gts[1] * gts[5]) / 1e6)
    out = {}

    def add(name, m):
        m = m & valid
        if m.any():
            h, _ = np.histogram(np.clip(db[m], BINS[0], BINS[-1] - 1e-6), BINS, weights=area[m])
            out[name] = out.get(name, 0) + h

    add('global', np.ones_like(valid))
    for name, g in geoms:
        if g.is_empty:
            continue
        m = features.rasterize([(g, 1)], out_shape=(ny, nx), transform=tr, dtype='uint8').astype(bool)
        add(name, m)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', help='runGeomosaicTiles run directory')
    ap.add_argument('--product', default='gamma0')
    ap.add_argument('--stride', type=int, default=8, help='sample every Nth pixel [8]')
    ap.add_argument('--processes', type=int, default=16)
    ap.add_argument('--southIceLat', type=float, default=55., help='sea-ice stand-in: ocean south of this [55]')
    ap.add_argument('--northIceLat', type=float, default=60., help='sea-ice stand-in: ocean north of this [60]')
    ap.add_argument('--iceMask', default=None, help='sea-ice extent polygons (GeoJSON, lon/lat)')
    args = ap.parse_args()
    import glob
    import pyproj
    vrtDir = f'{args.work}/vrt/{args.product}'
    vrts = [v for v in sorted(glob.glob(f'{vrtDir}/*.vrt'))
            if not os.path.basename(v).startswith(('band_', 'cap_', 'global'))]
    print(f'{len(vrts)} tiles; regions from Natural Earth', flush=True)
    regions, land, seaIce = regionGeometries(args)
    named = list(regions.items()) + [('land', land)] + list(seaIce.items())
    jobs = []
    for v in vrts:
        ds = gdal.Open(v)
        gt = ds.GetGeoTransform()
        x0, y1 = gt[0], gt[3]
        x1, y0 = x0 + ds.RasterXSize * gt[1], y1 + ds.RasterYSize * gt[5]
        epsg = int(ds.GetSpatialRef().GetAuthorityCode(None))
        if epsg == 4326:
            tb = box(x0, y0, x1, y1)
            geoms = [(n, g.intersection(tb)) for n, g in named]
            geoms.append(('ocean', tb.difference(land)))
            geoms = [(n, g.difference(land) if n in seaIce else g) for n, g in geoms]
        else:
            fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
            polar = box(-180, -90, 180, -60) if epsg == 3031 else box(-180, 60, 180, 90)
            tb = box(x0, y0, x1, y1)
            proj = lambda g: transform(fwd, g.intersection(polar).segmentize(0.1)).buffer(0).intersection(tb)
            landP = proj(land)
            geoms = [(n, proj(g)) for n, g in named]
            geoms.append(('ocean', tb.difference(landP)))
            geoms = [(n, g.difference(landP) if n in seaIce else g) for n, g in geoms]
        jobs.append((v, args.stride, geoms))
    total = {}
    with concurrent.futures.ProcessPoolExecutor(args.processes) as pool:
        for k, h in enumerate(pool.map(tileHist, jobs), 1):
            for n, c in h.items():
                total[n] = total.get(n, 0) + c
            if k % 100 == 0 or k == len(jobs):
                print(f'  {k}/{len(jobs)} tiles', flush=True)
    order = ['global', 'land', 'ocean'] + CONTINENTS + ['Greenland'] + list(seaIce)
    order = [n for n in order if n in total]
    out = f'{args.work}/stats'
    os.makedirs(out, exist_ok=True)
    centres = 0.5 * (BINS[:-1] + BINS[1:])
    with open(f'{out}/histograms.csv', 'w') as fp:
        fp.write('dB,' + ','.join(f'"{n} (km2)"' for n in order) + '\n')
        for i, c in enumerate(centres):
            fp.write(f'{c:.2f},' + ','.join(f'{total[n][i]:.1f}' for n in order) + '\n')
    rows = []
    for n in order:
        h = total[n]
        cdf = np.cumsum(h) / h.sum()
        p = [float(np.interp(q / 100., cdf, BINS[1:])) for q in PCTS]
        rows.append((n, h.sum(), p))
    with open(f'{out}/percentiles.csv', 'w') as fp:
        fp.write('region,area_km2,' + ','.join(f'p{q}' for q in PCTS) + '\n')
        for n, a, p in rows:
            fp.write(f'"{n}",{a:.0f},' + ','.join(f'{v:.2f}' for v in p) + '\n')
    print(f'{"region":<22}{"area M km2":>11}  ' + ' '.join(f'{"p" + str(q):>6}' for q in PCTS))
    for n, a, p in rows:
        print(f'{n:<22}{a / 1e6:>11.2f}  ' + ' '.join(f'{v:>6.1f}' for v in p))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ncol = 4
    nrow = math.ceil(len(order) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.8 * nrow), sharex=True, constrained_layout=True)
    for ax, (n, a, p) in zip(axs.flat, rows):
        h = total[n]
        ax.fill_between(centres, h / h.sum() / 0.1, step='mid', color='#1f6f8b', alpha=0.8, lw=0)
        for q, v in zip(PCTS, p):
            if q in (2, 50, 98):
                ax.axvline(v, color='#c77c1e', lw=1, ls='--' if q != 50 else '-')
        ax.axvspan(-24, -1, color='#999999', alpha=0.12, lw=0)
        ax.set_title(f'{n}  ({a / 1e6:.1f} M km2)\n2/50/98 %: {p[1]:.1f} / {p[5]:.1f} / {p[9]:.1f} dB', fontsize=9)
        ax.set_xlim(-30, 5)
    for ax in axs.flat[len(rows):]:
        ax.axis('off')
    for ax in axs[-1]:
        ax.set_xlabel('gamma0 (dB)')
    fig.suptitle(f'{os.path.basename(os.path.abspath(args.work))}: area-weighted backscatter (per dB). '
                 'Grey band: current stretch -24..-1 dB; orange: 2 / 50 / 98 %', fontsize=11)
    fig.savefig(f'{out}/histograms.png', dpi=110)
    print(f'{out}/histograms.png, histograms.csv, percentiles.csv')
    return 0


if __name__ == '__main__':
    sys.exit(main())
