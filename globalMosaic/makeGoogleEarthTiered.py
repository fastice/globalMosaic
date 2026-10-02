#!/usr/bin/env python3
'''
Google Earth product of a runGeomosaicTiles run with a zoom limit that follows the data:

  makeGoogleEarthTiered /scratch/ianj/mosaics/cycle30ascending_z11 [--processes 16]

Three KML superoverlays under <out> (default <work>/googleEarth), one doc.kml loading all three;
the deeper layers draw over the shallower ones where they have data:

  base/    zoom 0-9   everything (land, ocean, sea ice) from the mosaic, as makeGoogleEarth
  land/    zoom 5-10  land, lakes and ice shelves only (ocean and sea ice stop at 9): each tile
                      averaged 2 x 2 (2.4" -> 4.8", about the zoom-10 pixel) and then --smooth N x N,
                      no-data left out of every average
  detail/  zoom 7-11  full resolution where the data or the science call for it:
                      - land covered by 40/77 MHz granules (10 m data), except the Antarctic and
                        Greenland ice-sheet interiors
                      - the ice-sheet margins: within --marginKm of the Antarctic and Greenland coasts
                      - glaciers and ice caps (Natural Earth glaciated areas)
                      20 MHz and 5 MHz land stops at 10.
land/ and detail/ carry images only at their finest zoom (10, 11); their coarser levels are link-only
(linkOnly: KMLs without overlays, no PNGs) -- the index down to the fine tiles, while base shows
those zooms. So base plus any combination of land and detail works; land or detail alone is blank
until its finest zoom.

Sources for land/ and detail/ are written once to <out>/src10, <out>/src11 (int16 dB x 100,
DEFLATE); gdal2tiles is resumable (-e) but the sources are rebuilt on every run.
'''
import argparse
import concurrent.futures
import csv
import glob
import json
import os
import sys

import numpy as np
import yaml
from osgeo import gdal
from rasterio import features
from rasterio.transform import Affine
from scipy.ndimage import uniform_filter
from shapely.geometry import box, mapping, shape
from shapely.ops import transform, unary_union

from .globalGCOVTiles import unwrap
from .makeGoogleEarth import superoverlay

gdal.UseExceptions()
NODATA = -3000
POLAR = {3031: box(-180, -90, 180, -60), 3413: box(-75, 59, -10, 84)}   # Antarctica, Greenland
PSRES = 1000.                                                           # ice-class raster, m


def landGeometry():
    import cartopy.feature as cf
    return unary_union(list(cf.LAND.with_scale('50m').geometries()) +
                       list(cf.LAKES.with_scale('50m').geometries()) +
                       list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m').geometries()) +
                       [box(46, 36, 55.5, 47.5)])                         # the Caspian


def glacierGeometry():
    import cartopy.feature as cf
    return unary_union(list(cf.NaturalEarthFeature('physical', 'glaciated_areas', '50m').geometries()))


def iceClassRaster(land, epsg, marginKm, path):
    ''' 1 km polar stereographic raster: 1 = ice-sheet margin (land within marginKm of the coast),
    2 = interior. Antarctica (3031) and Greenland (3413). '''
    import pyproj
    fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
    # lat/lon polygons reaching the pole leave a slit along 180 once projected; close it (and any
    # sliver seams) so only true coast counts as the edge
    ice = land.intersection(POLAR[epsg])
    if epsg == 3413:
        # Greenland only: the box also holds Iceland and Canada's Arctic islands
        ice = max(getattr(ice, 'geoms', [ice]), key=lambda g: g.area)
    ice = transform(fwd, ice.segmentize(0.1)).buffer(5000).buffer(-5000)
    # the margin is measured from the GROUNDED coast (grounding lines and rock coast), and the ice
    # shelves are margin too, so grounding zones get full resolution, not only the shelf fronts
    import cartopy.feature as cf
    shelves = unary_union(list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys',
                                                      '50m').geometries())).intersection(POLAR[epsg])
    shelves = transform(fwd, shelves.segmentize(0.1)).buffer(0) if not shelves.is_empty else box(0, 0, 0, 0)
    grounded = ice.difference(shelves.buffer(2000))
    margin = grounded.intersection(grounded.boundary.buffer(marginKm * 1000.)).union(shelves.intersection(ice))
    x0, y0, x1, y1 = ice.bounds
    nx, ny = int((x1 - x0) / PSRES) + 1, int((y1 - y0) / PSRES) + 1
    tr = Affine(PSRES, 0, x0, 0, -PSRES, y1)
    r = features.rasterize([(ice, 2)], out_shape=(ny, nx), transform=tr, dtype='uint8')
    m = features.rasterize([(margin, 1)], out_shape=(ny, nx), transform=tr, dtype='uint8')
    r[m == 1] = 1
    np.savez_compressed(path, r=r, gt=np.array(tr.to_gdal()), epsg=epsg)
    return path


def iceClass(paths, epsg, X, Y):
    ''' Ice class (0 none, 1 margin, 2 interior) at points X, Y given in `epsg` (4326 lon/lat or a
    polar stereographic CRS), from the rasters in `paths`. '''
    import pyproj
    out = np.zeros(X.shape, np.uint8)
    for p in paths:
        z = np.load(p)
        e = int(z['epsg'])
        gt = z['gt']
        r = z['r']
        if epsg == e:
            x, y = X, Y
        else:
            x, y = pyproj.Transformer.from_crs(epsg, e, always_xy=True).transform(X, Y)
        j = np.floor((x - gt[0]) / gt[1]).astype(int)
        i = np.floor((y - gt[3]) / gt[5]).astype(int)
        ok = (i >= 0) & (i < r.shape[0]) & (j >= 0) & (j < r.shape[1])
        out[ok] = np.maximum(out[ok], r[i[ok], j[ok]])
    return out


def blockAverage(a, valid, k):
    ''' k x k block mean ignoring no-data; returns (mean, any-valid mask). '''
    ny, nx = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    s = np.where(valid, a, 0.)[:ny, :nx].reshape(ny // k, k, nx // k, k).sum((1, 3))
    w = valid[:ny, :nx].reshape(ny // k, k, nx // k, k).sum((1, 3))
    return np.where(w > 0, s / np.maximum(w, 1), 0.), w > 0


def writeTif(path, a, gt, srs):
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(path, a.shape[1], a.shape[0], 1, gdal.GDT_Int16,
                    options=['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=2', 'BIGTIFF=IF_SAFER'])
    ds.SetGeoTransform(gt)
    ds.SetProjection(srs)
    b = ds.GetRasterBand(1)
    b.SetNoDataValue(NODATA)
    b.WriteArray(a)
    ds = None


def tileSources(job):
    ''' Level-10 (land, averaged, smoothed) and level-11 (detail mask, full resolution) sources
    of one tile. '''
    name, vrt, epsg, landG, detailG, icePaths, out10, out11, smooth = job
    ds = gdal.Open(vrt)
    gt = ds.GetGeoTransform()
    srs = ds.GetProjection()
    a = ds.GetRasterBand(1).ReadAsArray().astype('f4')
    valid = a != NODATA
    notes = []
    # ---- level 10: 2 x 2 block mean, then an N x N no-data-aware box smooth, land only
    m, v = blockAverage(a, valid, 2)
    if smooth > 1:
        s = uniform_filter(np.where(v, m, 0.), smooth, mode='constant')
        w = uniform_filter(v.astype('f4'), smooth, mode='constant')
        m = np.where(w > 0, s / np.maximum(w, 1e-6), 0.)
    gt10 = (gt[0], gt[1] * 2, 0., gt[3], 0., gt[5] * 2)
    land = features.rasterize([(landG, 1)], out_shape=m.shape, transform=Affine.from_gdal(*gt10),
                              dtype='uint8').astype(bool) if not landG.is_empty else np.zeros(m.shape, bool)
    keep = v & land
    if keep.any():
        writeTif(f'{out10}/{name}.tif', np.where(keep, np.round(m), NODATA).astype('i2'), gt10, srs)
        notes.append('land10')
    # ---- level 11: full resolution inside the detail mask
    det = features.rasterize([(detailG, 1)], out_shape=a.shape, transform=Affine.from_gdal(*gt),
                             dtype='uint8').astype(bool) if not detailG.is_empty else np.zeros(a.shape, bool)
    if icePaths:
        # ice-sheet classes on this grid: interior leaves the detail set, margin joins it (on land)
        ny, nx = a.shape
        step = 4                                   # sample the 1 km classes every 4 pixels
        jj = np.arange(0, nx, step) + 0.5
        ii = np.arange(0, ny, step) + 0.5
        X = gt[0] + jj[None, :] * gt[1] + 0 * ii[:, None]
        Y = gt[3] + ii[:, None] * gt[5] + 0 * jj[None, :]
        c = iceClass(icePaths, epsg, X, Y)
        c = np.repeat(np.repeat(c, step, 0), step, 1)[:ny, :nx]
        landFull = features.rasterize([(landG, 1)], out_shape=a.shape, transform=Affine.from_gdal(*gt),
                                      dtype='uint8').astype(bool) if not landG.is_empty else np.zeros(a.shape, bool)
        det = (det & (c != 2)) | ((c == 1) & landFull)
    det &= valid
    if det.any():
        rows, cols = np.nonzero(det.any(1))[0], np.nonzero(det.any(0))[0]
        r0, r1, c0, c1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
        sub = np.where(det[r0:r1, c0:c1], a[r0:r1, c0:c1], NODATA).astype('i2')
        writeTif(f'{out11}/{name}.tif', sub,
                 (gt[0] + c0 * gt[1], gt[1], 0., gt[3] + r0 * gt[5], 0., gt[5]), srs)
        notes.append(f'detail11 {det.mean():.0%}')
    return name, notes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', help='runGeomosaicTiles run directory')
    ap.add_argument('--out', default=None, help='output directory [<work>/googleEarth]')
    ap.add_argument('--product', default='gamma0')
    ap.add_argument('--smooth', type=int, default=3, help='box smooth of the zoom-10 source, pixels [3]')
    ap.add_argument('--marginKm', type=float, default=200., help='ice-sheet margin width, km [200]')
    ap.add_argument('--processes', type=int, default=16)
    ap.add_argument('--tileSize', type=int, default=256, choices=[256, 512],
                    help='tile size, px; 512 writes a quarter as many files, same detail (levels then '
                    '0-8 / 4-9 / 6-10) [256]')
    ap.add_argument('--skipTiles', action='store_true', help='reuse existing src10/src11 (tiling only)')
    ap.add_argument('--title', default=None, help='name shown in Google Earth [NISAR <product> <run name>]')
    args = ap.parse_args()
    work = os.path.abspath(args.work)
    out = os.path.abspath(args.out or f'{work}/googleEarth')
    src10, src11, stage = f'{out}/src10', f'{out}/src11', f'{out}/stage'
    for d in (src10, src11, stage):
        os.makedirs(d, exist_ok=True)
    cfg = yaml.safe_load(open(f'{work}/run.yaml'))
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tileRun = cfg['tileRun'] if os.path.isabs(cfg['tileRun']) else os.path.join(repo, cfg['tileRun'])
    vrtDir = f'{work}/vrt/{args.product}'
    res = abs(gdal.Open(f'{vrtDir}/global.vrt').GetGeoTransform()[5])

    if not args.skipTiles:
        import pyproj
        print('masks: land, glaciers, ice-sheet margins, 40/77 MHz coverage', flush=True)
        land = landGeometry()
        glaciers = glacierGeometry()
        icePaths = {e: iceClassRaster(land, e, args.marginKm, f'{stage}/ice{e}.npz') for e in POLAR}
        # footprints (valid ones where scanned) of every catalogue next to the tiling
        cycDir = os.path.dirname(tileRun)
        geoms = {}
        for p in sorted(glob.glob(f'{cycDir}/catalogue*.geojson')):
            for f in json.load(open(p))['features']:
                geoms[f['properties']['name']] = f['geometry']
        for p in sorted(glob.glob(f'{cycDir}/validFootprints*.geojson')):
            for f in json.load(open(p))['features']:
                geoms[f['properties']['name']] = f['geometry']
        jobs = []
        for t in json.load(open(f'{tileRun}/tiles.geojson'))['features']:
            name, epsg = t['properties']['name'], t['properties']['epsg']
            vrt = f'{vrtDir}/{name}.vrt'
            if not os.path.exists(vrt):
                continue
            rows = list(csv.DictReader(open(f'{tileRun}/tiles/{name}.csv')))
            hi = [geoms[r['name']] for r in rows if int(r['bwMHz']) >= 40 and geoms.get(r['name'])]
            if epsg == 4326:
                tb = shape(t['geometry'])
                landG = land.intersection(tb)
                hiG = unary_union([unwrap(shape(g)) for g in hi]).intersection(tb) if hi else box(0, 0, 0, 0)
                detailG = unary_union([hiG, glaciers.intersection(tb)]).intersection(landG)
                ice = [icePaths[e] for e, b in POLAR.items() if b.intersects(tb)]
            else:
                fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
                h = t['properties']['halfWidthM'] + 20000.
                tb = box(-h, -h, h, h)
                capLL = POLAR[3031] if epsg == 3031 else box(-180, 60, 180, 90)
                landG = transform(fwd, land.intersection(capLL).segmentize(0.1)).buffer(0).intersection(tb)
                hiG = unary_union([transform(fwd, shape(g)).buffer(0) for g in hi]) if hi else box(0, 0, 0, 0)
                glG = transform(fwd, glaciers.intersection(capLL).segmentize(0.1)).buffer(0)
                detailG = unary_union([hiG, glG]).intersection(landG)
                ice = [icePaths[epsg]] if epsg in icePaths else []
            jobs.append((name, vrt, epsg, landG, detailG, ice, src10, src11, args.smooth))
        print(f'{len(jobs)} tiles -> {src10}, {src11}', flush=True)
        with concurrent.futures.ProcessPoolExecutor(args.processes) as pool:
            for k, (name, notes) in enumerate(pool.map(tileSources, jobs), 1):
                if k % 50 == 0 or k == len(jobs):
                    print(f'  {k}/{len(jobs)} {name}: {", ".join(notes) or "no land"}', flush=True)

    caps = sorted(glob.glob(f'{vrtDir}/cap_*.vrt'))
    runName = args.title or f'NISAR {args.product} {os.path.basename(work)}'
    layers = []
    for title, sub, srcs, r, zoom in (
            ('base (zoom 0-9, everything)', 'base', caps + [f'{vrtDir}/global.vrt'], res, '0-9'),
            ('land (zoom 5-10, smoothed)', 'land', sorted(glob.glob(f'{src10}/*.tif')), 2 * res, '5-10'),
            ('detail (zoom 7-11: 40/77 MHz land, ice margins, glaciers)', 'detail',
             sorted(glob.glob(f'{src11}/*.tif')), res, '7-11')):
        if not srcs:
            print(f'{sub}: nothing to tile')
            continue
        rc = superoverlay(srcs, r, f'{out}/{sub}', zoom, processes=args.processes,
                          title=f'{runName} - {title}', tileSize=args.tileSize)
        if rc:
            return rc
        if sub != 'base':               # base shows the coarse zooms: keep only this layer's finest images
            from .linkOnly import stripCoarse
            stripCoarse(f'{out}/{sub}', args.processes)
        layers.append((title, sub))
    links = '\n'.join(f'    <NetworkLink><name>{n}</name><Link><href>{d}/doc.kml</href></Link></NetworkLink>'
                      for n, d in layers)
    with open(f'{out}/doc.kml', 'w') as fp:
        fp.write('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n'
                 f'  <Document>\n    <name>{runName}</name>\n{links}\n  </Document>\n</kml>\n')
    print(f'done: {out}/doc.kml')
    return 0


if __name__ == '__main__':
    sys.exit(main())
