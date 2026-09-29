#!/usr/bin/env python3
'''
Global coherence mosaic from NISAR GUNW interferograms, on the lat/lon tiles of a GUNW tiling
(globalGCOVTiles on GUNW catalogues), then quick looks and a Google Earth superoverlay.

  mosaicCoherence cycle30gunw/ascending --work /scratch/ianj/mosaics/coherence30asc [--nProc 16]

For each tile, every listed GUNW's frequency-A unwrappedInterferogram/HH/coherenceMagnitude (80 m,
UTM or polar stereographic) is read over https (window covering the tile only), warped onto the
tile grid (average resampling) and averaged with the others where they overlap. NaN (no valid
data) is left out. Tiles are int16 coherence x 10000, no data -3000, in <work>/tiles/coherence;
VRTs in <work>/vrt/coherence (per tile, global.vrt for lat/lon, cap_<name>.vrt per polar cap).
Resumable: finished tiles are skipped. Then (unless --noGoogleEarth) quick looks
(<work>/quicklooks) and Google Earth zoom 0-9 (<work>/googleEarth/doc.kml, grey 0..1).

--res default 9" (0.0025 deg): Google Earth zoom 9 has 9.9" pixels; 9" divides the 6-deg tiles.
'''
import argparse
import concurrent.futures
import csv
import glob
import json
import os
import sys
import time

import numpy as np
from osgeo import gdal, osr

from .runGeomosaicTiles import tileGrid, log

gdal.UseExceptions()
NODATA = -3000
LAYER = '/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram'


def httpSetup(cookieDir):
    cookies = f'{cookieDir}/cookies.{os.getpid()}'
    for k, v in (('GDAL_HTTP_NETRC', 'YES'), ('GDAL_HTTP_COOKIEFILE', cookies),
                 ('GDAL_HTTP_COOKIEJAR', cookies), ('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR'),
                 ('GDAL_HTTP_MAX_RETRY', '5'), ('GDAL_HTTP_RETRY_DELAY', '10')):
        gdal.SetConfigOption(k, v)


def readWindow(url, dstSrs, bounds):
    ''' Coherence of one GUNW over the window covering `bounds` (in dstSrs), as a MEM dataset in the
    granule's own CRS, or None if they do not overlap. '''
    ds = gdal.OpenEx('/vsicurl/' + url, gdal.OF_MULTIDIM_RASTER)
    g = ds.GetRootGroup().OpenGroupFromFullname(LAYER)
    xs = g.OpenMDArray('xCoordinates').ReadAsArray()
    ys = g.OpenMDArray('yCoordinates').ReadAsArray()
    epsg = int(g.OpenMDArray('projection').GetAttribute('epsg_code').Read())
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(dstSrs, srs)
    x0, y0, x1, y1 = bounds
    n = 21
    edge = [(x0 + (x1 - x0) * i / (n - 1), y) for i in range(n) for y in (y0, y1)] + \
           [(x, y0 + (y1 - y0) * i / (n - 1)) for i in range(n) for x in (x0, x1)]
    pts = []
    for x, y in edge:
        try:
            pts.append(ct.TransformPoint(x, y)[:2])
        except Exception:
            pass
    if not pts:
        return None
    px, py = [p[0] for p in pts], [p[1] for p in pts]
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]
    c0 = max(int(np.floor((min(px) - xs[0]) / dx)) - 2, 0)
    c1 = min(int(np.ceil((max(px) - xs[0]) / dx)) + 3, len(xs))
    r0 = max(int(np.floor((max(py) - ys[0]) / dy)) - 2, 0)
    r1 = min(int(np.ceil((min(py) - ys[0]) / dy)) + 3, len(ys))
    if c1 - c0 < 2 or r1 - r0 < 2:
        return None
    coh = g.OpenGroup('HH').OpenMDArray('coherenceMagnitude').GetView(f'[{r0}:{r1},{c0}:{c1}]').ReadAsArray()
    if not np.isfinite(coh).any():
        return None
    mem = gdal.GetDriverByName('MEM').Create('', coh.shape[1], coh.shape[0], 1, gdal.GDT_Float32)
    mem.SetGeoTransform((xs[c0] - dx / 2, dx, 0., ys[r0] - dy / 2, 0., dy))
    mem.SetProjection(srs.ExportToWkt())
    b = mem.GetRasterBand(1)
    b.WriteArray(coh.astype('f4'))
    b.SetNoDataValue(float('nan'))
    return mem


def tileJob(job):
    name, grid, urls, out, cookieDir = job
    if os.path.exists(out):
        return name, 0, 0, 'skipped (done)'
    httpSetup(cookieDir)
    t0 = time.time()
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(grid['epsg'])
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    x0, y0 = grid['x0'], grid['y0']
    x1, y1 = x0 + grid['nx'] * grid['dx'], y0 + grid['ny'] * grid['dy']
    total = np.zeros((grid['ny'], grid['nx']), 'f8')
    count = np.zeros((grid['ny'], grid['nx']), 'i4')
    used, failed = 0, []
    for url in urls:
        mem = None
        for attempt in range(3):
            try:
                mem = readWindow(url, dst, (x0, y0, x1, y1))
                break
            except Exception as e:
                err = str(e)[:200]
                time.sleep(20)
        else:
            failed.append(f'{os.path.basename(url)[:60]}: {err}')
            continue
        if mem is None:
            continue
        w = gdal.Warp('', mem, format='MEM', dstSRS=dst.ExportToWkt(), outputBounds=(x0, y0, x1, y1),
                      width=grid['nx'], height=grid['ny'], resampleAlg='average',
                      srcNodata=float('nan'), dstNodata=float('nan'))
        a = w.GetRasterBand(1).ReadAsArray()
        ok = np.isfinite(a)
        total[ok] += a[ok]
        count[ok] += 1
        used += 1
    coh = np.where(count > 0, np.round(total / np.maximum(count, 1) * 10000), NODATA).astype('i2')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + '.tmp.tif'
    d = gdal.GetDriverByName('GTiff').Create(tmp, grid['nx'], grid['ny'], 1, gdal.GDT_Int16,
                                             options=['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=2'])
    d.SetGeoTransform((x0, grid['dx'], 0., y1, 0., -grid['dy']))
    d.SetProjection(dst.ExportToWkt())
    d.GetRasterBand(1).SetNoDataValue(NODATA)
    d.GetRasterBand(1).WriteArray(coh)                # warped north-up, like the geotransform
    d = None
    os.replace(tmp, out)
    note = f'{used}/{len(urls)} granules' + (f', {len(failed)} FAILED: ' + '; '.join(failed[:2]) if failed else '')
    return name, time.time() - t0, len(failed), note


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', help='GUNW tiling (globalGCOVTiles output)')
    ap.add_argument('--work', required=True, help='output directory')
    ap.add_argument('--res', type=float, default=9. / 3600, help='latitude spacing, deg [9 arcsec]')
    ap.add_argument('--psResM', type=float, default=300., help='polar cap spacing, m [300]')
    ap.add_argument('--nProc', type=int, default=16)
    ap.add_argument('--tiles', nargs='+', default=None, help='only these tiles')
    ap.add_argument('--zoom', default='0-9', help='Google Earth zoom levels [0-9]')
    ap.add_argument('--noGoogleEarth', action='store_true')
    args = ap.parse_args()
    work = os.path.abspath(args.work)
    os.makedirs(f'{work}/jobs', exist_ok=True)
    summary = open(f'{work}/summary.log', 'a')
    feats = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    if args.tiles:
        feats = [f for f in feats if f['properties']['name'] in args.tiles]
    jobs = []
    for f in feats:
        name, props = f['properties']['name'], f['properties']
        g = tileGrid(props, f['geometry'], args.res, args.psResM, 0.)
        # tileGrid: edges x0, y0 (south-west), sizes nx, ny, spacings dx, dy (> 0); no margin
        grid = dict(epsg=g['epsg'], x0=g['x0'], y0=g['y0'], nx=g['nx'], ny=g['ny'], dx=g['dx'], dy=g['dy'])
        urls = [r['url'] for r in csv.DictReader(open(f'{args.tileRun}/tiles/{name}.csv'))]
        jobs.append((name, grid, urls, f'{work}/tiles/coherence/{name}.tif', f'{work}/jobs'))
    jobs.sort(key=lambda j: -len(j[2]))                       # biggest first
    log(f'{len(jobs)} tiles, {args.nProc} processes, res {args.res * 3600:.1f} arcsec / caps '
        f'{args.psResM:.0f} m; GUNW coherence', summary)
    done = failedTiles = 0
    with concurrent.futures.ProcessPoolExecutor(args.nProc) as pool:
        for name, sec, nFail, note in pool.map(tileJob, jobs):
            done += 1
            failedTiles += bool(nFail)
            log(f'[{done}/{len(jobs)}, {failedTiles} with failures] {name}: {sec / 60:.1f} min {note}', summary)
    # VRTs: per tile, global lat/lon, one per cap
    vd = f'{work}/vrt/coherence'
    os.makedirs(vd, exist_ok=True)
    latlon, caps = [], []
    for name, grid, _, out, _ in jobs:
        if not os.path.exists(out):
            continue
        gdal.BuildVRT(f'{vd}/{name}.vrt', [out], srcNodata=NODATA, VRTNodata=NODATA)
        (latlon if grid['epsg'] == 4326 else caps).append((name, out))
    if latlon:
        gdal.BuildVRT(f'{vd}/global.vrt', [o for _, o in latlon], resolution='user', xRes=args.res,
                      yRes=args.res, srcNodata=NODATA, VRTNodata=NODATA)
    for name, out in caps:
        gdal.BuildVRT(f'{vd}/cap_{name}.vrt', [out], srcNodata=NODATA, VRTNodata=NODATA)
    log(f'VRTs: {len(latlon)} lat/lon tiles, {len(caps)} caps -> {vd}', summary)
    if args.noGoogleEarth or not latlon:
        return 0
    import subprocess
    subprocess.run([sys.executable, '-m', 'globalMosaic.makeQuickLook', work, '--product', 'coherence',
                    '--dbMin', '0', '--dbMax', '100', '--processes', str(args.nProc)])
    from .makeGoogleEarth import superoverlay
    rc = superoverlay([f'{vd}/cap_{n}.vrt' for n, _ in caps] + [f'{vd}/global.vrt'], args.res,
                      f'{work}/googleEarth', args.zoom, dbMin=0., dbMax=100., processes=args.nProc)
    log(f'Google Earth: {work}/googleEarth/doc.kml (rc {rc})', summary)
    return rc


if __name__ == '__main__':
    sys.exit(main())
