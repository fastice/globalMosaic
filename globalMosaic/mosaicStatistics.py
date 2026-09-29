#!/usr/bin/env python3
'''
Global mosaic of per-frame temporal statistics (frameStatistics) with geomosaic -geo, on the tiles
of a backscatter tiling.

  mosaicStatistics cycle30/ascending --frames /scratch/ianj/mosaics/stats30asc/frames \\
      --work /scratch/ianj/mosaics/stats30asc [--nProc 16] [--tiles ...]

One geomosaic run per tile and layer: the tile's frame-statistics files (the frames the tiling
lists for that tile) with `band:` = the layer and `weightBand:` 5 (n, cycles) for layers 1-4, so
overlapping frames are averaged weighted by their sample count; feathered with -fl as usual.
Linear Float32 output (NaN = no data), cropped to the tile:
  <work>/tiles/<layer>/<tile>.tif, layers mean (gamma0, linear), sigma, cv, cvc (speckle floor
  removed), n. VRTs per layer: <work>/vrt/<layer>/{<tile>, global, cap_<name>}.vrt.
Progress: <work>/quicklooks/progress.png (60S-60N and both polar views, of cv) every 5 minutes.
Resumable: finished tiles are skipped. The frame files are kept (inputs for reruns).
'''
import argparse
import concurrent.futures
import csv
import glob
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from osgeo import gdal

from .runGeomosaicTiles import tileGrid, log, Progress

gdal.UseExceptions()
LAYERS = {'mean': 1, 'sigma': 2, 'cv': 3, 'cvc': 4, 'n': 5}
NODATA = -3000


def frameKey(name):
    f = name.split('_')
    return f'{f[5]}_{f[7]}_{f[8]}'


def runJob(job):
    ''' One tile, all layers: geomosaic -geo per layer, then crop to the tile. '''
    name, grid, files, work, geomosaic, env = job
    outs = {L: f'{work}/tiles/{L}/{name}.tif' for L in LAYERS}
    if all(os.path.exists(o) for o in outs.values()):
        return name, 0., 'skipped (done)'
    t0 = time.time()
    jd = f'{work}/jobs/{name}'
    os.makedirs(jd, exist_ok=True)
    u = grid['unit']
    with open(f'{jd}/inputFile', 'w') as fp:
        fp.write(f'{(grid["x0"] + grid["dx"] / 2) / u:.9f} {(grid["y0"] + grid["dy"] / 2) / u:.9f} '
                 f'{grid["nx"] * grid["dx"] / u:.9f} {grid["ny"] * grid["dy"] / u:.9f} '
                 f'{grid["dx"] / u:.12f} {grid["dy"] / u:.12f}\n;\n0\n;\n')
    notes = []
    for L, band in LAYERS.items():
        if os.path.exists(outs[L]):
            continue
        with open(f'{jd}/{L}.yaml', 'w') as fp:
            fp.write(f'band: {band}\n' + ('weightBand: 5\n' if band != 5 else '') + 'files:\n')
            for f in files:
                fp.write(f'  - {f}\n')
        cmd = [geomosaic, '-GTiff', '-ompThreads', '1', '-fl', str(grid['fl']), '-epsg', str(grid['epsg']),
               '-geo', f'{L}.yaml', 'inputFile', 'none', f'out.{L}']
        with open(f'{jd}/{L}.log', 'w') as fp:
            rc = subprocess.run(cmd, cwd=jd, env=env, stdout=fp, stderr=subprocess.STDOUT).returncode
        raw = sorted(glob.glob(f'{jd}/out.{L}*.tif'))
        if rc != 0 or not raw:
            notes.append(f'{L} FAILED rc={rc}')
            continue
        os.makedirs(os.path.dirname(outs[L]), exist_ok=True)
        cropFloat(raw[0], outs[L], grid)
        os.remove(raw[0])
    return name, time.time() - t0, ', '.join(notes) or f'{len(files)} frames'


def cropFloat(src, dst, grid):
    ''' The tile's exact window (as runGeomosaicTiles.cropToTile), kept Float32 with NaN no data. '''
    w, s, e, n = grid['tile']
    ds = gdal.Open(src)
    gt = ds.GetGeoTransform()
    xoff, yoff = round((w - gt[0]) / gt[1]), round((gt[3] - n) / -gt[5])
    xs, ys = round((e - w) / gt[1]), round((n - s) / -gt[5])
    gdal.Translate(dst, ds, srcWin=[xoff, yoff, xs, ys], noData=float('nan'),
                   creationOptions=['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=3', 'BIGTIFF=IF_SAFER'])


def previewPiece(tif, out):
    ''' An int16 copy of a cv tile (x 10000, -3000 = no data) for the progress picture. '''
    d = gdal.Open(tif)
    a = d.GetRasterBand(1).ReadAsArray()
    b = np.where(np.isfinite(a), np.round(np.clip(a, 0, 3) * 10000), NODATA).astype('i2')
    m = gdal.GetDriverByName('GTiff').Create(out, b.shape[1], b.shape[0], 1, gdal.GDT_Int16,
                                             options=['COMPRESS=DEFLATE'])
    m.SetGeoTransform(d.GetGeoTransform())
    m.SetProjection(d.GetProjection())
    m.GetRasterBand(1).SetNoDataValue(NODATA)
    m.GetRasterBand(1).WriteArray(b)
    m = None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', help='backscatter tiling (tiles.geojson + tiles/<tile>.csv)')
    ap.add_argument('--frames', required=True, help='directory of frameStatistics outputs')
    ap.add_argument('--work', required=True)
    ap.add_argument('--res', type=float, default=2.4 / 3600, help='lat/lon spacing, deg [2.4 arcsec]')
    ap.add_argument('--psResM', type=float, default=80., help='polar cap spacing, m [80]')
    ap.add_argument('--featherKm', type=float, default=10.)
    ap.add_argument('--nProc', type=int, default=16)
    ap.add_argument('--tiles', nargs='+', default=None)
    ap.add_argument('--geomosaic', default='geomosaic')
    ap.add_argument('--googleEarth', default=None,
                    help='also build Google Earth of these layers, e.g. "cv cvc" (zoom --zoom) [none]')
    ap.add_argument('--zoom', default='0-9')
    ap.add_argument('--title', default=None, help='Google Earth name prefix [NISAR temporal statistics <work name>]')
    args = ap.parse_args()
    work = os.path.abspath(args.work)
    frames = os.path.abspath(args.frames)
    os.makedirs(work, exist_ok=True)
    summary = open(f'{work}/summary.log', 'a')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    feats = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    if args.tiles:
        feats = [f for f in feats if f['properties']['name'] in args.tiles]
    jobs, grids = [], {}
    for f in feats:
        name = f['properties']['name']
        keys = {frameKey(r['name']) for r in csv.DictReader(open(f'{args.tileRun}/tiles/{name}.csv'))}
        files = sorted(p for p in (f'{frames}/{k}.tif' for k in keys) if os.path.exists(p))
        if not files:
            continue
        grid = tileGrid(f['properties'], f['geometry'], args.res, args.psResM, args.featherKm)
        grids[name] = grid
        jobs.append((name, grid, files, work, args.geomosaic, env))
    jobs.sort(key=lambda j: -len(j[2]))
    log(f'{len(jobs)} tiles with frame statistics, {args.nProc} processes, res {args.res * 3600:.1f} '
        f'arcsec / caps {args.psResM:.0f} m; layers {", ".join(LAYERS)}', summary)
    progress = Progress(work, feats, [(j[0], 'cv', j[0]) for j in jobs], grids, 'cv', dbMin=0., dbMax=60.)
    lastDraw = time.time()
    done = failed = 0
    with concurrent.futures.ProcessPoolExecutor(args.nProc) as pool:
        for name, sec, note in pool.map(runJob, jobs):
            done += 1
            failed += 'FAILED' in note
            log(f'[{done}/{len(jobs)}, {failed} failed] {name}: {sec / 60:.1f} min {note}', summary)
            cv = f'{work}/tiles/cv/{name}.tif'
            if os.path.exists(cv):
                os.makedirs(f'{work}/quicklooks/pieces', exist_ok=True)
                progress.addJob(name, 'cv', previewPiece(cv, f'{work}/quicklooks/pieces/{name}.tif'))
            progress.count('FAILED' not in note, name)
            if time.time() - lastDraw > 300:
                progress.draw()
                lastDraw = time.time()
    progress.draw()
    for L in LAYERS:
        vd = f'{work}/vrt/{L}'
        os.makedirs(vd, exist_ok=True)
        latlon, caps = [], []
        for name, grid, *_ in jobs:
            t = f'{work}/tiles/{L}/{name}.tif'
            if os.path.exists(t):
                gdal.BuildVRT(f'{vd}/{name}.vrt', [t])
                (latlon if grid['epsg'] == 4326 else caps).append((name, t))
        if latlon:
            gdal.BuildVRT(f'{vd}/global.vrt', [t for _, t in latlon], resolution='user', xRes=args.res, yRes=args.res)
        for name, t in caps:
            gdal.BuildVRT(f'{vd}/cap_{name}.vrt', [t])
        log(f'{L}: {len(latlon)} lat/lon tiles, {len(caps)} caps -> {vd}', summary)
    # Google Earth: grey stretch per layer, linear (CV 0-0.6, n 0-10)
    from .makeGoogleEarth import superoverlay
    stretches = {'cv': (0., 0.6, 1.), 'cvc': (0., 0.6, 1.), 'n': (0., 10., 1.)}
    title = args.title or f'NISAR temporal statistics {os.path.basename(work)}'
    for L in (args.googleEarth or '').split():
        vd = f'{work}/vrt/{L}'
        srcs = sorted(glob.glob(f'{vd}/cap_*.vrt')) + ([f'{vd}/global.vrt'] if os.path.exists(f'{vd}/global.vrt') else [])
        if not srcs:
            continue
        if L not in stretches:
            log(f'Google Earth {L}: no linear stretch defined (cv, cvc, n); skipped', summary)
            continue
        lo, hi, sc = stretches[L]
        rc = superoverlay(srcs, args.res, f'{work}/googleEarth/{L}', args.zoom, dbMin=lo, dbMax=hi,
                          processes=args.nProc, title=f'{title} - {L}', nodata=float('nan'), scale=sc)
        log(f'Google Earth {L}: {work}/googleEarth/{L}/doc.kml (rc {rc})', summary)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
