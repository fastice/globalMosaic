#!/usr/bin/env python3
'''
Coherence mosaic from NISAR GUNW products, on ANY tiling runGeomosaicTiles understands --
lat/lon, EASE-Grid, Equal Earth or Wagner VII.

  gunwCoherence cycle30gunw/easeAsc --work /scratch/ianj/mosaics/coh30asc --nProc 16
  gunwCoherence <tiling> --work <dir> --tiles X+01Y+03        # one tile, for a test

Two stages, both resumable:

  1. describe  a ~1.5 kB VRT per granule (<cache>/<name>.coh.vrt) that points at its
               frequency-A coherence (HH, else VV) WHERE IT SITS at ASF, carrying the
               geotransform and CRS the HDF5 subdataset does not report. No pixel is copied.
  2. mosaic    one `geomosaic -geo` per tile over those VRTs, then crop to the tile.
               geomosaic reads windows over the network, exactly as it does for GCOV.

Why geomosaic rather than the gdal.Warp loop in mosaicCoherence: warping INTO the tile CRS
needs that CRS's inverse, and Wagner VII has none in PROJ. geomosaic carries its own Wagner
inverse (mosaicSource grimpProj.c), so the projected case works there and only there.
mosaicCoherence remains the tool for the lat/lon product -- it also builds the Google Earth
superoverlay, which this does not.

Output: <work>/tiles/coherence/<tile>.tif, Float32 coherence 0..1 with NaN no data, plus
<work>/vrt/coherence/<tile>.vrt and a run summary in <work>/summary.log.
'''
import argparse
import concurrent.futures
import csv
import glob
import json
import os
import subprocess
import sys
import time

import numpy as np
from osgeo import gdal, osr

from .mosaicStatistics import cropFloat
from .runGeomosaicTiles import Progress, log, tileGrid

gdal.UseExceptions()
LAYER = '/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram'


def httpSetup(cookieDir, chunk=None):
    '''
    Earthdata login over https: the redirect needs a cookie jar.

    `chunk` sets CPL_VSIL_CURL_CHUNK_SIZE, and the right value depends entirely on the access
    pattern. 4 MB is best for the large sequential window reads geomosaic does (measured
    8.37 MB/s per stream against ASF, bench/benchmarkASFRead.py). It is badly wrong for
    describing a granule: HDF5 metadata is scattered small reads, and each one then drags a
    whole 4 MB chunk -- 22.1 MB per granule against 0.5 MB at the 16 kB default, measured on
    one granule step by step. Over a cycle that is 92 GB instead of 2 GB, which is what 150
    MB/s during a metadata-only stage turned out to be. So: default here, 4 MB in the
    geomosaic child.
    '''
    os.makedirs(cookieDir, exist_ok=True)
    cookies = os.path.join(cookieDir, 'cookies')
    opts = [('GDAL_HTTP_NETRC', 'YES'), ('GDAL_HTTP_COOKIEFILE', cookies),
            ('GDAL_HTTP_COOKIEJAR', cookies), ('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR'),
            ('GDAL_HTTP_MAX_RETRY', '5'), ('GDAL_HTTP_RETRY_DELAY', '10')]
    if chunk:
        opts.append(('CPL_VSIL_CURL_CHUNK_SIZE', str(chunk)))
    for k, v in opts:
        gdal.SetConfigOption(k, v)


def makeVrt(job):
    """
    A VRT that points at the granule's coherence where it sits, at ASF.

    The array itself is never copied: GDAL can open
    HDF5:"/vsicurl/<url>"://<layer>/<pol>/coherenceMagnitude directly, but a NISAR GUNW carries
    its geolocation in sibling xCoordinates/yCoordinates arrays, so that subdataset comes back
    with an identity geotransform and no CRS -- unusable as a geomosaic -geo input. This reads
    only those coordinate arrays and the epsg (a few kB) and writes a ~1.5 kB VRT carrying the
    geotransform and projection. geomosaic then reads windows over the network exactly as it
    does for GCOV, instead of 20 MB per granule being downloaded and stored (89 GB for one
    cycle). Verified bit-identical to a locally extracted GeoTIFF over a test window.
    """
    name, url, out, cookieDir = job
    if os.path.exists(out):
        return name, 'cached'
    httpSetup(cookieDir)
    t0 = time.time()
    # Retry: with many processes opening HDF5 over https at once, ASF refuses some of them,
    # and a granule that fails here is silently absent from the mosaic. The same granule
    # opens perfectly seconds later, so the failures are the server's rate limit rather than
    # bad data -- 2,900 into one cycle 24 run, hundreds failed and every one retried cleanly.
    err = ''
    for attempt in range(4):
        try:
            ds = gdal.OpenEx('/vsicurl/' + url, gdal.OF_MULTIDIM_RASTER)
            g = ds.GetRootGroup().OpenGroupFromFullname(LAYER)
            # HH if there, else VV: a GUNW is single-pol per frequency and either can occur
            pols = [p for p in ('HH', 'VV', 'HV', 'VH') if p in g.GetGroupNames()]
            if not pols:
                return name, 'no polarization group'
            arr = g.OpenGroup(pols[0]).OpenMDArray('coherenceMagnitude')
            ny, nx = [d.GetSize() for d in arr.GetDimensions()]
            xs = g.OpenMDArray('xCoordinates').ReadAsArray()
            ys = g.OpenMDArray('yCoordinates').ReadAsArray()
            epsg = int(g.OpenMDArray('projection').GetAttribute('epsg_code').Read())
            break
        except Exception as e:
            err = str(e).replace('\n', ' ')
            if attempt < 3:
                time.sleep(10 * (attempt + 1))
    else:
        # the message is long and ends with the useful part, so keep the tail, not the head
        return name, f'read failed after 4 attempts: ...{err[-160:]}'
    dx, dy = float(xs[1] - xs[0]), float(ys[1] - ys[0])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    # LAYER has a leading slash and the subdataset syntax supplies its own: "://" + "science/..."
    sub = f'HDF5:"/vsicurl/{url}"://{LAYER.lstrip("/")}/{pols[0]}/coherenceMagnitude'
    vrt = (f'<VRTDataset rasterXSize="{nx}" rasterYSize="{ny}">\n'
           f'  <SRS>{srs.ExportToWkt()}</SRS>\n'
           f'  <GeoTransform>{xs[0] - dx / 2.}, {dx}, 0, {ys[0] - dy / 2.}, 0, {dy}</GeoTransform>\n'
           f'  <VRTRasterBand dataType="Float32" band="1">\n'
           f'    <NoDataValue>nan</NoDataValue>\n'
           f'    <SimpleSource>\n'
           f'      <SourceFilename relativeToVRT="0">{sub}</SourceFilename>\n'
           f'      <SourceBand>1</SourceBand>\n'
           f'      <SrcRect xOff="0" yOff="0" xSize="{nx}" ySize="{ny}"/>\n'
           f'      <DstRect xOff="0" yOff="0" xSize="{nx}" ySize="{ny}"/>\n'
           f'    </SimpleSource>\n'
           f'  </VRTRasterBand>\n</VRTDataset>\n')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + '.tmp'
    with open(tmp, 'w') as fp:
        fp.write(vrt)
    os.replace(tmp, out)
    return name, f'{nx}x{ny} EPSG:{epsg} {pols[0]} {time.time() - t0:.0f}s"""'.replace('"""', '')


def mosaicTile(job):
    ''' One geomosaic -geo run over a tile's coherence files. '''
    name, grid, files, jobDir, out, geomosaic, crs = job
    if os.path.exists(out):
        return name, 0, 0., 'skipped (done)'
    os.makedirs(jobDir, exist_ok=True)
    with open(f'{jobDir}/coh.yaml', 'w') as fp:
        fp.write('band: 1\nfiles:\n')
        for f in files:
            fp.write(f'  - {f}\n')
    # exactly as runGeomosaicTiles.writeJob: first-pixel CENTRE, extent and spacing, in the
    # grid's own unit (km for a projected tile, degrees for EPSG 4326), then no SAR inputs
    u = grid['unit']
    with open(f'{jobDir}/inputFile', 'w') as fp:
        fp.write(f'{(grid["x0"] + grid["dx"] / 2) / u:.9f} {(grid["y0"] + grid["dy"] / 2) / u:.9f} '
                 f'{grid["nx"] * grid["dx"] / u:.9f} {grid["ny"] * grid["dy"] / u:.9f} '
                 f'{grid["dx"] / u:.12f} {grid["dy"] / u:.12f}\n;\n0\n;\n')
    cmd = [geomosaic, '-GTiff', '-ompThreads', '1', '-fl', str(grid['fl'])] + crs + \
          ['-geo', 'coh.yaml', 'inputFile', 'none', 'out']
    t0 = time.time()
    cookies = os.path.abspath(os.path.join(os.path.dirname(jobDir), '..', 'cookies', 'cookies'))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               # the VRTs resolve to /vsicurl/, so the child needs the Earthdata login,
               # the cookie jar for its redirect, and the 4 MB chunk size the ASF benchmark
               # settled on (bench/benchmarkASFRead.py)
               GDAL_HTTP_NETRC='YES', GDAL_HTTP_COOKIEFILE=cookies, GDAL_HTTP_COOKIEJAR=cookies,
               GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_HTTP_MAX_RETRY='5',
               GDAL_HTTP_RETRY_DELAY='10', CPL_VSIL_CURL_CHUNK_SIZE='4194304')
    with open(f'{jobDir}/log', 'w') as fp:
        rc = subprocess.run(cmd, cwd=jobDir, env=env, stdout=fp,
                            stderr=subprocess.STDOUT).returncode
    raw = sorted(glob.glob(f'{jobDir}/out*.tif'))
    if rc != 0 or not raw:
        return name, rc or 1, time.time() - t0, f'FAILED rc={rc}'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cropFloat(raw[0], out, grid)
    return name, 0, time.time() - t0, f'{len(files)} granule(s)'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', help='a tiling directory (easeTiles or globalGCOVTiles)')
    ap.add_argument('--work', required=True, help='output directory')
    ap.add_argument('--cache', default=None,
                    help='extracted coherence GeoTIFFs [<work>/coherenceFrames]')
    ap.add_argument('--res', type=float, default=0.0025, help='lat/lon tile spacing, deg')
    ap.add_argument('--psResM', type=float, default=1000., help='projected tile spacing, m [1000]')
    ap.add_argument('--featherKm', type=float, default=10., help='feather width, km [10]')
    ap.add_argument('--nProc', type=int, default=16, help='concurrent jobs [16]')
    ap.add_argument('--tiles', nargs='+', default=None, help='only these tiles')
    ap.add_argument('--extractOnly', action='store_true', help='stop after stage 1')
    ap.add_argument('--progressMinutes', type=float, default=5.,
                    help='redraw quicklooks/progress.png at most this often [5]')
    ap.add_argument('--geomosaic', default='geomosaic', help='geomosaic executable')
    args = ap.parse_args()

    cache = args.cache or f'{args.work}/coherenceFrames'
    os.makedirs(args.work, exist_ok=True)
    summary = open(f'{args.work}/summary.log', 'a')     # log() wants a file object
    feats = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    if args.tiles:
        feats = [f for f in feats if f['properties']['name'] in args.tiles]
    if not feats:
        sys.exit('gunwCoherence: no tiles selected')

    # which granules each tile needs, and the union to extract
    want, urls = {}, {}
    for f in feats:
        name = f['properties']['name']
        path = f'{args.tileRun}/tiles/{name}.csv'
        if not os.path.exists(path):
            continue
        rows = list(csv.DictReader(open(path)))
        want[name] = [r['name'] for r in rows]
        for r in rows:
            urls[r['name']] = r['url']
    log(f'{len(feats)} tile(s), {len(urls)} distinct granule(s)', summary)

    todo = [(n, u, f'{cache}/{n}.coh.vrt', f'{args.work}/cookies')
            for n, u in sorted(urls.items()) if not os.path.exists(f'{cache}/{n}.coh.vrt')]
    log(f'stage 1: {len(todo)} granule(s) to describe, {len(urls) - len(todo)} cached', summary)
    bad = 0
    if todo:
        # PROCESSES, not threads. GDAL serialises every HDF5 access behind one global lock
        # (libhdf5 is not built thread-safe), so N threads read one granule at a time and the
        # whole stage runs at single-stream speed -- measured at 24 MB/s against ASF, where
        # the GCOV mosaic reaches ~1 GB/s precisely because it forks 32 geomosaic processes.
        with concurrent.futures.ProcessPoolExecutor(args.nProc) as pool:
            for k, (name, note) in enumerate(pool.map(makeVrt, todo), 1):
                if not note.endswith('s') and note not in ('cached',):
                    bad += 1
                    log(f'  [{k}/{len(todo)}] {name}: {note}', summary)
                elif k % 25 == 0 or k == len(todo):
                    log(f'  [{k}/{len(todo)}] {name}: {note}', summary)
    log(f'stage 1 done, {bad} granule(s) unusable', summary)
    if args.extractOnly:
        return 0

    # the tile grids, and the CRS argument geomosaic needs
    jobs = []
    for f in feats:
        p = f['properties']
        name = p['name']
        files = [f'{cache}/{n}.coh.vrt' for n in want.get(name, [])]
        files = [x for x in files if os.path.exists(x)]
        if not files:
            continue
        g = tileGrid(p, f['geometry'], args.res, args.psResM, args.featherKm)
        crs = (['-wkt', os.path.abspath(f'{args.tileRun}/{g["wkt"]}')] if g.get('wkt')
               else ['-epsg', str(g['epsg'])])
        jobs.append((name, g, files, f'{args.work}/jobs/{name}',
                     f'{args.work}/tiles/coherence/{name}.tif', args.geomosaic, crs))
    log(f'stage 2: {len(jobs)} tile(s) to mosaic, {args.nProc} processes', summary)
    # the same live picture the backscatter run draws, but 0..1 linear with NaN no data
    g0 = jobs[0][1] if jobs else {}
    projCRS = (open(os.path.abspath(f'{args.tileRun}/{g0["wkt"]}')).read() if g0.get('wkt')
               else (f'EPSG:{g0["epsg"]}' if g0.get('epsg') else None))
    progress = Progress(args.work, feats, [(j[0], 'coherence', j[0]) for j in jobs],
                        {j[0]: j[1] for j in jobs}, 'coherence', dbMin=0., dbMax=1.,
                        projCRS=projCRS, scale=1., nodata=float('nan'))
    lastDraw = time.time()
    failed = 0
    with concurrent.futures.ThreadPoolExecutor(args.nProc) as pool:
        futures = {pool.submit(mosaicTile, j): j[0] for j in jobs}
        for k, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            name, rc, sec, note = fut.result()
            failed += 1 if rc else 0
            progress.count(rc == 0, name)
            if rc == 0:
                progress.addJob(name, 'coherence', f'{args.work}/tiles/coherence/{name}.tif')
            log(f'[{k}/{len(jobs)}, {failed} failed] {name}: {sec / 60:.1f} min {note}', summary)
            if time.time() - lastDraw > args.progressMinutes * 60:
                progress.draw()
                lastDraw = time.time()
    progress.draw()

    tifs = sorted(glob.glob(f'{args.work}/tiles/coherence/*.tif'))
    if tifs:
        os.makedirs(f'{args.work}/vrt/coherence', exist_ok=True)
        gdal.BuildVRT(f'{args.work}/vrt/coherence/coherence.vrt', tifs)
        log(f'{len(tifs)} tile(s) -> {args.work}/vrt/coherence/coherence.vrt', summary)
    log(f'done; {failed} failed tile(s)', summary)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
