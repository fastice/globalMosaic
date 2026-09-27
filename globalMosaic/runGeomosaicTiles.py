#!/usr/bin/env python3
'''
Run geomosaic over the lat/lon tiles made by globalGCOVTiles.py and assemble a VRT hierarchy.
Standalone: the only program it calls is geomosaic; GDAL's Python API crops and builds VRTs.

  runGeomosaicTiles cycle30/both --work run_c030 --granules vsicurl --nProc 48 [--dryRun]
  runGeomosaicTiles --config run.yaml [--nProc 24 ...]

Run file (--config)
  A yaml whose keys are the option names below (tileRun, work, granules, nProc, date1, tiles,
  excludeTiles, ...). Options given on the command line override it. The resolved settings are
  written to <work>/run.yaml, as the record of what the run used.

Jobs
  One geomosaic run per tile and per frequency/polarization group (tiles/<tile>.<F><PP>.yaml from
  the tiler: frequency A = the science band, B = the 5 MHz band, which geomosaic reads on its own
  grid). Each job runs single-threaded -- OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1 and
  -ompThreads 1 -- and --nProc jobs run at once. A job is skipped when its outputs exist, so a
  rerun resumes.
    geomosaic -GTiff -removePad 5 -S1Cal -int16 -calOutput gamma0 -ompThreads 1 -fl <fl>
              -epsg <4326|3031|3413> [-date1 .. -date2 ..] -gcov <yaml> <inputFile> none <out>
  The input-file header origin is ALWAYS the centre of the first pixel: geomosaic's GeoTIFF writer
  subtracts half a pixel unconditionally (gdalIO tiffWriteCode.c). Its -center flag is a no-op;
  don't tie the convention to it.
  Needs a geomosaic with geographic output (-epsg 4326, header in degrees) and `dem none`.

Grid
  Latitude spacing --res (degrees, default 3 arcsec); longitude spacing --res x the band's lonMult
  (from tiles.geojson). Tile edges are exact multiples of both, so tiles abut exactly. Each job
  covers the tile plus a --featherKm margin (whole pixels) and is then cropped to the tile.

Granule paths (--granules)
  vsicurl        /vsicurl/<ASF url> from the tile CSV (Earthdata login in ~/.netrc); works anywhere
  s3             /vsis3/<ASF bucket>/<same path>: direct S3 reads, for runs in AWS us-west-2 only.
                 Temporary keys come from ASF's s3credentials endpoint (Earthdata login in
                 ~/.netrc); they last 1 h, so each job starts with keys that have >= 50 min left,
                 and a job that dies on an expired key is retried with fresh ones.
  DIR            <DIR>/<name>.h5 (local copies)
  TEMPLATE       any string with {name} or {url}, e.g. '/vsis3/my-bucket/gcov/{name}.h5'
  --factorFrom   passed through to the yaml for slim granules (shared RTC factor directory)

Outputs (--work DIR)
  jobs/<tile>.<grp>/{inputFile,gcov.yaml,out.*,log}   raw geomosaic runs (with margin)
  tiles/<product>/<tile>.<grp>.tif                     cropped to the tile, int16, nodata -3000
  vrt/<product>/<tile>.vrt       tile: frequency-B layer under frequency-A layer
  vrt/<product>/band_<lat0>_<lat1>.vrt                 one latitude band at its native spacing
  vrt/<product>/global.vrt       all lat/lon bands (resampled to the finest spacing on read)
  vrt/<product>/cap_<name>.vrt   polar cap tiles (EPSG 3031/3413, own CRS, not in global.vrt)
  products: gamma0 by default (-calOutput gamma0); --calOutput both adds sigma0
Polar caps
  Poleward of 84 deg the tiler makes one polar stereographic tile per pole (lat/lon degenerates
  there): header in km, --psResM spacing, -epsg 3031 / 3413.
'''
import argparse
import concurrent.futures
import datetime
import glob
import json
import math
import os
import subprocess
import sys
import threading
import time

import requests
import yaml
from osgeo import gdal

gdal.UseExceptions()
NODATA = -3000


def log(msg, fp=None):
    line = f'{time.strftime("%m-%d %H:%M:%S")} {msg}'
    print(line, flush=True)
    if fp:
        fp.write(line + '\n')
        fp.flush()


ASF_HTTPS = 'https://nisar.asf.earthdatacloud.nasa.gov/NISAR/'
ASF_BUCKET = 'sds-n-cumulus-prod-nisar-products'
ASF_S3CREDS = 'https://nisar.asf.earthdatacloud.nasa.gov/s3credentials'


def granulePath(spec, name, url):
    if spec == 'vsicurl':
        return f'/vsicurl/{url}'
    if spec == 's3':
        if not url.startswith(ASF_HTTPS):
            raise ValueError(f'not an ASF NISAR url: {url}')
        return f'/vsis3/{ASF_BUCKET}/{url[len(ASF_HTTPS):]}'
    if '{' in spec:
        return spec.format(name=name, url=url)
    return os.path.join(spec, f'{name}.h5')


def tileGrid(props, geom, res, psResM, featherKm):
    ''' Output grid of one job: the tile plus the feather margin, on the tile's pixel lattice.
    Lat/lon tiles (EPSG 4326) are in degrees; polar cap tiles (kind 'ps') are squares about the
    pole in metres. Edges (x0, y0) are pixel EDGES; the header adds half a pixel, because geomosaic
    always reads the header origin as the first pixel's centre. '''
    if props.get('kind') == 'ps':
        half = math.ceil(props['halfWidthM'] / psResM) * psResM
        pad = math.ceil(featherKm * 1000. / psResM)
        x0 = -half - pad * psResM
        n = int(round(2 * half / psResM)) + 2 * pad
        return dict(epsg=props['epsg'], unit=1000., tile=(-half, -half, half, half), dx=psResM,
                    dy=psResM, x0=x0, y0=x0, nx=n, ny=n, fl=round(featherKm * 1000. / psResM))
    lons = [c[0] for c in geom['coordinates'][0]]
    lats = [c[1] for c in geom['coordinates'][0]]
    w, e, s, n = min(lons), max(lons), min(lats), max(lats)
    dLat, dLon = res, res * props['lonMult']
    latMid = 0.5 * (s + n)
    padY = math.ceil(featherKm / 111.32 / dLat)
    # a tile spanning all longitudes already wraps: no east-west margin
    padX = 0 if e - w >= 360 else \
        math.ceil(featherKm / (111.32 * max(math.cos(math.radians(latMid)), 0.05)) / dLon)
    x0, y0 = w - padX * dLon, max(s - padY * dLat, -90.)
    nx = round((e - w) / dLon) + 2 * padX
    ny = round((min(n + padY * dLat, 90.) - y0) / dLat)
    return dict(epsg=4326, unit=1., tile=(w, s, e, n), dx=dLon, dy=dLat, x0=x0, y0=y0, nx=nx, ny=ny,
                fl=round(featherKm / (res * 111.32)))


def writeJob(jobDir, grid, yamlIn, csvRows, spec, factorFrom):
    os.makedirs(jobDir, exist_ok=True)
    # geomosaic input file: first-pixel CENTRE (always -- the writer subtracts half a pixel), size
    # and spacing -- degrees for EPSG 4326, km for a projected tile -- then no range/Doppler inputs
    u = grid['unit']
    with open(f'{jobDir}/inputFile', 'w') as fp:
        fp.write(f'{(grid["x0"] + grid["dx"] / 2) / u:.9f} {(grid["y0"] + grid["dy"] / 2) / u:.9f} '
                 f'{grid["nx"] * grid["dx"] / u:.9f} {grid["ny"] * grid["dy"] / u:.9f} '
                 f'{grid["dx"] / u:.12f} {grid["dy"] / u:.12f}\n;\n0\n;\n')
    urls = {r['name']: r['url'] for r in csvRows}
    out = []
    for line in open(yamlIn):
        s = line.strip()
        if s.startswith('- ') and s.endswith('.h5'):
            name = s[2:-3]
            out.append(f'  - {granulePath(spec, name, urls.get(name, ""))}\n')
        else:
            out.append(line)
    if factorFrom:
        out.insert(1, f'factorFrom: {factorFrom}\n')
    with open(f'{jobDir}/gcov.yaml', 'w') as fp:
        fp.writelines(out)


REMOTE_ERRORS = ('could not open remote input', 'CURL error', 'Could not resolve host',
                 'HTTP response code', 'ExpiredToken', 'AccessDenied', 'InvalidAccessKeyId')


class S3Keys:
    ''' ASF temporary S3 keys (1 h), shared by all jobs and refreshed when < minLeft remain. '''

    def __init__(self, minLeft=3000.):
        self.lock, self.keys, self.expires, self.minLeft = threading.Lock(), None, 0., minLeft

    def env(self):
        with self.lock:
            if self.keys is None or self.expires - time.time() < self.minLeft:
                r = requests.get(ASF_S3CREDS, timeout=60)     # Earthdata login from ~/.netrc
                r.raise_for_status()
                k = r.json()
                self.keys = dict(AWS_ACCESS_KEY_ID=k['accessKeyId'],
                                 AWS_SECRET_ACCESS_KEY=k['secretAccessKey'],
                                 AWS_SESSION_TOKEN=k['sessionToken'], AWS_REGION='us-west-2')
                exp = datetime.datetime.strptime(k['expiration'], '%Y-%m-%d %H:%M:%S%z')
                self.expires = exp.timestamp()
            return dict(self.keys)


def runJob(jobDir, cmdBase, products, retries=2, s3Keys=None):
    ''' Run one geomosaic job; retry when it died on a remote (network) error -- geomosaic makes a
    failed /vsi open fatal, and over hundreds of remote jobs a DNS/network blip is certain. '''
    if all(os.path.exists(f'{jobDir}/out.{p}.tif') for p in products):
        return jobDir, 0, 0., 'skipped (done)'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               GOTO_NUM_THREADS='1',
               # remote granules: Earthdata login from ~/.netrc; GDAL follows the signing redirect on
               # each open, so the (unsigned) ASF URLs in the list do not expire mid-job
               GDAL_HTTP_NETRC='YES', GDAL_HTTP_COOKIEFILE=f'{jobDir}/cookies',
               GDAL_HTTP_COOKIEJAR=f'{jobDir}/cookies', GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
               GDAL_HTTP_MAX_RETRY='5', GDAL_HTTP_RETRY_DELAY='10')
    t0 = time.time()
    for attempt in range(retries + 1):
        if s3Keys is not None:
            env.update(s3Keys.env())
        with open(f'{jobDir}/log', 'w') as fp:
            rc = subprocess.run(cmdBase + ['-gcov', 'gcov.yaml', 'inputFile', 'none', 'out'],
                                cwd=jobDir, env=env, stdout=fp, stderr=subprocess.STDOUT).returncode
        missing = [p for p in products if not os.path.exists(f'{jobDir}/out.{p}.tif')]
        if rc == 0 and not missing:
            break
        logText = open(f'{jobDir}/log', errors='ignore').read()
        if attempt < retries and any(e in logText for e in REMOTE_ERRORS):
            os.replace(f'{jobDir}/log', f'{jobDir}/log.attempt{attempt + 1}')
            time.sleep(60)
            continue
        break
    note = ('missing ' + ','.join(missing)) if missing else ''
    if attempt:
        note += f' ({attempt} remote retr{"y" if attempt == 1 else "ies"})'
    return jobDir, rc, time.time() - t0, note


def cropToTile(src, dst, grid):
    ''' Cut the margin off: the tile's exact pixel window, from the output's own geotransform. '''
    w, s, e, n = grid['tile']
    ds = gdal.Open(src)
    gt = ds.GetGeoTransform()
    xoff, yoff = round((w - gt[0]) / gt[1]), round((gt[3] - n) / -gt[5])
    xs, ys = round((e - w) / gt[1]), round((n - s) / -gt[5])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    gdal.Translate(dst, ds, srcWin=[xoff, yoff, xs, ys], noData=NODATA,
                   creationOptions=['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=2', 'BIGTIFF=IF_SAFER'])


def parseWithConfig(ap):
    '''
    Parse the command line over the --config yaml: yaml keys become defaults (unknown keys are an
    error), so anything given on the command line wins. tileRun, work and granules are required
    from one or the other.
    '''
    pre, _ = ap.parse_known_args()
    if pre.config:
        with open(pre.config) as fp:
            cfg = yaml.safe_load(fp) or {}
        known = {a.dest for a in ap._actions}
        bad = sorted(set(cfg) - known)
        if bad:
            ap.error(f'{pre.config}: unknown keys {bad}')
        ap.set_defaults(**cfg)
    args = ap.parse_args()
    missing = [k for k in ('tileRun', 'work', 'granules') if getattr(args, k) is None]
    if missing:
        ap.error(f'missing {missing} (command line or --config)')
    return args


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', nargs='?', default=None, help='output directory of globalGCOVTiles.py')
    ap.add_argument('--config', default=None, help='run yaml (keys = option names); command line overrides')
    ap.add_argument('--work', default=None, help='work/output directory')
    ap.add_argument('--granules', default=None, help='vsicurl | s3 | DIR | TEMPLATE with {name}/{url}')
    ap.add_argument('--factorFrom', default=None, help='shared RTC factor directory (slim granules)')
    ap.add_argument('--res', type=float, default=3. / 3600, help='latitude spacing, deg [3 arcsec]')
    ap.add_argument('--psResM', type=float, default=100., help='polar cap tile spacing, m [100]')
    ap.add_argument('--calOutput', choices=['gamma0', 'sigma0', 'both'], default='gamma0',
                    help='geomosaic -calOutput [gamma0: skips the RTC factor band, ~3x faster]')
    ap.add_argument('--featherKm', type=float, default=10., help='feather length and margin [10]')
    ap.add_argument('--removePad', type=int, default=5)
    ap.add_argument('--date1', default=None, help='MM-DD-YYYY (geomosaic -date1)')
    ap.add_argument('--date2', default=None, help='MM-DD-YYYY (geomosaic -date2)')
    ap.add_argument('--nProc', type=int, default=16, help='concurrent single-threaded jobs [16]')
    ap.add_argument('--tiles', nargs='+', default=None, help='only these tile names')
    ap.add_argument('--excludeTiles', nargs='+', default=None, help='skip these tile names')
    ap.add_argument('--geomosaic', default='geomosaic', help='geomosaic executable [on PATH]')
    ap.add_argument('--retries', type=int, default=2,
                    help='re-run a job that died on a remote/network error [2]')
    ap.add_argument('--dryRun', action='store_true', help='write job dirs, print commands, run nothing')
    args = parseWithConfig(ap)

    os.makedirs(args.work, exist_ok=True)
    with open(f'{args.work}/run.yaml', 'w') as fp:
        yaml.safe_dump({k: v for k, v in vars(args).items() if k != 'config'}, fp, sort_keys=False)
    summary = open(f'{args.work}/summary.log', 'a')
    feats = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    if args.tiles:
        feats = [f for f in feats if f['properties']['name'] in args.tiles]
    if args.excludeTiles:
        feats = [f for f in feats if f['properties']['name'] not in args.excludeTiles]
    import csv
    jobs, grids, bands, caps = [], {}, {}, []
    for f in feats:
        name, props = f['properties']['name'], f['properties']
        grid = tileGrid(props, f['geometry'], args.res, args.psResM, args.featherKm)
        rows = list(csv.DictReader(open(f'{args.tileRun}/tiles/{name}.csv')))
        for y in sorted(glob.glob(f'{args.tileRun}/tiles/{name}.*.yaml')):
            grp = y.split('.')[-2]                       # e.g. AHH, BHH, BVV
            jobDir = f'{args.work}/jobs/{name}.{grp}'
            writeJob(jobDir, grid, y, rows, args.granules, args.factorFrom)
            jobs.append((name, grp, jobDir))
            grids[jobDir] = grid
        if grid['epsg'] == 4326:
            s = grid['tile'][1]
            latBand = (math.floor(s / 6) * 6, math.floor(s / 6) * 6 + props['height'])
            bands.setdefault(latBand, dict(dLon=grid['dx'], tiles=[]))['tiles'].append(name)
        else:
            caps.append(name)
    products = ['sigma0', 'gamma0'] if args.calOutput == 'both' else [args.calOutput]
    cmdBase = [args.geomosaic, '-GTiff', '-removePad', str(args.removePad), '-S1Cal',
               '-int16', '-calOutput', args.calOutput, '-ompThreads', '1']
    if args.date1 and args.date2:
        cmdBase += ['-date1', args.date1, '-date2', args.date2]
    jobCmd = lambda jobDir: cmdBase + ['-fl', str(grids[jobDir]['fl']), '-epsg', str(grids[jobDir]['epsg'])]
    log(f'{len(feats)} tiles ({len(caps)} polar stereographic caps), {len(jobs)} jobs, {args.nProc} '
        f'processes, res {args.res * 3600:.2f} arcsec / caps {args.psResM:.0f} m, {args.calOutput}; '
        f'granules {args.granules}', summary)
    if args.dryRun:
        for name, grp, jobDir in jobs[:3]:
            log(f'  {jobDir}: {" ".join(jobCmd(jobDir))} -gcov gcov.yaml inputFile none out', summary)
            log(f'    inputFile: {open(f"{jobDir}/inputFile").readline().strip()}', summary)
        return 0

    failed = []
    s3Keys = S3Keys() if args.granules == 's3' else None
    with concurrent.futures.ThreadPoolExecutor(args.nProc) as pool:
        for jobDir, rc, sec, note in pool.map(lambda j: runJob(j[2], jobCmd(j[2]), products, args.retries, s3Keys), jobs):
            ok = rc == 0 and not note.startswith('missing')
            status = 'OK' if ok else f'FAILED rc={rc}'
            log(f'{os.path.basename(jobDir)}: {sec / 60:.1f} min {status} {note}', summary)
            if status != 'OK':
                failed.append(jobDir)
                # the reason: last error-looking lines of the job's geomosaic log
                errs = [l.strip() for l in open(f'{jobDir}/log', errors='ignore')
                        if any(k in l.lower() for k in ('error', 'fail', 'denied', 'http'))]
                for l in errs[-3:]:
                    log(f'    {l[:300]}', summary)
    if len(failed) == len(jobs):
        log(f'all {len(jobs)} jobs failed; no VRTs built (see {args.work}/jobs/*/log)', summary)
        return 1

    # crop, then the VRT hierarchy: tile (B under A) -> latitude band -> global
    for product in products:
        tileVrts = {}
        for name, grp, jobDir in jobs:
            if jobDir in failed:
                continue
            dst = f'{args.work}/tiles/{product}/{name}.{grp}.tif'
            if not os.path.exists(dst):
                cropToTile(f'{jobDir}/out.{product}.tif', dst, grids[jobDir])
            tileVrts.setdefault(name, []).append((grp, dst))
        os.makedirs(f'{args.work}/vrt/{product}', exist_ok=True)
        for name, parts in tileVrts.items():
            # VRT sources are painted in order: frequency B first, frequency A on top
            parts = sorted(parts, key=lambda p: p[0][0] != 'B')
            gdal.BuildVRT(f'{args.work}/vrt/{product}/{name}.vrt', [p for _, p in parts],
                          srcNodata=NODATA, VRTNodata=NODATA)
        bandVrts = []
        for (lat0, lat1), b in sorted(bands.items()):
            srcs = [f'{args.work}/vrt/{product}/{t}.vrt' for t in b['tiles'] if t in tileVrts]
            if not srcs:
                continue
            out = f'{args.work}/vrt/{product}/band_{lat0:+03d}_{lat1:+03d}.vrt'
            gdal.BuildVRT(out, srcs, resolution='user', xRes=b['dLon'], yRes=args.res,
                          srcNodata=NODATA, VRTNodata=NODATA)
            bandVrts.append(out)
        # polar caps are in their own CRS: separate VRTs, not part of the lat/lon global VRT
        for name in caps:
            if name in tileVrts:
                gdal.BuildVRT(f'{args.work}/vrt/{product}/cap_{name}.vrt',
                              [f'{args.work}/vrt/{product}/{name}.vrt'], srcNodata=NODATA, VRTNodata=NODATA)
        gdal.BuildVRT(f'{args.work}/vrt/{product}/global.vrt', bandVrts, resolution='user',
                      xRes=args.res, yRes=args.res, resampleAlg='nearest',
                      srcNodata=NODATA, VRTNodata=NODATA)
        log(f'{product}: {len(tileVrts)} tile VRTs, {len(bandVrts)} band VRTs, global.vrt', summary)
    log(f'done; {len(failed)} failed job(s){": " + ", ".join(os.path.basename(f) for f in failed) if failed else ""}',
        summary)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
