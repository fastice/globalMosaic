"""
Measure the fastest way to read NISAR GCOV HDF5 from ASF with GDAL.

Run this IN us-west-2. Numbers from outside the region are dominated by round-trip latency
(~2 s per request against ~1 ms in region) and rank the options differently -- measured from
Seattle a 4096x4096 window took 11.1 s (1.5 Mpx/s) against 0.33 s (51.5 Mpx/s) off local disk,
and essentially all of that gap was latency, not bandwidth.

Every variant runs in a FRESH SUBPROCESS: /vsicurl keeps a per-process block cache, so a second
read in the same process measures the cache, not the network.
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import time

BUCKET = 's3://sds-n-cumulus-prod-nisar-products'
PREFIX = 'NISAR_L2_GCOV_PROVISIONAL_V1'
BAND = 'science/LSAR/GCOV/grids/frequencyA/HHHH'

# name -> (extra GDAL config, strips). One row per thing worth knowing.
VARIANTS = {
    'default':            ({}, 1),
    'chunk1M':            ({'CPL_VSIL_CURL_CHUNK_SIZE': '1048576'}, 1),
    'chunk4M':            ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304'}, 1),
    'chunk16M':           ({'CPL_VSIL_CURL_CHUNK_SIZE': '16777216'}, 1),
    'chunk4M+http2':      ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304',
                            'GDAL_HTTP_VERSION': '2', 'GDAL_HTTP_MULTIPLEX': 'YES'}, 1),
    'chunk4M+vsicache':   ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304',
                            'VSI_CACHE': 'TRUE', 'VSI_CACHE_SIZE': '536870912'}, 1),
    'chunk4M+threads':    ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304',
                            'GDAL_NUM_THREADS': 'ALL_CPUS'}, 1),
    'chunk4M+8strips':    ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304'}, 8),
    'chunk4M+64strips':   ({'CPL_VSIL_CURL_CHUNK_SIZE': '4194304'}, 64),
}

WORKER = r'''
import os, sys, time, json
from osgeo import gdal
gdal.UseExceptions()
path, xoff, yoff, w, h, strips = sys.argv[1], *map(int, sys.argv[2:7])
t0 = time.time()
ds = gdal.Open(path)
tOpen = time.time() - t0
band = ds.GetRasterBand(1)
t1 = time.time()
n = 0
rows = max(1, h // strips)
y = yoff
while y < yoff + h:
    nr = min(rows, yoff + h - y)
    n += band.ReadAsArray(xoff, y, w, nr).size
    y += nr
tRead = time.time() - t1
print(json.dumps({'open': tOpen, 'read': tRead, 'mpx': n / 1e6}))
'''


def gdalPath(granule, mode):
    """Build the GDAL subdataset path for a granule stem, over s3 or https."""
    if mode == 's3':
        loc = '/vsis3/%s/%s/%s/%s.h5' % (BUCKET.replace('s3://', ''), PREFIX, granule, granule)
    else:
        loc = '/vsicurl/%s' % resolveSigned(granule)
    return 'HDF5:"%s"://%s' % (loc, BAND)


def resolveSigned(granule):
    """Follow the ASF -> Earthdata -> CloudFront chain and keep the final signed URL."""
    import requests
    url = ('https://nisar.asf.earthdatacloud.nasa.gov/NISAR/%s/%s/%s.h5'
           % (PREFIX, granule, granule))
    r = requests.get(url, allow_redirects=True, stream=True,
                     headers={'Range': 'bytes=0-0'}, timeout=120)
    r.close()
    return r.url


def s3Credentials():
    """ASF in-region credentials; they last about an hour."""
    import earthaccess
    earthaccess.login(strategy='netrc')
    c = earthaccess.get_s3_credentials(daac='ASF')
    return {'AWS_ACCESS_KEY_ID': c['accessKeyId'],
            'AWS_SECRET_ACCESS_KEY': c['secretAccessKey'],
            'AWS_SESSION_TOKEN': c['sessionToken'],
            'AWS_REGION': 'us-west-2'}


def runOne(path, cfg, strips, win, creds):
    env = dict(os.environ)
    env['GDAL_DISABLE_READDIR_ON_OPEN'] = 'EMPTY_DIR'
    env.update(creds)
    env.update(cfg)
    cmd = [sys.executable, '-c', WORKER, path] + [str(v) for v in win] + [str(strips)]
    t0 = time.time()
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if p.returncode != 0:
        return None, (p.stderr.strip().splitlines() or ['failed'])[-1]
    return json.loads(p.stdout.strip().splitlines()[-1]), None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('granule', help='granule stem (no .h5)')
    ap.add_argument('--mode', choices=['s3', 'https', 'both'], default='both',
                    help='access path; s3 only works in us-west-2 [both]')
    ap.add_argument('--window', type=int, nargs=4, default=[8000, 12000, 4096, 4096],
                    metavar=('XOFF', 'YOFF', 'W', 'H'))
    ap.add_argument('--repeat', type=int, default=2, help='runs per variant, best is kept [2]')
    ap.add_argument('--variants', nargs='+', default=None, help='subset of variant names')
    ap.add_argument('--nProc', type=int, default=0,
                    help='also measure N concurrent readers with the best variant [off]')
    args = ap.parse_args()

    names = args.variants or list(VARIANTS)
    modes = ['s3', 'https'] if args.mode == 'both' else [args.mode]
    mpx = args.window[2] * args.window[3] / 1e6

    for mode in modes:
        creds = {}
        if mode == 's3':
            try:
                creds = s3Credentials()
            except Exception as e:
                print('s3: cannot get credentials (%s) -- skipping' % type(e).__name__)
                continue
        try:
            path = gdalPath(args.granule, mode)
        except Exception as e:
            print('%s: cannot build path (%s) -- skipping' % (mode, e))
            continue

        print('\n=== %s ===  window %s = %.1f Mpx' % (mode, args.window, mpx))
        # warm the CDN/S3 once so the first variant is not penalised
        runOne(path, {'CPL_VSIL_CURL_CHUNK_SIZE': '4194304'}, 1, args.window, creds)
        rows = []
        for name in names:
            cfg, strips = VARIANTS[name]
            best, err = None, None
            for _ in range(args.repeat):
                r, e = runOne(path, cfg, strips, args.window, creds)
                if r is None:
                    err = e
                    break
                if best is None or r['read'] < best['read']:
                    best = r
            if best is None:
                print('  %-20s FAILED: %s' % (name, err))
                continue
            rows.append((best['read'], name, best['open']))
            print('  %-20s open %6.2f s  read %7.2f s  %8.2f Mpx/s'
                  % (name, best['open'], best['read'], mpx / best['read']))
        if rows:
            rows.sort()
            print('  fastest: %s (%.2f Mpx/s)' % (rows[0][1], mpx / rows[0][0]))
            if args.nProc:
                cfg, strips = VARIANTS[rows[0][1]]
                print('  --- %d concurrent readers, variant %s ---' % (args.nProc, rows[0][1]))
                t0 = time.time()
                procs = [subprocess.Popen(
                    [sys.executable, '-c', WORKER, path] +
                    [str(v) for v in args.window] + [str(strips)],
                    env={**os.environ, 'GDAL_DISABLE_READDIR_ON_OPEN': 'EMPTY_DIR',
                         **creds, **cfg},
                    stdout=subprocess.DEVNULL) for _ in range(args.nProc)]
                for p in procs:
                    p.wait()
                el = time.time() - t0
                print('  aggregate %.2f Mpx/s over %d processes (%.2f s wall)'
                      % (mpx * args.nProc / el, args.nProc, el))
    return 0


if __name__ == '__main__':
    sys.exit(main())
