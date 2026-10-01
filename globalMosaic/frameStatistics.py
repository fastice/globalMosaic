#!/usr/bin/env python3
'''
Temporal backscatter statistics per (track, frame) from a stack of NISAR GCOV cycles.

  frameStatistics cycle30/ascending --catalogues cycle30/catalogue0{23..31}.geojson ... \\
      --out /scratch/ianj/mosaics/stats30asc/frames [--frames 166_138_4005 ...] [--nProc 8]

Frames: every (track, frame, mode) among the granules the tiling selected, so the statistics cover
what the backscatter mosaic covers; direction and channel as the tiling chose them (file names).
For each frame, every cycle's granule(s) of that track/frame/mode is read over https -- only the
chunks of the power layer used by the tiling (HHHH, or VVVV) and of the mask (rangeReader: a 40 MHz
granule costs 1.7 GB instead of the 5.4 GB GDAL's read-ahead fetched). Valid samples
(mask 1..254, finite, > 0) are averaged in LINEAR power over --resM blocks laid on multiples of
--resM in the frame's CRS, so every cycle's blocks coincide (GCOV origins are snapped to the
posting). Granules of one cycle (a partial plus a full frame) are merged into that cycle's one
sample. Across cycles: sum, sum of squares, count per cell.

Output <out>/<track>_<frame>_<mode>.tif (co-pol: HH, or VV) and, for dual- and quad-pol frames,
<out>/<track>_<frame>_<mode>.cross.tif (cross-pol: HV for dual-H and quad, VH for dual-V; one read
of each granule serves both, sharing the mask), the frame's CRS at --resM, Float32, NaN = no data:
  1 mean gamma0 (linear)   2 sigma (sample std, linear)   3 CV = sigma / mean (raw: it includes
  the speckle floor, about 1/sqrt(looks per cell))   4 n, the number of cycles with data
Bands 1-3 are NaN where n < --minCount; metadata 'polarization' (e.g. HVHV). Per-cycle arrays live in
memory only. Resumable per layer (a frame with its co-pol file only reads the cross-pol layer); <out>/frames.log records each frame's cycles and time.
'''
import argparse
import collections
import concurrent.futures
import csv
import glob
import os
import sys
import time

import numpy as np
from osgeo import gdal, osr

from .globalGCOVTiles import loadGranules

gdal.UseExceptions()
BANDS = ['mean gamma0 (linear)', 'sigma (linear)', 'CV', 'n (cycles)']


def httpSetup(cookieDir):
    cookies = f'{cookieDir}/cookies.{os.getpid()}'
    for k, v in (('GDAL_HTTP_NETRC', 'YES'), ('GDAL_HTTP_COOKIEFILE', cookies),
                 ('GDAL_HTTP_COOKIEJAR', cookies), ('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR'),
                 ('GDAL_HTTP_MAX_RETRY', '5'), ('GDAL_HTTP_RETRY_DELAY', '10'),
                 # big reads: bandwidth-bound instead of request-latency-bound
                 ('CPL_VSIL_CURL_CHUNK_SIZE', str(10 * 1024 * 1024)),     # GDAL maximum; larger falls back to 16 KB
                 ('CPL_VSIL_CURL_CACHE_SIZE', str(1024 * 1024 * 1024)),
                 ('GDAL_HTTP_MULTIRANGE', 'YES'), ('GDAL_HTTP_MERGE_CONSECUTIVE_RANGES', 'YES')):
        gdal.SetConfigOption(k, v)


def frameKey(name):
    f = name.split('_')
    return f'{f[5]}_{f[7]}_{f[8]}'             # track, frame, mode


def reduceGranule(url, freq, chans, resM):
    ''' One granule reduced to --resM blocks on the lattice of multiples of resM, for each power
    layer in chans (e.g. ['HH', 'HV'], read with one shared mask; a layer the granule lacks is left
    out): (column of the first block on the lattice, row of the first block, epsg,
    {chan: (mean power sums, valid native pixels per block)}, native spacing). Reads only those
    layers' and the mask's chunks (rangeReader), not the rest of the file. '''
    import h5py
    from .rangeReader import RangeFile, chunkIndex, rangesFor
    rf = RangeFile(url)
    g = h5py.File(rf, 'r')[f'/science/LSAR/GCOV/grids/frequency{freq}']
    xs, ys = g['xCoordinates'][()], g['yCoordinates'][()]
    epsg = int(np.asarray(g['projection'].attrs['epsg_code']).ravel()[0])
    dx, dy = float(xs[1] - xs[0]), float(ys[1] - ys[0])      # dy < 0
    k = int(round(resM / abs(dx)))
    if abs(k * abs(dx) - resM) > 1e-6 * resM:
        raise ValueError(f'{resM} m is not a whole number of {abs(dx)} m pixels')
    left, top = xs[0] - dx / 2, ys[0] - dy / 2                # pixel edges
    # first block edge on the lattice at or after the first pixel edge
    bx = int(np.ceil(round(left / resM, 6)))
    by = int(np.floor(round(top / resM, 6)))
    c0 = int(round((bx * resM - left) / dx))
    r0 = int(round((top - by * resM) / abs(dy)))
    nbx, nby = (len(xs) - c0) // k, (len(ys) - r0) // k
    if nbx < 1 or nby < 1:
        return None
    P = {c: g[f'{c}{c}'] for c in chans if f'{c}{c}' in g}
    if not P:
        return None
    M = g['mask']
    idx = {c: chunkIndex(a) for c, a in P.items()}
    iM = chunkIndex(M)
    out = {c: (np.zeros((nby, nbx)), np.zeros((nby, nbx), 'i4')) for c in P}
    step = max(1, 2048 // k)                                  # block rows per full-width read
    for b0 in range(0, nby, step):
        b1 = min(b0 + step, nby)
        ra, rb, ca, cb = r0 + b0 * k, r0 + b1 * k, c0, c0 + nbx * k
        ranges = rangesFor(iM, M.chunks, ra, rb, ca, cb)
        for c, a in P.items():
            ranges += rangesFor(idx[c], a.chunks, ra, rb, ca, cb)
        rf.prefetch(ranges)
        m = M[ra:rb, ca:cb]
        inMask = (m > 0) & (m < 255)
        nb = b1 - b0
        for c, a in P.items():
            p = a[ra:rb, ca:cb]
            ok = inMask & np.isfinite(p) & (p > 0)
            s, w = out[c]
            s[b0:b1] = np.where(ok, p, 0.).reshape(nb, k, nbx, k).sum((1, 3), dtype='f8')
            w[b0:b1] = ok.reshape(nb, k, nbx, k).sum((1, 3))
        rf.drop()
    return bx, by, epsg, out, abs(dx)


def writeStats(out, samples, resM, minCount, chan):
    ''' Temporal statistics of one layer from its per-cycle samples [(cycle, [(bx, by, (s, w))])]
    to <out> (4 bands, see the module doc). Returns (cycles, cells with n >= minCount), or None when
    no cell has valid data. '''
    epsg0, pix = samples[0][2], samples[0][3]
    # frame lattice extent = union of every part
    bx0 = min(p[0] for _, ps, _, _ in samples for p in ps)
    by0 = max(p[1] for _, ps, _, _ in samples for p in ps)
    bx1 = max(p[0] + p[2][0].shape[1] for _, ps, _, _ in samples for p in ps)
    by1 = min(p[1] - p[2][0].shape[0] for _, ps, _, _ in samples for p in ps)
    nx, ny = bx1 - bx0, by0 - by1
    S = np.zeros((ny, nx))
    S2 = np.zeros((ny, nx))
    N = np.zeros((ny, nx), 'i2')
    for cycle, parts, _, _ in samples:
        s = np.zeros((ny, nx))
        w = np.zeros((ny, nx))
        for bx, by, (ps, pw) in parts:
            i, j = by0 - by, bx - bx0
            s[i:i + ps.shape[0], j:j + ps.shape[1]] += ps
            w[i:i + ps.shape[0], j:j + ps.shape[1]] += pw
        ok = w > 0
        v = np.where(ok, s / np.maximum(w, 1), 0.)
        S += v
        S2 += v * v
        N += ok
    if not N.any():
        return None
    good = N >= minCount
    n = np.maximum(N, 1)
    mean = S / n
    var = np.maximum(S2 - n * mean * mean, 0.) / np.maximum(n - 1, 1)
    sigma = np.sqrt(var)
    cv = np.where(mean > 0, sigma / np.where(mean > 0, mean, 1), np.nan)
    stack = [np.where(good, a, np.nan).astype('f4') for a in (mean, sigma, cv)]
    stack.append(np.where(N > 0, N, np.nan).astype('f4'))
    # crop to the cells with data (cycles' grids differ in extent; the union is mostly empty)
    rr, cc = np.nonzero((N > 0).any(1))[0], np.nonzero((N > 0).any(0))[0]
    stack = [a[rr[0]:rr[-1] + 1, cc[0]:cc[-1] + 1] for a in stack]
    by0, bx0 = by0 - rr[0], bx0 + cc[0]
    ny, nx = stack[0].shape
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg0)
    tmp = out + '.tmp.tif'
    d = gdal.GetDriverByName('GTiff').Create(tmp, nx, ny, len(BANDS), gdal.GDT_Float32,
                                             options=['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=3'])
    d.SetGeoTransform((bx0 * resM, resM, 0., by0 * resM, 0., -resM))
    d.SetProjection(srs.ExportToWkt())
    for b, (a, desc) in enumerate(zip(stack, BANDS), 1):
        band = d.GetRasterBand(b)
        band.SetNoDataValue(float('nan'))
        band.SetDescription(desc)
        band.WriteArray(a)
    cycles = [c for c, _, _, _ in samples]
    d.SetMetadata({'cycles': ' '.join(cycles), 'native_pixel_m': str(pix), 'minCount': str(minCount),
                   'polarization': f'{chan}{chan}'})
    d = None
    os.replace(tmp, out)
    return cycles, int(good.sum())


def frameJob(job):
    ''' One frame: per layer (co-pol HH/VV -> <key>.tif, cross-pol HV/VH -> <key>.cross.tif) the
    statistics over cycles; only layers without an output yet are read. '''
    key, byCycle, outs, chans, resM, minCount, cookieDir = job
    todo = {L: c for L, c in chans.items() if c and not os.path.exists(outs[L])}
    if not todo:
        return key, 0., 'skipped (done)'
    httpSetup(cookieDir)
    t0 = time.time()
    samples = {L: [] for L in todo}
    failed, epsg0 = [], None
    for cycle in sorted(byCycle):
        # one sample per cycle: merge that cycle's granules (block sums and counts add up)
        parts = []
        for name, url, freq in byCycle[cycle]:
            for attempt in range(3):
                try:
                    r = reduceGranule(url, freq, list(todo.values()), resM)
                    break
                except Exception as e:
                    err = str(e)[:150]
                    time.sleep(15)
            else:
                failed.append(f'{name[17:47]}: {err}')
                continue
            if r is not None:
                parts.append(r)
        if not parts:
            continue
        epsg0 = epsg0 or parts[0][2]
        parts = [p for p in parts if p[2] == epsg0]
        for L, c in todo.items():
            ps = [(bx, by, layers[c]) for bx, by, _, layers, _ in parts if c in layers]
            if ps:
                samples[L].append((cycle, ps, epsg0, parts[0][4]))
    notes = []
    for L, c in todo.items():
        if not samples[L]:
            notes.append(f'{L} {c}{c}: no data')
            continue
        r = writeStats(outs[L], samples[L], resM, minCount, c)
        if r is None:
            notes.append(f'{L} {c}{c}: no valid data ({len(samples[L])} cycles read, all masked)')
        else:
            notes.append(f'{L} {c}{c}: {len(r[0])} cycles ({" ".join(r[0])}), {r[1]} cells n>={minCount}')
    note = '; '.join(notes)
    if failed:
        note += f'; {len(failed)} FAILED: ' + '; '.join(failed[:2])
    return key, time.time() - t0, note


def safeFrameJob(job):
    ''' frameJob, but one bad frame is logged as FAILED instead of stopping the whole run. '''
    try:
        return frameJob(job)
    except Exception as e:
        return job[0], 0., f'FAILED {type(e).__name__}: {str(e)[:150]}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', help='the backscatter tiling whose frames to use (e.g. cycle30/ascending)')
    ap.add_argument('--catalogues', nargs='+', required=True, help='GCOV catalogues of the cycles to stack')
    ap.add_argument('--out', required=True, help='directory for the per-frame statistics files')
    ap.add_argument('--resM', type=float, default=80., help='statistics cell size, m [80]')
    ap.add_argument('--minCount', type=int, default=3, help='minimum cycles for statistics [3]')
    ap.add_argument('--frames', nargs='+', default=None, help='only these track_frame_mode keys')
    ap.add_argument('--nProc', type=int, default=8)
    ap.add_argument('--layers', nargs='+', choices=['co', 'cross'], default=['co', 'cross'],
                    help='co-pol (HH, or VV) and/or cross-pol (HV for dual-H/quad, VH for dual-V) [both]')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    # frames (and the channel) the tiling uses
    want = {}
    for p in glob.glob(f'{args.tileRun}/tiles/*.csv'):
        for r in csv.DictReader(open(p)):
            want.setdefault(frameKey(r['name']), (r['name'].split('_')[6], r['freq'], r['pol']))
    if args.frames:
        want = {k: v for k, v in want.items() if k in args.frames}
    # every cycle's granules of those frames, same direction, latest version per scene
    stacks = collections.defaultdict(lambda: collections.defaultdict(list))
    for cat in args.catalogues:
        for g in loadGranules(cat, allowV=True):
            k = frameKey(g['name'])
            if k in want and g['direction'] == want[k][0]:
                stacks[k][g['cycle']].append((g['name'], g['url'], want[k][1]))
    jobs = []
    for k in sorted(stacks):
        pol = want[k][2]
        chans = {'co': 'HH' if pol in ('SH', 'DH', 'QP') else 'VV',
                 'cross': {'DH': 'HV', 'QP': 'HV', 'DV': 'VH'}.get(pol) if 'cross' in args.layers else None}
        if 'co' not in args.layers:
            chans['co'] = None
        jobs.append((k, dict(stacks[k]), {'co': f'{args.out}/{k}.tif', 'cross': f'{args.out}/{k}.cross.tif'},
                     chans, args.resM, args.minCount, args.out))
    jobs.sort(key=lambda j: -sum(len(v) for v in j[1].values()))
    print(f'{len(jobs)} frames, {sum(sum(len(v) for v in j[1].values()) for j in jobs)} granules, '
          f'{args.nProc} processes -> {args.out}', flush=True)
    log = open(f'{args.out}/frames.log', 'a')
    with concurrent.futures.ProcessPoolExecutor(args.nProc) as pool:
        for i, (k, sec, note) in enumerate(pool.map(safeFrameJob, jobs), 1):
            line = f'{time.strftime("%m-%d %H:%M:%S")} [{i}/{len(jobs)}] {k}: {sec / 60:.1f} min {note}'
            print(line, flush=True)
            log.write(line + '\n')
            log.flush()
    for c in glob.glob(f'{args.out}/cookies.*'):
        os.remove(c)
    return 0


if __name__ == '__main__':
    sys.exit(main())
