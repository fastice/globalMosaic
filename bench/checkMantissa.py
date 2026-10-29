#!/usr/bin/env python3
"""
Across every NISAR product type: is `truncate_mantissa` applied, and what would it save?

ISCE3 exposes `data_type: complex64_zero_mantissa` only in the `focus` (RSLC) and `gslc`
workflows, even though `isce3.core.types.truncate_mantissa()` takes real OR complex arrays and
reads the mantissa width from the dtype. This checks, from the data, which products actually
carry truncated mantissas and measures what the rest would gain.

Method per product: resolve one granule, open its header over HTTP (no download), list the
float32/complex64 raster layers with their on-disk sizes, then on the largest layer read a
window that actually contains data and (1) test whether the low 13 mantissa bits are in use,
(2) recompress that window as-is, zero-mantissa, and true fp16 under the product's own filter
settings.

    checkMantissa.py                 # all products
    checkMantissa.py --levels GCOV RUNW --json out.json
"""
import argparse
import json
import os
import sys
import time

LEVELS = ['RSLC', 'GSLC', 'GCOV', 'RIFG', 'RUNW', 'ROFF', 'GUNW', 'GOFF']
MASK32 = 0xFFFFE000          # keep sign + exp8 + top 10 mantissa bits


def resolve(url, tries=5):
    import requests
    for i in range(tries):
        try:
            s = requests.Session()
            r = s.get(url, allow_redirects=True, stream=True,
                      headers={'Range': 'bytes=0-0'}, timeout=180)
            u = r.url
            r.close()
            return u
        except Exception:
            time.sleep(6)
    return None


def openRemote(signed, tries=4):
    import h5py
    import fsspec
    for i in range(tries):
        try:
            return h5py.File(fsspec.open(signed, 'rb').open(), 'r')
        except Exception:
            time.sleep(6)
    return None


def rasterLayers(f, minPx=1e6):
    """2-D float32 or complex64 layers, biggest first, with on-disk bytes."""
    import h5py
    out = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Dataset) or obj.ndim != 2 or obj.size < minPx:
            return
        t = obj.id.get_type()
        try:
            cplx = t.get_nmembers() == 2
        except Exception:
            cplx = False
        if obj.dtype.kind == 'f' or cplx:
            out.append(dict(name=name, shape=obj.shape, cplx=cplx,
                            disk=obj.id.get_storage_size(), px=obj.size,
                            filt=obj.compression, lvl=obj.compression_opts,
                            shuffle=obj.shuffle))
    f.visititems(visit)
    out.sort(key=lambda d: -d['disk'])
    return out


def findWindow(d, cplx, n=1024, probes=8):
    """A window with real data in it -- most of these grids are mostly nodata."""
    import numpy as np
    H, W = d.shape
    best = None
    for fy in (0.5, 0.45, 0.55, 0.35, 0.65, 0.3, 0.7, 0.5):
        for fx in (0.5, 0.4, 0.6, 0.45):
            y = min(max(0, int(H * fy)), max(0, H - n))
            x = min(max(0, int(W * fx)), max(0, W - n))
            a = d[y:y + n, x:x + n]
            v = a.real if cplx else a
            good = np.isfinite(v) & (v != 0)
            if good.mean() > 0.5:
                return a, y, x, good.mean()
            if best is None or good.mean() > best[3]:
                best = (a, y, x, good.mean())
            probes -= 1
            if probes <= 0:
                return best
    return best


def store(arr, lvl, shuffle, tmp='/tmp/claude-3707/cm.h5'):
    import h5py
    if os.path.exists(tmp):
        os.unlink(tmp)
    ch = (min(512, arr.shape[0]), min(512, arr.shape[1]))
    with h5py.File(tmp, 'w') as g:
        g.create_dataset('d', data=arr, chunks=ch,
                         compression=('gzip' if lvl else None),
                         compression_opts=(lvl or None), shuffle=bool(shuffle))
    b = os.path.getsize(tmp)
    os.unlink(tmp)
    return b


def analyse(level, verbose=True):
    import asf_search as asf
    import numpy as np
    res = dict(level=level)
    try:
        r = asf.search(dataset='NISAR', processingLevel=level, maxResults=10)
        gl = list(r)
    except Exception as e:
        res['error'] = 'search failed: %s' % type(e).__name__
        return res
    if not gl:
        res['error'] = 'no granules'
        return res
    g = gl[0]
    res['granule'] = g.properties['fileName']
    signed = resolve(g.properties['url'])
    if signed is None:
        res['error'] = 'could not resolve URL'
        return res
    f = openRemote(signed)
    if f is None:
        res['error'] = 'could not open remotely'
        return res
    try:
        lay = rasterLayers(f)
        if not lay:
            res['error'] = 'no raster layers'
            return res
        res['totalDisk'] = sum(l['disk'] for l in lay)
        top = lay[0]
        res['layer'] = top['name'].split('/')[-1]
        res['layerDisk'] = top['disk']
        res['cplx'] = top['cplx']
        res['bitsPerPx'] = top['disk'] * 8 / top['px'] * (1 if not top['cplx'] else 1)
        res['filt'] = '%s-%s shuffle=%s' % (top['filt'], top['lvl'], top['shuffle'])
        d = f[top['name']]
        win = findWindow(d, top['cplx'])
        if win is None:
            res['error'] = 'no window'
            return res
        a, y, x, frac = win
        res['validFrac'] = float(frac)
        a = a.astype('complex64') if top['cplx'] else a.astype('float32')
        comp = (np.concatenate([a.real.ravel(), a.imag.ravel()]) if top['cplx']
                else a.ravel()).astype('float32')
        nz = comp[np.isfinite(comp) & (comp != 0)]
        if nz.size < 1000:
            res['error'] = 'window has no data'
            return res
        live = (((np.abs(nz).view('uint32')) & 0x7FFFFF) & 0x1FFF) != 0
        res['liveFrac'] = float(live.mean())
        res['applied'] = bool(live.mean() < 0.5)
        lvl = top['lvl'] if top['filt'] == 'gzip' else 1
        b0 = store(a, lvl, top['shuffle'])
        if top['cplx']:
            zr = (a.real.view(np.uint32) & np.uint32(MASK32)).view(np.float32)
            zi = (a.imag.view(np.uint32) & np.uint32(MASK32)).view(np.float32)
            zm = (zr + 1j * zi).astype('complex64')
            f16 = np.empty(a.shape, dtype=np.dtype([('r', '<f2'), ('i', '<f2')]))
            f16['r'] = a.real.astype('float16')
            f16['i'] = a.imag.astype('float16')
        else:
            zm = (a.view(np.uint32) & np.uint32(MASK32)).view(np.float32)
            f16 = a.astype('float16')
        bz = store(zm, lvl, top['shuffle'])
        bh = store(f16, lvl, top['shuffle'])
        res['ratioZero'] = b0 / bz
        res['ratioF16'] = b0 / bh
        res['saveBytes'] = int(top['disk'] * (1 - bz / b0)) if not res['applied'] else 0
    except Exception as e:
        res['error'] = '%s: %s' % (type(e).__name__, str(e)[:80])
    finally:
        try:
            f.close()
        except Exception:
            pass
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--levels', nargs='+', default=LEVELS)
    ap.add_argument('--json', default=None)
    args = ap.parse_args()
    rows = []
    for lvl in args.levels:
        print('--- %s ---' % lvl, flush=True)
        r = analyse(lvl)
        rows.append(r)
        if 'error' in r:
            print('    %s' % r['error'], flush=True)
        else:
            print('    %-26s %8.1f MB  live %6.2f%%  %s  zero-mant %.2fx  fp16 %.2fx'
                  % (r['layer'], r['layerDisk'] / 1e6, 100 * r['liveFrac'],
                     'APPLIED ' if r['applied'] else 'not appl', r['ratioZero'], r['ratioF16']),
                  flush=True)
    print('\n%-7s %-24s %-9s %10s %8s %9s %9s %10s'
          % ('LEVEL', 'largest layer', 'status', 'layer MB', 'live %', 'zero-m', 'fp16', 'save MB'))
    for r in rows:
        if 'error' in r:
            print('%-7s %s' % (r['level'], r['error']))
            continue
        print('%-7s %-24s %-9s %10.1f %7.2f%% %8.2fx %8.2fx %10.1f'
              % (r['level'], r['layer'][:24],
                 'APPLIED' if r['applied'] else 'NOT appl',
                 r['layerDisk'] / 1e6, 100 * r['liveFrac'],
                 r['ratioZero'], r['ratioF16'], r['saveBytes'] / 1e6))
    if args.json:
        json.dump(rows, open(args.json, 'w'), indent=2)
        print('\nwrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
