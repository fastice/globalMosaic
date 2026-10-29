#!/usr/bin/env python3
"""
Savings from zero-mantissa truncation across one full orbit frame's product set.

Format-preserving only: every raster keeps its declared dtype (float32 / complex64) and its
existing filters. The only change is zeroing low mantissa bits before compression, which is
what ISCE3's `truncate_mantissa()` already does for RSLC and GSLC.

For each product, every float32/complex64 raster is examined: its on-disk bytes are read from
HDF5, a window with real data in it is recompressed under the layer's own filter settings at
each mantissa width, and the ratio is projected onto the whole layer.

    frameSavings.py --track 131 --dir D --frame 062 --bits 6 10
"""
import argparse
import json
import os
import re
import sys
import time

MANT = 23


def mask(nbits):
    import numpy as np
    return np.uint32(((0xFFFFFFFF << (MANT - nbits)) & 0xFFFFFFFF))


def resolve(url, tries=4):
    import requests
    for _ in range(tries):
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


def openRemote(signed, tries=3):
    import h5py
    import fsspec
    for _ in range(tries):
        try:
            return h5py.File(fsspec.open(signed, 'rb').open(), 'r')
        except Exception:
            time.sleep(6)
    return None


def store(arr, lvl, shuffle):
    import h5py
    p = '/tmp/claude-3707/fs.h5'
    if os.path.exists(p):
        os.unlink(p)
    ch = (min(512, arr.shape[0]), min(512, arr.shape[1]))
    with h5py.File(p, 'w') as g:
        g.create_dataset('d', data=arr, chunks=ch,
                         compression=('gzip' if lvl else None),
                         compression_opts=(lvl or None), shuffle=bool(shuffle))
    b = os.path.getsize(p)
    os.unlink(p)
    return b


def nodataWindow(d, cplx, n=1024):
    """A window that is mostly NOT data -- needed to price the nodata separately."""
    import numpy as np
    H, W = d.shape
    n = min(n, H, W)
    worst = None
    for (fy, fx) in ((0.02, 0.02), (0.02, 0.9), (0.9, 0.02), (0.9, 0.9), (0.05, 0.5)):
        y = min(max(0, int(H * fy)), H - n)
        x = min(max(0, int(W * fx)), W - n)
        a = d[y:y + n, x:x + n]
        v = a.real if cplx else a
        good = np.isfinite(v) & (v != 0)
        if good.mean() < 0.02:
            return a, good.mean()
        if worst is None or good.mean() < worst[1]:
            worst = (a, good.mean())
    return worst


def window(d, cplx, n=1024):
    import numpy as np
    H, W = d.shape
    n = min(n, H, W)
    best = None
    for fy in (0.5, 0.4, 0.6, 0.3):
        for fx in (0.5, 0.4, 0.6):
            y = min(max(0, int(H * fy)), H - n)
            x = min(max(0, int(W * fx)), W - n)
            a = d[y:y + n, x:x + n]
            v = a.real if cplx else a
            good = np.isfinite(v) & (v != 0)
            if good.mean() > 0.5:
                return a, good.mean()
            if best is None or good.mean() > best[1]:
                best = (a, good.mean())
    return best


def layerRows(f, bits, minPx=1e6):
    import h5py
    import numpy as np
    rows = []
    cands = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Dataset) or obj.ndim != 2 or obj.size < minPx:
            return
        t = obj.id.get_type()
        try:
            cplx = t.get_nmembers() == 2
        except Exception:
            cplx = False
        if obj.dtype.kind == 'f' or cplx:
            cands.append((name, cplx))
    f.visititems(visit)
    for name, cplx in cands:
        d = f[name]
        disk = d.id.get_storage_size()
        if disk < 1e6:
            continue
        w = window(d, cplx)
        if w is None:
            continue
        a, frac = w
        a = a.astype('complex64') if cplx else a.astype('float32')
        comp = (np.concatenate([a.real.ravel(), a.imag.ravel()]) if cplx
                else a.ravel()).astype('float32')
        nz = comp[np.isfinite(comp) & (comp != 0)]
        if nz.size < 1000:
            continue
        live = float(((((np.abs(nz).view('uint32')) & 0x7FFFFF) & 0x1FFF) != 0).mean())
        lvl = d.compression_opts if d.compression == 'gzip' else 0
        b0 = store(a, lvl, d.shuffle)
        # Geocoded products are largely nodata, and nodata already compresses to almost
        # nothing. Projecting a valid-window ratio onto the whole layer would claim a saving
        # on bytes that are not there. Price nodata separately and infer how much of the
        # layer's on-disk bytes are actually attributable to valid samples:
        #     disk = nValid*bitsValid/8 + nNodata*bitsNodata/8
        # then apply the truncation saving only to the valid part.
        nw = nodataWindow(d, cplx)
        bitsNodata = 0.0
        if nw is not None:
            na = nw[0].astype('complex64') if cplx else nw[0].astype('float32')
            bitsNodata = store(na, lvl, d.shuffle) * 8.0 / na.size
        bitsValid = b0 * 8.0 / a.size
        npx = float(d.size)
        denom = (bitsValid - bitsNodata)
        vf = ((disk * 8.0 / npx) - bitsNodata) / denom if denom > 1e-9 else frac
        vf = min(max(vf, 0.0), 1.0)
        row = dict(name=name.split('/')[-1], disk=disk, live=live, cplx=cplx,
                   applied=bool(live < 0.5), bitsValid=bitsValid, bitsNodata=bitsNodata,
                   validFrac=vf, npx=npx)
        for nb in bits:
            m = mask(nb)
            if cplx:
                zr = (a.real.view(np.uint32) & m).view(np.float32)
                zi = (a.imag.view(np.uint32) & m).view(np.float32)
                z = (zr + 1j * zi).astype('complex64')
            else:
                z = (a.view(np.uint32) & m).view(np.float32)
            bz = store(z, lvl, d.shuffle)
            row['r%d' % nb] = b0 / bz
            row['save%d' % nb] = (0 if row['applied'] else
                                  max(0.0, npx * vf * (bitsValid - bz * 8.0 / a.size) / 8.0))
        rows.append(row)
    return rows


def main():
    import asf_search as asf
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--track', default='131')
    ap.add_argument('--dir', default='D')
    ap.add_argument('--frame', default='062')
    ap.add_argument('--bits', type=int, nargs='+', default=[6, 10])
    ap.add_argument('--json', default=None)
    args = ap.parse_args()
    pat = '_%s_%s_%s_' % (args.track, args.dir, args.frame)
    LEV = ['RSLC', 'GSLC', 'GCOV', 'RIFG', 'RUNW', 'ROFF', 'GUNW', 'GOFF']
    out = {}
    for lv in LEV:
        print('=== %s ===' % lv, flush=True)
        pick = None
        try:
            for g in asf.search(dataset='NISAR', processingLevel=lv, maxResults=400):
                if pat in g.properties.get('fileName', ''):
                    pick = g
                    break
        except Exception as e:
            print('   search failed: %s' % type(e).__name__, flush=True)
            continue
        if pick is None:
            print('   no granule for this frame', flush=True)
            continue
        signed = resolve(pick.properties['url'])
        f = openRemote(signed) if signed else None
        if f is None:
            print('   could not open', flush=True)
            continue
        try:
            rows = layerRows(f, args.bits)
        except Exception as e:
            print('   %s: %s' % (type(e).__name__, str(e)[:70]), flush=True)
            rows = []
        finally:
            try:
                f.close()
            except Exception:
                pass
        out[lv] = dict(granule=pick.properties['fileName'], rows=rows)
        for r in rows:
            print('   %-24s %8.1f MB live %5.1f%% valid %5.1f%% %s %s' %
                  (r['name'][:24], r['disk'] / 1e6, 100 * r['live'], 100 * r['validFrac'],
                   'APPLIED' if r['applied'] else 'not-app',
                   ' '.join('%db %.2fx sv %6.1fMB' % (nb, r['r%d' % nb], r['save%d' % nb] / 1e6)
                            for nb in args.bits)), flush=True)
    if args.json:
        json.dump(out, open(args.json, 'w'), indent=2)
        print('\nwrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
