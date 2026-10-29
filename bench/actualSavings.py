#!/usr/bin/env python3
"""
Download a frame's products and MEASURE the saving by rewriting them, not by extrapolating.

Earlier estimates sampled a 1024x1024 window per layer and projected the compression ratio onto
the whole layer. That assumes the window is representative and that our write path matches the
producer's. This does the real thing: download each product, rewrite every float32/complex64
raster with its low mantissa bits zeroed, keep every dtype/chunk/filter/attribute identical,
and compare actual file sizes.

RSLC and GSLC are skipped -- ISCE3 already applies `complex64_zero_mantissa` to those
(verified: 0.00% of their mantissa low bits are in use).

    actualSavings.py --track 131 --dir D --frame 062 --work /fast/disk/dir
"""
import argparse
import json
import os
import subprocess
import sys
import time

MANT = 23
DEFAULT_BITS = 10
PER_LAYER_BITS = {'numberOfLooks': 6}      # needs ~0.01 looks; 6 bits gives 0.06
LEVELS = ['GCOV', 'RIFG', 'RUNW', 'ROFF', 'GUNW', 'GOFF']


def mask(nbits):
    import numpy as np
    return np.uint32(((0xFFFFFFFF << (MANT - nbits)) & 0xFFFFFFFF))


def download(url, dest, limit='40M'):
    if os.path.exists(dest) and os.path.getsize(dest) > 1e6:
        return True
    cmd = ['aria2c', '-x8', '-s8', '--max-overall-download-limit=' + limit,
           '--netrc-path=' + os.path.expanduser('~/.netrc'), '--auto-file-renaming=false',
           '--allow-overwrite=true', '--console-log-level=error', '--summary-interval=0',
           '-d', os.path.dirname(dest), '-o', os.path.basename(dest), url]
    return subprocess.run(cmd).returncode == 0


# HDF5 object-reference attributes address objects in the SOURCE file. Writing them into a
# different file segfaults -- no Python traceback, so it presents as an external kill. These
# four are the dimension-scale bookkeeping; HDF5 rebuilds them via make_scale/attach_scale.
REF_ATTRS = ('DIMENSION_LIST', 'REFERENCE_LIST', 'CLASS', 'NAME')


def copyAttrs(src, dst):
    for k, v in src.attrs.items():
        if k in REF_ATTRS:
            continue
        try:
            dst.attrs[k] = v
        except Exception:
            pass


def truncateFile(src, dst):
    """Rewrite src to dst with mantissas truncated. Everything else byte-for-byte equivalent."""
    import h5py
    import numpy as np
    if os.path.exists(dst):
        os.unlink(dst)
    touched = []
    with h5py.File(src, 'r') as fi, h5py.File(dst, 'w') as fo:
        def isRaster(d):
            if d.ndim != 2 or d.size < 1e6:
                return False
            if d.dtype.kind == 'f':
                return True
            return d.dtype.kind == 'c'

        def visit(name, obj):
            # visititems also yields committed (named) Datatype objects, which are neither
            # Group nor Dataset and have no .ndim -- RIFG and GUNW carry them.
            if isinstance(obj, h5py.Datatype):
                try:
                    fi.copy(name, fo, name=name)
                except Exception:
                    pass
                return
            if isinstance(obj, h5py.Group):
                g = fo.require_group(name)
                copyAttrs(obj, g)
                return
            if not isRaster(obj):
                # h5py's Group.copy() SEGFAULTS on compound (complex) datasets -- reproduced
                # deterministically on GCOV's calibrationInformation/crosstalk/* arrays, which
                # are complex64. It dies with no Python traceback, so it looks like a kill.
                # Complex types are a known sharp edge in h5py (see the ComplexFloat16Decoder
                # comment in isce3/core/types.py citing h5py issue 2156). Write those out
                # explicitly instead; everything else copies verbatim.
                if obj.dtype.kind == 'c':
                    d2 = fo.create_dataset(name, shape=obj.shape, dtype=obj.dtype,
                                           chunks=obj.chunks, compression=obj.compression,
                                           compression_opts=obj.compression_opts,
                                           shuffle=obj.shuffle)
                    d2[...] = obj[...]
                    copyAttrs(obj, d2)
                else:
                    # fi.copy() also carries DIMENSION_LIST across, so write these out too
                    # rather than relying on it.
                    try:
                        d2 = fo.create_dataset(name, data=obj[()], chunks=obj.chunks,
                                               compression=obj.compression,
                                               compression_opts=obj.compression_opts,
                                               shuffle=obj.shuffle)
                    except Exception:
                        d2 = fo.create_dataset(name, data=obj[()])
                    copyAttrs(obj, d2)
                return
            nb = PER_LAYER_BITS.get(name.split('/')[-1], DEFAULT_BITS)
            m = mask(nb)
            out = fo.create_dataset(name, shape=obj.shape, dtype=obj.dtype,
                                    chunks=obj.chunks, compression=obj.compression,
                                    compression_opts=obj.compression_opts,
                                    shuffle=obj.shuffle, fletcher32=obj.fletcher32)
            step = max(1, (obj.chunks or (512, 512))[0]) * 8
            t0 = time.time()
            print('      %-26s %2d bits  %.0f Mpx ...' % (name.split('/')[-1][:26], nb,
                  obj.size / 1e6), flush=True)
            for y in range(0, obj.shape[0], step):
                a = obj[y:y + step]
                if a.dtype.kind == 'c':
                    r = np.ascontiguousarray(a.real).view(np.uint32)
                    i = np.ascontiguousarray(a.imag).view(np.uint32)
                    r = (r & m).view(np.float32)
                    i = (i & m).view(np.float32)
                    a = (r + 1j * i).astype(a.dtype)
                else:
                    a = (np.ascontiguousarray(a).view(np.uint32) & m).view(a.dtype)
                out[y:y + a.shape[0]] = a
            copyAttrs(obj, out)
            print('      %-26s done in %.0f s' % (name.split('/')[-1][:26], time.time() - t0),
                  flush=True)
            touched.append((name.split('/')[-1], nb, obj.size))
        fi.visititems(visit)
        copyAttrs(fi, fo)
    return touched


def main():
    import asf_search as asf
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--track', default='131')
    ap.add_argument('--dir', default='D')
    ap.add_argument('--frame', default='062')
    ap.add_argument('--work', required=True)
    ap.add_argument('--levels', nargs='+', default=LEVELS)
    ap.add_argument('--json', default=None)
    ap.add_argument('--keep', action='store_true', help='keep the rewritten copies')
    args = ap.parse_args()
    os.makedirs(args.work, exist_ok=True)
    pat = '_%s_%s_%s_' % (args.track, args.dir, args.frame)
    out = {}
    for lv in args.levels:
        print('=== %s ===' % lv, flush=True)
        pick = None
        for g in asf.search(dataset='NISAR', processingLevel=lv, maxResults=400):
            if pat in g.properties.get('fileName', ''):
                pick = g
                break
        if pick is None:
            print('   no granule for this frame', flush=True)
            continue
        name = pick.properties['fileName']
        src = os.path.join(args.work, name)
        t0 = time.time()
        if not download(pick.properties['url'], src):
            print('   download failed', flush=True)
            continue
        sz0 = os.path.getsize(src)
        print('   downloaded %.2f GB in %.0f s' % (sz0 / 1e9, time.time() - t0), flush=True)
        dst = src.replace('.h5', '.trunc.h5')
        t1 = time.time()
        try:
            touched = truncateFile(src, dst)
        except Exception as e:
            print('   rewrite FAILED: %s: %s' % (type(e).__name__, str(e)[:90]), flush=True)
            continue
        sz1 = os.path.getsize(dst)
        print('   rewrote in %.0f s : %.2f GB -> %.2f GB  (%.2fx, saves %.2f GB / %.1f%%)'
              % (time.time() - t1, sz0 / 1e9, sz1 / 1e9, sz0 / sz1,
                 (sz0 - sz1) / 1e9, 100 * (sz0 - sz1) / sz0), flush=True)
        for n, nb, px in touched:
            print('      %-26s %2d bits  %.1f Mpx' % (n[:26], nb, px / 1e6), flush=True)
        out[lv] = dict(granule=name, before=sz0, after=sz1,
                       layers=[dict(name=n, bits=nb, px=px) for n, nb, px in touched])
        if not args.keep:
            os.unlink(dst)
    print('\n%-6s %10s %10s %9s %9s' % ('LEVEL', 'before GB', 'after GB', 'ratio', 'saved %'))
    b = a = 0
    for lv in args.levels:
        r = out.get(lv)
        if not r:
            continue
        b += r['before']
        a += r['after']
        print('%-6s %10.2f %10.2f %8.2fx %8.1f%%'
              % (lv, r['before'] / 1e9, r['after'] / 1e9, r['before'] / r['after'],
                 100 * (r['before'] - r['after']) / r['before']))
    if b:
        print('%-6s %10.2f %10.2f %8.2fx %8.1f%%'
              % ('TOTAL', b / 1e9, a / 1e9, b / a, 100 * (b - a) / b))
    if args.json:
        json.dump(out, open(args.json, 'w'), indent=2)
        print('wrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
