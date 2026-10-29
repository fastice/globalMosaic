#!/usr/bin/env python3
"""
fp16 vs fp32: read and write, across C/GDAL/h5py and two GDAL+HDF5 stacks.

Two separate questions, because NISAR products differ:

  GCOV   archive bands are float32 and slim products are float16, so this is a NATIVE
         comparison and the one the geomosaic fast path exists for.
  RSLC   ASF ships complex64 -- an HDF5 COMPOUND of two 32-bit floats, not a native complex
         type, and NOT complex32 (checked on routine and urgent products, 20 and 40 MHz).
         So fp16 here is a conversion YOU would make; this measures what it costs and saves.

Cache handling: every read is preceded by posix_fadvise(POSIX_FADV_DONTNEED) on the file, which
needs no root, unlike drop_caches. Without it the second read of a file is served from page
cache and every number after the first is meaningless. A "cached" row is reported alongside so
the gap is visible rather than assumed.

    fp16bench.py --source <product.h5> --outdir <dir on fast local disk>
    fp16bench.py --outdir <dir> --synthetic          # no product to hand
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STACKS = {
    'system': dict(bin=os.path.join(HERE, 'fp16read.sys'),
                   py=None),
    'conda':  dict(bin=os.path.join(HERE, 'fp16read.conda'),
                   py=sys.executable),
}


def dropCache(path):
    """Evict this file from the page cache. No root needed, unlike /proc/sys/vm/drop_caches."""
    try:
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)
    except OSError:
        pass


def timeIt(fn, path, cached=False):
    if not cached:
        dropCache(path)
    t0 = time.time()
    n = fn()
    return time.time() - t0, n


def runC(exe, args, path, cached=False):
    if not os.path.exists(exe):
        return None
    if not cached:
        dropCache(path)
    p = subprocess.run([exe] + [str(a) for a in args], capture_output=True, text=True)
    if p.returncode != 0:
        return None
    try:
        dt, mpxs = p.stdout.strip().split()
        return float(dt)
    except ValueError:
        return None


def buildSources(args):
    """Make float32 / float16 (and complex) test files from a real product, or synthetically."""
    import numpy as np
    import h5py
    N = args.size
    out = {}
    if args.synthetic or not args.source:
        rng = np.random.default_rng(0)
        a = np.exp(rng.normal(-3.0, 1.0, size=(N, N))).astype('float32')
        src = 'synthetic log-normal'
    else:
        with h5py.File(args.source, 'r') as f:
            ds = args.dataset or findBand(f)
            d = f[ds]
            y = min(args.yoff, max(0, d.shape[0] - N))
            x = min(args.xoff, max(0, d.shape[1] - N))
            blk = d[y:y + N, x:x + N]
            src = '%s[%d:%d, %d:%d] %s' % (ds, y, y + N, x, x + N, d.dtype)
            if np.iscomplexobj(blk):
                out['c64'] = blk.astype('complex64')
                a = np.abs(blk).astype('float32')
            else:
                a = blk.astype('float32')
        # scale so the values sit where float16 has resolution, as gamma0 in linear power does
        a = (a / max(np.percentile(a, 99), 1e-30) * 0.1).astype('float32')
    out['f32'] = a
    out['f16'] = a.astype('float16')
    if 'c64' in out:
        # NISAR stores complex as an HDF5 COMPOUND {r, i}, not a native complex type -- h5py
        # writes numpy complex64 that way automatically. "complex32" is the same layout with
        # float16 members, so build it as a structured dtype rather than an (N, N, 2) array:
        # same bytes, but a 2-D dataset with the real chunking and the real type.
        c = out['c64']
        c32 = np.empty(c.shape, dtype=np.dtype([('r', '<f2'), ('i', '<f2')]))
        c32['r'] = c.real.astype('float16')
        c32['i'] = c.imag.astype('float16')
        out['c32'] = c32
    return out, src


def findBand(f):
    import h5py
    hits = []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and obj.ndim == 2 and obj.size > 4e6:
            hits.append((obj.size, name))
    f.visititems(visit)
    if not hits:
        raise SystemExit('no 2-D dataset larger than 4 Mpx found; pass --dataset')
    hits.sort()
    return hits[-1][1]


def writeH5(path, name, arr, deflate, shuffle):
    import h5py
    if os.path.exists(path):
        os.unlink(path)
    t0 = time.time()
    with h5py.File(path, 'w') as f:
        f.create_dataset(name, data=arr, chunks=(min(512, arr.shape[0]), min(512, arr.shape[1])),
                         compression=('gzip' if deflate else None),
                         compression_opts=(deflate or None), shuffle=shuffle)
    return time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True, help='work dir, on FAST LOCAL disk not NFS')
    ap.add_argument('--source', default=None, help='a real NISAR product to take data from')
    ap.add_argument('--dataset', default=None, help='dataset path inside it [largest 2-D]')
    ap.add_argument('--size', type=int, default=8192, help='square window [8192]')
    ap.add_argument('--xoff', type=int, default=4096)
    ap.add_argument('--yoff', type=int, default=4096)
    ap.add_argument('--deflate', type=int, default=4, help='gzip level, 0 for none [4]')
    ap.add_argument('--shuffle', action='store_true', help='HDF5 shuffle filter (NISAR uses it)')
    ap.add_argument('--synthetic', action='store_true')
    ap.add_argument('--json', default=None, help='also write results here')
    args = ap.parse_args()

    import numpy as np
    os.makedirs(args.outdir, exist_ok=True)
    res = {'host': os.uname().nodename, 'outdir': args.outdir,
           'deflate': args.deflate, 'shuffle': args.shuffle, 'size': args.size}

    print('building test data ...', flush=True)
    arrays, srcDesc = buildSources(args)
    res['source'] = srcDesc
    print('  source: %s' % srcDesc)
    N = args.size
    mpx = N * N / 1e6

    # ---------------- write ----------------
    print('\n=== WRITE (h5py, gzip %d, shuffle %s) ===' % (args.deflate, args.shuffle))
    paths, writes = {}, {}
    for tag in ('f32', 'f16'):
        p = os.path.join(args.outdir, tag + '.h5')
        dt = writeH5(p, 'data', arrays[tag], args.deflate, args.shuffle)
        paths[tag] = p
        writes[tag] = dt
        print('  %-5s %7.2f s  %7.2f Mpx/s  file %8.1f MB'
              % (tag, dt, mpx / dt, os.path.getsize(p) / 1e6))
    if 'c64' in arrays:
        for tag in ('c64', 'c32'):
            p = os.path.join(args.outdir, tag + '.h5')
            dt = writeH5(p, 'data', arrays[tag], args.deflate, args.shuffle)
            paths[tag] = p
            writes[tag] = dt
            print('  %-5s %7.2f s  %7.2f Mpx/s  file %8.1f MB'
                  % (tag, dt, mpx / dt, os.path.getsize(p) / 1e6))
    res['write'] = {k: {'sec': v, 'bytes': os.path.getsize(paths[k])} for k, v in writes.items()}

    # ---------------- read ----------------
    print('\n=== READ (cold: posix_fadvise DONTNEED before each) ===')
    print('  %-28s %8s %8s %10s' % ('reader', 'f32 s', 'f16 s', 'f16 vs f32'))
    rows = {}

    def show(label, a, b):
        if a is None or b is None:
            print('  %-28s %8s %8s %10s'
                  % (label, '-' if a is None else '%.2f' % a,
                     '-' if b is None else '%.2f' % b, 'n/a'))
            return
        print('  %-28s %8.2f %8.2f %9.2fx' % (label, a, b, a / b))
        rows[label] = {'f32': a, 'f16': b}

    for stack, cfg in STACKS.items():
        exe = cfg['bin']
        if not os.path.exists(exe):
            continue
        for mode, lab in (('gdal', 'C GDAL'), ('h5', 'C HDF5 convert'), ('lut', 'C HDF5 + LUT')):
            a = runC(exe, [mode, paths['f32'], 'data', 0, 0, N, N], paths['f32'])
            b = runC(exe, [mode, paths['f16'], 'data', 0, 0, N, N], paths['f16'])
            show('%s / %s' % (lab, stack), a, b)

    from osgeo import gdal
    gdal.UseExceptions()

    def gdalRead(p):
        def go():
            ds = gdal.Open('HDF5:"%s"://data' % p)
            return ds.GetRasterBand(1).ReadAsArray().size
        return timeIt(go, p)[0]
    show('python GDAL %s' % gdal.__version__.split('.')[0] + '.' + gdal.__version__.split('.')[1],
         gdalRead(paths['f32']), gdalRead(paths['f16']))

    import h5py

    def h5Read(p, native):
        def go():
            with h5py.File(p, 'r') as f:
                d = f['data']
                if native:
                    return d.id.read_direct_chunk if False else d[...].size
                out = np.empty(d.shape, 'float32')
                d.read_direct(out)
                return out.size
        return timeIt(go, p)[0]
    show('h5py -> native dtype', h5Read(paths['f32'], True), h5Read(paths['f16'], True))
    show('h5py -> float32 (HDF5 converts)', h5Read(paths['f32'], False), h5Read(paths['f16'], False))

    # cached, to show what the page cache hides
    def cachedGdal(p):
        def go():
            ds = gdal.Open('HDF5:"%s"://data' % p)
            return ds.GetRasterBand(1).ReadAsArray().size
        return timeIt(go, p, cached=True)[0]
    cachedGdal(paths['f32']); cachedGdal(paths['f16'])   # prime
    show('python GDAL (WARM cache)', cachedGdal(paths['f32']), cachedGdal(paths['f16']))
    res['read'] = rows

    # ---------------- conversion: numpy vs HDF5 ----------------
    print('\n=== CONVERT fp16 -> fp32 (same bytes, different converter) ===')
    with h5py.File(paths['f16'], 'r') as f:
        raw = f['data'][...]                 # native float16, in memory
    t0 = time.time(); _ = raw.astype('float32'); tNumpy = time.time() - t0
    dropCache(paths['f16'])
    t0 = time.time()
    with h5py.File(paths['f16'], 'r') as f:
        out = np.empty(f['data'].shape, 'float32'); f['data'].read_direct(out)
    tH5 = time.time() - t0
    print('  numpy astype (in memory)      %7.3f s  %8.1f Mpx/s' % (tNumpy, mpx / tNumpy))
    print('  HDF5 convert on read (cold)   %7.3f s  %8.1f Mpx/s' % (tH5, mpx / tH5))
    res['convert'] = {'numpyAstype': tNumpy, 'hdf5OnRead': tH5}

    if 'c64' in arrays:
        print('\n=== complex64 vs complex32-as-2xfp16 ===')
        for tag in ('c64', 'c32'):
            print('  %-5s file %8.1f MB  write %6.2f s'
                  % (tag, os.path.getsize(paths[tag]) / 1e6, writes[tag]))
        print('  storage ratio  %.2fx' % (os.path.getsize(paths['c64'])
                                          / max(1, os.path.getsize(paths['c32']))))

    if args.json:
        with open(args.json, 'w') as f:
            json.dump(res, f, indent=2)
        print('\nwrote %s' % args.json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
