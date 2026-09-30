#!/usr/bin/env python3
'''
Valid-data footprints of GCOV granules from their masks, for the tiler.

  scanMasks cycle30/catalogue.geojson --out cycle30/validFootprints.geojson [--all] [--threads 16]

A GCOV footprint is the acquisition extent, but only mask values 1..254 (the subswath number) are
valid: 0 means partially focused, 255 outside the swath, and geomosaic (useMask) drops both. A
partial frame at the start of a datatake can be entirely mask 0, so it covers nothing although its
footprint says otherwise. For each granule (partial frames `_P_` only, unless --all) this reads
frequency A's mask (or B's, when A is absent) every --stride pixels over https, polygonizes the
valid part and writes it in EPSG:4326 with validFrac (valid share of the in-swath samples).
Resumable: granules already in --out are skipped.
'''
import argparse
import concurrent.futures
import json
import os
import sys

import numpy as np
import pyproj
from osgeo import gdal
from rasterio import features
from rasterio.transform import Affine
from shapely.geometry import mapping, shape
from shapely.ops import transform, unary_union

gdal.UseExceptions()


def validFootprint(url, stride):
    # per-process Earthdata cookie jar: workers sharing one jar break each other's login redirect
    cookies = os.path.expanduser(f'~/.scanMasks.cookies.{os.getpid()}')
    for k, v in (('GDAL_HTTP_NETRC', 'YES'), ('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR'),
                 ('GDAL_HTTP_COOKIEFILE', cookies), ('GDAL_HTTP_COOKIEJAR', cookies),
                 ('GDAL_HTTP_MAX_RETRY', '5'), ('GDAL_HTTP_RETRY_DELAY', '10')):
        gdal.SetThreadLocalConfigOption(k, v)
    ds = gdal.OpenEx('/vsicurl/' + url, gdal.OF_MULTIDIM_RASTER)
    root = ds.GetRootGroup()
    g = None
    for freq in ('A', 'B'):
        try:
            g = root.OpenGroupFromFullname(f'/science/LSAR/GCOV/grids/frequency{freq}')
            break
        except Exception:
            continue
    if g is None:
        return None, 0.
    xs = g.OpenMDArray('xCoordinates').ReadAsArray()
    ys = g.OpenMDArray('yCoordinates').ReadAsArray()
    epsg = int(g.OpenMDArray('projection').GetAttribute('epsg_code').Read())
    m = g.OpenMDArray('mask').GetView(f'[::{stride},::{stride}]').ReadAsArray()
    inSwath = m != 255
    valid = inSwath & (m != 0)
    frac = float(valid.sum() / max(inSwath.sum(), 1))
    if not valid.any():
        return None, frac
    dx, dy = (xs[1] - xs[0]) * stride, (ys[1] - ys[0]) * stride
    aff = Affine(dx, 0, xs[0] - dx / 2, 0, dy, ys[0] - dy / 2)
    polys = [shape(p) for p, v in features.shapes(valid.astype(np.uint8), mask=valid, transform=aff) if v]
    geom = unary_union(polys).buffer(abs(dx)).buffer(-abs(dx))          # close pinholes
    toLL = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True).transform
    return transform(toLL, geom.simplify(abs(dx))), frac


def scanOne(name, url, stride):
    try:
        geom, frac = validFootprint(url, stride)
        return name, geom, frac, None
    except Exception as e:
        return name, None, None, str(e)[:200] or type(e).__name__


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('catalogue', help='catalogue.geojson from globalGCOVTiles')
    ap.add_argument('--out', required=True, help='output geojson (resumable)')
    ap.add_argument('--all', action='store_true', help='scan every granule, not only partial frames')
    ap.add_argument('--direction', choices=['A', 'D', 'both'], default='both',
                    help='only granules of this orbit direction, by file name [both]')
    ap.add_argument('--stride', type=int, default=32, help='mask subsampling [32 (640 m at 20 m)]')
    ap.add_argument('--threads', type=int, default=16, help='worker processes [16]')
    args = ap.parse_args()

    cat = json.load(open(args.catalogue))['features']
    todo = [f['properties'] for f in cat
            if (args.all or f['properties']['name'].split('_')[15] == 'P')
            and args.direction in ('both', f['properties']['name'].split('_')[6])]
    done = {}
    if os.path.exists(args.out):
        done = {f['properties']['name']: f for f in json.load(open(args.out))['features']}
    todo = [p for p in todo if p['name'] not in done]
    print(f'{len(todo)} granules to scan ({len(done)} already done)', flush=True)


    def save():
        with open(args.out + '.tmp', 'w') as fp:
            json.dump({'type': 'FeatureCollection', 'features': list(done.values())}, fp)
        os.replace(args.out + '.tmp', args.out)

    # processes, not threads: GDAL's HDF5 layer is not thread-safe (threads hang)
    with concurrent.futures.ProcessPoolExecutor(args.threads) as pool:
        futs = [pool.submit(scanOne, p['name'], p['url'], args.stride) for p in todo]
        for k, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            name, geom, frac, err = fut.result()
            if err is not None:
                print(f'  {name}: ERROR {err}', flush=True)
                continue
            done[name] = {'type': 'Feature', 'properties': {'name': name, 'validFrac': round(frac, 4)},
                          'geometry': mapping(geom) if geom is not None else None}
            if k % 50 == 0 or k == len(todo):
                save()
                print(f'{k}/{len(todo)} scanned', flush=True)
    save()
    fr = [f['properties']['validFrac'] for f in done.values()]
    print(f'done: {len(fr)} granules; all-invalid {sum(v == 0 for v in fr)}, '
          f'< 50% valid {sum(v < 0.5 for v in fr)}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
