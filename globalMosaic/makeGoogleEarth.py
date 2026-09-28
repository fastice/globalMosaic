#!/usr/bin/env python3
'''
Google Earth product from a runGeomosaicTiles run: an EPSG:4326 KML superoverlay of compressed
8-bit PNG tiles (grey + transparency), built with gdal2tiles from the run's VRTs.

  makeGoogleEarth run_c030 --out run_c030/googleEarth [--zoom 0-9] [--processes 16]
  makeGoogleEarth --sources ant_gamma0_100m.tif --out antGE        # any rasters, any CRS

  1. mosaic: vrt/<product>/global.vrt (lat/lon bands) plus every polar cap VRT warped to
     EPSG:4326 at the same spacing, in one VRT (caps painted first, bands on top). With
     --sources, the given rasters instead (int16 dB x 100, nodata -3000, as geomosaic writes);
     any not in EPSG:4326 are warped (virtually) at --res.
  2. stretch: int16 dB x 100 -> byte 1..255 over --dbMin..--dbMax (a virtual layer; nothing
     is written); nodata becomes a mask, so it is transparent in the PNGs.
  3. tiles: gdal2tiles --profile=geodetic -k (KML superoverlay), --zoom levels, PNG tiles,
     resumable (-e). Open <out>/doc.kml in Google Earth.

Zoom (gdal2tiles geodetic profile, 256-pixel tiles; a level-z tile spans 360 / 2**z deg, as its
KML shows): pixel = 1.40625 / 2**z deg -- 8: 19.8" (~610 m), 9: 9.9" (~305 m), 10: 4.9" (~153 m),
11: 2.5" (~76 m, full resolution for a 3" mosaic), 12: 1.2" (~38 m). Size grows ~4x per level.
'''
import argparse
import glob
import os
import re
import subprocess
import sys

from osgeo import gdal

gdal.UseExceptions()
NODATA = -3000


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', nargs='?', default=None, help='runGeomosaicTiles --work directory')
    ap.add_argument('--sources', nargs='+', default=None,
                    help='rasters to use instead of a run directory (any CRS)')
    ap.add_argument('--res', type=float, default=3. / 3600,
                    help='lat/lon spacing for warped sources, deg [3 arcsec]')
    ap.add_argument('--out', required=True, help='output directory (doc.kml + tile pyramid)')
    ap.add_argument('--product', default='gamma0', help='gamma0 or sigma0 [gamma0]')
    ap.add_argument('--zoom', default='0-9', help="gdal2tiles zoom levels, e.g. '0-9' or '0-10' [0-9]")
    ap.add_argument('--dbMin', type=float, default=-24., help='stretch minimum, dB [-24]')
    ap.add_argument('--dbMax', type=float, default=-1., help='stretch maximum, dB [-1]')
    ap.add_argument('--resampling', default='average', help='gdal2tiles resampling [average]')
    ap.add_argument('--processes', type=int, default=16, help='gdal2tiles processes [16]')
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    stage = f'{args.out}/stage'
    os.makedirs(stage, exist_ok=True)
    if args.sources:
        inputs, res = [os.path.abspath(s) for s in args.sources], args.res
    else:
        if args.work is None:
            ap.error('give a run directory or --sources')
        vrtDir = os.path.abspath(f"{args.work}/vrt/{args.product}")
        ds = gdal.Open(f'{vrtDir}/global.vrt')
        res = abs(ds.GetGeoTransform()[5])
        ds = None
        # polar caps first (own CRS), bands on top
        inputs = sorted(glob.glob(f'{vrtDir}/cap_*.vrt')) + [f'{vrtDir}/global.vrt']
    # 1. anything not already lat/lon is warped to EPSG:4326 at the mosaic spacing (virtual)
    srcs = []
    for src in inputs:
        srs = gdal.Open(src).GetSpatialRef()
        if srs is not None and srs.GetAuthorityCode(None) == '4326':
            srcs.append(src)
            continue
        warped = f'{stage}/{os.path.splitext(os.path.basename(src))[0]}.4326.vrt'
        gdal.Warp(warped, src, format='VRT', dstSRS='EPSG:4326', xRes=res, yRes=res,
                  srcNodata=NODATA, dstNodata=NODATA, resampleAlg='average')
        srcs.append(warped)
    mosaic = f'{stage}/mosaic.vrt'
    gdal.BuildVRT(mosaic, srcs, resolution='user', xRes=res, yRes=res, srcNodata=NODATA,
                  VRTNodata=NODATA)
    # 2. 8-bit stretch; the source nodata becomes a mask band (transparent in the PNGs)
    byte = f'{stage}/mosaic.byte.vrt'
    gdal.Translate(byte, mosaic, format='VRT', outputType=gdal.GDT_Byte, noData=0,
                   scaleParams=[[args.dbMin * 100, args.dbMax * 100, 1, 255]], maskBand='auto')
    # VRT scaling does not clamp to 1..255: data darker than dbMin would land on 0 = nodata
    # (transparent). A LUT clamps at its end points; source nodata is applied before it.
    with open(byte) as fp:
        xml = fp.read()
    xml = re.sub(r'\s*<ScaleOffset>.*?</ScaleOffset>\s*<ScaleRatio>.*?</ScaleRatio>',
                 f'\n      <LUT>{args.dbMin * 100:g}:1,{args.dbMax * 100:g}:255</LUT>', xml)
    xml = re.sub(r'\s*<Scale>.*?</Scale>', '', xml)
    with open(byte, 'w') as fp:
        fp.write(xml)
    # 3. KML superoverlay of PNG tiles
    cmd = ['gdal2tiles.py', '--profile=geodetic', '-k', '-z', args.zoom, '-r', args.resampling,
           f'--processes={args.processes}', '-w', 'none', '-e', '--tiledriver=PNG', byte, args.out]
    print(' '.join(cmd), flush=True)
    rc = subprocess.run(cmd).returncode
    if rc == 0:
        n = sum(len(f) for _, _, f in os.walk(args.out) if True)
        print(f'done: {args.out}/doc.kml; files under {args.out}: {n}')
    return rc


if __name__ == '__main__':
    sys.exit(main())
