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
import math
import os
import re
import subprocess
import sys
from xml.sax.saxutils import escape

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
    ap.add_argument('--zoom', default='0-9', help="zoom levels as 256-px equivalents, e.g. '0-9' or '0-10' [0-9]")
    ap.add_argument('--tileSize', type=int, default=256, choices=[256, 512],
                    help='tile size, px; 512 writes a quarter as many files, same detail [256]')
    ap.add_argument('--dbMin', type=float, default=-24., help='stretch minimum, dB [-24]')
    ap.add_argument('--dbMax', type=float, default=-1., help='stretch maximum, dB [-1]')
    ap.add_argument('--resampling', default='average', help='gdal2tiles resampling [average]')
    ap.add_argument('--processes', type=int, default=16, help='gdal2tiles processes [16]')
    ap.add_argument('--title', default=None,
                    help='name shown in Google Earth [NISAR <product> <run directory name>]')
    ap.add_argument('--retitle', action='store_true',
                    help='only rename an existing product (<out>/doc.kml) to --title, no tiling')
    args = ap.parse_args()
    title = args.title or f'NISAR {args.product} {os.path.basename(os.path.abspath(args.work or args.out))}'
    if args.retitle:
        setTitle(f'{args.out}/doc.kml', title)
        print(f'{args.out}/doc.kml -> "{title}"')
        return 0

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
    return superoverlay(inputs, res, args.out, args.zoom, args.dbMin, args.dbMax,
                        args.resampling, args.processes, title, tileSize=args.tileSize)


def setTitle(kml, title):
    ''' Rename a superoverlay: the first <name> of its doc.kml (what Google Earth lists). '''
    s = open(kml).read()
    s = re.sub(r'<name>.*?</name>', f'<name>{escape(title)}</name>', s, count=1)
    open(kml, 'w').write(s)


def tileZoom(zoom, tileSize=256):
    ''' gdal2tiles zoom levels for a tile size, from levels given as 256-px equivalents: a 512-px tile
    at level z has the pixel of a 256-px tile at z + 1, so '7-11' -> '6-10' (never below 0). '''
    k = int(round(math.log2(tileSize / 256)))
    lo, hi = (int(v) for v in zoom.split('-'))
    return f'{max(lo - k, 0)}-{hi - k}'


def superoverlay(inputs, res, out, zoom, dbMin=-24., dbMax=-1., resampling='average', processes=16,
                 title=None, nodata=NODATA, scale=100., tileSize=256):
    ''' nodata / scale: int16 dB x 100 with -3000 by default; Float32 linear layers (e.g. CV) pass
    nodata=float('nan'), scale=1 and the stretch in their own units as dbMin/dbMax.
    zoom: levels as 256-px equivalents (same pixel size whatever tileSize); tileSize 512 writes a
    quarter as many files (much faster on network file systems such as EFS). '''
    ''' KML superoverlay of 8-bit PNG tiles of int16 dB x 100 rasters (any CRS) at <out>/doc.kml.
    Returns gdal2tiles' exit code. '''
    os.makedirs(out, exist_ok=True)
    stage = f'{out}/stage'
    os.makedirs(stage, exist_ok=True)
    # 1. anything not already lat/lon is warped to EPSG:4326 at the mosaic spacing (virtual)
    srcs = []
    for src in inputs:
        srs = gdal.Open(src).GetSpatialRef()
        if srs is not None and srs.GetAuthorityCode(None) == '4326':
            srcs.append(src)
            continue
        warped = f'{stage}/{os.path.splitext(os.path.basename(src))[0]}.4326.vrt'
        gdal.Warp(warped, src, format='VRT', dstSRS='EPSG:4326', xRes=res, yRes=res,
                  srcNodata=nodata, dstNodata=nodata, resampleAlg='average')
        srcs.append(warped)
    mosaic = f'{stage}/mosaic.vrt'
    gdal.BuildVRT(mosaic, srcs, resolution='user', xRes=res, yRes=res, srcNodata=nodata,
                  VRTNodata=nodata)
    # 2. 8-bit stretch; the source nodata becomes a mask band (transparent in the PNGs)
    byte = f'{stage}/mosaic.byte.vrt'
    gdal.Translate(byte, mosaic, format='VRT', outputType=gdal.GDT_Byte, noData=0,
                   scaleParams=[[dbMin * scale, dbMax * scale, 1, 255]], maskBand='auto')
    # VRT scaling does not clamp to 1..255: data darker than dbMin would land on 0 = nodata
    # (transparent). A LUT clamps at its end points; source nodata is applied before it.
    with open(byte) as fp:
        xml = fp.read()
    xml = re.sub(r'\s*<ScaleOffset>.*?</ScaleOffset>\s*<ScaleRatio>.*?</ScaleRatio>',
                 f'\n      <LUT>{dbMin * scale:g}:1,{dbMax * scale:g}:255</LUT>', xml)
    xml = re.sub(r'\s*<Scale>.*?</Scale>', '', xml)
    with open(byte, 'w') as fp:
        fp.write(xml)
    # 3. KML superoverlay of PNG tiles
    cmd = ['gdal2tiles.py', '--profile=geodetic', '-k', '-z', tileZoom(zoom, tileSize), '-r', resampling,
           f'--processes={processes}', '-w', 'none', '-e', '--tiledriver=PNG', f'--tilesize={tileSize}'] + \
        (['-t', title] if title else []) + [byte, out]
    print(' '.join(cmd), flush=True)
    rc = subprocess.run(cmd).returncode
    if rc == 0 and title:
        setTitle(f'{out}/doc.kml', title)        # also when gdal2tiles resumed an older product
    if rc == 0:                                  # ocean overlays at sea level, not under the water
        from .oceanAltitude import liftOcean
        liftOcean(out, processes=processes)
    if rc == 0:
        n = sum(len(f) for _, _, f in os.walk(out))
        print(f'done: {out}/doc.kml; files under {out}: {n}')
    return rc

if __name__ == '__main__':
    sys.exit(main())
