#!/usr/bin/env python3
'''
Quick-look PNGs of a runGeomosaicTiles run, for viewing in a browser or Jupyter.

  makeQuickLook run_c030 [--tileRes 0.01] [--globalRes 0.05] [--processes 24]

Writes <work>/quicklooks/ (next to the progress images runGeomosaicTiles keeps there):
  tiles/<tile>.png   each tile (frequency B under A) at --tileRes deg (0.01 deg ~ 1 km)
  global.png         all lat/lon tiles at --globalRes deg (0.05 deg ~ 5 km), from the tile
                     quick looks, so the full-resolution mosaic is read only once
  cap_<name>.png     each polar cap in its own polar stereographic CRS at --capResM
Grey scale: --dbMin..--dbMax dB -> 1..255, block-averaged in dB; no data is transparent.
'''
import argparse
import concurrent.futures
import glob
import os
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()
NODATA = -3000


def stretch(a, dbMin=-24., dbMax=-1.):
    ''' int16 dB x 100 -> byte 1..255 over dbMin..dbMax; nodata -> 0. '''
    b = np.clip((a / 100. - dbMin) / (dbMax - dbMin) * 254. + 1., 1, 255).astype(np.uint8)
    b[a == NODATA] = 0
    return b


def writePng(ds, png, dbMin, dbMax):
    ''' int16 dB x 100 dataset -> 8-bit grey PNG, nodata (0) transparent. '''
    b = stretch(ds.GetRasterBand(1).ReadAsArray(), dbMin, dbMax)
    mem = gdal.GetDriverByName('MEM').Create('', b.shape[1], b.shape[0], 1, gdal.GDT_Byte)
    mem.GetRasterBand(1).WriteArray(b)
    mem.GetRasterBand(1).SetNoDataValue(0)
    gdal.GetDriverByName('PNG').CreateCopy(png, mem, options=['ZLEVEL=9'])


def reduce(src, dst, xRes, yRes, fmt='MEM'):
    return gdal.Translate(dst, src, format=fmt, xRes=xRes, yRes=yRes, resampleAlg='average',
                          noData=NODATA, outputType=gdal.GDT_Int16)


def tileLook(src, qlDir, res, dbMin, dbMax):
    name = os.path.basename(src)[:-4]
    ds = reduce(src, f'{qlDir}/tiles/{name}.tif', res, res, fmt='GTiff')
    writePng(ds, f'{qlDir}/tiles/{name}.png', dbMin, dbMax)
    return f'{qlDir}/tiles/{name}.tif'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', help='runGeomosaicTiles work directory')
    ap.add_argument('--product', default='gamma0', help='gamma0 or sigma0 [gamma0]')
    ap.add_argument('--tileRes', type=float, default=0.01, help='tile quick-look spacing, deg [0.01]')
    ap.add_argument('--globalRes', type=float, default=0.05, help='global quick-look spacing, deg [0.05]')
    ap.add_argument('--capResM', type=float, default=1000., help='polar cap quick-look spacing, m [1000]')
    ap.add_argument('--dbMin', type=float, default=-24., help='stretch minimum, dB [-24]')
    ap.add_argument('--dbMax', type=float, default=-1., help='stretch maximum, dB [-1]')
    ap.add_argument('--processes', type=int, default=24, help='parallel tile quick looks [24]')
    args = ap.parse_args()

    vrtDir = f'{args.work}/vrt/{args.product}'
    qlDir = f'{args.work}/quicklooks'
    os.makedirs(f'{qlDir}/tiles', exist_ok=True)
    everything = sorted(glob.glob(f'{vrtDir}/*.vrt'))
    caps = [v for v in everything if os.path.basename(v).startswith('cap_')]
    tiles = [v for v in everything if v not in caps and not os.path.basename(v).startswith('band_')
             and os.path.basename(v) != 'global.vrt']
    with concurrent.futures.ProcessPoolExecutor(args.processes) as pool:
        small = list(pool.map(tileLook, tiles, [qlDir] * len(tiles), [args.tileRes] * len(tiles),
                              [args.dbMin] * len(tiles), [args.dbMax] * len(tiles)))
    print(f'{len(small)} tile quick looks in {qlDir}/tiles')
    if small:
        gdal.BuildVRT(f'{qlDir}/tiles.vrt', small, srcNodata=NODATA, VRTNodata=NODATA)
        writePng(reduce(f'{qlDir}/tiles.vrt', '', args.globalRes, args.globalRes),
                 f'{qlDir}/global.png', args.dbMin, args.dbMax)
        print(f'{qlDir}/global.png')
    for cap in caps:
        png = f'{qlDir}/{os.path.basename(cap)[:-4]}.png'
        writePng(reduce(cap, '', args.capResM, args.capResM), png, args.dbMin, args.dbMax)
        print(png)
    return 0


if __name__ == '__main__':
    sys.exit(main())
