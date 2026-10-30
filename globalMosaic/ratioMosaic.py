#!/usr/bin/env python3
'''
Dual-pol ratio of a like-pol and a cross-pol mosaic (openSEPPO's -dpratio).

  ratioMosaic like/vrt/gamma0/global.vrt cross/vrt/gamma0/global.vrt --out ratio.tif
  ratioMosaic like.vrt cross.vrt --out ratio.tif --mode amp

geomosaic writes gamma0 as int16 dB x 100 with no data -3000, so the ratio reduces to
arithmetic on two dB rasters:

  dB    likeDb - crossDb                                  float32, nan no data
  pwr   10 ** ((likeDb - crossDb) / 10)                    float32, nan no data
  amp   10 ** ((likeDb - crossDb) / 20) * 1000             uint16, clipped [1, 65535], 0 no data
  DN    (likeDb - crossDb) * dnSlope + 127                 uint8,  clipped [1, 255],   0 no data

The first three are exact: a power ratio is a dB difference, and an amplitude ratio is half of
one in dB. DN is not -- openSEPPO defines it on its own 8-bit rendering, whose dB-to-byte
scaling is not given in the recipe, so `--dnSlope` is explicit here and defaults to 1 DN per dB
(the reading under which the published formula is self-consistent). Use dB or amp if you want a
number you can defend.

No data follows the recipe: a pixel is no data where EITHER input is no data, and for amp/DN
that is decided BEFORE clipping, so a valid-but-saturated pixel clips to the end of the range
rather than being thrown away. The two inputs must be on the same grid -- which they are if
both mosaics were built from one `easeTiles` run, the point of its --crossPol flag.

Reads and writes in blocks, so a global 1 km mosaic does not have to fit in memory.
'''
import argparse
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()

DBSCALE = 100.0        # geomosaic's int16 gamma0 is dB x 100
DBNODATA = -3000       # ... with this sentinel
BLOCK = 2048


def openPair(likePath, crossPath):
    '''Both rasters, checked for a common grid -- a silent mismatch would ratio one place
    against another.'''
    like = gdal.Open(likePath)
    cross = gdal.Open(crossPath)
    if (like.RasterXSize, like.RasterYSize) != (cross.RasterXSize, cross.RasterYSize):
        sys.exit(f'ratioMosaic: size mismatch, like {like.RasterXSize}x{like.RasterYSize} '
                 f'vs cross {cross.RasterXSize}x{cross.RasterYSize}; the two mosaics must come '
                 f'from one tiling (easeTiles --crossPol writes it)')
    gtL, gtC = like.GetGeoTransform(), cross.GetGeoTransform()
    if max(abs(a - b) for a, b in zip(gtL, gtC)) > 1e-6:
        sys.exit(f'ratioMosaic: geotransform mismatch\n  like  {gtL}\n  cross {gtC}')
    return like, cross


def readDb(band, x, y, w, h, nodata):
    '''One block as dB float32, with no data as NaN. int16 dB x 100 is rescaled; a float band is
    taken as dB already.'''
    a = band.ReadAsArray(x, y, w, h)
    if a is None:
        sys.exit('ratioMosaic: read failed')
    if np.issubdtype(a.dtype, np.integer):
        bad = (a == (DBNODATA if nodata is None else int(nodata)))
        out = a.astype('f4') / DBSCALE
    else:
        out = a.astype('f4')
        bad = ~np.isfinite(out) if nodata is None else (out == np.float32(nodata))
    out[bad] = np.nan
    return out


MODES = {
    #        dtype   nodata  clip
    'dB':   ('f4',   np.nan, None),
    'pwr':  ('f4',   np.nan, None),
    'amp':  ('u2',   0,      (1, 65535)),
    'DN':   ('u1',   0,      (1, 255)),
}


def combine(diffDb, mode, dnSlope):
    '''The recipe, from the dB difference. Returns (values, validMask).'''
    valid = np.isfinite(diffDb)
    if mode == 'dB':
        return diffDb, valid
    if mode == 'pwr':
        return np.power(10.0, diffDb / 10.0, dtype='f4'), valid
    if mode == 'amp':
        return np.power(10.0, diffDb / 20.0, dtype='f4') * 1000.0, valid
    return diffDb * dnSlope + 127.0, valid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('like', help='like-pol mosaic (HHHH or VVVV)')
    ap.add_argument('cross', help='cross-pol mosaic (HVHV or VHVH), same grid')
    ap.add_argument('--out', required=True, help='output GeoTIFF')
    ap.add_argument('--mode', choices=list(MODES), default='dB',
                    help='pwr/amp/DN/dB as in the recipe [dB]')
    ap.add_argument('--dnSlope', type=float, default=1.0,
                    help='DN per dB for --mode DN [1.0]; see the note above')
    ap.add_argument('--likeNodata', type=float, default=None,
                    help='override the no-data value read from the like raster')
    ap.add_argument('--crossNodata', type=float, default=None,
                    help='override the no-data value read from the cross raster')
    args = ap.parse_args()

    like, cross = openPair(args.like, args.cross)
    bL, bC = like.GetRasterBand(1), cross.GetRasterBand(1)
    ndL = args.likeNodata if args.likeNodata is not None else bL.GetNoDataValue()
    ndC = args.crossNodata if args.crossNodata is not None else bC.GetNoDataValue()
    dtype, outNodata, clip = MODES[args.mode]
    nx, ny = like.RasterXSize, like.RasterYSize

    drv = gdal.GetDriverByName('GTiff')
    gdalType = {'f4': gdal.GDT_Float32, 'u2': gdal.GDT_UInt16, 'u1': gdal.GDT_Byte}[dtype]
    out = drv.Create(args.out, nx, ny, 1, gdalType,
                     options=['COMPRESS=DEFLATE', 'TILED=YES', 'BIGTIFF=IF_SAFER'])
    out.SetGeoTransform(like.GetGeoTransform())
    out.SetProjection(like.GetProjection())
    ob = out.GetRasterBand(1)
    ob.SetNoDataValue(float(outNodata) if np.isfinite(outNodata) else float('nan'))
    ob.SetDescription(f'{args.mode} ratio: like {args.like} / cross {args.cross}')

    nValid = nTotal = nClipped = 0
    for y in range(0, ny, BLOCK):
        h = min(BLOCK, ny - y)
        for x in range(0, nx, BLOCK):
            w = min(BLOCK, nx - x)
            diff = readDb(bL, x, y, w, h, ndL) - readDb(bC, x, y, w, h, ndC)
            vals, valid = combine(diff, args.mode, args.dnSlope)
            # the recipe decides no data BEFORE clamping, so a saturated pixel saturates
            # rather than disappearing
            if clip is not None:
                lo, hi = clip
                nClipped += int(np.count_nonzero(valid & ((vals < lo) | (vals > hi))))
                vals = np.clip(np.nan_to_num(vals, nan=0.0), lo, hi)
            vals = np.where(valid, vals, outNodata).astype(dtype)
            ob.WriteArray(vals, x, y)
            nValid += int(np.count_nonzero(valid))
            nTotal += valid.size
    ob.FlushCache()
    out = None
    print(f'{args.out}: {args.mode}, {nValid}/{nTotal} valid ({100 * nValid / nTotal:.1f}%)'
          + (f', {nClipped} clipped to {clip}' if clip else ''))


if __name__ == '__main__':
    main()
