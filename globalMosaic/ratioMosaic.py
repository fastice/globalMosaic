#!/usr/bin/env python3
'''
Dual-pol ratio and RGB composite from a like-pol and a cross-pol mosaic (openSEPPO -dpratio).

  ratioMosaic like/vrt/gamma0/ease_proj.vrt cross/vrt/gamma0/ease_proj.vrt --out ratio.tif
  ratioMosaic LIKE CROSS --mode rgb --out browse.tif       # the purple/yellow false colour

geomosaic writes gamma0 as int16 dB x 100 with no data -3000, so everything here is arithmetic
on two dB rasters.

  dB    likeDb - crossDb                           float32, nan no data
  pwr   10 ** ((likeDb - crossDb) / 10)            float32, nan no data
  amp   10 ** ((likeDb - crossDb) / 20) * 1000     uint16, clipped [1, 65535], 0 no data
  DN    DN(like) - DN(cross) + 127                 uint8,  clipped [1, 255],   0 no data
  rgb   R = DN(like), G = DN(cross), B = the DN ratio above

`rgb` is the published rendering, and it is a composite rather than a colour map: bare ground
comes out deep indigo (cross-pol at the noise floor, so B dominates), vegetated ground greener
and less blue as the ratio falls. WHERE THERE IS NO CROSS-POL -- single-pol SH/SV
acquisitions, which is how the ice sheets are taken -- all three bands are set to DN(like), so
those areas render as GREYSCALE backscatter instead of dropping out.

Note what the blue channel can do. Cross-pol is always weaker than like-pol, so
DN(like) >= DN(cross) and B = DN(like) - DN(cross) + 127 is confined to [127, 255]: blue never
drops below mid scale. So a raw DN composite cannot render forest as pale yellow -- that needs
B dimmer than R and G, which this arithmetic forbids. The published map's pastel look must
therefore come from a per-band DISPLAY stretch applied afterwards, which openSEPPO does not
document. This writes the DN bands faithfully and leaves the stretch to the renderer.

DN follows openSEPPO's power_to_dn_uint8 exactly: a linear dB scale,

    dB = -31.15 + DN * 0.15        i.e. minDb -31.0 -> DN 1, maxDb 7.1 -> DN 255

clamped at both ends, with DN 0 reserved for no data. DN and rgb take the difference of the
CLAMPED DNs, which is not the same as scaling the dB difference: cross-pol below -31 dB pins
to DN 1, which is routine over desert and water, and the two then disagree. dB, pwr and amp
are plain dB arithmetic and exact.

The two rasters must share a CRS and pixel size. They need NOT share an extent -- a ratio
cover has fewer tiles than a full-cover like-pol mosaic, so the cross-pol mosaic is usually
smaller. The smaller one is aligned onto the other's grid with a nearest-neighbour VRT, which
is exact because both come from one tiling and sit on the same pixel lattice.

Reads and writes in blocks, so a global 1 km mosaic need not fit in memory.
'''
import argparse
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()

DBSCALE = 100.0        # geomosaic's int16 gamma0 is dB x 100
DBNODATA = -3000       # ... with this sentinel
BLOCK = 2048
# openSEPPO power_to_dn_uint8 defaults
MINDB, MAXDB, MINDN, MAXDN = -31.0, 7.1, 1, 255

#        bands dtype nodata  clip
MODES = {
    'dB':  (1, 'f4', np.nan, None),
    'pwr': (1, 'f4', np.nan, None),
    'amp': (1, 'u2', 0,      (1, 65535)),
    'DN':  (1, 'u1', 0,      (1, 255)),
    'rgb': (3, 'u1', 0,      (1, 255)),
}


def align(ref, other):
    '''`other` resampled onto `ref`'s exact grid, if it is not already on it. Nearest
    neighbour, which is lossless here: both mosaics come from one tiling, so they share a
    pixel lattice and this only trims or pads the extent.'''
    gtR, gtO = ref.GetGeoTransform(), other.GetGeoTransform()
    if (ref.RasterXSize, ref.RasterYSize) == (other.RasterXSize, other.RasterYSize) \
            and max(abs(a - b) for a, b in zip(gtR, gtO)) < 1e-6:
        return other
    if abs(gtR[1] - gtO[1]) > 1e-6 or abs(gtR[5] - gtO[5]) > 1e-6:
        sys.exit(f'ratioMosaic: pixel sizes differ, {gtR[1]}x{gtR[5]} vs {gtO[1]}x{gtO[5]}; '
                 'the two mosaics must be built at the same --psResM')
    x0, y0 = gtR[0], gtR[3]
    x1 = x0 + gtR[1] * ref.RasterXSize
    y1 = y0 + gtR[5] * ref.RasterYSize
    nd = other.GetRasterBand(1).GetNoDataValue()
    print(f'  aligning {other.RasterXSize}x{other.RasterYSize} onto '
          f'{ref.RasterXSize}x{ref.RasterYSize}', flush=True)
    return gdal.Warp('', other, format='VRT', outputBounds=(min(x0, x1), min(y0, y1),
                                                            max(x0, x1), max(y0, y1)),
                     xRes=abs(gtR[1]), yRes=abs(gtR[5]), resampleAlg='near',
                     srcNodata=nd, dstNodata=nd)


def readDb(band, x, y, w, h, nodata):
    '''One block as dB float32, no data as NaN. int16 dB x 100 is rescaled; a float band is
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


def dbToDn(db, minDb=MINDB, maxDb=MAXDB):
    '''openSEPPO's power_to_dn_uint8, from dB instead of linear power. DN 0 is no data; values
    outside [minDb, maxDb] clamp to the ends rather than dropping out.'''
    scale = (MAXDN - MINDN) / (maxDb - minDb)
    offset = MINDN - minDb * scale
    valid = np.isfinite(db)
    raw = db * scale + offset
    dn = np.zeros(db.shape, dtype='u1')
    inRange = valid & (raw >= MINDN) & (raw <= MAXDN)
    dn[inRange] = np.round(raw[inRange]).astype('u1')
    dn[valid & (raw < MINDN)] = MINDN
    dn[valid & (raw > MAXDN)] = MAXDN
    return dn


def samplePercentiles(bL, bC, ndL, ndC, mode, minDb, maxDb, lo, hi, maxPx=4000000):
    """Per-band stretch limits, from a decimated read of the whole mosaic.

    Decimated rather than full: the percentiles of a 1-in-N sample of a global raster are
    indistinguishable from the real ones for this purpose, and it avoids a second full pass.
    """
    nx, ny = bL.XSize, bL.YSize
    step = max(1, int(round((nx * ny / maxPx) ** 0.5)))
    w, h = max(1, nx // step), max(1, ny // step)
    def readSmall(b, nd):
        a = b.ReadAsArray(0, 0, nx, ny, buf_xsize=w, buf_ysize=h)
        if np.issubdtype(a.dtype, np.integer):
            bad = (a == (DBNODATA if nd is None else int(nd)))
            out = a.astype('f4') / DBSCALE
        else:
            out = a.astype('f4')
            bad = ~np.isfinite(out) if nd is None else (out == np.float32(nd))
        out[bad] = np.nan
        return out
    bands, valid = combine(readSmall(bL, ndL), readSmall(bC, ndC), mode, minDb, maxDb)
    lims = []
    for v in bands:
        s = v[valid & np.isfinite(v) & (v > 0)]
        lims.append((float(np.percentile(s, lo)), float(np.percentile(s, hi)))
                    if s.size else (1.0, 255.0))
    print(f'  stretch from a {w}x{h} sample ({lo}-{hi} percentile): '
          + ', '.join(f'{a:.0f}-{b:.0f}' for a, b in lims), flush=True)
    return lims


def applyStretch(v, lim, clip):
    """Linear stretch of one band onto the output range."""
    lo, hi = lim
    if hi <= lo:
        return v
    return (v - lo) / (hi - lo) * (clip[1] - clip[0]) + clip[0]


def combine(likeDb, crossDb, mode, minDb, maxDb):
    '''The recipe. Returns (list of bands, validMask).'''
    if mode in ('dB', 'pwr', 'amp'):
        d = likeDb - crossDb
        valid = np.isfinite(d)
        if mode == 'dB':
            return [d], valid
        if mode == 'pwr':
            return [np.power(10.0, d / 10.0, dtype='f4')], valid
        return [np.power(10.0, d / 20.0, dtype='f4') * 1000.0], valid
    # DN and rgb: difference the CLAMPED DNs, as openSEPPO does
    dnL = dbToDn(likeDb, minDb, maxDb)
    dnC = dbToDn(crossDb, minDb, maxDb)
    ratio = dnL.astype('f4') - dnC.astype('f4') + 127.0
    both = (dnL != 0) & (dnC != 0)
    if mode == 'DN':
        return [ratio], both
    # rgb: grey where there is no cross-pol, so single-pol ground still renders
    greyOnly = (dnL != 0) & (dnC == 0)
    r = dnL.astype('f4')
    g = np.where(both, dnC.astype('f4'), np.where(greyOnly, r, 0.0))
    b = np.where(both, ratio, np.where(greyOnly, r, 0.0))
    return [r, g, b], (dnL != 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('like', help='like-pol mosaic (HHHH or VVVV). For --mode rgb use the '
                                 'FULL-cover mosaic, so single-pol ground renders grey')
    ap.add_argument('cross', help='cross-pol mosaic (HVHV or VHVH)')
    ap.add_argument('--out', required=True, help='output GeoTIFF')
    ap.add_argument('--mode', choices=list(MODES), default='dB',
                    help='pwr/amp/DN/dB, or rgb for the false-colour composite [dB]')
    ap.add_argument('--minDb', type=float, default=MINDB,
                    help=f'dB mapped to DN {MINDN} [{MINDB}]')
    ap.add_argument('--maxDb', type=float, default=MAXDB,
                    help=f'dB mapped to DN {MAXDN} [{MAXDB}]')
    ap.add_argument('--stretch', type=float, nargs=2, metavar=('LO', 'HI'), default=None,
                    help='per-band percentile stretch for --mode rgb, e.g. --stretch 2 98. '
                         'OFF by default, which writes raw DN. The published map applies a '
                         'stretch from SEPPO (a separate package from openSEPPO), which is '
                         'not documented, so there is nothing to reproduce -- this is for '
                         'dialling in an appearance. Note it makes the output a BROWSE image, '
                         'not DN data.')
    ap.add_argument('--likeNodata', type=float, default=None)
    ap.add_argument('--crossNodata', type=float, default=None)
    args = ap.parse_args()

    like = gdal.Open(args.like)
    cross = align(like, gdal.Open(args.cross))
    bL, bC = like.GetRasterBand(1), cross.GetRasterBand(1)
    ndL = args.likeNodata if args.likeNodata is not None else bL.GetNoDataValue()
    ndC = args.crossNodata if args.crossNodata is not None else bC.GetNoDataValue()
    nBands, dtype, outNodata, clip = MODES[args.mode]
    nx, ny = like.RasterXSize, like.RasterYSize

    gdalType = {'f4': gdal.GDT_Float32, 'u2': gdal.GDT_UInt16, 'u1': gdal.GDT_Byte}[dtype]
    out = gdal.GetDriverByName('GTiff').Create(
        args.out, nx, ny, nBands, gdalType,
        options=['COMPRESS=DEFLATE', 'TILED=YES', 'BIGTIFF=IF_SAFER',
                 'PHOTOMETRIC=RGB'] if nBands == 3 else
                ['COMPRESS=DEFLATE', 'TILED=YES', 'BIGTIFF=IF_SAFER'])
    out.SetGeoTransform(like.GetGeoTransform())
    out.SetProjection(like.GetProjection())
    names = ('like (R)', 'cross (G)', 'ratio (B)') if nBands == 3 else (f'{args.mode} ratio',)
    for i in range(nBands):
        ob = out.GetRasterBand(i + 1)
        ob.SetNoDataValue(float(outNodata) if np.isfinite(outNodata) else float('nan'))
        ob.SetDescription(names[i])

    lims = None
    if args.stretch is not None:
        if args.mode != 'rgb':
            sys.exit('ratioMosaic: --stretch applies to --mode rgb only')
        lims = samplePercentiles(bL, bC, ndL, ndC, args.mode, args.minDb, args.maxDb,
                                 args.stretch[0], args.stretch[1])

    nValid = nTotal = nGrey = 0
    for y in range(0, ny, BLOCK):
        h = min(BLOCK, ny - y)
        for x in range(0, nx, BLOCK):
            w = min(BLOCK, nx - x)
            lk = readDb(bL, x, y, w, h, ndL)
            cr = readDb(bC, x, y, w, h, ndC)
            vals, valid = combine(lk, cr, args.mode, args.minDb, args.maxDb)
            grey = None
            if args.mode == 'rgb':
                grey = np.isfinite(lk) & ~np.isfinite(cr)
                nGrey += int(np.count_nonzero(grey))
            for i, v in enumerate(vals):
                if lims is not None:
                    # Single-pol ground has R == G == B, and stretching each band by its own
                    # limits would tint it -- the greyscale fallback would come out coloured.
                    # So those pixels take the RED band's stretch in all three.
                    v = np.where(grey, applyStretch(v, lims[0], clip),
                                 applyStretch(v, lims[i], clip))
                if clip is not None:
                    v = np.clip(np.nan_to_num(v, nan=0.0), *clip)
                out.GetRasterBand(i + 1).WriteArray(
                    np.where(valid, v, outNodata).astype(dtype), x, y)
            nValid += int(np.count_nonzero(valid))
            nTotal += valid.size
    out.FlushCache()
    out = None
    msg = f'{args.out}: {args.mode}, {nValid}/{nTotal} valid ({100 * nValid / nTotal:.1f}%)'
    if args.mode == 'rgb':
        msg += f', {nGrey} single-pol pixels rendered grey'
    print(msg)


if __name__ == '__main__':
    main()
