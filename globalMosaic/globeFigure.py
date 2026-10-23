#!/usr/bin/env python3
'''
Draw a tiled mosaic as a globe on its own projection: the whole Earth as a filled ocean, with
the backscatter painted on top wherever there is data.

  globeFigure /scratch/ianj/mosaics/ease/global/vrt/gamma0/ease_proj.vrt --out globe.png
  globeFigure .../global/tiles/gamma0 --out globe.png --width 4000 --coast
  globeFigure .../ease_6933.vrt --out g.png --ocean '#0b3d66' --dbMin -26 --dbMax -2

Why a separate figure from the driver's quicklook: the quicklook is a bare raster, so
everywhere without data is transparent and the map has no shape -- a NISAR mosaic covers land
and sea ice, so the result reads as scattered fragments rather than a world. Here the
projection's own boundary is drawn first and filled, which puts the fragments on a globe.

The boundary is the image of the graticule edge -- lon = +/-180 from pole to pole, plus the
two pole lines -- pushed through the projection's FORWARD transform. Forward is all that is
needed, which matters because Wagner VII (ESRI:54076) has no inverse in PROJ.

The raster can be a VRT, a single GeoTIFF, or a directory of tiles (a VRT is built from
*.tif in it). It is read once, decimated to the output width with averaging, and stretched
from the stored int16 (dB x 100) to grey.
'''
import argparse
import glob
import os

import numpy as np

NODATA = -3000


def mapBoundary(crs, nLat=721, nLon=1441):
    ''' The projection's outline in projected metres: the globe's edge, forward-projected.
    Returned closed, so it can be filled. '''
    import pyproj
    fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
    lat = np.linspace(-90., 90., nLat)
    lon = np.linspace(-180., 180., nLon)
    # down the dateline, along the north pole, back up the antimeridian, along the south pole
    lons = np.concatenate([np.full(nLat, -180.), lon, np.full(nLat, 180.), lon[::-1]])
    lats = np.concatenate([lat, np.full(nLon, 90.), lat[::-1], np.full(nLon, -90.)])
    x, y = fwd(lons, lats)
    x, y = np.asarray(x), np.asarray(y)
    good = np.isfinite(x) & np.isfinite(y)
    return x[good], y[good]


def readRaster(src, width):
    ''' (grey array, extent) from a VRT, a GeoTIFF or a directory of tiles, decimated to
    `width` columns by AVERAGING -- these mosaics are narrow strips of coverage, and nearest
    sampling at 20x throws most of them away. '''
    from osgeo import gdal
    gdal.UseExceptions()
    if os.path.isdir(src):
        tifs = sorted(glob.glob(f'{src}/*.tif'))
        if not tifs:
            raise SystemExit(f'globeFigure: no .tif in {src}')
        src = gdal.BuildVRT('', tifs, VRTNodata=NODATA)
    ds = gdal.Open(src) if isinstance(src, str) else src
    scale = width / ds.RasterXSize
    # Translate, not Warp: the CRS is unchanged, and Warp builds a coordinate transformation
    # even for a plain resize -- which fails on Wagner VII, whose inverse PROJ does not have
    # ("Inverse operation for wag7 is not available").
    out = gdal.Translate('', ds, format='MEM', resampleAlg='average',
                         width=width, height=max(int(round(ds.RasterYSize * scale)), 1),
                         noData=NODATA)
    a = out.ReadAsArray().astype('f4')
    if a.ndim == 3:
        a = a[0]
    gt = out.GetGeoTransform()
    crs = out.GetSpatialRef()
    ny, nx = a.shape
    ext = [gt[0], gt[0] + nx * gt[1], gt[3] + ny * gt[5], gt[3]]
    del out, ds
    return np.where(a == NODATA, np.nan, a), ext, crs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('raster', help='mosaic VRT, GeoTIFF, or a directory of tiles')
    ap.add_argument('--out', default='globe.png', help='output image [globe.png]')
    ap.add_argument('--width', type=int, default=3000, help='output width, pixels [3000]')
    ap.add_argument('--dbMin', type=float, default=-24., help='black at this dB [-24]')
    ap.add_argument('--dbMax', type=float, default=-1., help='white at this dB [-1]')
    ap.add_argument('--ocean', default='#1b4f72', help="ocean fill colour ['#1b4f72']")
    ap.add_argument('--edge', default='#0d2b3e', help='map outline colour')
    ap.add_argument('--background', default='white', help='colour outside the map [white]')
    ap.add_argument('--coast', action='store_true', help='draw coastlines over the ocean')
    ap.add_argument('--graticule', type=float, default=0.,
                    help='graticule spacing in degrees, 0 for none [0]')
    ap.add_argument('--crs', default=None,
                    help="CRS for the outline; default is the raster's own (a .wkt file is "
                         'read if given)')
    ap.add_argument('--title', default=None, help='title drawn above the globe')
    ap.add_argument('--dpi', type=int, default=100)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    grey, ext, rasterCrs = readRaster(args.raster, args.width)
    crs = args.crs
    if crs is None:
        crs = rasterCrs.ExportToWkt()
    elif os.path.isfile(crs):
        crs = open(crs).read()

    bx, by = mapBoundary(crs)
    fig, ax = plt.subplots(figsize=(args.width / args.dpi,
                                    args.width / args.dpi * (np.ptp(by) / np.ptp(bx))))
    fig.patch.set_facecolor(args.background)
    ax.set_facecolor(args.background)
    # ocean first: the globe's shape, so land and sea ice sit ON something
    ax.fill(bx, by, facecolor=args.ocean, edgecolor='none', zorder=1)
    if args.coast or args.graticule:
        overlays(ax, crs, args)
    with np.errstate(invalid='ignore'):
        g = np.clip((grey / 100. - args.dbMin) / (args.dbMax - args.dbMin), 0., 1.)
    rgba = np.dstack([g, g, g, np.where(np.isfinite(g), 1., 0.)])
    im = ax.imshow(rgba, extent=ext, interpolation='nearest', zorder=3)
    # Clip to the map: tiles built before the Wagner VII inverse rejected |lon| > 180 carry a
    # mirrored sliver of data outside the lens, and nothing outside it is meaningful anyway.
    from matplotlib.patches import Polygon
    im.set_clip_path(Polygon(np.column_stack([bx, by]), closed=True,
                             transform=ax.transData, facecolor='none', edgecolor='none'))
    ax.plot(bx, by, color=args.edge, linewidth=0.8, zorder=4)
    ax.set_xlim(bx.min(), bx.max())
    ax.set_ylim(by.min(), by.max())
    ax.set_aspect('equal')
    ax.axis('off')
    if args.title:
        ax.set_title(args.title, fontsize=13)
    fig.tight_layout(pad=0.2)
    fig.savefig(args.out, dpi=args.dpi, facecolor=fig.get_facecolor())
    valid = np.isfinite(grey).mean() * 100.
    print(f'{args.out}: {grey.shape[1]} x {grey.shape[0]} px, {valid:.1f}% of the frame has data')


def overlays(ax, crs, args):
    ''' Coastlines and graticule, forward-projected. Segments that jump the dateline are
    broken, or they draw a line straight across the map. '''
    import pyproj
    fwd = pyproj.Transformer.from_crs(4326, crs, always_xy=True).transform
    span = 0.25 * (np.abs(mapBoundary(crs)[0]).max())

    def draw(lon, lat, **kw):
        x, y = fwd(np.asarray(lon, float), np.asarray(lat, float))
        x, y = np.asarray(x, float).copy(), np.asarray(y, float).copy()
        bad = ~np.isfinite(x) | ~np.isfinite(y)
        bad[1:] |= np.abs(np.diff(x)) > span
        x[bad], y[bad] = np.nan, np.nan
        ax.plot(x, y, **kw)

    if args.coast:
        try:
            import cartopy.feature as cf
            for geom in cf.COASTLINE.with_scale('110m').geometries():
                for line in (geom.geoms if hasattr(geom, 'geoms') else [geom]):
                    lon, lat = np.array(line.coords).T
                    draw(lon, lat, color='#9fb8c8', linewidth=0.3, zorder=2)
        except ImportError:
            print('globeFigure: cartopy not installed, skipping coastlines')
    if args.graticule:
        step = args.graticule
        for lon in np.arange(-180., 180.0001, step):
            draw(np.full(181, lon), np.linspace(-90., 90., 181),
                 color='#ffffff', linewidth=0.25, alpha=0.35, zorder=2)
        for lat in np.arange(-90. + step, 90., step):
            draw(np.linspace(-180., 180., 361), np.full(361, lat),
                 color='#ffffff', linewidth=0.25, alpha=0.35, zorder=2)


if __name__ == '__main__':
    main()
