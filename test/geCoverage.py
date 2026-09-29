#!/usr/bin/env python3
'''
Opaque (data) area per zoom level of a Google Earth tile tree, in a latitude band. Every level
should show about the same area; a level with less is where the data fades out on zooming in.
  python3 test/geCoverage.py /scratch/ianj/mosaics/cycle30ascending/googleEarth [--south -80 --north -55]
Also counts KML tiles with no PNG and unreadable PNGs (truncated writes).
'''
import argparse
import glob
import math
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from osgeo import gdal

gdal.UseExceptions()


def tileArea(job):
    ''' Opaque area (km2) of one geodetic tile, or None if the PNG is unreadable. '''
    png, z, y, south, north = job
    size = 360. / 2 ** z                          # degrees per tile, gdal2tiles geodetic (not TMS-compatible)
    try:
        ds = gdal.Open(png)
        alpha = ds.GetRasterBand(ds.RasterCount).ReadAsArray()
    except RuntimeError:
        return None
    n = alpha.shape[0]
    lat = -90 + y * size + (n - 0.5 - np.arange(n)) * size / n      # row centres, north up
    inBand = (lat >= south) & (lat <= north)
    dA = (size / n * 111.32) ** 2 * np.cos(np.radians(lat))
    return float(((alpha > 0).sum(axis=1) * dA * inBand).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tree', help='googleEarth directory (<z>/<x>/<y>.png)')
    ap.add_argument('--south', type=float, default=-80.)
    ap.add_argument('--north', type=float, default=-55.)
    ap.add_argument('--processes', type=int, default=16)
    args = ap.parse_args()
    zooms = sorted(int(d) for d in os.listdir(args.tree) if d.isdigit())
    for z in zooms:
        size = 360. / 2 ** z
        y0, y1 = math.floor((args.south + 90) / size), math.ceil((args.north + 90) / size)
        jobs, noPng = [], 0
        for kml in glob.glob(f'{args.tree}/{z}/*/*.kml'):
            y = int(os.path.basename(kml)[:-4])
            if not y0 <= y < y1:
                continue
            png = kml[:-4] + '.png'
            if not os.path.exists(png):
                noPng += 1
                continue
            jobs.append((png, z, y, args.south, args.north))
        with ProcessPoolExecutor(args.processes) as pool:
            areas = list(pool.map(tileArea, jobs, chunksize=64))
        bad = sum(a is None for a in areas)
        total = sum(a for a in areas if a)
        print(f'zoom {z:2d}: {len(jobs):7d} PNGs, {noPng:5d} KML without PNG, {bad:4d} unreadable, '
              f'data {total / 1e6:8.3f} M km2', flush=True)


if __name__ == '__main__':
    main()
