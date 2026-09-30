#!/usr/bin/env python3
'''
Land that no granule of the given catalogues covers in one orbit direction -- for the acquisition
planners: these gaps cannot be filled from the archive, only by acquiring.

  notAcquired --catalogues cycle30/catalogue0*.geojson cycle30/catalogue.geojson \\
      --direction A --out cycle30/notAcquiredAscending

Acquisition footprints (not valid-data ones: the question is what was acquired) of every granule
with a usable channel are rasterised at --cellDeg; land is Natural Earth 10 m land + Antarctic ice
shelves, south of --north (the tiling's limit). Connected uncovered areas of at least --minKm2 are
listed with country, nearest named place and the share the other direction covers:
  <out>/notAcquired.csv, notAcquired.geojson (the areas), notAcquired.png (map)
'''
import argparse
import csv
import json
import os
import sys

import numpy as np
from rasterio import features
from rasterio.transform import from_origin
from scipy import ndimage
from shapely.geometry import Point, mapping, shape
from shapely.ops import unary_union
from shapely.strtree import STRtree

from .globalGCOVTiles import loadGranules


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalogues', nargs='+', required=True)
    ap.add_argument('--direction', choices=['A', 'D'], default='A')
    ap.add_argument('--out', required=True)
    ap.add_argument('--cellDeg', type=float, default=0.02, help='raster cell, deg [0.02, ~2 km]')
    ap.add_argument('--minKm2', type=float, default=4., help='smallest area listed, km2 [4]')
    ap.add_argument('--north', type=float, default=77.5, help='northern limit, deg [77.5, the tiling]')
    args = ap.parse_args()
    import cartopy.feature as cf
    import cartopy.io.shapereader as shpreader
    D = args.cellDeg
    nx, ny = int(round(360 / D)), int(round(180 / D))
    tr = from_origin(-180, 90, D, D)
    grans = []
    for c in args.catalogues:
        grans += loadGranules(c, allowV=True)
    other = 'D' if args.direction == 'A' else 'A'

    def cover(direction):
        return features.rasterize(((g['geom'], 1) for g in grans if g['direction'] == direction),
                                  out_shape=(ny, nx), transform=tr, dtype='uint8').astype(bool)
    have, otherHave = cover(args.direction), cover(other)
    landGeoms = list(cf.LAND.with_scale('10m').geometries()) + \
        list(cf.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '10m').geometries())
    land = features.rasterize(((g, 1) for g in landGeoms), out_shape=(ny, nx), transform=tr,
                              dtype='uint8').astype(bool)
    lat = 90 - (np.arange(ny) + 0.5) * D
    cellKm2 = np.broadcast_to((D * 111.32) ** 2 * np.cos(np.radians(lat))[:, None], (ny, nx))
    gap = land & ~have
    gap[lat > args.north, :] = False
    lab, n = ndimage.label(gap, structure=np.ones((3, 3)))
    idx = np.arange(1, n + 1)
    areas = ndimage.sum(cellKm2, lab, idx)
    otherKm2 = ndimage.sum(cellKm2 * otherHave, lab, idx)
    objs = ndimage.find_objects(lab)
    countries = [(r.geometry, r.attributes['NAME_LONG']) for r in shpreader.Reader(
        shpreader.natural_earth('10m', 'cultural', 'admin_0_countries')).records()]
    cTree = STRtree([g for g, _ in countries])
    places = [(Point(r.geometry.x, r.geometry.y), r.attributes['NAME']) for r in shpreader.Reader(
        shpreader.natural_earth('10m', 'cultural', 'populated_places')).records()]
    pTree = STRtree([p for p, _ in places])
    rows, feats = [], []
    for i in np.argsort(-areas):
        a = float(areas[i])
        if a < args.minKm2:
            break
        sl = objs[i]
        m = lab[sl] == i + 1
        geom = unary_union([shape(g) for g, v in features.shapes(
            m.astype('uint8'), mask=m, transform=from_origin(-180 + sl[1].start * D, 90 - sl[0].start * D, D, D)) if v])
        c = geom.representative_point()
        names = sorted({countries[j][1] for j in cTree.query(geom, predicate='intersects')})
        if not names:
            names = [countries[cTree.nearest(c)][1]]
        k = pTree.nearest(c)
        x0, y0, x1, y1 = geom.bounds
        rows.append(dict(areaKm2=round(a), country=', '.join(names),
                         nearestPlace=f'{places[k][1]} ({places[k][0].distance(c) * 111.32:.0f} km)',
                         lon=round(c.x, 3), lat=round(c.y, 3), west=round(x0, 2), south=round(y0, 2),
                         east=round(x1, 2), north=round(y1, 2),
                         otherDirectionPct=round(100 * float(otherKm2[i]) / a)))
        feats.append({'type': 'Feature', 'geometry': mapping(geom.simplify(D / 2)), 'properties': rows[-1]})
    os.makedirs(args.out, exist_ok=True)
    json.dump({'type': 'FeatureCollection', 'features': feats}, open(f'{args.out}/notAcquired.geojson', 'w'))
    with open(f'{args.out}/notAcquired.csv', 'w') as fp:
        w = csv.DictWriter(fp, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    total = sum(r['areaKm2'] for r in rows)
    print(f'{len(rows)} land areas >= {args.minKm2:g} km2 without {args.direction} acquisitions, '
          f'{total:,} km2 -> {args.out}/notAcquired.csv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    fig = plt.figure(figsize=(18, 9))
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.add_feature(cf.LAND.with_scale('50m'), color='#e6e6e6')
    ax.coastlines('50m', lw=0.3)
    ax.add_geometries([shape(f['geometry']) for f in feats], ccrs.PlateCarree(), facecolor='#d62728',
                      edgecolor='#d62728', lw=0.8)
    big = [r for r in rows if r['areaKm2'] < 200000][:40]
    ax.scatter([r['lon'] for r in rows], [r['lat'] for r in rows], s=12, facecolor='none',
               edgecolor='#d62728', lw=0.6, transform=ccrs.PlateCarree())
    for r in big:
        ax.annotate(r['nearestPlace'].split(' (')[0], (r['lon'], r['lat']), fontsize=6, xytext=(3, 3),
                    textcoords='offset points', transform=ccrs.PlateCarree())
    ax.set_title(f'Land with no {"ascending" if args.direction == "A" else "descending"} acquisition in '
                 f'{len(args.catalogues)} catalogues: {len(rows)} areas, {total:,} km2 (circles mark each area)')
    fig.savefig(f'{args.out}/notAcquired.png', dpi=110, bbox_inches='tight')
    return 0


if __name__ == '__main__':
    sys.exit(main())
