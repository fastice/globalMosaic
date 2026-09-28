#!/usr/bin/env python3
'''
Coverage map for a globalGCOVTiles.py run: the footprints of the granules it selected, coloured by
the bandwidth of the channel used, with land that no selected granule covers in red. Three
panels: mid-latitudes (plate carree), Arctic and Antarctic (polar stereographic).

  coverageMaps RUNDIR [--title TEXT] [--latMax 77.5]   -> RUNDIR/coverage.png, coverage.txt
'''
import argparse
import csv
import glob
import json
import math
import os

import matplotlib
matplotlib.use('Agg')
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from shapely.geometry import box, shape
from shapely.validation import make_valid
from shapely.ops import unary_union

from globalMosaic import globalGCOVTiles as gtiles

# opaque colours: region layers meet at 60 deg and a transparent fill would darken the overlap
TIERS = [('77 / 40 MHz', lambda bw: bw >= 40, '#6f93c9'), ('20 MHz', lambda bw: bw == 20, '#a9cfe9'),
         ('5 MHz', lambda bw: bw == 5, '#f1c98a')]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runDir')
    ap.add_argument('--title', default=None)
    ap.add_argument('--latMax', type=float, default=77.5, help='northern limit of the tiling [77.5]')
    ap.add_argument('--validFootprints', nargs='+', default=None,
                    help='scanMasks outputs [validFootprints*.geojson in the run dir or its parent]')
    args = ap.parse_args()
    # the run's own catalogue, else the shared one in its parent (e.g. cycle30/catalogue.geojson)
    here = f'{args.runDir}/catalogue.geojson'
    cat = json.load(open(here if os.path.exists(here) else
                         f'{os.path.dirname(os.path.abspath(args.runDir))}/catalogue.geojson'))
    geoms = {f['properties']['name']: f['geometry'] for f in cat['features']}
    # fill-cycle catalogues (catalogueNNN.geojson) next to the main one
    catDir = args.runDir if os.path.exists(here) else os.path.dirname(os.path.abspath(args.runDir))
    for extra in sorted(glob.glob(f'{catDir}/catalogue[0-9][0-9][0-9].geojson')):
        for f in json.load(open(extra))['features']:
            geoms.setdefault(f['properties']['name'], f['geometry'])
    # valid-data footprints (scanMasks) replace acquisition footprints where available, so the
    # map shows what geomosaic will actually fill (mask 0 = partially focused is dropped)
    vfPaths = args.validFootprints or (sorted(glob.glob(f'{args.runDir}/validFootprints*.geojson')) or
        sorted(glob.glob(f'{os.path.dirname(os.path.abspath(args.runDir))}/validFootprints*.geojson')))
    for vfPath in vfPaths:
        for f in json.load(open(vfPath))['features']:
            geoms[f['properties']['name']] = f['geometry']
        print(f'valid-data footprints from {vfPath}')
    geoms = {k: v for k, v in geoms.items() if v is not None}
    sel = {}
    for f in glob.glob(f'{args.runDir}/tiles/*.csv'):
        for row in csv.DictReader(open(f)):
            sel[row['name']] = int(row['bwMHz'])
    # Coverage and gaps per region, each in its own projection: polar regions in the polar
    # Lambert equal-area grids (EASE 2.0, EPSG 6932/6931), so footprints near the pole and across
    # the antimeridian are not distorted by lat/lon; mid-latitudes in lat/lon.
    import pyproj
    from shapely.ops import transform
    land = unary_union(list(cf.LAND.with_scale('50m').geometries()) + list(cf.NaturalEarthFeature(
        'physical', 'antarctic_ice_shelves_polys', '50m').geometries())).buffer(0)
    regions = [('Antarctic (<60S)', -90, -60, 6932, ccrs.LambertAzimuthalEqualArea(central_latitude=-90)),
               ('mid-latitude', -60, 60, None, ccrs.PlateCarree()),
               ('Arctic (>60N)', 60, args.latMax, 6931, ccrs.LambertAzimuthalEqualArea(central_latitude=90))]
    eqArea = pyproj.Transformer.from_crs(4326, ccrs.EckertIV().proj4_init, always_xy=True).transform
    result = {}
    for rname, s0, n0, epsg, crs in regions:
        clip = box(-180, s0, 180, n0)
        if epsg:
            fwd = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform
            toGeom = lambda g: make_valid(transform(fwd, shape(g)))     # corners projected directly
            # a latitude band is a disk or annulus about the pole in a polar projection
            from shapely.geometry import Point
            radius = lambda lat: math.hypot(*fwd(0., lat))
            poleward, equatorward = (s0, n0) if s0 < 0 else (n0, s0)
            clipR = Point(0, 0).buffer(radius(equatorward), 720)
            if abs(poleward) < 90:
                clipR = clipR.difference(Point(0, 0).buffer(radius(poleward), 720))
            # densify first: clipped edges along a parallel would become straight chords
            landR = transform(fwd, land.intersection(clip).segmentize(0.25)).buffer(200).buffer(-200)
            landR = landR.intersection(clipR)
            eps, areaOf = 2000., (lambda g: g.area / 1e6)
        else:
            toGeom = lambda g: gtiles.unwrapValid(shape(g))
            landR, clipR = land.intersection(clip), clip
            eps, areaOf = 0.02, (lambda g: transform(eqArea, g).area / 1e6)
        tiers = {}
        for label, test, _ in TIERS:
            polys = [toGeom(geoms[n]) for n, bw in sel.items() if test(bw) and n in geoms]
            tiers[label] = unary_union(polys).buffer(0).intersection(clipR) if polys else None
        cov = unary_union([g for g in tiers.values() if g is not None]).buffer(0)
        # drop zero-width slivers narrower than ~2 km; real gaps survive the opening
        gaps = landR.difference(cov).buffer(-eps).buffer(eps).intersection(landR)
        result[rname] = dict(crs=crs, tiers=tiers, gaps=gaps, land=landR, aLand=areaOf(landR),
                             aGap=areaOf(gaps))
    aLand = sum(r['aLand'] for r in result.values())
    aGap = sum(r['aGap'] for r in result.values())
    stats = [f'{len(sel)} granules selected; land (to {args.latMax}N, incl. Antarctic ice shelves) '
             f'{aLand / 1e6:.1f} M km2; uncovered land {aGap / 1e6:.2f} M km2 ({100 * aGap / aLand:.1f}%)']
    for rname, r in result.items():
        stats.append(f'  {rname:18s} uncovered {r["aGap"] / 1e6:6.2f} of {r["aLand"] / 1e6:6.1f} M km2 '
                     f'({100 * r["aGap"] / max(r["aLand"], 1):.1f}%)')
    open(f'{args.runDir}/coverage.txt', 'w').write('\n'.join(stats) + '\n')

    fig = plt.figure(figsize=(16, 13), dpi=110)
    panels = [(fig.add_axes([0.02, 0.50, 0.96, 0.41], projection=ccrs.PlateCarree()),
               [-180, 180, -60, args.latMax], 'Mid-latitudes (60S - %.1fN)' % args.latMax),
              (fig.add_axes([0.04, 0.03, 0.44, 0.44], projection=ccrs.NorthPolarStereo()),
               [-180, 180, 55, 90], 'Arctic'),
              (fig.add_axes([0.52, 0.03, 0.44, 0.44], projection=ccrs.SouthPolarStereo()),
               [-180, 180, -90, -55], 'Antarctic')]
    for ax, ext, name in panels:
        ax.set_extent(ext, crs=ccrs.PlateCarree())
        ax.add_feature(cf.OCEAN.with_scale('110m'), facecolor='#eef2f5', zorder=0)
        for r in result.values():
            ax.add_geometries([r['land']], r['crs'], facecolor='#d9d9d9', edgecolor='none', zorder=1)
            for label, _, color in TIERS:
                if r['tiers'][label] is not None and not r['tiers'][label].is_empty:
                    ax.add_geometries([r['tiers'][label]], r['crs'], facecolor=color,
                                      edgecolor='none', zorder=2)
            if not r['gaps'].is_empty:
                ax.add_geometries([r['gaps']], r['crs'], facecolor='#d62728', edgecolor='#d62728',
                                  linewidth=0.3, zorder=3)
        ax.coastlines('50m', linewidth=0.3, color='#333')
        if name == 'Arctic':
            ax.add_geometries([box(-180, args.latMax, 180, 90)], ccrs.PlateCarree(), facecolor='none',
                              edgecolor='#888', hatch='//', linewidth=0, zorder=4)
        ax.gridlines(linewidth=0.3, alpha=0.4)
        ax.set_title(name, fontsize=11)
    fig.legend(handles=[Patch(color=c, label=l) for l, _, c in TIERS] +
               [Patch(color='#d62728', label='land not covered'),
                Patch(facecolor='none', edgecolor='#888', hatch='//', label=f'north of {args.latMax}N (not tiled)')],
               loc='upper center', ncol=5, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.text(0.5, 0.96, (args.title or os.path.basename(os.path.abspath(args.runDir))) + '   |   ' + stats[0],
             ha='center', fontsize=10)
    fig.savefig(f'{args.runDir}/coverage.png')
    print('\n'.join(stats))


if __name__ == '__main__':
    main()
