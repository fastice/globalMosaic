#!/usr/bin/env python3
'''
Tile names of a tiling that cover an area, to hand to --tiles.

  tilesFor run/tiles --bbox -75 -60 -30 15            # South America
  tilesFor run/tiles --point 25 -19.5                 # the tile holding one place
  tilesFor run/tiles --bbox 10 35 45 72 --names       # just the names, space separated

A tiling covers the globe, so a regional mosaic is the global tiling with only some tiles
built. That is deliberate: the grid never depends on the data or on the area asked for, so
two regions, two cycles and two products all land on the same pixels and can be compared
directly. Building a different grid per request would throw that away.

  easeTiles --catalogue cat.geojson --hemisphere global --crs ESRI:54076 --out run/tiles
  runGeomosaicTiles run/tiles --work run --tiles $(tilesFor run/tiles --bbox -75 -60 -30 15 --names)
'''
import argparse
import json
import sys

from shapely.geometry import box, shape


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tileRun', help='a tiling directory (its tiles.geojson)')
    ap.add_argument('--bbox', type=float, nargs=4, metavar=('W', 'S', 'E', 'N'),
                    help='area in degrees')
    ap.add_argument('--point', type=float, nargs=2, metavar=('LON', 'LAT'),
                    help='a single place')
    ap.add_argument('--names', action='store_true',
                    help='print only the names, space separated, for --tiles')
    ap.add_argument('--withData', action='store_true',
                    help='only tiles the tiler found granules for (the default lists those '
                         'too; this drops any with selected == 0)')
    args = ap.parse_args()
    if not args.bbox and not args.point:
        sys.exit('tilesFor: give --bbox or --point')

    area = (box(*args.bbox) if args.bbox
            else shape(dict(type='Point', coordinates=args.point)).buffer(1e-6))
    feats = json.load(open(f'{args.tileRun}/tiles.geojson'))['features']
    hit = []
    for f in feats:
        p = f['properties']
        if args.withData and not p.get('selected'):
            continue
        try:
            g = shape(f['geometry']).buffer(0)
        except Exception:
            continue
        # a tiling of a projected CRS stores its outline in lat/lon; one that spans the
        # dateline is tested against both images of the area
        if g.intersects(area) or g.intersects(area.buffer(0)):
            hit.append((p['name'], p.get('selected', 0), p.get('coverage', 0.)))
    if not hit:
        sys.exit('tilesFor: no tile covers that area')
    if args.names:
        print(' '.join(n for n, _, _ in hit))
        return
    print(f'{len(hit)} tile(s):')
    print('%-10s %9s %9s' % ('tile', 'granules', 'coverage'))
    for n, s, c in sorted(hit):
        print('%-10s %9d %8.1f%%' % (n, s, c * 100))


if __name__ == '__main__':
    main()
