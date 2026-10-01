#!/usr/bin/env python3
'''
Make the coarse levels of a tiered Google Earth layer link-only: their tile KMLs keep their
Regions and NetworkLinks (the index Google Earth walks down to the fine tiles) but lose their
images, and those PNGs are deleted. The base layer shows the same imagery at those zooms.

  python3 -m globalMosaic.linkOnly /scratch/ianj/mosaics/cycle30ascending_z11/googleEarth/land \\
      /scratch/ianj/mosaics/cycle30ascending_z11/googleEarth/detail [--processes 16]

Every zoom level below a tree's finest is stripped (land: 5-9 of 5-10, detail: 7-10 of 7-11), which
removes about a quarter of each layer's PNGs. Opened without base, such a layer is blank until its
finest zoom. Run after oceanAltitude (split tiles' .land.png/.sea.png go too) and before tarring;
rerunning is harmless. Not reversible except by rebuilding the layer -- never run it on base.
'''
import argparse
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

OVERLAY = re.compile(r'    <GroundOverlay>.*?</GroundOverlay>\n', re.S)


def stripCoarse(tree, processes=16):
    zooms = sorted(int(d) for d in os.listdir(tree) if d.isdigit())
    if len(zooms) < 2:
        return
    coarse = [f'{tree}/{z}' for z in zooms[:-1]]
    with ProcessPoolExecutor(processes) as pool:
        # big levels split by column directory so the work spreads over the processes
        jobs = [c for z in coarse for c in sorted(glob.glob(f'{z}/*'))]
        res = list(pool.map(stripCol, jobs, chunksize=16))
    print(f'{tree}: zooms {zooms[0]}-{zooms[-2]} link-only ({zooms[-1]} keeps its images): '
          f'{sum(r[0] for r in res)} KMLs stripped, {sum(r[1] for r in res)} PNGs removed', flush=True)


def stripCol(col):
    ''' One column directory <z>/<x>/ (same as stripDir, one level down). '''
    kmls = pngs = 0
    for kml in glob.glob(f'{col}/*.kml'):
        s = open(kml).read()
        t = OVERLAY.sub('', s)
        if t != s:
            tmp = kml + '.tmp'
            with open(tmp, 'w') as fp:
                fp.write(t)
            os.replace(tmp, kml)
            kmls += 1
    for png in glob.glob(f'{col}/*.png'):
        os.remove(png)
        pngs += 1
    return kmls, pngs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tree', nargs='+', help='layer directories (land/, detail/) -- not base/')
    ap.add_argument('--processes', type=int, default=16)
    args = ap.parse_args()
    for t in args.tree:
        if os.path.basename(os.path.normpath(t)) == 'base':
            sys.exit('refusing to strip base/: it is the layer that shows the coarse zooms')
        stripCoarse(t, args.processes)
    return 0


if __name__ == '__main__':
    sys.exit(main())
