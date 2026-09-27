#!/usr/bin/env python3
'''
Progress picture of a runGeomosaicTiles run, drawn from what is on disk -- run it any time, from
any shell, without touching the run (e.g. one started before the driver drew its own).

  monitorRun /scratch/ianj/mosaics/cycle30ascending      # or: monitorRun --config test/fullAsc.yaml

Reads <work>/run.yaml (the run's settings), the tiling it names, the finished tiles in
<work>/tiles/<product>/ and the job lines in <work>/summary.log, and writes
<work>/quicklooks/monitor.png (60S-60N in lat/lon, polar stereographic poleward of 60N/60S;
finished grey, planned blue, failed red) and monitor.txt. Low-res pieces of finished tiles are
kept in quicklooks/progress/ (shared with the driver's own progress picture), so reruns are quick.
'''
import argparse
import datetime
import glob
import json
import os
import re
import sys
import time

import yaml
from osgeo import gdal

from .runGeomosaicTiles import Progress, tileGrid

gdal.UseExceptions()
JOBLINE = re.compile(r'^(\d\d-\d\d \d\d:\d\d:\d\d) (?:\[\d+/\d+, \d+ failed\] )?(\S+?)\.([AB][HV][HV]): '
                     r'[\d.]+ min (OK|FAILED)(.*)')
START = re.compile(r'^(\d\d-\d\d \d\d:\d\d:\d\d) \d+ tiles \(')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('work', nargs='?', default=None, help='run (work) directory')
    ap.add_argument('--config', default=None, help='run yaml instead of a directory (uses its work:)')
    ap.add_argument('--stem', default='monitor', help='output name in <work>/quicklooks [monitor]')
    args = ap.parse_args()
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    work = args.work
    if work is None:
        if args.config is None:
            ap.error('give a run directory or --config')
        cfgPath = args.config if os.path.exists(args.config) else os.path.join(repo, args.config)
        work = yaml.safe_load(open(cfgPath))['work']
        work = work if os.path.isabs(work) else os.path.join(repo, work)
    cfg = yaml.safe_load(open(f'{work}/run.yaml'))
    tileRun = cfg['tileRun'] if os.path.isabs(cfg['tileRun']) else os.path.join(repo, cfg['tileRun'])
    product = 'gamma0' if cfg.get('calOutput', 'gamma0') != 'sigma0' else 'sigma0'

    # the run's jobs, as the driver enumerates them
    feats = json.load(open(f'{tileRun}/tiles.geojson'))['features']
    if cfg.get('tiles'):
        feats = [f for f in feats if f['properties']['name'] in cfg['tiles']]
    if cfg.get('excludeTiles'):
        feats = [f for f in feats if f['properties']['name'] not in cfg['excludeTiles']]
    jobs, grids = [], {}
    for f in feats:
        name = f['properties']['name']
        grid = tileGrid(f['properties'], f['geometry'], cfg['res'], cfg['psResM'], cfg['featherKm'])
        for y in sorted(glob.glob(f'{tileRun}/tiles/{name}.*.yaml')):
            jobDir = f'{work}/jobs/{name}.{y.split(".")[-2]}'
            jobs.append((name, y.split('.')[-2], jobDir))
            grids[jobDir] = grid

    # outcomes from the log since the last (re)start; a finished tile on disk always counts as done
    t0, last, ran = None, {}, 0
    for line in open(f'{work}/summary.log', errors='ignore'):
        m = START.match(line)
        if m:
            t0, last, ran = m.group(1), {}, 0
            continue
        m = JOBLINE.match(line)
        if m:
            last[(m.group(2), m.group(3))] = m.group(4)
            ran += 'skipped' not in m.group(5)          # jobs actually run since the (re)start
    p = Progress(work, feats, jobs, grids, product)
    if t0:
        now = datetime.datetime.now()
        start = datetime.datetime.strptime(f'{now.year}-{t0}', '%Y-%m-%d %H:%M:%S')
        p.t0 = time.time() - (now - start).total_seconds()
    for name, grp, _ in jobs:
        tif = f'{work}/tiles/{product}/{name}.{grp}.tif'
        if os.path.exists(tif):
            piece = f'{p.dir}/progress/{name}.{grp}.tif'
            if os.path.exists(piece) and os.path.getmtime(piece) >= os.path.getmtime(tif):
                p.pieces.insert(0, piece) if grp.startswith('B') else p.pieces.append(piece)
                if grids[f'{work}/jobs/{name}.{grp}']['epsg'] != 4326:
                    p.pieces.remove(piece)
                    p.caps.setdefault(name, []).insert(0 if grp.startswith('B') else 99, piece)
            else:
                p.addJob(name, grp, tif)
            p.count(True, name)
        elif last.get((name, grp)) == 'FAILED':
            p.count(False, name)
    p.etaDone = ran        # ETA from the rate of jobs actually run since the (re)start
    p.draw(args.stem)
    print(open(f'{p.dir}/{args.stem}.txt').read().strip())
    print(f'{p.dir}/{args.stem}.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
