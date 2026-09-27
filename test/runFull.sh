#!/bin/bash
# Full cycle-30 run: 930 jobs (all tiles plus the south cap). Resumable: rerunning skips
# finished jobs. Jobs are single-threaded and mostly waiting on downloads, so NPROC can exceed
# the core count; 16 cores / 256 GB here. Log in $WORK.log.
#   NPROC=24 test/runFull.sh              # TILING=cycle30/ascending for one direction
#   NOCAP=1 test/runFull.sh               # leave out the ~13 h south-cap job
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
TILING=${TILING:-cycle30/both}
WORK=${WORK:-run_c030}
TILES=()
if [ -n "${NOCAP:-}" ]; then
    TILES=(--tiles $(python3 -c "import json; print(' '.join(f['properties']['name'] for f in json.load(open('$TILING/tiles.geojson'))['features'] if 'PS' not in f['properties']['name']))"))
fi
nohup runGeomosaicTiles $TILING --work $WORK --granules vsicurl --nProc ${NPROC:-16} "${TILES[@]}" \
    > $WORK.log 2>&1 &
echo "started pid $!; log: $(pwd)/$WORK.log"
