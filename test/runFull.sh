#!/bin/bash
# Full cycle-30 run in the background: mosaic (test/full.yaml), then the Google Earth tiles if
# the mosaic finished cleanly. Everything goes in the yaml's work directory: <work>/run.log,
# <work>/vrt/, <work>/googleEarth/doc.kml. Extra options go to the mosaicker and override the
# yaml, e.g.  test/runFull.sh --nProc 32 ;  another run file:  CONFIG=my.yaml test/runFull.sh
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
CONFIG=${CONFIG:-test/full.yaml}
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")
# RUNROOT (e.g. a scratch drive) holds the run directory instead of the repo:
#   RUNROOT=/scratch/me bash test/runFull.sh   ->  /scratch/me/cycle30both
if [ -n "${RUNROOT:-}" ]; then WORK=$RUNROOT/$(basename $WORK); fi
mkdir -p $WORK
nohup bash -c '
    set -e
    runGeomosaicTiles --config "$0" --work "$1" "${@:2}"
    echo "=== mosaic done $(date); Google Earth tiles"
    makeGoogleEarth "$1" --out "$1/googleEarth" --processes ${NPROC:-24} --zoom ${ZOOM:-0-9}
    echo "=== all done $(date): open $1/googleEarth/doc.kml"
' "$CONFIG" "$WORK" "$@" > $WORK/run.log 2>&1 &
echo "started pid $!; log: $(realpath $WORK/run.log)"
