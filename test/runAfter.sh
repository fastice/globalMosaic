#!/bin/bash
# Wait (in the background) until the current run -- mosaic, quick looks, Google Earth -- has
# finished, then start test/runFull.sh with the same arguments. Survives logging out.
#   bash test/runAfter.sh                    # default run (test/fullAsc.yaml)
#   bash test/runAfter.sh test/fullDesc.yaml
# Log: ~/runAfter.log
cd "$(dirname "$0")/.."
nohup bash -c '
    while pgrep -f "[r]unGeomosaicTiles|[m]akeQuickLook|[m]akeGoogleEarth|[g]dal2tiles" > /dev/null; do
        sleep 300
    done
    echo "=== previous run finished $(date); starting test/runFull.sh $*"
    bash test/runFull.sh "$@"
' runAfter "$@" > ~/runAfter.log 2>&1 &
echo "waiting in the background (pid $!); log: ~/runAfter.log"
