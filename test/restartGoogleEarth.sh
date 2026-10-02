#!/bin/bash
# Restart a run's Google Earth step without waiting for runFull's `rm -rf googleEarth` of the old
# product (millions of files): stops runFull's wrapper and that rm, moves the old googleEarth aside
# (instant), deletes it in the background at low priority, and builds the new product, then the
# tars (tarGoogleEarth.sh; one per layer for a tiered run). Log: <work>/ge.log.
#   bash test/restartGoogleEarth.sh                       # test/fullAsc11.yaml
#   bash test/restartGoogleEarth.sh test/fullAsc.yaml
#   NPROC=16 bash test/restartGoogleEarth.sh
# Only after the mosaic itself is done (run.log: '=== quick looks ... Google Earth tiles').
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
CONFIG=${1:-test/fullAsc11.yaml}
eval "$(python3 - "$CONFIG" <<'PY'
import os, shlex, sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
w = c['work'].rstrip('/')
print(f"W={shlex.quote(w)}")
print(f"ZOOM={shlex.quote(str(c.get('zoom') or '0-9'))}")
print(f"TITLE={shlex.quote(c.get('title') or 'NISAR gamma0 ' + os.path.basename(w))}")
print(f"TILESIZE={int(c.get('tileSize') or 256)}")
PY
)"
grep -q 'Google Earth tiles' "$W/run.log" || { echo "$W/run.log: the mosaic is not done yet -- not restarting"; exit 1; }
mine() { pgrep -f "$1" | grep -vx -e $$ -e $PPID || true; }
# 1. runFull's wrapper for this run (it would go on to tile and tar), then its rm
pids=$(mine "runGeomosaicTiles --config.* $W( |$)")
[ -n "$pids" ] && { echo "stopping runFull wrapper: $pids"; kill $pids; }
pids=$(mine "rm -rf $W/googleEarth $W/googleEarth.tar")
[ -n "$pids" ] && { echo "stopping rm: $pids"; kill $pids; }
sleep 2
# 2. old product aside, deleted in the background
if [ -e "$W/googleEarth" ]; then
    OLD=$W/googleEarth.old.$(date +%s)
    mv "$W/googleEarth" "$OLD"
    nohup nice -n 19 ionice -c 3 rm -rf "$OLD" > /dev/null 2>&1 &
    echo "old product moved to $OLD; deleting it in the background (pid $!)"
fi
rm -f "$W/googleEarth.tar"
# 3. new product, then the tars
export W ZOOM TITLE TILESIZE NPROC=${NPROC:-16}
nohup bash -c '
    set -e
    if [ "$ZOOM" = tiered ]; then
        makeGoogleEarthTiered "$W" --out "$W/googleEarth" --processes $NPROC --title "$TITLE" --tileSize $TILESIZE
        bash test/tarGoogleEarth.sh "$W" split
    else
        makeGoogleEarth "$W" --out "$W/googleEarth" --processes $NPROC --zoom $ZOOM --title "$TITLE" --tileSize $TILESIZE
        bash test/tarGoogleEarth.sh "$W"
    fi
    echo "=== all done $(date): $W/googleEarth/doc.kml"
' > "$W/ge.log" 2>&1 &
echo "Google Earth build started (pid $!): tail -f $W/ge.log"
