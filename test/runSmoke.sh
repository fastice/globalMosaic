#!/bin/bash
# Smoke test in the foreground: mosaic (test/smoke.yaml, about 10 min on petermann, nearly all
# download from ASF), then its Google Earth tiles. Everything goes in the yaml's work directory:
# <work>/run.log, <work>/vrt/, <work>/googleEarth/doc.kml.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
CONFIG=${CONFIG:-test/smoke.yaml}
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")
mkdir -p $WORK
{
    time runGeomosaicTiles --config $CONFIG "$@"
    grep -h 'read+average' $WORK/jobs/*/log
    time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-4} --zoom ${ZOOM:-0-9}
    echo "=== all done $(date): open $WORK/googleEarth/doc.kml"
} 2>&1 | tee $WORK/run.log
