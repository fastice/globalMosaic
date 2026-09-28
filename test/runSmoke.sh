#!/bin/bash
# Smoke test in the foreground: mosaic (test/smoke.yaml, about 10 min on petermann, nearly all
# download from ASF), then quick-look PNGs (<work>/quicklooks/) and Google Earth tiles. Everything goes in the yaml's work directory:
# <work>/run.log, <work>/vrt/, <work>/googleEarth/doc.kml.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
CONFIG=${CONFIG:-test/smokeAsc.yaml}
# the run yaml can also be given first:  bash test/runSmoke.sh test/smokeAsc.yaml
if [[ "${1:-}" == *.yaml ]]; then CONFIG=$1; shift; fi
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")
mkdir -p $WORK
{
    time runGeomosaicTiles --config $CONFIG "$@"
    grep -h 'read+average' $WORK/jobs/*/log
    makeQuickLook $WORK --processes ${NPROC:-24}
    time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-24} --zoom ${ZOOM:-0-9}
    tar -cf $WORK/googleEarth.tar -C $WORK --exclude=googleEarth/stage googleEarth
    echo "=== all done $(date): open $WORK/googleEarth/doc.kml; archive $WORK/googleEarth.tar"
} 2>&1 | tee $WORK/run.log
