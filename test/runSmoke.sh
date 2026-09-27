#!/bin/bash
# Smoke test in the foreground: mosaic (test/smoke.yaml, about 10 min on petermann, nearly all
# download from ASF), then Google Earth tiles for it (a minute or two).
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('test/smoke.yaml'))['work'])")
time runGeomosaicTiles --config test/smoke.yaml "$@"
grep -h 'read+average' $WORK/jobs/*/log
ls $WORK/vrt/gamma0/
time makeGoogleEarth $WORK --out $WORK/googleEarth --processes 4 --zoom 0-9
echo "open $WORK/googleEarth/doc.kml"
