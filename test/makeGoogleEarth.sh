#!/bin/bash
# Google Earth tiles alone (runFull.sh already does this after the mosaic), e.g. to redo them
# with another stretch or zoom:  ZOOM=0-10 test/makeGoogleEarth.sh test/fullAsc.yaml --dbMin -28
set -euo pipefail
cd "$(dirname "$0")/.."
CONFIG=test/fullAsc.yaml
if [[ "${1:-}" == *.yaml ]]; then CONFIG=$1; shift; fi
WORK=${WORK:-$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")}
time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-24} --zoom ${ZOOM:-0-9} "$@"
