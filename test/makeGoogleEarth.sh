#!/bin/bash
# Google Earth tiles alone (runFull.sh already does this after the mosaic), e.g. to redo them
# with another stretch or zoom:  WORK=runs/cycle30both ZOOM=0-10 test/makeGoogleEarth.sh --dbMin -28
set -euo pipefail
cd "$(dirname "$0")/.."
WORK=${WORK:-${RUNROOT:+$RUNROOT/cycle30both}}
WORK=${WORK:-runs/cycle30both}
time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-24} --zoom ${ZOOM:-0-9} "$@"
