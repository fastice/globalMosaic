#!/bin/bash
# Google Earth tiles alone (runFull.sh already does this after the mosaic), e.g. to redo them
# with another stretch or zoom:  WORK=run_c030 ZOOM=0-10 test/makeGoogleEarth.sh --dbMin -28
set -euo pipefail
cd "$(dirname "$0")/.."
WORK=${WORK:-run_c030}
time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-24} --zoom ${ZOOM:-0-9} "$@"
