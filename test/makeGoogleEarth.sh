#!/bin/bash
# Google Earth KML superoverlay from a finished run; open $WORK/googleEarth/doc.kml.
set -euo pipefail
cd "$(dirname "$0")/.."
WORK=${WORK:-run_c030}
time makeGoogleEarth $WORK --out $WORK/googleEarth --processes ${NPROC:-16} --zoom ${ZOOM:-0-9}
