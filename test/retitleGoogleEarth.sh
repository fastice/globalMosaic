#!/bin/bash
# Rename existing Google Earth products (the name Google Earth lists) without rebuilding them.
#   bash test/retitleGoogleEarth.sh
set -euo pipefail
cd "$(dirname "$0")/.."
M=/scratch/ianj/mosaics
r() { [ -f "$1/doc.kml" ] && python3 -m globalMosaic.makeGoogleEarth --out "$1" --retitle --title "$2" || true; }
r $M/cycle30ascending/googleEarth      "NISAR gamma0, cycle 30 ascending, 3 arcsec"
r $M/cycle30ascending_z11/googleEarth  "NISAR gamma0, cycle 30 ascending, 2.4 arcsec (zoom to 11 where it matters)"
for L in base land detail; do r $M/cycle30ascending_z11/googleEarth/$L "NISAR gamma0 2.4 arcsec - $L layer"; done
r $M/coherence30asc/googleEarth         "NISAR 12-day coherence, cycle 30 ascending"
r $M/smokeAsc/googleEarth               "NISAR gamma0 smoke test (ascending)"
