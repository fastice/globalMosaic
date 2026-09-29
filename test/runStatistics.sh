#!/bin/bash
# Temporal backscatter statistics (mean, sigma, CV, n) of every ascending
# cycle from 2026-06-17 (cycles 023-031), per (track, frame) of the ascending tiling, then the
# global mosaic with geomosaic -geo (80 m caps, 2.4" lat/lon), then Google Earth zoom 0-9 of cv. In the background; log: <work>/run.log. Resumable: finished frames and tiles are skipped.
# The per-frame statistics (<work>/frames) are kept, for rerunning the mosaic.
#   bash test/runStatistics.sh                      # ascending  -> /scratch/ianj/mosaics/stats30asc
#   DIR=descending bash test/runStatistics.sh       # descending -> /scratch/ianj/mosaics/stats30desc
# Needs a geomosaic with -geo (installGeomosaic once the helheim change is pushed).
set -euo pipefail
cd "$(dirname "$0")/.."
DIR=${DIR:-ascending}
case $DIR in ascending) TAG=asc ;; descending) TAG=desc ;; *) echo "DIR must be ascending or descending"; exit 1 ;; esac
WORK=${WORK:-/scratch/ianj/mosaics/stats30$TAG}
mkdir -p $WORK
CATS="cycle30/catalogue023.geojson cycle30/catalogue024.geojson cycle30/catalogue025.geojson cycle30/catalogue026.geojson cycle30/catalogue027.geojson cycle30/catalogue028.geojson cycle30/catalogue029.geojson cycle30/catalogue.geojson cycle30/catalogue031.geojson"
nohup bash -c "
    set -e
    python3 -m globalMosaic.frameStatistics cycle30/$DIR --catalogues $CATS --out $WORK/frames --nProc ${NPROC:-16}
    echo \"=== frame statistics done \$(date)\"
    python3 -m globalMosaic.mosaicStatistics cycle30/$DIR --frames $WORK/frames --work $WORK --nProc ${NPROC:-16} \
        --googleEarth 'cv' --zoom ${ZOOM:-0-9} --title "NISAR temporal statistics, $DIR, 2026-06-17 to cycle 31"
    tar -cf $WORK/googleEarth.tar -C $WORK --exclude='googleEarth/*/stage' googleEarth
    echo \"=== all done \$(date): $WORK/googleEarth/cv/doc.kml, archive $WORK/googleEarth.tar\"
" > $WORK/run.log 2>&1 &
echo "started pid $!; log: $WORK/run.log"
