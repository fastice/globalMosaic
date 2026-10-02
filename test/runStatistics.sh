#!/bin/bash
# Temporal backscatter statistics (mean, sigma, CV, n), co-pol and cross-pol, of every ascending
# cycle from 2026-06-17 (cycles 023-031), per (track, frame) of the ascending tiling, then the
# global mosaic with geomosaic -geo (80 m caps, 2.4" lat/lon), then Google Earth zoom 0-9 of cv. In the background; log: <work>/run.log. Resumable: finished frames and tiles are skipped.
# The per-frame statistics (<work>/frames: <frame>.tif co-pol, <frame>.cross.tif HV/VH) are kept,
# for rerunning the mosaic. Mosaics: co-pol in <work>, cross-pol in <work>_cross.
#   bash test/runStatistics.sh                      # ascending  -> /scratch/ianj/mosaics/stats30asc
#   DIR=descending bash test/runStatistics.sh       # descending -> /scratch/ianj/mosaics/stats30desc
#   bash test/runStatistics.sh -kill               # stop a running one (restart with the plain command)
# Needs a geomosaic with -geo (installGeomosaic once the helheim change is pushed).
set -euo pipefail
if [[ "${1:-}" == -kill || "${1:-}" == --kill ]]; then
    # the background wrapper, the python drivers and their workers, and geomosaic -geo jobs
    # (not gdal2tiles: another run's Google Earth build may be using it)
    pat='globalMosaic\.(frameStatistics|mosaicStatistics)|geomosaic .*-geo '
    pids=$(pgrep -f "$pat" | grep -vx -e $$ -e $PPID || true)
    [ -n "$pids" ] && kill $pids
    sleep 3
    left=$( (pgrep -f "$pat" || true) | (grep -vx -e $$ -e $PPID || true) | wc -l)
    n=$(echo $pids | wc -w)
    echo "stopped $n processes; $left still running. Finished frames and tiles are kept: rerun to resume."
    exit 0
fi
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
    for L in co cross; do
        W=$WORK; T='co-pol (HH, or VV)'
        [ \$L = cross ] && { W=${WORK}_cross; T='cross-pol (HV, or VH)'; }
        python3 -m globalMosaic.mosaicStatistics cycle30/$DIR --frames $WORK/frames --work \$W --layer \$L \
            --nProc ${NPROC:-16} --googleEarth 'cv' --zoom ${ZOOM:-0-9} --tileSize ${TILESIZE:-512} \
            --title \"NISAR temporal statistics, $DIR, \$T, 2026-06-17 to cycle 31\"
        tar -cf \$W/googleEarth.tar -C \$W --exclude='googleEarth/*/stage' googleEarth
        echo \"=== \$L done \$(date): \$W/googleEarth/cv/doc.kml\"
    done
    echo \"=== all done \$(date): $WORK/googleEarth/cv/doc.kml (co-pol), ${WORK}_cross/googleEarth/cv/doc.kml (cross-pol)\"
" > $WORK/run.log 2>&1 &
echo "started pid $!; log: $WORK/run.log"
