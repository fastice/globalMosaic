#!/bin/bash
# Timing test of the temporal statistics on AWS: two frames (087_024_2005 = 20 MHz, Iberia;
# 055_020_4005 = 40 MHz, Pakistan), all ascending cycles 023-031, then their 4 tiles through
# geomosaic -geo. No Google Earth. Foreground, a few minutes to an hour.
# On petermann (over the internet) the frames took 7 and 60 min.
#   bash test/statsTest.sh                  # -> /scratch/ianj/mosaics/statsTest
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
WORK=${WORK:-/scratch/ianj/mosaics/statsTest}
mkdir -p $WORK
geomosaic -geo 2>&1 | grep -q 'unknown\|Unknown' && { echo "geomosaic has no -geo: run installGeomosaic"; exit 1; } || true
CATS="cycle30/catalogue023.geojson cycle30/catalogue024.geojson cycle30/catalogue025.geojson cycle30/catalogue026.geojson cycle30/catalogue027.geojson cycle30/catalogue028.geojson cycle30/catalogue029.geojson cycle30/catalogue.geojson cycle30/catalogue031.geojson"
{
    echo "=== frames $(date)"
    time python3 -m globalMosaic.frameStatistics cycle30/ascending --catalogues $CATS --out $WORK/frames \
        --frames 087_024_2005 055_020_4005 --nProc 2
    echo "=== mosaic $(date)"
    time python3 -m globalMosaic.mosaicStatistics cycle30/ascending --frames $WORK/frames --work $WORK \
        --tiles N42W006_6x6 N42W012_6x6 N30E072_6x6 N36E072_6x6 --nProc 4
    echo "=== done $(date): per-frame times in $WORK/frames/frames.log; tiles in $WORK/tiles/{mean,sigma,cv,n}"
} 2>&1 | tee $WORK/run.log
