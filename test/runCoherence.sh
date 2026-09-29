#!/bin/bash
# Global 12-day coherence mosaic from GUNW (cycle-030 pairs, 029/028/027 fill), then quick looks and
# Google Earth zoom 0-9, in the background. Log: <work>/run.log. Resumable (finished tiles skipped).
#   bash test/runCoherence.sh                 # ascending -> /scratch/ianj/mosaics/coherence30asc
#   DIR=descending bash test/runCoherence.sh  # descending -> /scratch/ianj/mosaics/coherence30desc
set -euo pipefail
cd "$(dirname "$0")/.."
DIR=${DIR:-ascending}
WORK=${WORK:-/scratch/ianj/mosaics/coherence30${DIR:0:3}}
mkdir -p $WORK
nohup python3 -m globalMosaic.mosaicCoherence cycle30gunw/$DIR --work $WORK --nProc ${NPROC:-16} \
    > $WORK/run.log 2>&1 &
echo "started pid $!; log: $WORK/run.log"
