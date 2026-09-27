#!/bin/bash
# Smoke test: two 6x6 deg tiles tested on petermann (N42W012 freq A, N24E006 freq A + B).
# About 10 min on petermann; nearly all of it is download from ASF.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
TILING=${TILING:-cycle30/both}          # or cycle30/ascending, cycle30/descending
WORK=${WORK:-runTest}
time runGeomosaicTiles $TILING --work $WORK --granules vsicurl --nProc ${NPROC:-4} \
    --tiles N42W012_6x6 N24E006_6x6
grep -h 'read+average' $WORK/jobs/*/log
ls $WORK/vrt/gamma0/
