#!/bin/bash
# Full cycle-30 run (test/full.yaml) in the background; log in run_c030.log.
# Extra options override the yaml, e.g.  test/runFull.sh --nProc 32
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
nohup runGeomosaicTiles --config test/full.yaml "$@" > run_c030.log 2>&1 &
echo "started pid $!; log: $(pwd)/run_c030.log"
