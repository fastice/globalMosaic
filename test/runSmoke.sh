#!/bin/bash
# Smoke test (test/smoke.yaml): about 10 min on petermann, nearly all download from ASF.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
time runGeomosaicTiles --config test/smoke.yaml "$@"
grep -h 'read+average' runTest/jobs/*/log
ls runTest/vrt/gamma0/
