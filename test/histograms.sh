#!/bin/bash
# Area-weighted backscatter histograms of a run: global, land, ocean, continents, Greenland, sea-ice
# stand-ins. Safe while a run goes (uses finished tiles). Output in <work>/stats/.
#   bash test/histograms.sh                    # the default full run (test/fullAsc.yaml)
#   bash test/histograms.sh test/fullAsc11.yaml
set -euo pipefail
cd "$(dirname "$0")/.."
CONFIG=${1:-test/fullAsc.yaml}
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")
python3 -m globalMosaic.backscatterHistograms $WORK --stride ${STRIDE:-8} --processes ${NPROC:-16}
