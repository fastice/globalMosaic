#!/bin/bash
# Land without data in a run's finished tiles, worst tiles first (safe while a run goes):
#   bash test/gaps.sh                     # the default full run (test/fullAsc.yaml)
#   bash test/gaps.sh test/smokeAsc.yaml
# Full list in <work>/quicklooks/gaps.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
CONFIG=${1:-test/fullAsc.yaml}
WORK=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['work'])")
python3 -m globalMosaic.gapReport $WORK --top ${TOP:-40}
