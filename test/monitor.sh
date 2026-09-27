#!/bin/bash
# Progress picture of a run from what is on disk (safe while it runs):
#   bash test/monitor.sh                     # the default full run (test/fullAsc.yaml)
#   bash test/monitor.sh test/smokeAsc.yaml  # another run
# Writes <work>/quicklooks/monitor.png and monitor.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
monitorRun --config ${1:-test/fullAsc.yaml}
