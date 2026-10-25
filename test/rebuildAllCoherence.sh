#!/bin/bash
# Delete the existing per-cycle coherence products and rebuild cycles 24..32 from scratch.
#
#   bash test/rebuildAllCoherence.sh            # cycles 024..032
#   CYCLES="024 025" bash test/rebuildAllCoherence.sh
#   SCAN=1 bash test/rebuildAllCoherence.sh     # also scan valid-data footprints (slow)
#
# From scratch, not resumed: both the antimeridian fix (a footprint crossing 180 used to
# project across the whole map and falsely cover every tile at its latitude) and --allGranules
# (every pair over a tile, so each cell is the mean coherence rather than one picked pair)
# change which granules a tile lists, and an existing tiles.geojson would be reused as-is.
#
# The catalogues in cycle*gunw/ are NOT deleted: cycle 24's is committed, the rest are cached
# from their ASF search and cost ~10 minutes each to rebuild.
set -u
ROOT=${ROOT:-/scratch/ianj/mosaics/coherence}
CYCLES=${CYCLES:-"024 025 026 027 028 029 030 031 032"}
export SCAN=${SCAN:-0}
export NPROC=${NPROC:-96}
export DIRECTION=${DIRECTION:-A}
export ROOT
# Descending must not land on top of ascending: they share the catalogues but nothing else.
if [ "$DIRECTION" != A ] && [ "$ROOT" = /scratch/ianj/mosaics/coherence ]; then
    echo "DIRECTION=$DIRECTION with the default ROOT would overwrite the ascending product." >&2
    echo "Give it its own, e.g. ROOT=/scratch/ianj/mosaics/coherence$DIRECTION" >&2
    exit 2
fi

cd "$(dirname "$0")/.." || exit 1

P='[g]unwCoherence|[r]unCycleCoherence|[g]eomosaic|[e]aseTiles|[r]unGeomosaicTiles|[s]canMasks'
if pgrep -f "$P" >/dev/null; then
    echo "=== $(date) stopping what is running ==="
    pgrep -af "$P" | cut -c1-100
    pkill -f "$P"; sleep 5; pkill -9 -f "$P"; sleep 2
fi
if pgrep -f "$P" >/dev/null; then
    echo "SOMETHING IS STILL RUNNING -- stopping here rather than deleting under it" >&2
    pgrep -af "$P" | cut -c1-100 >&2
    exit 1
fi

echo "=== $(date) removing the previous coherence products in $ROOT ==="
rm -rf "$ROOT"/cycle0* "$ROOT"/coherence_cycle0*.png
mkdir -p "$ROOT"

for c in $CYCLES; do
    echo "######## $(date) cycle $c ########"
    bash test/runCycleCoherence.sh "$c" || echo "  cycle $c FAILED; continuing"
done
echo "######## $(date) all cycles done ########"
ls -l "$ROOT"/coherence_cycle0*.png 2>/dev/null
