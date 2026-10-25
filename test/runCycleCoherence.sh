#!/bin/bash
# Coherence map for ONE cycle, on the same Wagner VII grid as the backscatter product.
#
#   bash test/runCycleCoherence.sh 024
#   for c in 024 025 026 027 028 029 030 031 032; do
#       bash test/runCycleCoherence.sh $c
#   done
#   RES=2000 NPROC=16 bash test/runCycleCoherence.sh 025
#   CRS=EPSG:6933 bash test/runCycleCoherence.sh 026        # no rebuilt geomosaic needed
#
# One cycle per run, deliberately: a cycle's GUNWs are its own 12-day pairs, so each map is a
# snapshot of that cycle and the series shows how coherence changes. No --fillCatalogues --
# borrowing another cycle's pairs would blur exactly the thing being looked at. A track missed
# in a cycle is therefore simply absent from that cycle's map.
#
# The tile grid depends only on CRS and TILEKM, not on the data, so every cycle lands on the
# same pixels and the maps can be differenced directly.
#
# Resumable: tiling is skipped when tiles.geojson exists, extracted coherence is cached per
# granule, and finished tiles are skipped.
set -u
CYCLE=${1:-}
if [ -z "$CYCLE" ]; then echo "usage: $0 <cycle, e.g. 024>" >&2; exit 2; fi

cd "$(dirname "$0")/.." || exit 1
ROOT=${ROOT:-/scratch/ianj/mosaics/coherence}
CRS=${CRS:-ESRI:54076}
TILEKM=${TILEKM:-600}
RES=${RES:-1000}                 # output spacing, m
NPROC=${NPROC:-32}
DIRECTION=${DIRECTION:-A}
# The catalogue directory is named either way round in practice -- cycle30gunw beside
# cycle024gunw -- so accept both rather than silently re-searching ASF for one that is there.
SHORT=$((10#$CYCLE))
SCAN=${SCAN:-1}                  # 0 to skip the valid-footprint scan (and keep the banding)
SCANTHREADS=${SCANTHREADS:-24}
CAT=${CAT:-}
if [ -z "$CAT" ]; then
    for d in "cycle${CYCLE}gunw" "cycle${SHORT}gunw"; do
        if [ -f "$d/catalogue.geojson" ]; then CAT=$d/catalogue.geojson; break; fi
    done
    CAT=${CAT:-cycle${CYCLE}gunw/catalogue.geojson}
fi
W=$ROOT/cycle$CYCLE
T=$W/tiles-global

mkdir -p "$W"

if [ ! -f "$CAT" ]; then
    echo "=== $(date) searching ASF for cycle $CYCLE GUNW ==="
    globalGCOVTiles --cycle "$CYCLE" --product GUNW --direction "$DIRECTION" \
        --bbox -180 -90 180 90 --out "$(dirname "$CAT")" || exit 1
fi
echo "=== $(date) catalogue: $CAT ($(python3 -c "
import json;print(len(json.load(open('$CAT'))['features']))") granules) ==="

# Valid-data footprints. Without them the tiler selects on the ACQUISITION footprint, which
# overstates: a frame's data can stop ~100 km short of its footprint, the cover counts that
# ground as covered, drops the granule that really holds it, and leaves a band along the frame
# boundary -- at the same latitudes on every track. Measured on
# NISAR_L2_PR_GUNW_024_071_A_168: footprint to lat -17.75, data only to -18.64.
VF=$(dirname "$CAT")/validFootprints.geojson
if [ "$SCAN" = 1 ] && [ ! -f "$VF" ]; then
    echo "=== $(date) scanning valid-data footprints (the slow step, ~1 h) ==="
    scanMasks "$CAT" --all --direction "$DIRECTION" --threads "$SCANTHREADS" --out "$VF" \
        || echo "  WARNING: the scan failed; tiling will fall back to acquisition footprints"
fi
vfArg=()
[ -f "$VF" ] && vfArg=(--validFootprints "$VF")

if [ -f "$T/tiles.geojson" ]; then
    echo "=== $(date) tiling exists, skipping ==="
else
    echo "=== $(date) tiling cycle $CYCLE ($TILEKM km, $CRS) ==="
    # --allGranules: every pair over a tile, not the fewest that cover it. geomosaic averages
    # overlapping inputs, so each cell is the MEAN coherence of every available pair rather
    # than whichever one a minimal cover picked. Costs more reads; the backscatter product
    # keeps the minimal cover.
    easeTiles --catalogue "$CAT" --hemisphere global --crs "$CRS" --tileKm "$TILEKM" \
        --direction "$DIRECTION" --allowV --allGranules "${vfArg[@]}" --out "$T" || exit 1
fi

echo "=== $(date) coherence mosaic at ${RES} m, $NPROC jobs ==="
gunwCoherence "$T" --work "$W" --psResM "$RES" --nProc "$NPROC" \
    || echo "  WARNING: some tiles failed (see $W/summary.log); rerun to retry them"

echo "=== $(date) figure ==="
# 0..1, never dB: coherence is a ratio. The same stretch for every cycle, so the series is
# comparable by eye.
globeFigure "$W/tiles/coherence" --out "$ROOT/coherence_cycle$CYCLE.png" \
    --width "${QLPX:-6000}" --coast --dbMin 0 --dbMax 1 \
    --title "NISAR GUNW coherence, cycle $CYCLE"

echo "=== $(date) done: $ROOT/coherence_cycle$CYCLE.png ==="
