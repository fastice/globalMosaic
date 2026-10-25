#!/bin/bash
# Backscatter mosaic for ONE cycle, on the same Wagner VII grid as the coherence series.
#
#   bash test/runCycleGCOV.sh 024
#   RES=500 NPROC=16 bash test/runCycleGCOV.sh 025
#
# One cycle per run and NO --fillCatalogues: each map is that cycle's own acquisitions, so the
# series shows change rather than a blend. That is the difference from runEaseMosaic.sh, which
# builds a single best-coverage mosaic from every cycle at once.
#
# Minimal cover, not --allGranules: for backscatter one good observation per place is what is
# wanted, with the better bandwidth preferred. (The coherence series uses --allGranules so
# that overlapping pairs average.)
#
# The grid depends only on CRS and TILEKM, never on the data, so every cycle -- and the
# coherence series -- lands on the same pixels and the maps difference directly.
set -u
CYCLE=${1:-}
if [ -z "$CYCLE" ]; then echo "usage: $0 <cycle, e.g. 024>" >&2; exit 2; fi

cd "$(dirname "$0")/.." || exit 1
ROOT=${ROOT:-/scratch/ianj/mosaics/gcov}
CRS=${CRS:-ESRI:54076}
TILEKM=${TILEKM:-600}
RES=${RES:-1000}
NPROC=${NPROC:-32}
DIRECTION=${DIRECTION:-A}
GRANULES=${GRANULES:-vsicurl}
SHORT=$((10#$CYCLE))
W=$ROOT/cycle$CYCLE
T=$W/tiles-global

mkdir -p "$W"

# The cycle catalogues live in cycle30/ from the Google Earth work: catalogue0NN.geojson, and
# catalogue.geojson for cycle 030. Anything missing is searched from ASF once.
CAT=${CAT:-}
if [ -z "$CAT" ]; then
    for c in "cycle30/catalogue$CYCLE.geojson" "cycle30/catalogue0$SHORT.geojson"; do
        [ -f "$c" ] && { CAT=$c; break; }
    done
    if [ -z "$CAT" ] && [ "$SHORT" = 30 ]; then CAT=cycle30/catalogue.geojson; fi
fi
if [ -z "$CAT" ] || [ ! -f "$CAT" ]; then
    CAT=cycle${CYCLE}gcov/catalogue.geojson
    if [ ! -f "$CAT" ]; then
        echo "=== $(date) searching ASF for cycle $CYCLE GCOV ==="
        globalGCOVTiles --cycle "$CYCLE" --product GCOV --direction "$DIRECTION" \
            --bbox -180 -90 180 90 --out "cycle${CYCLE}gcov" || exit 1
    fi
fi
echo "=== $(date) cycle $CYCLE catalogue: $CAT ==="

# Valid-data footprints, where they exist for these granules. They are not what fixed the
# banding (that was the antimeridian split) but they keep frames that hold no data from being
# chosen over frames that do.
vfArg=()
ls cycle30/validFootprints*.geojson >/dev/null 2>&1 && \
    vfArg=(--validFootprints $(ls cycle30/validFootprints*.geojson))

if [ -f "$T/tiles.geojson" ]; then
    echo "=== $(date) tiling exists, skipping ==="
else
    echo "=== $(date) tiling cycle $CYCLE ($TILEKM km, $CRS) ==="
    easeTiles --catalogue "$CAT" --hemisphere global --crs "$CRS" --tileKm "$TILEKM" \
        --direction "$DIRECTION" --allowV "${vfArg[@]}" --out "$T" || exit 1
fi

echo "=== $(date) mosaicking at ${RES} m, $NPROC jobs ==="
runGeomosaicTiles "$T" --work "$W" --granules "$GRANULES" --psResM "$RES" --nProc "$NPROC" \
    || echo "  WARNING: some tiles failed (see $W/summary.log); rerun to retry them"

echo "=== $(date) figure ==="
SFX=""
[ "$DIRECTION" != A ] && SFX="_$DIRECTION"
globeFigure "$W/tiles/gamma0" --out "$ROOT/gamma0_cycle$CYCLE$SFX.png" \
    --width "${QLPX:-6000}" --coast \
    --title "NISAR GCOV gamma0, cycle $CYCLE, direction $DIRECTION"

echo "=== $(date) done: $ROOT/gamma0_cycle$CYCLE$SFX.png ==="
