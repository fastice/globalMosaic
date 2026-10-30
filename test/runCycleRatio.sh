#!/bin/bash
# Dual-pol ratio mosaic for ONE cycle: like-pol and cross-pol on one grid, then the ratio.
#
#   bash test/runCycleRatio.sh 033
#   MODE=amp RES=500 NPROC=40 bash test/runCycleRatio.sh 033
#
# Follows openSEPPO's -dpratio (doc/ratio.md): like/cross, with -pwr/-amp/-DN/-dB. geomosaic
# writes gamma0 as int16 dB x 100, so the ratio is arithmetic on two dB mosaics -- see
# ratioMosaic.
#
# TWO tilings, not one, and this is the point of the whole script: the cover is restricted to
# granules that HAVE a cross-pol term (--ratioCover / --crossPol drop SH and SV), so the same
# granule wins each pixel in both passes and the ratio compares one acquisition with itself.
# The selection is deterministic given the same candidates, so the two easeTiles runs agree.
#
# An existing like-pol mosaic from runCycleGCOV.sh CANNOT be substituted here. Its cover
# includes single-pol granules; the cross-pol pass would skip those and fill the pixels from a
# different acquisition, and the ratio would be meaningless in exactly the places it looked
# fine. The like-pol built here is a slightly sparser product than that one, by design.
#
# Both mosaics skip tiles that are already done (runGeomosaicTiles returns "skipped (done)"
# when out.<product>.tif exists), so a re-run costs only what is missing.
set -u
CYCLE=${1:-}
if [ -z "$CYCLE" ]; then echo "usage: $0 <cycle, e.g. 033>" >&2; exit 2; fi

cd "$(dirname "$0")/.." || exit 1
ROOT=${ROOT:-/scratch/ianj/mosaics/ratio}
CRS=${CRS:-ESRI:54076}
TILEKM=${TILEKM:-600}
RES=${RES:-1000}
NPROC=${NPROC:-80}
DIRECTION=${DIRECTION:-A}
GRANULES=${GRANULES:-vsicurl}
MODE=${MODE:-dB}
SHORT=$((10#$CYCLE))
W=$ROOT/cycle$CYCLE
TL=$W/tiles-like
TC=$W/tiles-cross
mkdir -p "$W"

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

vfArg=()
ls cycle30/validFootprints*.geojson >/dev/null 2>&1 && \
    vfArg=(--validFootprints $(ls cycle30/validFootprints*.geojson))

# One cover, two terms. --allowV so DV granules (VV/VH) reach it at all.
for pair in "like:$TL:--ratioCover" "cross:$TC:--crossPol"; do
    lab=${pair%%:*}; rest=${pair#*:}; dir=${rest%%:*}; flag=${rest#*:}
    if [ -f "$dir/tiles.geojson" ]; then
        echo "=== $(date) $lab tiling exists, skipping ==="
    else
        echo "=== $(date) tiling $lab ($TILEKM km, $CRS, cross-pol-capable granules only) ==="
        easeTiles --catalogue "$CAT" --hemisphere global --crs "$CRS" --tileKm "$TILEKM" \
            --direction "$DIRECTION" --allowV "$flag" "${vfArg[@]}" --out "$dir" || exit 1
    fi
done

for pair in "like:$TL" "cross:$TC"; do
    lab=${pair%%:*}; dir=${pair#*:}
    echo "=== $(date) mosaicking $lab at ${RES} m, $NPROC jobs ==="
    runGeomosaicTiles "$dir" --work "$W/$lab" --granules "$GRANULES" --psResM "$RES" \
        --nProc "$NPROC" \
        || echo "  WARNING: some $lab tiles failed (see $W/$lab/summary.log); rerun to retry"
done

L=$W/like/vrt/gamma0/global.vrt
C=$W/cross/vrt/gamma0/global.vrt
for f in "$L" "$C"; do
    [ -f "$f" ] || { echo "missing $f -- not building the ratio" >&2; exit 1; }
done

SFX=""
[ "$DIRECTION" != A ] && SFX="_$DIRECTION"
echo "=== $(date) ratio ($MODE) ==="
ratioMosaic "$L" "$C" --mode "$MODE" --out "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif" || exit 1

echo "=== $(date) figure ==="
globeFigure "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif" \
    --out "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.png" --width "${QLPX:-6000}" --coast \
    --title "NISAR dual-pol ratio ($MODE), cycle $CYCLE, direction $DIRECTION"

echo "=== $(date) done: $ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif ==="
