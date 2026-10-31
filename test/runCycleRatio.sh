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

# The assembled mosaic has different names by CRS: a PROJECTED global run (Wagner VII, EASE,
# Equal Earth) writes vrt/<product>/ease_<epsg>.vrt and deliberately NO global.vrt, because it
# has no latitude bands; a lat/lon run writes global.vrt. Look for either, or the ratio step
# aborts after both mosaics have already been built.
pickVrt() {
    local d=$1 f
    for f in "$d"/ease_*.vrt "$d"/global.vrt; do
        [ -f "$f" ] && { echo "$f"; return 0; }
    done
    return 1
}
L=$(pickVrt "$W/like/vrt/gamma0") || {
    echo "no assembled like-pol VRT in $W/like/vrt/gamma0 -- not building the ratio" >&2
    echo "    see $W/like/summary.log" >&2; exit 1; }
C=$(pickVrt "$W/cross/vrt/gamma0") || {
    echo "no assembled cross-pol VRT in $W/cross/vrt/gamma0 -- not building the ratio" >&2
    echo "    see $W/cross/summary.log" >&2; exit 1; }
echo "=== $(date) like  $L"
echo "=== $(date) cross $C"

SFX=""
[ "$DIRECTION" != A ] && SFX="_$DIRECTION"
echo "=== $(date) ratio ($MODE) ==="
# MINCROSS drops pixels whose cross-pol is at the noise floor, where the ratio is a lower
# bound rather than a measurement -- open water above all, which otherwise renders as the
# brightest thing in the image. Off unless asked for: MINCROSS=-28 bash test/runCycleRatio.sh
mcArg=()
[ -n "${MINCROSS:-}" ] && mcArg=(--minCross "$MINCROSS")
ratioMosaic "$L" "$C" --mode "$MODE" "${mcArg[@]}" \
    --out "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif" || exit 1

# globeFigure's defaults are wrong for a ratio and fail SILENTLY: for a Float32 raster it
# stretches 0..1 (the coherence range), and a dB ratio runs 0..25, so every pixel saturates
# and the figure comes out ALL WHITE. Give it a range per mode. rgb is skipped entirely --
# globeFigure renders band 1 only, so it would show the red channel alone, not the composite.
case "$MODE" in
    dB)  VMIN=${VMIN:-0};    VMAX=${VMAX:-15} ;;
    pwr) VMIN=${VMIN:-1};    VMAX=${VMAX:-30} ;;
    amp) VMIN=${VMIN:-1000}; VMAX=${VMAX:-10000} ;;
    DN)  VMIN=${VMIN:-127};  VMAX=${VMAX:-255} ;;
    rgb) VMIN=; VMAX= ;;
esac
if [ "$MODE" = rgb ]; then
    echo "=== $(date) no figure for rgb: it is a 3-band composite, open the tif directly ==="
else
    echo "=== $(date) figure (stretch $VMIN .. $VMAX) ==="
    globeFigure "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif" \
        --out "$ROOT/ratio${MODE}_cycle$CYCLE$SFX.png" --width "${QLPX:-6000}" --coast \
        --dbMin "$VMIN" --dbMax "$VMAX" \
        --title "NISAR dual-pol ratio ($MODE), cycle $CYCLE, direction $DIRECTION"
fi

echo "=== $(date) done: $ROOT/ratio${MODE}_cycle$CYCLE$SFX.tif ==="
