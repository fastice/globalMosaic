#!/bin/bash
# False-colour dual-pol composite TILES from a like-pol and a cross-pol tile directory.
#
#   bash test/runRatioTiles.sh LIKE_TILES CROSS_TILES OUT_DIR
#   bash test/runRatioTiles.sh /scratch/ianj/mosaics/cycle30ascending_z11/tiles/gamma0 \
#                              /scratch/ianj/mosaics/cycle30ascendingCross_z11/tiles/gamma0 \
#                              /scratch/ianj/mosaics/ratio30_z11
#   MODE=dB bash test/runRatioTiles.sh ...        # the dB ratio instead of the composite
#
# R = like-pol DN, G = cross-pol DN, B = the DN ratio (openSEPPO -dpratio). Where there is no
# cross-pol the three bands are equal, so single-pol ground -- the ice sheets -- renders grey.
#
# THE NUMBERS BELOW ARE A LOOK, NOT A STANDARD. SEPPO's scaling (which the published Earth Big
# Data map uses) is not documented, so these were settled by eye on the cycle-30 ascending z11
# product, 2026-11-01, and each one does something specific:
#
#   minCross -28   drop the ratio where cross-pol is at the noise floor (NESZ ~ -25 dB).
#                  Without it open water is the brightest thing in the image, because smooth
#                  water barely depolarises so HV sits at the floor and the ratio is huge.
#   G1 = 155       the balance point for green. 165 leaves an orange cast (R sits above G for
#                  typical vegetation); 144 overshoots into chartreuse.
#   B0 = 122       near the ratio band's true floor of 127 -- B = DN(like) - DN(cross) + 127
#                  and cross is never stronger than like, so nothing below 127 is physical.
#                  Values under it come from RFI, where HV exceeds HH.
#   grey 20 205    grey is a DIFFERENT population from dual-pol land and much darker. Sharing
#                  the red band's limits clipped Sahara sand (DN ~60) and blue ice to black;
#                  205 rather than 185 stops Greenland's interior saturating.
#
# To derive a new set for another product, do NOT stretch per tile -- each tile would get its
# own limits and the mosaic would come out patchwork. Instead:
#
#   gdalbuildvrt -te -180 -60 180 70 -tr 0.05 0.05 -r nearest -overwrite Ls.vrt <like>/global.vrt
#   gdalbuildvrt -te -180 -60 180 70 -tr 0.05 0.05 -r nearest -overwrite Cs.vrt <cross>/global.vrt
#   ratioMosaic Ls.vrt Cs.vrt --base Ls.vrt --mode rgb --minCross -28 --stretch 2 98 --out /dev/shm/x.tif
#
# and take the two printed lines. The latitude window keeps ice and sea ice out of the colour
# statistics; -r nearest rather than average matters because averaging smooths speckle and
# yields limits the full-resolution tiles do not match.
set -u
LIKE=${1:-}; CROSS=${2:-}; OUT=${3:-}
if [ -z "$OUT" ]; then
    echo "usage: $0 <likeTiles> <crossTiles> <outDir>" >&2; exit 2
fi
MODE=${MODE:-rgb}
NPROC=${NPROC:-8}
MINCROSS=${MINCROSS:--28}
LIMS=${LIMS:-"65 185 28 155 122 214"}
GREY=${GREY:-"20 205"}

mkdir -p "$OUT/$MODE"
echo "=== $(date) $MODE tiles: $(ls "$CROSS"/*.tif 2>/dev/null | wc -l) cross-pol tiles, $NPROC at a time"
[ "$MODE" = rgb ] && echo "    colour $LIMS   grey $GREY   minCross $MINCROSS"

# The like-pol tile of a pair is named for the LIKE term: .AHV. -> .AHH., .BVH. -> .BVV.
# A cross-pol tile with no match is skipped: the two covers can pick different frequency or
# polarization groups for a tile, and pairing an H-basis like-pol with a V-basis cross-pol
# would be meaningless rather than merely wrong-looking.
export LIKE OUT MODE MINCROSS LIMS GREY
ls "$CROSS"/*.tif | xargs -P "$NPROC" -I{} bash -c '
    f={}; n=$(basename "$f")
    ln=$(echo "$n" | sed -E "s/\.([AB])HV\./.\1HH./; s/\.([AB])VH\./.\1VV./")
    [ -f "$LIKE/$ln" ] || { echo "  no like-pol for $n"; exit 0; }
    if [ "$MODE" = rgb ]; then
        ratioMosaic "$LIKE/$ln" "$f" --mode rgb --minCross "$MINCROSS" --base "$LIKE/$ln" \
            --stretchLimits $LIMS --greyLimits $GREY --out "$OUT/$MODE/$n" > /dev/null
    else
        ratioMosaic "$LIKE/$ln" "$f" --mode "$MODE" --minCross "$MINCROSS" \
            --out "$OUT/$MODE/$n" > /dev/null
    fi'

N=$(ls "$OUT/$MODE" 2>/dev/null | wc -l)
echo "=== $(date) $N tiles -> $OUT/$MODE"
# Same guard as runCycleGCOV.sh: with no tiles, gdalbuildvrt writes nothing and globeFigure
# dies with "No such file or directory" -- after which the script used to print a cheerful
# "done:" line naming a PNG that was never written.
if [ "$N" -eq 0 ]; then
    echo "=== $(date) no tiles were built -- not building a VRT or a figure ===" >&2
    echo "    check that ratioMosaic is on PATH (pip install -e .) and that the tile names pair up" >&2
    exit 1
fi

# Lay the output out the way runGeomosaicTiles does -- vrt/<product>/<tile>.vrt per tile, a
# global.vrt, and a run.yaml naming the tiling -- so makeGoogleEarthTiered can consume it
# directly (it reads cfg['tileRun'] to find the per-tile CSVs and the catalogues beside them).
# TILERUN is the tiling the CROSS-pol mosaic was built from; without it the Google Earth step
# cannot run, but the tiles and the quicklook are still produced.
mkdir -p "$OUT/vrt/$MODE"
for f in "$OUT/$MODE"/*.tif; do
    b=$(basename "$f" .tif); echo "${b%%.*}"
done | sort -u | while read -r tile; do
    gdalbuildvrt -vrtnodata 0 -overwrite "$OUT/vrt/$MODE/$tile.vrt" \
        "$OUT/$MODE/$tile".*.tif > /dev/null 2>&1
done
gdalbuildvrt -vrtnodata 0 -overwrite "$OUT/vrt/$MODE/global.vrt" \
    "$OUT/vrt/$MODE"/*.vrt > /dev/null 2>&1
ln -sf "vrt/$MODE/global.vrt" "$OUT/$MODE.vrt" 2>/dev/null || \
    gdalbuildvrt -vrtnodata 0 -overwrite "$OUT/$MODE.vrt" "$OUT/$MODE"/*.tif > /dev/null
if [ -n "${TILERUN:-}" ]; then
    printf 'tileRun: %s\nwork: %s\ngranules: vsicurl\n' "$TILERUN" "$OUT" > "$OUT/run.yaml"
    echo "=== $(date) run.yaml -> tileRun $TILERUN   (makeGoogleEarthTiered $OUT --product $MODE)"
else
    echo "=== $(date) no TILERUN set: no run.yaml, so makeGoogleEarthTiered cannot run on this" >&2
fi
echo "=== $(date) $(ls "$OUT/vrt/$MODE"/*.vrt | wc -l) tile VRTs, $OUT/vrt/$MODE/global.vrt"

if [ "$MODE" = rgb ]; then
    globeFigure "$OUT/$MODE.vrt" --out "$OUT/${MODE}_quicklook.png" --width "${QLPX:-6000}" \
        --coast --title "${TITLE:-NISAR dual-pol composite}"
else
    globeFigure "$OUT/$MODE.vrt" --out "$OUT/${MODE}_quicklook.png" --width "${QLPX:-6000}" \
        --coast --dbMin "${VMIN:-2}" --dbMax "${VMAX:-11}" \
        --title "${TITLE:-NISAR dual-pol ratio (dB)}"
fi
echo "=== $(date) done: $OUT/${MODE}_quicklook.png"
