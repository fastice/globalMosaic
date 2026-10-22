#!/bin/bash
# Equal-area global mosaic on EASE-Grid 2.0 (EPSG:6931 north, 6932 south) at 1 km, both
# hemispheres, then screen-size quick looks. The two hemispheres meet at the equator, so
# together they cover the globe.
#
#   bash test/runEaseMosaic.sh                          # 1 km, vsicurl granules
#   RES=500 NPROC=24 bash test/runEaseMosaic.sh         # finer, more jobs
#   GRANULES=/Volumes/insar4/ian/Data/NISAR/Antarctica-GCOV bash test/runEaseMosaic.sh
#   HEMIS=south bash test/runEaseMosaic.sh              # one hemisphere only
#   DIRECTION=both bash test/runEaseMosaic.sh            # both orbit directions [A]
#   QLONLY=1 bash test/runEaseMosaic.sh                 # redo the quick looks, skip mosaicking
#
# Resumable: tiling is skipped when tiles.geojson exists, and runGeomosaicTiles skips finished
# tiles, so rerunning continues where it stopped.
set -euo pipefail
cd "$(dirname "$0")/.."

CAT=${CAT:-cycle30/catalogue.geojson}
ROOT=${ROOT:-/scratch/ianj/mosaics/ease}
RES=${RES:-1000}                 # output spacing, m
TILEKM=${TILEKM:-600}            # tile size, km (~5.4 deg, like the 6x6 deg lat/lon tiles)
NPROC=${NPROC:-16}
GRANULES=${GRANULES:-vsicurl}
PRODUCT=${PRODUCT:-gamma0}
QLPX=${QLPX:-1800}               # quick-look width in pixels (about screen size)
DBMIN=${DBMIN:--24}              # stretch, dB
DBMAX=${DBMAX:--1}
HEMIS=${HEMIS:-south north}
QLONLY=${QLONLY:-0}
DIRECTION=${DIRECTION:-A}         # A, D or both; ascending, as the Google Earth product is

mkdir -p "$ROOT/quicklooks"

for H in $HEMIS; do
    T=$ROOT/tiles-$H
    W=$ROOT/$H
    if [ "$QLONLY" = 0 ]; then
        if [ -f "$T/tiles.geojson" ]; then
            echo "=== $(date) $H: tiling exists, skipping ==="
        else
            echo "=== $(date) $H: tiling ($TILEKM km tiles) ==="
            easeTiles --catalogue "$CAT" --hemisphere "$H" --tileKm "$TILEKM" \
                      --direction "$DIRECTION" --out "$T"
        fi
        echo "=== $(date) $H: mosaicking at ${RES} m, $NPROC jobs ==="
        # runGeomosaicTiles exits 1 if ANY tile failed. Over hundreds of remote jobs a transient
        # failure is likely, and set -e would then abandon the quick looks and the other
        # hemisphere. Failed tiles are retried on the next run, so carry on and say so.
        runGeomosaicTiles "$T" --work "$W" --granules "$GRANULES" --psResM "$RES" --nProc "$NPROC" \
            || echo "  WARNING: $H had failed tiles (see $W/summary.log); continuing -- rerun to retry them"
    fi

    echo "=== $(date) $H: quick look ==="
    # every EASE tile of a hemisphere shares one CRS and grid, so a single VRT covers it
    shopt -s nullglob
    tifs=("$W/tiles/$PRODUCT"/*.tif)
    shopt -u nullglob
    if [ ${#tifs[@]} -eq 0 ]; then
        echo "  no tiles in $W/tiles/$PRODUCT -- skipped"
        continue
    fi
    vrt=$ROOT/quicklooks/$H.vrt
    gdalbuildvrt -q -overwrite "$vrt" "${tifs[@]}"
    # int16 dB x 100 -> byte; nodata -3000 becomes 0 (black) because JPEG has no alpha
    gdal_translate -q -of JPEG -ot Byte -outsize "$QLPX" 0 \
        -scale $((DBMIN * 100)) $((DBMAX * 100)) 1 255 \
        "$vrt" "$ROOT/quicklooks/$H.jpg"
    # ... and a PNG, where no-data stays transparent
    gdal_translate -q -of PNG -ot Byte -outsize "$QLPX" 0 -a_nodata 0 \
        -scale $((DBMIN * 100)) $((DBMAX * 100)) 1 255 \
        "$vrt" "$ROOT/quicklooks/$H.png"
    echo "  ${#tifs[@]} tiles -> $ROOT/quicklooks/$H.jpg (+ .png, $QLPX px wide)"
done

echo "=== $(date) done ==="
ls -l "$ROOT/quicklooks"/*.jpg "$ROOT/quicklooks"/*.png 2>/dev/null || true
