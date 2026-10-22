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
#   HEMIS=global bash test/runEaseMosaic.sh              # one global map, Equal Earth
#   HEMIS=global CRS=EPSG:6933 bash test/runEaseMosaic.sh  # ... in EASE-Grid 2.0 Global instead
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
NPROC=${NPROC:-32}       # geomosaic sits near 30% CPU: it waits on the network, so oversubscribe
GRANULES=${GRANULES:-vsicurl}
PRODUCT=${PRODUCT:-gamma0}
QLPX=${QLPX:-1800}               # quick-look width in pixels (about screen size)
DBMIN=${DBMIN:--24}              # stretch, dB
DBMAX=${DBMAX:--1}
HEMIS=${HEMIS:-south north}
QLONLY=${QLONLY:-0}
DIRECTION=${DIRECTION:-A}         # A, D or both; ascending, as the Google Earth product is
# --hemisphere global only. Empty = easeTiles' default, Equal Earth, which has no EPSG code and
# so reaches geomosaic as -wkt: that needs the GP_GENERIC SRS support (installGeomosaic builds
# it). EPSG:6933 (EASE-Grid 2.0 Global, cylindrical equal-area) is also equal-area and goes in
# as -epsg, which every geomosaic understands -- no rebuild.
CRS=${CRS:-}
# Granules from other cycles, used ONLY where the main catalogue leaves a tile uncovered --
# the same mechanism the Google Earth product used (globalGCOVTiles --fillCatalogues). One
# cycle of ascending NISAR does not cover the Earth: cycle30 alone left the Sahara and much
# of the sea ice empty. Set to "" for a single-cycle run.
FILLCATS=${FILLCATS:-cycle30/catalogue029.geojson cycle30/catalogue028.geojson cycle30/catalogue025.geojson cycle30/catalogue024.geojson cycle30/catalogue023.geojson cycle30/catalogue031.geojson cycle30/catalogue026.geojson cycle30/catalogue027.geojson}

mkdir -p "$ROOT/quicklooks"

for H in $HEMIS; do
    T=$ROOT/tiles-$H
    W=$ROOT/$H
    if [ "$QLONLY" = 0 ]; then
        if [ -f "$T/tiles.geojson" ]; then
            echo "=== $(date) $H: tiling exists, skipping ==="
        else
            echo "=== $(date) $H: tiling ($TILEKM km tiles) ==="
            crsArg=()
            [ "$H" = global ] && [ -n "$CRS" ] && crsArg=(--crs "$CRS")
            fillArg=()
            [ -n "$FILLCATS" ] && fillArg=(--fillCatalogues $FILLCATS)
            easeTiles --catalogue "$CAT" --hemisphere "$H" --tileKm "$TILEKM" \
                      --direction "$DIRECTION" "${crsArg[@]}" "${fillArg[@]}" --out "$T"
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
    # int16 dB x 100 -> byte; nodata -3000 becomes 0 (black) because JPEG has no alpha.
    # -r average, not the default nearest: the quicklook decimates ~20x and NISAR coverage is
    # narrow strips, so nearest keeps one pixel in 400 and thins the map for no reason.
    gdal_translate -q -of JPEG -ot Byte -outsize "$QLPX" 0 -r average \
        -scale $((DBMIN * 100)) $((DBMAX * 100)) 1 255 \
        "$vrt" "$ROOT/quicklooks/$H.jpg"
    # ... and a PNG, where no-data stays transparent
    gdal_translate -q -of PNG -ot Byte -outsize "$QLPX" 0 -a_nodata 0 -r average \
        -scale $((DBMIN * 100)) $((DBMAX * 100)) 1 255 \
        "$vrt" "$ROOT/quicklooks/$H.png"
    echo "  ${#tifs[@]} tiles -> $ROOT/quicklooks/$H.jpg (+ .png, $QLPX px wide)"
done

echo "=== $(date) done ==="
ls -l "$ROOT/quicklooks"/*.jpg "$ROOT/quicklooks"/*.png 2>/dev/null || true
