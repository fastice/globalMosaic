#!/bin/bash
# Tar a run's Google Earth product without what Google Earth does not read: every stage/, the tiered
# builder's src10/src11 GeoTIFFs, and the original PNG of each tile oceanAltitude split into
# <y>.land.png + <y>.sea.png (kept on disk only for --undo).
#   bash test/tarGoogleEarth.sh /scratch/ianj/mosaics/cycle30ascending_z11          # one googleEarth.tar
#   bash test/tarGoogleEarth.sh /scratch/ianj/mosaics/cycle30ascending_z11 split    # ge_base.tar (with the
#        top doc.kml), ge_land.tar, ge_detail.tar -- unpack all into one place for the tiered product
set -euo pipefail
W=$(realpath "$1")
MODE=${2:-one}
cd "$W"
X=$(mktemp)
find googleEarth -name '*.sea.png' | sed 's/\.sea\.png$/.png/' > "$X"
echo "$(wc -l < "$X") split tiles: their original PNGs are left out"
EXC=(--exclude='*/stage' --exclude='googleEarth/src1[01]' -X "$X")
if [ "$MODE" = split ]; then
    for L in base land detail; do
        [ -d googleEarth/$L ] || continue
        extra=(); [ $L = base ] && extra=(googleEarth/doc.kml)
        tar -cf ge_$L.tar "${EXC[@]}" "${extra[@]}" googleEarth/$L
        echo "$W/ge_$L.tar $(du -h --apparent-size ge_$L.tar | cut -f1)"
    done
else
    tar -cf googleEarth.tar "${EXC[@]}" googleEarth
    echo "$W/googleEarth.tar $(du -h --apparent-size googleEarth.tar | cut -f1)"
fi
rm -f "$X"
