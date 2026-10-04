#!/bin/bash
# Tar a run's Google Earth product without what Google Earth does not read: every stage/, the tiered
# builder's src10/src11/src20 GeoTIFFs, and the original PNG of each tile oceanAltitude split into
# <y>.land.png + <y>.sea.png (kept on disk only for --undo).
#   bash test/tarGoogleEarth.sh /scratch/ianj/mosaics/cycle30ascending_z11          # one googleEarth.tar
#   bash test/tarGoogleEarth.sh /scratch/ianj/mosaics/cycle30ascending_z11 split    # ge_<layer>.tar for every
#        layer present (base, land, detail, detail20), each with the top doc.kml -- unpack into one place
#   bash test/tarGoogleEarth.sh /scratch/ianj/mosaics/cycle30ascending_z11 split detail20   # only these layers
set -euo pipefail
W=$(realpath "$1")
MODE=${2:-one}
cd "$W"
X=$(mktemp)
find googleEarth -name '*.sea.png' | sed 's/\.sea\.png$/.png/' > "$X"
echo "$(wc -l < "$X") split tiles: their original PNGs are left out"
EXC=(--exclude='*/stage' --exclude='googleEarth/src[0-9]*' -X "$X")
if [ "$MODE" = split ]; then
    shift 2 || true
    LAYERS=${*:-base land detail detail20}
    for L in $LAYERS; do
        [ -d googleEarth/$L ] || continue
        tar -cf ge_$L.tar "${EXC[@]}" googleEarth/doc.kml googleEarth/$L
        echo "$W/ge_$L.tar $(du -h --apparent-size ge_$L.tar | cut -f1)"
    done
else
    tar -cf googleEarth.tar "${EXC[@]}" googleEarth
    echo "$W/googleEarth.tar $(du -h --apparent-size googleEarth.tar | cut -f1)"
fi
rm -f "$X"
