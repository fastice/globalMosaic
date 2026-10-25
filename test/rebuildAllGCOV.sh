#!/bin/bash
# Wait, then delete the existing backscatter products and rebuild cycles 24..32 from scratch.
#
#   SLEEPHOURS=4.5 nohup bash test/rebuildAllGCOV.sh > /scratch/ianj/mosaics/gcov.log 2>&1 &
#   CYCLES="024 025" SLEEPHOURS=0 bash test/rebuildAllGCOV.sh
#
# The wait is so this can be started before leaving while another job finishes. Nothing is
# deleted and nothing is killed until the wait is over.
set -u
ROOT=${ROOT:-/scratch/ianj/mosaics/gcov}
LEGACY=${LEGACY:-/scratch/ianj/mosaics/ease}
CYCLES=${CYCLES:-"024 025 026 027 028 029 030 031 032"}
SLEEPHOURS=${SLEEPHOURS:-4.5}
export NPROC=${NPROC:-32}
export RES=${RES:-1000}
export DIRECTION=${DIRECTION:-A}
export ROOT
REMOVELEGACY=${REMOVELEGACY:-1}

cd "$(dirname "$0")/.." || exit 1
SECS=$(python3 -c "print(int(float('$SLEEPHOURS') * 3600))")
echo "=== $(date) waiting ${SLEEPHOURS} h (${SECS} s) before touching anything ==="
echo "    then: stop running jobs, delete $ROOT/cycle0*"
[ "$REMOVELEGACY" = 1 ] && echo "          and the single-mosaic product $LEGACY"
echo "    then: build cycles $CYCLES at ${RES} m, $NPROC jobs, direction $DIRECTION"
sleep "$SECS"

P='[g]unwCoherence|[r]unCycleCoherence|[r]unCycleGCOV|[g]eomosaic|[e]aseTiles|[r]unGeomosaicTiles|[s]canMasks|[r]unEaseMosaic'
if pgrep -f "$P" >/dev/null; then
    echo "=== $(date) stopping what is still running ==="
    pgrep -af "$P" | cut -c1-100
    pkill -f "$P"; sleep 5; pkill -9 -f "$P"; sleep 2
fi
if pgrep -f "$P" >/dev/null; then
    echo "SOMETHING IS STILL RUNNING -- stopping rather than deleting under it" >&2
    pgrep -af "$P" | cut -c1-100 >&2
    exit 1
fi

echo "=== $(date) removing previous backscatter products ==="
rm -rf "$ROOT"/cycle0* "$ROOT"/gamma0_cycle0*.png
if [ "$REMOVELEGACY" = 1 ] && [ -d "$LEGACY" ]; then
    echo "    removing $LEGACY (the old single global mosaic)"
    rm -rf "$LEGACY"
fi
mkdir -p "$ROOT"

for c in $CYCLES; do
    echo "######## $(date) cycle $c ########"
    bash test/runCycleGCOV.sh "$c" || echo "  cycle $c FAILED; continuing"
done
echo "######## $(date) all cycles done ########"
ls -l "$ROOT"/gamma0_cycle0*.png 2>/dev/null
