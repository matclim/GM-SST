#!/usr/bin/env bash
# tests/run_all.sh — validate the straw-drift reco across directions & momenta.
# Runnable from anywhere; resolves build dir and field map relative to itself.
#
#   ./tests/run_all.sh
#   FIELD=/path/map.root OUTDIR=/tmp/val ./tests/run_all.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$HERE/.." && pwd)"
BUILD="${BUILD:-$ROOT_DIR/build}"
FIELD="${FIELD:-$ROOT_DIR/../fieldmaps/NEWMAPS/2026_07_02_MainSpectrometerField_V21_2455.root}"
OUTDIR="${OUTDIR:-$BUILD}"
STRAWS="${STRAWS:-$OUTDIR/straws.root}"
ORIGIN_Z=500
NEV=200
POS="0 0 24000"

[ -x "$BUILD/run_StrawTracker" ] || { echo "build run_StrawTracker first"; exit 1; }
[ -x "$BUILD/run_reco" ]         || { echo "build run_reco first"; exit 1; }
[ -f "$FIELD" ]                  || { echo "no field map: $FIELD"; exit 1; }
mkdir -p "$OUTDIR"

# Straw table (geometry) — dump once if missing.
[ -f "$STRAWS" ] || { echo ">>> dumping straw table"; "$BUILD/run_StrawTracker" --dump-straws "$STRAWS"; }

# name         dir(x y z)        energyMeV  seed
CONFIGS=(
  "straight_10   0     0    1     10000   1"
  "tiltx_10      0.02  0    1     10000   2"
  "tilty_10      0     0.02 1     10000   3"
  "ntilty_10     0    -0.02 1     10000   4"
  "diag_10       0.02  0.02 1     10000   5"
  "bigtilty_10   0     0.05 1     10000   6"
  "straight_5    0     0    1      5000   7"
  "straight_20   0     0    1     20000   8"
)

for cfg in "${CONFIGS[@]}"; do
  read -r NAME DX DY DZ E SEED <<< "$cfg"
  echo "=============================================================="
  echo ">>> $NAME  dir=($DX $DY $DZ)  E=$E MeV  seed=$SEED"
  HITS="$OUTDIR/muons_${NAME}.root"
  RECO="$OUTDIR/reco_${NAME}.root"

  "$BUILD/run_StrawTracker" \
    --n-events "$NEV" --seed "$SEED" \
    --particle mu- --energy-MeV "$E" \
    --pos-mm $POS --dir "$DX" "$DY" "$DZ" \
    --field-map "$FIELD" --output "$HITS"

  IN="$HITS"; [ -f "$IN" ] || IN="${HITS%.root}_t-1.root"
  [ -f "$IN" ] || { echo "!! no hits for $NAME"; exit 1; }

  "$BUILD/run_reco" \
    --hits "$IN" --field "$FIELD" --straws "$STRAWS" \
    --p-guess "$E" --p-true "$E" --field-origin-z "$ORIGIN_Z" \
    --reco-out "$RECO"
done

echo "=============================================================="
echo "Done. Summarise with:  root -l -b -q $HERE/summary.C"
