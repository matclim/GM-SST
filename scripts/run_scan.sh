#!/usr/bin/env bash
# run_scan.sh ─────────────────────────────────────────────────────────────────
# Grid-scan driver for the straw-tracker simulation.
#
# Sweeps three axes:
#   Momentum: 5, 10, 20, 50, 100 GeV/c
#   Angle:    0, 5, 15 degrees, tilted in the y–z (bending) plane
#   Position: (x, y) over the full tracker aperture, step STEP_CM (default 10)
#
# For every (p, θ, x, y) the simulation is run with NEVENTS events. Output
# ROOT files are organised as:
#   $OUTDIR/p<P>GeV_a<A>deg/hits_x<X>cm_y<Y>cm.root
#
# Existing output files are skipped, so the script is resumable: re-run it
# after a crash/interruption and it picks up where it left off.
#
# Usage:
#   ./run_scan.sh              # run with defaults (asks for confirmation)
#   YES=1 ./run_scan.sh        # skip confirmation
#   DRYRUN=1 ./run_scan.sh     # print the config and the first few jobs, exit
#   STEP_CM=25 ./run_scan.sh   # coarser grid
#   JOBS=8 ./run_scan.sh       # 8 parallel workers
#
# Env vars:
#   OUTDIR     output root directory               (default: scan_output)
#   STEP_CM    grid step in cm                     (default: 10)
#   NEVENTS    events per configuration            (default: 100)
#   JOBS       parallel workers                    (default: 4)
#   EXEC       path to run_StrawTracker            (default: ./run_StrawTracker)
#   FIELDMAP   path to the field-map file          (default: ../../fieldmaps/V21_2000A_B_full_fieldmap_SHiP.txt)
#   POSZ_MM    z of the muon gun in lab frame      (default: 24000)
#   YES        1 to skip confirmation prompt       (default: 0)
#   DRYRUN     1 to print config and exit          (default: 0)
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

# ── Defaults (override via env) ───────────────────────────────────────────────
: "${OUTDIR:=scan_output}"
: "${STEP_CM:=10}"
: "${NEVENTS:=100}"
: "${JOBS:=4}"
: "${EXEC:=./run_StrawTracker}"
: "${FIELDMAP:=../../fieldmaps/V21_2000A_B_full_fieldmap_SHiP.txt}"
: "${POSZ_MM:=24000}"
: "${YES:=0}"
: "${DRYRUN:=0}"

# Grid extents (cm). 399 × 599 cm matches the stock tracker aperture.
XMIN_CM=-199
XMAX_CM=199
YMIN_CM=-299
YMAX_CM=299

MOMENTA=(5 10 20 50 100)     # GeV/c
ANGLES=(0 5 15)              # degrees, tilted in the y–z plane

# ── Helpers ──────────────────────────────────────────────────────────────────

# Muon kinetic energy (MeV) for a target lab momentum p (GeV/c).
# E_kin = sqrt(p² + m²) − m, with m = 105.66 MeV.
p_to_Ekin() {
    awk -v p="$1" 'BEGIN {
        p *= 1000           # GeV → MeV
        m  = 105.66
        printf "%.3f", sqrt(p*p + m*m) - m
    }'
}

# Direction components "dy dz" for polar angle θ (deg) from +z, tilted in y–z.
angle_to_dy_dz() {
    awk -v th="$1" 'BEGIN {
        pi = 3.141592653589793
        t  = th * pi / 180
        printf "%.10f %.10f", sin(t), cos(t)
    }'
}

# ── Sanity checks ────────────────────────────────────────────────────────────

if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: executable '$EXEC' not found or not executable." >&2
    exit 1
fi
if [[ ! -r "$FIELDMAP" ]]; then
    echo "WARNING: field-map file '$FIELDMAP' not readable. Jobs will fall back" >&2
    echo "         to the uniform dipole unless this path is fixed."           >&2
fi
if (( STEP_CM < 1 )); then
    echo "ERROR: STEP_CM must be >= 1 (got $STEP_CM)." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Job generator ────────────────────────────────────────────────────────────
# Emits one shell command per line on stdout. Commands are independent and
# can be run in any order / in parallel.

generate_jobs() {
    local p_gev e_mev ang dy dz x_cm y_cm x_mm y_mm subdir outfile
    for p_gev in "${MOMENTA[@]}"; do
        e_mev=$(p_to_Ekin "$p_gev")
        for ang in "${ANGLES[@]}"; do
            read -r dy dz <<< "$(angle_to_dy_dz "$ang")"
            subdir="$OUTDIR/p${p_gev}GeV_a${ang}deg"
            mkdir -p "$subdir"

            for (( x_cm = XMIN_CM ; x_cm <= XMAX_CM ; x_cm += STEP_CM )); do
                for (( y_cm = YMIN_CM ; y_cm <= YMAX_CM ; y_cm += STEP_CM )); do
                    x_mm=$(( x_cm * 10 ))
                    y_mm=$(( y_cm * 10 ))
                    outfile="$subdir/hits_x${x_cm}cm_y${y_cm}cm.root"

                    # Resume support: Geant4 appends a thread-id suffix
                    # (_t0 / _t-1 / ...) to the --output name, so glob.
                    if compgen -G "${outfile%.root}*.root" > /dev/null; then
                        continue
                    fi

                    printf '%s --n-events %d --output %s --particle mu- --energy-MeV %s --pos-mm %d %d %d --dir 0 %s %s --field-map %s\n' \
                        "$EXEC" "$NEVENTS" "$outfile" "$e_mev" \
                        "$x_mm" "$y_mm" "$POSZ_MM" "$dy" "$dz" "$FIELDMAP"
                done
            done
        done
    done
}

# ── Summary ──────────────────────────────────────────────────────────────────

JOBFILE=$(mktemp)
trap 'rm -f "$JOBFILE"' EXIT
generate_jobs > "$JOBFILE"
N_JOBS=$(wc -l < "$JOBFILE")

# Theoretical grid size (before skipping already-done outputs).
N_X=$(( (XMAX_CM - XMIN_CM) / STEP_CM + 1 ))
N_Y=$(( (YMAX_CM - YMIN_CM) / STEP_CM + 1 ))
N_GRID=$(( N_X * N_Y * ${#MOMENTA[@]} * ${#ANGLES[@]} ))
N_EVT=$(( N_JOBS * NEVENTS ))

echo "# ── Scan configuration ──────────────────────────────────────────────"
echo "#   momenta (GeV/c) : ${MOMENTA[*]}"
echo "#   angles  (deg)   : ${ANGLES[*]}"
echo "#   x range (cm)    : [$XMIN_CM, $XMAX_CM]  (${N_X} points)"
echo "#   y range (cm)    : [$YMIN_CM, $YMAX_CM]  (${N_Y} points)"
echo "#   grid step (cm)  : $STEP_CM"
echo "#   events/config   : $NEVENTS"
echo "#   gun z (mm)      : $POSZ_MM"
echo "# ── Workload ────────────────────────────────────────────────────────"
echo "#   total configs   : $N_GRID"
echo "#   jobs pending    : $N_JOBS   (others already have output)"
echo "#   total events    : $N_EVT"
echo "#   parallel workers: $JOBS"
echo "#   output dir      : $OUTDIR"
echo "# ────────────────────────────────────────────────────────────────────"

if (( N_JOBS == 0 )); then
    echo "Nothing to do — all outputs already exist."
    exit 0
fi

if [[ "$DRYRUN" == "1" ]]; then
    echo ""
    echo "Dry run — first 3 commands:"
    head -3 "$JOBFILE" | sed 's/^/  /'
    exit 0
fi

if [[ "$YES" != "1" ]]; then
    read -r -p "Proceed with $N_JOBS jobs? [y/N] " ans
    [[ "$ans" == "y" || "$ans" == "Y" ]] || { echo "Cancelled."; exit 0; }
fi

# ── Execution ────────────────────────────────────────────────────────────────

LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"

echo "Starting $N_JOBS jobs with $JOBS parallel workers."
echo "Per-job stdout/stderr → $LOGDIR/"
echo ""

# Run each command, redirecting its log to a numbered file so failures stay
# diagnosable but the console stays quiet. xargs -P gives us a worker pool.
JOB_COUNTER_FILE=$(mktemp)
echo 0 > "$JOB_COUNTER_FILE"

run_one() {
    local cmd="$1"
    local n
    # Atomic counter for a unique log-file id.
    n=$(flock "$JOB_COUNTER_FILE" bash -c '
        read -r x < "$0"
        echo $(( x + 1 )) > "$0"
        echo "$x"' "$JOB_COUNTER_FILE")
    local log="$LOGDIR/job_$(printf '%06d' "$n").log"
    echo "=== $(date -u +%FT%TZ) === $cmd" > "$log"
    eval "$cmd" >> "$log" 2>&1
    local rc=$?
    echo "=== exit $rc ===" >> "$log"
    if (( rc != 0 )); then
        echo "FAILED (rc=$rc): see $log" >&2
    fi
    return $rc
}
export -f run_one
export JOB_COUNTER_FILE LOGDIR

xargs -a "$JOBFILE" -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
XARGS_RC=$?

rm -f "$JOB_COUNTER_FILE"

N_FAIL=$(grep -l '^=== exit [1-9]' "$LOGDIR"/*.log 2>/dev/null | wc -l)
N_OK=$(( N_JOBS - N_FAIL ))

echo ""
echo "# ── Results ────────────────────────────────────────────────────────"
echo "#   completed : $N_OK / $N_JOBS"
echo "#   failures  : $N_FAIL"
echo "#   outputs   : $OUTDIR/"
echo "#   logs      : $LOGDIR/"
if (( N_FAIL > 0 )); then
    echo ""
    echo "To list failing jobs:"
    echo "  grep -l '^=== exit [1-9]' $LOGDIR/*.log"
fi

exit "$XARGS_RC"
