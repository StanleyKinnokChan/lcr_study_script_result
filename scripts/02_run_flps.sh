#!/bin/bash
# Run fLPS 2.0 on all downloaded proteomes (PARALLEL VERSION WITH PROGRESS).
# Parameters match Teekas et al. (2024) Open Biology:
#   - min length: 3 aa (-m 3)
#   - max unique residues filtered in 03_analyse.py (fLPS2 has no -u flag)
#   - exclude X residues (default in fLPS2)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
export FLPS="$PROJECT_DIR/flps2/flps"
export PROTEOME_DIR="$PROJECT_DIR/data/proteomes"
export OUT_DIR="$PROJECT_DIR/results/flps"

mkdir -p "$OUT_DIR"

if [ ! -f "$FLPS" ]; then
    echo "ERROR: fLPS binary not found at $FLPS"
    echo "Run scripts/00_setup.sh first."
    exit 1
fi

# Detect number of CPU cores
CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# 1. Pre-calculate the total number of files to process
TOTAL=$(find "$PROTEOME_DIR" -maxdepth 1 -name "*.longest.fa" 2>/dev/null | wc -l | tr -d ' ')

if [ "$TOTAL" -eq 0 ]; then
    echo "No .longest.fa files found in $PROTEOME_DIR"
    exit 0
fi

echo "fLPS binary: $FLPS"
echo "Input dir:   $PROTEOME_DIR"
echo "Output dir:  $OUT_DIR"
echo "Total files: $TOTAL"
echo "Running on:  $CORES concurrent threads"
echo "--------------------------------------------------"

# 2. Define the worker function
run_flps() {
    local fa="$1"
    [ -f "$fa" ] || return 0
    
    local species=$(basename "$fa" .longest.fa)
    local out="$OUT_DIR/${species}.flps.txt"

    if [ -f "$out" ]; then
        echo "[$species] Already done, skipping."
    else
        echo "[$species] Running fLPS..."
        "$FLPS" -m 3 "$fa" > "$out"
        echo "[$species] finished -> $(wc -l < "$out") lines"
    fi
    
    # Send a secret marker to stdout so the parent script knows a job finished
    echo "__JOB_DONE__"
}
export -f run_flps

# 3. Execute in parallel and pipe the output to a while loop to track progress safely
DONE=0

find "$PROTEOME_DIR" -maxdepth 1 -name "*.longest.fa" -print0 | \
    xargs -0 -P "$CORES" -I {} bash -c 'run_flps "{}"' | \
    while IFS= read -r line; do
        if [[ "$line" == "__JOB_DONE__" ]]; then
            DONE=$((DONE + 1))
            REMAINING=$((TOTAL - DONE))
            # Print a distinct progress indicator
            echo ">>> PROGRESS: $DONE / $TOTAL completed ($REMAINING left) <<<"
        else
            # Print the normal output from the worker
            echo "$line"
        fi
    done

echo "--------------------------------------------------"
echo "All tasks finished successfully!"
echo "Results in: $OUT_DIR"