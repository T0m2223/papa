#!/bin/bash
#SBATCH --job-name=gpsa-threads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=gpsa-threads-%j.out
#
# Strong-scaling sweep for the report graph (speedup vs number of threads).
#
# Run the whole script inside ONE Slurm step, e.g.:
#   srun -N1 -n1 -c 32 ./thread_scaling.sh X.txt Y.txt
# or submit it:
#   sbatch thread_scaling.sh
# or inside an allocation:
#   salloc -N1 -n1 -c 32
#   ./thread_scaling.sh X.txt Y.txt
#
# The binary must already exist as ./gpsa.
#
# Defaults: default dataset X.txt/Y.txt, best granularity 65536 (256x256),
# threads 1 2 4 8 16 32, best of 3 repeats.
#
# Override examples:
#   THREADS="2 4 8 16 32" REPEATS=5 ./thread_scaling.sh
#   GRAIN=65536 BLOCK=256 ./thread_scaling.sh X2.txt Y2.txt

set -u

BIN=./gpsa
XF="${1:-X.txt}"
YF="${2:-Y.txt}"
THREADS="${THREADS:-1 2 4 8 16 32}"
REPEATS="${REPEATS:-3}"
GRAIN="${GRAIN:-65536}"
BLOCK="${BLOCK:-256}"

export OMP_PROC_BIND=close
export OMP_PLACES=cores

if [[ ! -x "$BIN" ]]; then
    echo "Error: $BIN not found/executable. Compile first." >&2
    exit 1
fi
if [[ ! -f "$XF" || ! -f "$YF" ]]; then
    echo "Error: missing input file(s): $XF / $YF" >&2
    exit 1
fi

time_of() {
    grep -oE 'completed in [0-9.eE+-]+ seconds' "$1" | head -1 | awk '{print $3}'
}
lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a < b)}'; }

echo "host     : $(hostname)"
echo "binary   : $BIN"
echo "dataset  : $XF x $YF"
echo "grain    : $GRAIN (block ${BLOCK}x${BLOCK})"
echo "threads  : $THREADS"
echo "repeats  : $REPEATS (best of)"
echo

# ---------- sequential baseline ----------
best_seq=""
for ((r = 1; r <= REPEATS; r++)); do
    "$BIN" --x "$XF" --y "$YF" --exec-mode 1 --save-to aligned-sequential.txt \
        > "seq-$r.log" 2>&1
    t=$(time_of "seq-$r.log")
    [[ -z "$t" ]] && { echo "Sequential run failed; see seq-$r.log" >&2; exit 1; }
    if [[ -z "$best_seq" ]] || lt "$t" "$best_seq"; then best_seq="$t"; fi
done

# sequence lengths (inner-loop cells = N*M)
read -r N M < <(grep -oE 'sizes [0-9]+ and [0-9]+' seq-1.log | awk '{print $2, $4}')
if [[ -z "${N:-}" || -z "${M:-}" ]]; then
    echo "Could not parse sequence sizes from seq-1.log" >&2
    exit 1
fi
CELLS=$(( N * M ))
GFLOP=$(awk -v c="$CELLS" 'BEGIN{printf "%.9f", 5*c/1e9}')

echo "sequential baseline: ${best_seq}s  (${N} x ${M} cells, ${GFLOP} GFLOP)"
echo
echo "CSV (paste into the report charts):"
printf "version,threads,time_s,speedup,gflops\n"

# ---------- thread sweep ----------
for T in $THREADS; do
    for v in taskloop tasks; do
        if [[ "$v" == "taskloop" ]]; then mode=2; else mode=3; fi

        best=""
        for ((r = 1; r <= REPEATS; r++)); do
            log="run-${v}-t${T}-${r}.log"
            OMP_NUM_THREADS="$T" "$BIN" --x "$XF" --y "$YF" --exec-mode "$mode" \
                --grain-size "$GRAIN" --block-size-x "$BLOCK" --block-size-y "$BLOCK" \
                > "$log" 2>&1
            t=$(time_of "$log")
            [[ -z "$t" ]] && { echo "  $v T=$T run $r failed (see $log)" >&2; continue; }
            if [[ -z "$best" ]] || lt "$t" "$best"; then best="$t"; fi
        done

        if [[ -n "$best" ]]; then
            read -r sp gl < <(awk -v s="$best_seq" -v t="$best" -v g="$GFLOP" \
                'BEGIN{printf "%.3f %.3f", s/t, g/t}')
            printf "%s,%s,%s,%s,%s\n" "$v" "$T" "$best" "$sp" "$gl"
        else
            printf "%s,%s,FAIL,,,\n" "$v" "$T"
        fi
    done
done

echo
echo "Logs: seq-*.log and run-*.log"
