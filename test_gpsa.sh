#!/bin/bash
#SBATCH --job-name=gpsa-tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=gpsa-tests-%j.out
#
# Run the whole script inside one Slurm step, e.g.:
#   srun -N1 -n1 -c 32 ./test_gpsa.sh
# or submit it:
#   sbatch test_gpsa.sh
# or use an allocation:
#   salloc -N1 -n1 -c 32
#   ./test_gpsa.sh
#
# The script then calls ./gpsa directly; the gpsa runs inherit the step's
# CPUs, and OMP_NUM_THREADS controls their OpenMP threads.
#
# Optional args: ./test_gpsa.sh <x-file> <y-file>
# The binary must already exist as ./gpsa in the current directory.

set -u

BIN=./gpsa
XF="${1:-X3.txt}"
YF="${2:-Y3.txt}"
CORES="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-32}}"
REPEATS="${REPEATS:-3}"

export OMP_NUM_THREADS="$CORES"
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

# print runtime in seconds from a gpsa log
time_of() {
    grep -oE 'completed in [0-9.eE+-]+ seconds' "$1" | head -1 | awk '{print $3}'
}

# true if $1 < $2 (floats)
lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a < b)}'; }

echo "=============================================="
echo " gpsa tests"
echo " host      : $(hostname)"
echo " binary    : $BIN"
echo " cores     : $CORES (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo " dataset   : $XF x $YF"
echo " repeats   : $REPEATS (best of)"
echo "=============================================="

# ---------------- sequential baseline ----------------
best_seq=""
for ((r = 1; r <= REPEATS; r++)); do
    "$BIN" --x "$XF" --y "$YF" --exec-mode 1 --save-to aligned-sequential.txt \
        > "seq-$r.log" 2>&1
    t=$(time_of "seq-$r.log")
    if [[ -z "$t" ]]; then
        echo "Sequential run failed; see seq-$r.log" >&2
        exit 1
    fi
    if [[ -z "$best_seq" ]] || lt "$t" "$best_seq"; then
        best_seq="$t"
    fi
done

echo "sequential baseline: ${best_seq}s"
echo

printf "%-10s %-10s %-12s %-10s %-10s %-8s\n" \
    "version" "grain" "block" "time(s)" "speedup" "check"
printf "%-10s %-10s %-12s %-10s %-10s %-8s\n" \
    "--------" "--------" "--------" "--------" "--------" "--------"

# The three granularities required by the assignment.
# 4096 (64x64), 65536 (256x256), 16777216 (4096x4096).
for spec in "4096 64" "65536 256" "16777216 4096"; do
    set -- $spec
    G=$1
    B=$2

    for v in taskloop tasks; do
        if [[ "$v" == "taskloop" ]]; then
            mode=2
            afile=aligned-taskloop.txt
        else
            mode=3
            afile=aligned-tasks.txt
        fi

        best=""
        check="FAIL"
        for ((r = 1; r <= REPEATS; r++)); do
            log="run-${v}-${G}-${r}.log"
            "$BIN" --x "$XF" --y "$YF" --exec-mode "$mode" \
                --grain-size "$G" --block-size-x "$B" --block-size-y "$B" \
                > "$log" 2>&1
            t=$(time_of "$log")
            if [[ -z "$t" ]]; then
                echo "  $v grain=$G run $r failed (see $log)" >&2
                continue
            fi
            if [[ -z "$best" ]] || lt "$t" "$best"; then
                best="$t"
            fi
            if cmp -s aligned-sequential.txt "$afile"; then
                check="OK"
            else
                check="FAIL"
            fi
        done

        if [[ -n "$best" ]]; then
            speedup=$(awk -v s="$best_seq" -v t="$best" \
                'BEGIN{ printf "%.2fx", s/t }')
        else
            speedup="?"
        fi

        printf "%-10s %-10s %-12s %-10s %-10s %-8s\n" \
            "$v" "$G" "${B}x${B}" "${best:-FAIL}" "$speedup" "$check"
    done
done

echo
echo "done. Logs: seq-*.log and run-*.log"
