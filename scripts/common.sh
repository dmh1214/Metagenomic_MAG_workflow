#!/usr/bin/env bash
set -Eeo pipefail
SCRIPT_ROOT="${1:?Frozen script directory required}"
source "$SCRIPT_ROOT/run.env.sh"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_NUM_THREADS="$THREADS"
NATIVE=""
module purge
[[ -z "$PYTHON_MODULE" ]] || module load "$PYTHON_MODULE"
finish() {
    local status=$?
    [[ -z "$NATIVE" ]] || python3 "$SCRIPT_ROOT/helpers.py" logs "$NATIVE" "$LOG_DIR/native" "$RUN" || true
    return "$status"
}
trap finish EXIT
load() { [[ -z "$1" ]] || module load "$1"; }
fresh() { [[ ! -e "$1" ]] || { echo "Output exists: $1" >&2; exit 1; }; mkdir -p "$1"; }
group_row() {
    local index="${1:-${SLURM_ARRAY_TASK_ID:-1}}"
    GROUP=$(sed -n "${index}p" "$SCRIPT_ROOT/groups.txt")
    ASSEMBLY="$ASSEMBLY_ROOT/$GROUP/final.contigs.fa"
    GROUP_OUT="$RUN/groups/$GROUP"
}
all_reads() {
    SAMPLES=(); R1S=(); R2S=()
    while IFS=$'\t' read -r sample r1 r2 group; do
        [[ "$sample" == sample ]] && continue
        SAMPLES+=("$sample"); R1S+=("$r1"); R2S+=("$r2")
    done < "$SCRIPT_ROOT/samples.tsv"
}
group_bams() {
    BAMS=()
    for sample in "${SAMPLES[@]}"; do BAMS+=("$GROUP_OUT/bam/$sample.bam"); done
}
