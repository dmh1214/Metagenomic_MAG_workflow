#!/usr/bin/env bash
source "${1:?}/common.sh"
group_row
load "$MEGAHIT_MODULE"
r1=(); r2=()
while IFS=$'\t' read -r sample read1 read2 group; do
    [[ "$sample" == sample ]] && continue
    if [[ "$STRATEGY" == joint || "$group" == "$GROUP" ]]; then r1+=("$read1"); r2+=("$read2"); fi
done < "$SCRIPT_ROOT/samples.tsv"
NATIVE="$ASSEMBLY_ROOT/$GROUP"
if [[ "${2:-}" == continue ]]; then
    megahit --continue -o "$NATIVE"
else
    [[ ! -e "$NATIVE" ]]
    memory=$(( SLURM_MEM_PER_NODE * 1024 * 1024 * 85 / 100 ))
    megahit -1 "$(IFS=,; echo "${r1[*]}")" -2 "$(IFS=,; echo "${r2[*]}")" \
        -o "$NATIVE" -t "$THREADS" -m "$memory" --min-contig-len 1000 --presets meta-large
fi
