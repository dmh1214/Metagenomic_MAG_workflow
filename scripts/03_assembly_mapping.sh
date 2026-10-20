#!/usr/bin/env bash
source "${1:?}/common.sh"
all_reads
task=$((SLURM_ARRAY_TASK_ID-1))
group_row "$((task / SAMPLE_COUNT + 1))"
i=$((task % SAMPLE_COUNT))
load "$BOWTIE2_MODULE"; load "$SAMTOOLS_MODULE"
sample="${SAMPLES[$i]}"
mkdir -p "$GROUP_OUT/bam"
bam="$GROUP_OUT/bam/$sample.bam"
[[ ! -e "$bam" ]]
bowtie2 --very-sensitive -p "$((THREADS-4))" -x "$GROUP_OUT/index/assembly" \
    -1 "${R1S[$i]}" -2 "${R2S[$i]}" 2> "$LOG_DIR/${GROUP}_${sample}.bowtie2.log" \
    | samtools sort -@ 2 -m 2G -T "$bam.tmp" -o "$bam" -
samtools index -@ 2 "$bam"
