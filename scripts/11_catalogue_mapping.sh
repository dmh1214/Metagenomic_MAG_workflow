#!/usr/bin/env bash
source "${1:?}/common.sh"
all_reads
i=$((SLURM_ARRAY_TASK_ID-1))
load "$BOWTIE2_MODULE"; load "$SAMTOOLS_MODULE"
sample="${SAMPLES[$i]}"
mkdir -p "$CATALOGUE_BAMS"
bam="$CATALOGUE_BAMS/$sample.bam"
[[ ! -e "$bam" ]]
bowtie2 --very-sensitive -p "$((THREADS-4))" -x "$RUN/catalogue_index/reference" \
    -1 "${R1S[$i]}" -2 "${R2S[$i]}" 2> "$LOG_DIR/catalogue_${sample}.bowtie2.log" \
    | samtools sort -@ 2 -m 2G -T "$bam.tmp" -o "$bam" -
samtools index -@ 2 "$bam"
python3 "$SCRIPT_ROOT/read_recruitment.py" save-input \
    "$LOG_DIR/catalogue_${sample}.bowtie2.log" "$CATALOGUE_BAMS/$sample.input_reads.json"
