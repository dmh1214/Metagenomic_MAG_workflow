#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$SAMTOOLS_MODULE"
python3 "$SCRIPT_ROOT/read_recruitment.py" count \
    --samples "$SCRIPT_ROOT/samples.tsv" --bam-dir "$CATALOGUE_BAMS" \
    --fasta-dir "$ANALYSIS/filtered_MAGs" --out "$RUN/read_recruitment" \
    --threads "$THREADS"
