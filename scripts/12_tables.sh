#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$SAMTOOLS_MODULE"; load "$COVERM_MODULE"; load "$SUBREAD_MODULE"
python3 "$SCRIPT_ROOT/helpers.py" tables-prepare "$ANALYSIS" "$TABLES"
all_reads
bams=()
for sample in "${SAMPLES[@]}"; do bams+=("$CATALOGUE_BAMS/$sample.bam"); done
coverm genome --bam-files "${bams[@]}" --genome-fasta-files "$ANALYSIS/filtered_MAGs/"*.fa \
    --methods mean covered_fraction --min-covered-fraction 0.10 \
    --min-read-percent-identity 95 --min-read-aligned-percent 75 --contig-end-exclusion 0 \
    --threads "$THREADS" --output-format dense --output-file "$TABLES/coverm.raw.tsv"
featureCounts -T "$THREADS" -F SAF -a "$TABLES/genes.saf" -o "$TABLES/featureCounts.raw.tsv" \
    -p --countReadPairs -B -C -s 0 -Q 20 --primary -O --fraction "${bams[@]}"
python3 "$SCRIPT_ROOT/helpers.py" tables-finish "$TABLES" "$SCRIPT_ROOT/samples.tsv"
