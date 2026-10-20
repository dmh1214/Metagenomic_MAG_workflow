#!/usr/bin/env bash
source "${1:?}/common.sh"
task=$((SLURM_ARRAY_TASK_ID-1))
group_row "$((task / 3 + 1))"
all_reads; group_bams
case "$((task % 3))" in
0)
    load "$METABAT2_MODULE"
    NATIVE="$GROUP_OUT/metabat2"
    fresh "$NATIVE"
    mkdir "$NATIVE/bins"
    jgi_summarize_bam_contig_depths --outputDepth "$NATIVE/depth.tsv" --minContigLength 1000 "${BAMS[@]}"
    metabat2 -i "$ASSEMBLY" -a "$NATIVE/depth.tsv" -o "$NATIVE/bins/bin" -m 2500 -t "$THREADS" --seed 42
    ;;
1)
    load "$SEMIBIN2_MODULE"
    NATIVE="$GROUP_OUT/semibin2"
    [[ ! -e "$NATIVE" ]]
    SemiBin2 single_easy_bin -i "$ASSEMBLY" -b "${BAMS[@]}" -o "$NATIVE" \
        --self-supervised --sequencing-type short_read --engine cpu --threads "$THREADS" \
        --minfasta-kbs 200 --compression none
    ;;
2)
    load "$VAMB_MODULE"
    NATIVE="$GROUP_OUT/vamb"
    [[ ! -e "$NATIVE" ]]
    vamb bin default --outdir "$NATIVE" --fasta "$ASSEMBLY" \
        --bamdir "$GROUP_OUT/bam" --minfasta 200000 -p "$THREADS" --seed 42 -o
    ;;
esac
