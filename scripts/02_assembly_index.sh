#!/usr/bin/env bash
source "${1:?}/common.sh"
group_row
load "$BOWTIE2_MODULE"
fresh "$GROUP_OUT/index"
bowtie2-build --large-index --threads "$THREADS" "$ASSEMBLY" "$GROUP_OUT/index/assembly"
