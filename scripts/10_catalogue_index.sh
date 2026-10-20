#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$BOWTIE2_MODULE"
fresh "$RUN/catalogue_index"
python3 "$SCRIPT_ROOT/helpers.py" reference "$ANALYSIS/filtered_MAGs" "$RUN/catalogue_index/reference.fa"
bowtie2-build --large-index --threads "$THREADS" "$RUN/catalogue_index/reference.fa" "$RUN/catalogue_index/reference"
