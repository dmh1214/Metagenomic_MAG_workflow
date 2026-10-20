#!/usr/bin/env bash
source "${1:?}/common.sh"
if [[ "$STRATEGY" == grouped ]]; then
    module purge
    load "$DREP_MODULE"
    NATIVE="$RUN/drep"
    [[ ! -e "$NATIVE" ]]
    dRep dereplicate "$NATIVE" -g "$RUN/passing_MAGs/"*.fa --genomeInfo "$RUN/genome_info.csv" \
        -pa 0.90 -sa "$ANI" --S_algorithm fastANI -nc 0.10 -l 50000 \
        -comp "$MIN_COMPLETENESS" -con "$MAX_CONTAMINATION" -p "$THREADS"
    chosen="$NATIVE/dereplicated_genomes"
else
    chosen="$RUN/passing_MAGs"
fi
load "$PYTHON_MODULE"
python3 "$SCRIPT_ROOT/helpers.py" finalize "$RUN" "$chosen"
