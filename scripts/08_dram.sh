#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$DRAM_MODULE"
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
for table in 11 4; do
    shopt -s nullglob
    files=("$ANALYSIS/filtered_MAGs_table$table/"*.fa)
    (( ${#files[@]} > 0 )) || continue
    NATIVE="$ANALYSIS/DRAM_table$table"
    if [[ -f "$NATIVE/.done" ]]; then NATIVE=""; continue; fi
    [[ ! -e "$NATIVE" ]]
    DRAM.py annotate -i "$ANALYSIS/filtered_MAGs_table$table/*.fa" -o "$NATIVE" \
        --config_loc "$DRAM_CONFIG" --threads "$THREADS" --trans_table "$table"
    touch "$NATIVE/.done"
    python3 "$SCRIPT_ROOT/helpers.py" logs "$NATIVE" "$LOG_DIR/native" "$RUN"
    NATIVE=""
done
if [[ ! -s "$ANALYSIS/gene_to_MAG.tsv" ]]; then
    python3 "$SCRIPT_ROOT/collect_dram_annotations.py" --run-dir "$ANALYSIS" --output-dir "$ANALYSIS"
fi
