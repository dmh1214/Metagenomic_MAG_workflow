#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$DRAM_MODULE"
NATIVE="$RUN/distill"
mkdir -p "$NATIVE"
for table in 11 4; do
    input="$ANALYSIS/DRAM_table$table"
    [[ -s "$input/annotations.tsv" ]] || continue
    output="$RUN/distill/table$table"
    [[ ! -f "$output/.done" ]] || continue
    [[ ! -e "$output" ]]
    args=(-i "$input/annotations.tsv" -o "$output" --config_loc "$DRAM_CONFIG")
    [[ ! -s "$input/trnas.tsv" ]] || args+=(--trna_path "$input/trnas.tsv")
    [[ ! -s "$input/rrnas.tsv" ]] || args+=(--rrna_path "$input/rrnas.tsv")
    python3 "$SCRIPT_ROOT/distill.py" "$(command -v DRAM.py)" "${args[@]}"
    touch "$output/.done"
done
