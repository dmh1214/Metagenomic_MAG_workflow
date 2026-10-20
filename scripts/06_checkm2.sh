#!/usr/bin/env bash
source "${1:?}/common.sh"
python3 "$SCRIPT_ROOT/helpers.py" merge "$RUN" "$SCRIPT_ROOT/groups.txt"
load "$CHECKM2_MODULE"
NATIVE="$RUN/checkm2"
[[ ! -e "$NATIVE" ]]
checkm2 predict --input "$RUN/merged_MAGs" -x fa --output-directory "$NATIVE" \
    --threads "$THREADS" --database_path "$CHECKM2_DB"
python3 "$SCRIPT_ROOT/helpers.py" select "$RUN" "$MIN_COMPLETENESS" "$MAX_CONTAMINATION"
