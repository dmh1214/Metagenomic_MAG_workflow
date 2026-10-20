#!/usr/bin/env bash
source "${1:?}/common.sh"
load "$GTDB_MODULE"
export GTDBTK_DATA_PATH="$GTDB_DB"
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
NATIVE="$RUN/gtdb"
fresh "$NATIVE"
python "$SCRIPT_ROOT/helpers.py" patch-gtdb "$NATIVE/python_overlay"
export PYTHONPATH="$NATIVE/python_overlay"
python "$NATIVE/python_overlay/gtdbtk_entry.py" classify_wf --genome_dir "$ANALYSIS/filtered_MAGs" --extension fa \
    --out_dir "$NATIVE/taxonomy" --skip_ani_screen --cpus "$THREADS" --pplacer_cpus 1
python "$SCRIPT_ROOT/helpers.py" taxonomy "$NATIVE/taxonomy" "$ANALYSIS/all_MAG_taxonomy.tsv"
