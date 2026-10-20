#!/usr/bin/env bash
source "${1:?}/common.sh"
group_row
python3 "$SCRIPT_ROOT/helpers.py" das-input "$ASSEMBLY" "$GROUP_OUT"
set +u
source "$CONDA_SH"
conda activate "$DASTOOL_ENV"
NATIVE="$GROUP_OUT/dastool"
fresh "$NATIVE"
DAS_Tool -i "$GROUP_OUT/das_input/metabat2.tsv,$GROUP_OUT/das_input/semibin2.tsv,$GROUP_OUT/das_input/vamb.tsv" \
    -l metabat2,semibin2,vamb -c "$ASSEMBLY" -o "$NATIVE/final" \
    --search_engine diamond --score_threshold 0.5 --write_bins --threads "$THREADS"
