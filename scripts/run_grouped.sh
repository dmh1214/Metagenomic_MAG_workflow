#!/usr/bin/env bash
set -e
input="${1:?Usage: bash scripts/run_grouped.sh /path/to/grouped_clean_data [options]}"
shift
python3 "$(dirname -- "${BASH_SOURCE[0]}")/submit.py" --input "$input" --strategy grouped "$@"
