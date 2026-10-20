import os
import runpy
import sys
import pandas as pd

# Keep numeric MAG IDs as strings; numeric score columns retain their own types.
entrypoint, arguments = sys.argv[1], sys.argv[2:]
paths = {os.path.realpath(arguments[i+1]) for i, value in enumerate(arguments)
         if value in ("-i", "--trna_path", "--rrna_path")}
read_csv = pd.read_csv


def string_mag_ids(path, *args, **kwargs):
    if isinstance(path, (str, os.PathLike)) and os.path.realpath(path) in paths:
        kwargs["dtype"] = {"fasta": str}
        kwargs["low_memory"] = False
    return read_csv(path, *args, **kwargs)


pd.read_csv = string_mag_ids
sys.argv = [entrypoint, "distill", *arguments]
runpy.run_path(entrypoint, run_name="__main__")
