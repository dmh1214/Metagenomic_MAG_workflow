#!/usr/bin/env python3
"""Merge DRAM translation-table outputs by column name, preserving missing data.

Never reruns annotation and refuses to overwrite existing merged files.
"""
import argparse
import csv
import json
import tempfile
from pathlib import Path

csv.field_size_limit(10_000_000)


def collect(run_dir, output_dir):
    paths = [run_dir / f"DRAM_table{table}" / "annotations.tsv" for table in (11, 4)]
    paths = [path for path in paths if path.is_file() and path.stat().st_size]
    if not paths:
        raise ValueError("No DRAM annotations.tsv found")
    output_dir.mkdir(parents=True, exist_ok=True)
    names = ["all_gene_annotations.tsv", "gene_to_MAG.tsv", "annotation_merge_QC.json"]
    for name in names:
        if (output_dir / name).exists():
            raise ValueError(f"Output exists; select another --output-dir: {output_dir / name}")
    fields, headers = ["gene_id"], {}
    for path in paths:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"))
        if len(set(header)) != len(header) or "fasta" not in header:
            raise ValueError(f"Invalid DRAM header: {path}")
        if header[0] not in ("", "gene_id", "Unnamed: 0"):
            raise ValueError(f"Unexpected first gene-ID column: {path}: {header[0]}")
        header[0] = "gene_id"
        headers[path] = header
        for field in header:
            if field not in fields:
                fields.append(field)
    seen, stats = set(), {}
    # All parsing completes before final names are published.
    with tempfile.TemporaryDirectory(prefix="annotation_merge_", dir=output_dir) as temp:
        temporary = Path(temp)
        with open(temporary / names[0], "w", encoding="utf-8", newline="") as annotated, \
                open(temporary / names[1], "w", encoding="utf-8", newline="") as mapping:
            aw = csv.DictWriter(annotated, fieldnames=fields, delimiter="\t", lineterminator="\n")
            mw = csv.writer(mapping, delimiter="\t", lineterminator="\n")
            aw.writeheader()
            mw.writerow(["gene_id", "MAG_ID"])
            for path in paths:
                count = 0
                with open(path, encoding="utf-8-sig", newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    next(reader)
                    for values in reader:
                        if len(values) != len(headers[path]):
                            raise ValueError(f"Ragged annotation row: {path}, row {count+2}")
                        row = dict(zip(headers[path], values))
                        gene = row["gene_id"]
                        if not gene or gene in seen or not row["fasta"]:
                            raise ValueError(f"Missing/duplicate gene or missing MAG: {gene}")
                        seen.add(gene)
                        aw.writerow(row)
                        mw.writerow([gene, row["fasta"]])
                        count += 1
                stats[str(path)] = count
        (temporary / names[2]).write_text(json.dumps({"source_rows": stats, "total_genes": len(seen),
            "columns": fields, "missing_fields": "empty, not zero; no database hit inferred"}, indent=2), encoding="utf-8")
        for name in names:
            (temporary / name).rename(output_dir / name)
    print(f"Merged {len(seen)} genes by column name: {output_dir}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    collect(args.run_dir, args.output_dir)
