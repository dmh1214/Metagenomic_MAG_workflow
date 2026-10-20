"""DRAM coordinate extraction and abundance matrix conversion."""
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
csv.field_size_limit(10_000_000)

def rows(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")

def write_table(path, header, data):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(data)

def read_catalogue(run_dir):
    fasta_files = sorted((run_dir / "filtered_MAGs").glob("*.fa"))
    if not fasta_files:
        raise ValueError("No filtered_MAGs/*.fa found")
    contigs = {}
    for fasta in fasta_files:
        name, length = None, 0
        def save():
            if name is not None:
                if name in contigs:
                    raise ValueError(f"Contig belongs to multiple MAGs: {name}")
                contigs[name] = (fasta.stem, length)
        with open(fasta, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    save()
                    name, length = line[1:].split()[0], 0
                else:
                    length += len(line.strip())
            save()
    return fasta_files, contigs

def prepare(run_dir, out, fasta_files, contigs):
    write_table(out / "contig_to_MAG.tsv", ["MAG_ID", "contig", "length"],
                ((mag, contig, length) for contig, (mag, length) in contigs.items()))
    genes, copies, stats = set(), defaultdict(int), {}
    with open(out / "genes.saf", "w", encoding="utf-8", newline="") as saf, \
            open(out / "gene_annotation_compact.tsv", "w", encoding="utf-8", newline="") as compact:
        sw, cw = csv.writer(saf, delimiter="\t"), csv.writer(compact, delimiter="\t")
        sw.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        cw.writerow(["gene_id", "MAG_ID", "contig", "start", "end", "length", "KO", "translation_table"])
        for table in (11, 4):
            path = run_dir / f"DRAM_table{table}" / "annotations.tsv"
            if not path.exists():
                continue
            total, annotated = 0, 0
            for row in rows(path):
                gene = row.get("", row.get("gene_id"))
                if not gene or gene in genes:
                    raise ValueError(f"Missing/duplicate gene ID: {gene}")
                genes.add(gene)
                mag, contig = row["fasta"], row["scaffold"]
                if contig not in contigs or contigs[contig][0] != mag:
                    raise ValueError(f"Annotation/FASTA mismatch: {gene} {mag} {contig}")
                start, end = int(row["start_position"]), int(row["end_position"])
                if not 1 <= start <= end <= contigs[contig][1]:
                    raise ValueError(f"Invalid gene coordinates: {gene}")
                ko = sorted(set(re.findall(r"K\d{5}", row.get("ko_id", row.get("kegg_id", "")))))
                strand = "+" if row["strandedness"] in ("1", "+", "+1") else "-"
                sw.writerow([gene, contig, start, end, strand])
                cw.writerow([gene, mag, contig, start, end, end-start+1, ";".join(ko), table])
                for item in ko:
                    copies[(mag, item)] += 1
                total += 1
                annotated += bool(ko)
            stats[str(table)] = {"genes": total, "genes_with_KO": annotated}
    if not genes:
        raise ValueError("No annotation genes found")
    write_table(out / "MAG_KO_copy_number_long.tsv", ["MAG_ID", "KO", "gene_copies"],
                ((mag, ko, n) for (mag, ko), n in sorted(copies.items())))
    mags = [f.stem for f in fasta_files]
    kos = sorted({k for _, k in copies})
    write_table(out / "KO_by_MAG_copy_number.tsv", ["KO"] + mags,
                ([ko] + [copies.get((mag, ko), 0) for mag in mags] for ko in kos))
    (out / "annotation_QC.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    return stats

def parse_coverm(path, samples, out):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != 1 + 2*len(samples):
            raise ValueError(f"Unexpected CoverM columns: {header}")
        for i, sample in enumerate(samples):
            if sample not in header[1+2*i] or not header[1+2*i].lower().endswith("mean"):
                raise ValueError(f"Unexpected CoverM sample/method order: {header}")
            if "covered fraction" not in header[2+2*i].lower():
                raise ValueError(f"Unexpected CoverM covered-fraction column: {header}")
        coverage, breadth = {}, {}
        for row in reader:
            if row[0] == "unmapped":
                continue
            coverage[row[0]] = [float(row[1+2*i]) for i in range(len(samples))]
            breadth[row[0]] = [float(row[2+2*i]) for i in range(len(samples))]
    totals = [sum(v[i] for v in coverage.values()) for i in range(len(samples))]
    write_table(out / "MAG_mean_coverage.tsv", ["MAG_ID"] + samples,
                ([mag] + coverage[mag] for mag in sorted(coverage)))
    write_table(out / "MAG_covered_fraction.tsv", ["MAG_ID"] + samples,
                ([mag] + breadth[mag] for mag in sorted(breadth)))
    write_table(out / "MAG_relative_within_catalogue.tsv", ["MAG_ID"] + samples,
                ([mag] + [v/t if t else 0 for v, t in zip(coverage[mag], totals)] for mag in sorted(coverage)))

def count_rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
        header = next(reader)
        if header[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
            raise ValueError("Unexpected featureCounts format")
        yield header
        yield from reader

def transpose_numeric_table(source, target):
    """Export the user-facing sample x feature orientation as well as feature x sample."""
    with open(source, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        features, columns = [], [[] for _ in header[1:]]
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"Ragged matrix: {source}")
            features.append(row[0])
            for column, value in zip(columns, row[1:]):
                column.append(value)
    write_table(target, ["sample"]+features,
                ([sample]+values for sample, values in zip(header[1:], columns)))

def aggregate(out, samples):
    # Retain only compact MAG/KO joins; original huge Pfam columns never enter memory.
    mapping = {r["gene_id"]: (r["MAG_ID"], r["KO"].split(";") if r["KO"] else [])
               for r in rows(out / "gene_annotation_compact.tsv")}
    counts_path = out / "featureCounts.raw.tsv"
    iterator = count_rows(counts_path)
    header = next(iterator)
    observed = [Path(x).stem for x in header[6:]]
    if observed != samples:
        raise ValueError(f"featureCounts sample order mismatch: {observed} != {samples}")
    totals = [0.0] * len(samples)
    for row in iterator:
        length_kb = float(row[5])/1000
        for i, value in enumerate(row[6:]):
            totals[i] += float(value)/length_kb
    ko_counts, ko_tpm, mag_ko = {}, {}, {}
    annotated_tpm = [0.0]*len(samples)
    with open(out / "gene_counts.tsv", "w", encoding="utf-8", newline="") as cf, \
            open(out / "gene_TPM_within_catalogue.tsv", "w", encoding="utf-8", newline="") as tf:
        cw, tw = csv.writer(cf, delimiter="\t"), csv.writer(tf, delimiter="\t")
        cw.writerow(["gene_id"]+samples)
        tw.writerow(["gene_id"]+samples)
        iterator = count_rows(counts_path)
        next(iterator)
        seen = set()
        for row in iterator:
            gene = row[0]
            if gene in seen or gene not in mapping:
                raise ValueError(f"Unexpected/duplicate counted gene: {gene}")
            seen.add(gene)
            mag, kos = mapping[gene]
            values = list(map(float, row[6:]))
            tpm = [v/(float(row[5])/1000)/t*1e6 if t else 0 for v, t in zip(values, totals)]
            cw.writerow([gene]+values)
            tw.writerow([gene]+tpm)
            if kos:
                annotated_tpm = [a+b for a, b in zip(annotated_tpm, tpm)]
            # Equal split across distinct KO IDs to conserve gene abundance.
            for ko in kos:
                for target, key, contribution in ((ko_counts, ko, values), (ko_tpm, ko, tpm), (mag_ko, (mag, ko), tpm)):
                    accumulator = target.setdefault(key, [0.0]*len(samples))
                    for i, value in enumerate(contribution):
                        accumulator[i] += value/len(kos)
        if seen != set(mapping):
            raise ValueError("Not all annotation genes are present in featureCounts results")
    for name, values in (("KO_counts_fractional.tsv", ko_counts), ("KO_TPM_within_catalogue.tsv", ko_tpm)):
        write_table(out / name, ["KO"]+samples, ([ko]+values[ko] for ko in sorted(values)))
    write_table(out / "MAG_KO_TPM_long.tsv", ["MAG_ID", "KO", "sample", "TPM"],
                ((mag, ko, sample, value) for (mag, ko), values in sorted(mag_ko.items())
                 for sample, value in zip(samples, values)))
    write_table(out / "sample_function_QC.tsv", ["sample", "sum_gene_RPK", "KO_annotated_gene_TPM", "KO_annotated_fraction"],
                ((s, t, a, a/1e6) for s, t, a in zip(samples, totals, annotated_tpm)))
    for name in ("gene_counts", "gene_TPM_within_catalogue", "KO_counts_fractional", "KO_TPM_within_catalogue"):
        transpose_numeric_table(out / f"{name}.tsv", out / f"sample_x_{name}.tsv")
