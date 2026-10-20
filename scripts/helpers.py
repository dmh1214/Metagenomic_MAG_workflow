#!/usr/bin/env python3
"""File conversion helpers; analysis commands are in the numbered shell scripts."""
import csv
import gzip
import os
import re
import shutil
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)


def read_table(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, columns, rows):
    with open(path, "x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def fastas(folder):
    paths = sorted(p for p in Path(folder).iterdir() if p.is_file() and
                   re.search(r"\.(fa|fna|fasta)(\.gz)?$", p.name) and
                   not any(x in p.name.lower() for x in ("unbinned", "unassigned", "too_short", "lowdepth")))
    if not paths:
        raise ValueError(f"No bin FASTAs: {folder}")
    return paths


def fasta_lines(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        yield from handle


def link(source, destination):
    if os.name == "nt":
        shutil.copyfile(source, destination)
    else:
        destination.symlink_to(source.resolve())


def das_input(assembly, folder):
    folder = Path(folder)
    destination = folder / "das_input"
    destination.mkdir()
    for label, directory in (("metabat2", folder/"metabat2/bins"),
                             ("semibin2", folder/"semibin2/output_bins"), ("vamb", folder/"vamb/bins")):
        with open(destination/f"{label}.tsv", "x") as handle:
            for path in fastas(directory):
                name = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", path.name)
                for line in fasta_lines(path):
                    if line.startswith(">"):
                        handle.write(f"{line[1:].split()[0]}\t{label}_{name}\n")


def merge(run, groups_file):
    run = Path(run)
    output = run / "merged_MAGs"
    output.mkdir()
    provenance = []
    with open(run/"contig_provenance.tsv", "x") as mapping:
        mapping.write("contig\tMAG_ID\tassembly_group\toriginal_contig\n")
        for group in Path(groups_file).read_text().splitlines():
            for path in fastas(run/"groups"/group/"dastool/final_DASTool_bins"):
                mag = group + "__" + re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", path.name)
                number = 0
                with open(output/(mag+".fa"), "x") as handle:
                    for line in fasta_lines(path):
                        if line.startswith(">"):
                            number += 1
                            contig = mag + "__c" + str(number)
                            mapping.write(f"{contig}\t{mag}\t{group}\t{line[1:].split()[0]}\n")
                            handle.write(">"+contig+"\n")
                        else:
                            handle.write(line)
                provenance.append((mag, group, str(path)))
    write_table(run/"MAG_provenance.tsv", ["MAG_ID", "assembly_group", "source_fasta"], provenance)


def quality_name(value):
    return re.sub(r"\.(fa|fna|fasta)$", "", value)


def select(run, min_comp, max_cont):
    run = Path(run)
    selected = run/"passing_MAGs"
    selected.mkdir()
    quality = []
    for row in read_table(run/"checkm2/quality_report.tsv"):
        name = quality_name(row["Name"])
        comp, cont = float(row["Completeness"]), float(row["Contamination"])
        if comp >= float(min_comp) and cont <= float(max_cont):
            table = int(float(row["Translation_Table_Used"]))
            if table not in (4, 11):
                raise ValueError(f"Unsupported genetic code: {name}: {table}")
            link(run/"merged_MAGs"/(name+".fa"), selected/(name+".fa"))
            quality.append((name+".fa", comp, cont))
    if not quality:
        raise ValueError("No MAGs passed quality filtering")
    with open(run/"genome_info.csv", "x", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genome", "completeness", "contamination"])
        writer.writerows(quality)


def finalize(run, chosen):
    run = Path(run)
    analysis = run/"analysis"
    analysis.mkdir()
    for name in ("filtered_MAGs", "filtered_MAGs_table11", "filtered_MAGs_table4"):
        (analysis/name).mkdir()
    quality = {quality_name(row["Name"]): row for row in read_table(run/"checkm2/quality_report.tsv")}
    selected = []
    for path in fastas(chosen):
        name = path.stem
        table = int(float(quality[name]["Translation_Table_Used"]))
        link(path, analysis/"filtered_MAGs"/path.name)
        link(path, analysis/f"filtered_MAGs_table{table}"/path.name)
        selected.append((name, table, quality[name]["Completeness"], quality[name]["Contamination"]))
    write_table(analysis/"selected_MAGs.tsv", ["MAG_ID", "translation_table", "Completeness", "Contamination"], selected)
    shutil.copy2(run/"checkm2/quality_report.tsv", analysis/"CheckM2_quality_report.tsv")
    shutil.copy2(run/"MAG_provenance.tsv", analysis/"MAG_provenance.tsv")


def patch_gtdb(overlay):
    import gtdbtk
    overlay = Path(overlay)
    source = Path(gtdbtk.__file__).resolve().parent
    shutil.copytree(source, overlay/"gtdbtk", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    markers = overlay/"gtdbtk/markers.py"
    text = markers.read_text()
    markers.write_text(text.replace("list_masked_seq.tostring()", "list_masked_seq.tobytes()"))
    shutil.copy2(shutil.which("gtdbtk"), overlay/"gtdbtk_entry.py")


def taxonomy(folder, output):
    records = []
    for path in sorted(Path(folder).glob("gtdbtk.*.summary.tsv")):
        for row in read_table(path):
            records.append((row["user_genome"], row["classification"], row.get("classification_method", ""), row.get("warnings", "")))
    write_table(output, ["MAG_ID", "classification", "classification_method", "warnings"], records)


def reference(folder, output):
    with open(output, "x") as handle:
        for path in fastas(folder):
            for line in fasta_lines(path):
                handle.write(line)


def tables_prepare(analysis, output):
    import table_functions as tables
    analysis, output = Path(analysis), Path(output)
    output.mkdir()
    files, contigs = tables.read_catalogue(analysis)
    tables.prepare(analysis, output, files, contigs)
    for name in ("all_MAG_taxonomy.tsv", "CheckM2_quality_report.tsv", "selected_MAGs.tsv", "MAG_provenance.tsv"):
        source = analysis/name
        if source.exists():
            shutil.copy2(source, output/name)


def tables_finish(output, manifest):
    import table_functions as tables
    output = Path(output)
    samples = [row["sample"] for row in read_table(manifest)]
    tables.parse_coverm(output/"coverm.raw.tsv", samples, output)
    for name in ("MAG_mean_coverage", "MAG_covered_fraction", "MAG_relative_within_catalogue"):
        tables.transpose_numeric_table(output/f"{name}.tsv", output/f"sample_x_{name}.tsv")
    tables.aggregate(output, samples)
    print("Finished abundance tables:", output)


def logs(source, destination, run):
    source, destination, run = Path(source), Path(destination), Path(run)
    if not source.exists():
        return
    for path in sorted(source.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(run)
        if path.suffix in (".log", ".out", ".err") or path.name == "log.txt" or any(p in ("log", "logs") for p in path.relative_to(source).parts[:-1]):
            target = destination/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), target)
            path.symlink_to(target.resolve())


if __name__ == "__main__":
    action = sys.argv[1]
    actions = {"das-input": das_input, "merge": merge, "select": select, "finalize": finalize,
               "patch-gtdb": patch_gtdb, "taxonomy": taxonomy, "reference": reference,
               "tables-prepare": tables_prepare, "tables-finish": tables_finish, "logs": logs}
    actions[action](*sys.argv[2:])
