# MAG workflow

[中文说明](README.zh-CN.md)

A modular Slurm workflow from clean paired-end reads to metagenome-assembled genomes (MAGs), gene annotations, and sample-level abundance tables.

- **Joint coassembly:** assemble all samples together.
- **Grouped coassembly:** assemble predefined groups separately, refine bins within each assembly, and dereplicate the pooled MAGs.
- Both routes share downstream classification, DRAM annotation, read mapping, and MAG/gene/KO quantification.

![Workflow overview](assets/workflow.png)

[Editable SVG](assets/workflow.svg)

## Quick start

Configure your HPC environment once:

```bash
cp config.example.json config.json
# Edit module names, database paths, and resource settings in config.json.
```

Preview either route; previews do not submit analysis jobs:

```bash
bash scripts/run_joint.sh /path/to/clean_data
bash scripts/run_grouped.sh /path/to/grouped_clean_data
```

Add `--submit` to schedule the workflow:

```bash
bash scripts/run_joint.sh /path/to/clean_data --submit
```

For another dataset on the same cluster, normally only the input path changes.

## Input layout

Joint coassembly accepts paired reads in one directory or its subdirectories. Grouped coassembly uses the read directory structure as assembly groups:

```text
clean_data/
├── group_A/
│   ├── sample01_R1.fq.gz
│   └── sample01_R2.fq.gz
└── group_B/
    ├── sample02_R1.fq.gz
    └── sample02_R2.fq.gz
```

Supported paired names include `_R1/_R2`, `_R1_001/_R2_001`, `_good_1/_good_2`, and `_1/_2`, with `.fq` or `.fastq`, optionally gzip-compressed. Sample names must be unique across the dataset. Merge sequencing lanes per sample upstream.

For flat reads requiring explicit assembly groups, place a tab-delimited `samples.tsv` in the input directory, or pass its path as the input. Paths are relative to the sample sheet or absolute. See [the example](examples/samples.tsv). Grouping cannot be inferred biologically from filenames alone.

## Run in sections

```bash
bash scripts/run_joint.sh /path/to/clean_data --stage assembly --submit
# After successful assembly:
bash scripts/run_joint.sh /path/to/clean_data --stage binning --submit
```

| Stage | Work performed |
|---|---|
| `all` | Complete workflow; the default |
| `assembly` | MEGAHIT assembly |
| `assembly_continue` | Continue MEGAHIT in its original output directory |
| `binning` | Assembly indexing/mapping, MetaBAT2/SemiBin2/VAMB, DAS Tool, quality filtering, final catalogue |
| `quality` | Pool and namespace bins, CheckM2, quality filtering |
| `catalog` | Finalize the catalogue; grouped route also runs dRep |
| `annotate` | Submit GTDB-Tk and DRAM as separate jobs |
| `taxonomy` | GTDB-Tk only |
| `dram` | DRAM annotation and column-aligned annotation merging |
| `distill` | DRAM metabolic summaries |
| `quantify` | Index final MAGs, map all samples, produce abundance matrices |
| `tables` | Quantify from existing annotations and compatible BAMs |

Within a section, Slurm `afterok` dependencies connect jobs automatically. Between sections, wait for success or pass `--after JOB_ID`. Completed DRAM/distill subgroups with this workflow's `.done` markers are retained.

## Outputs

For input `/project/clean_data`, the default result directory is `/project/MAG_results/clean_data_joint/` or `clean_data_grouped/`. Override with `--run-dir`.

Main matrices in `tables/`:

- `sample_x_MAG_relative_within_catalogue.tsv`
- `sample_x_MAG_mean_coverage.tsv` and `sample_x_MAG_covered_fraction.tsv`
- `sample_x_gene_counts.tsv` and `sample_x_gene_TPM_within_catalogue.tsv`
- `sample_x_KO_counts_fractional.tsv` and `sample_x_KO_TPM_within_catalogue.tsv`
- `MAG_KO_TPM_long.tsv`: MAG-resolved functional contributions.

Feature-by-sample versions, gene coordinates, gene–KO–MAG links, quality information, and provenance are also retained. Logs are stored under `logs/<dataset>/<submission>/`; script/configuration snapshots are stored under `scripts/.runs/`.

Existing results can be supplied with `--assembly-root`, `--analysis-dir`, and `--bam-dir`. See the [Chinese guide](README.zh-CN.md#复用已有结果) for examples. BAMs must match the reference sequences and gene coordinates exactly; separate-assembly BAMs are not interchangeable with final-catalogue BAMs.

## Environment

Linux, Bash, Slurm, Environment Modules, and Python 3.8+ are required. The dispatcher and local tests use the Python standard library. DRAM supplies pandas for metabolic distillation.

Analysis tools: MEGAHIT, Bowtie2, samtools, MetaBAT2, SemiBin2, **VAMB 5.x**, DAS Tool, CheckM2, dRep, GTDB-Tk, **DRAM 1.5.x**, CoverM, and Subread/featureCounts.

The example contains module names from the source workflow, not a universal installation recipe. Adjust them and configure compatible CheckM2/GTDB/DRAM databases once per cluster. Leave partition, QOS, and account blank to use scheduler defaults. Resource requests are starting values, not measured requirements for every dataset.

## Analysis conventions

- Quality thresholds default to completeness ≥50% and contamination ≤10%.
- Grouped coassembly uses 95% ANI dereplication; joint coassembly retains the original non-dereplicated route.
- Final contigs are namespaced before annotation, preventing cross-assembly ID collisions. All samples are mapped to the same final catalogue.
- MAG relative abundance is normalized **within the retained catalogue**, not the entire community or soil mass.
- Gene counts are paired-fragment counts. Overlapping genes and multiple KO assignments receive fractional contributions; these counts are not necessarily integers.
- TPM uses all counted genes in the catalogue, including genes without KO annotations.
- DRAM translation tables 11/4 are handled separately and merged by column name.
- A private GTDB-Tk package copy applies the known `tostring()`/`tobytes()` compatibility fix. Shared installations are not edited.

Existing outputs are not automatically deleted or overwritten. Retry only failed sections after preserving incomplete outputs. A local test pass does not establish successful execution of biological tools on the target cluster.

## Development

```bash
python3 -m unittest discover -s tests -v
for script in scripts/*.sh; do bash -n "$script" || exit 1; done
```

GitHub Actions runs these lightweight tests without installing biological tools. Validation scope is recorded in [docs/VALIDATION.md](docs/VALIDATION.md).

No license has been selected yet.

## Tool documentation

[MEGAHIT](https://github.com/voutcn/megahit) · [SemiBin](https://github.com/BigDataBiology/SemiBin) · [VAMB](https://vamb.readthedocs.io/en/stable/how_to_run.html) · [DAS Tool](https://github.com/cmks/DAS_Tool) · [dRep](https://drep.readthedocs.io/en/latest/module_descriptions.html) · [DRAM](https://github.com/WrightonLabCSU/DRAM) · [CoverM](https://github.com/wwood/CoverM) · [featureCounts](https://subread.sourceforge.net/featureCounts.html)
