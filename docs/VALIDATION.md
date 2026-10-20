# Validation

Local tests: `python3 -m unittest discover -s tests -v`

Seven tests passed:

- Paired-file discovery, joint/grouped selection and missing-mate detection.
- Relative paths and generic groups in samples.tsv.
- Task-array sizes and ordered Slurm dependencies, including single-stage submission.
- Fake sbatch submission of the full 14-job/array plan, without running analysis.
- Namespacing colliding assembly contigs, CheckM2 selection, translation tables 11/4,
  column-aligned DRAM merging, fractional KO aggregation and sample × feature exports.
- Bowtie2 paired-end input parsing, read-end denominator conversion and mixed-input rejection.
- Simulated recruitment counts using saved input totals and the legacy FASTQ fallback,
  catalogue-contig selection, MAPQ20 output and standalone/automatic Slurm dependencies.

Python syntax and JSON configuration were parsed locally. The cluster modules, databases,
shell jobs, large-data resource requirements and full biological analyses have not been run
by this validation. Existing successful project scripts supplied the analysis commands.

The recruitment update passed all seven tests with simulated samtools counts; no server
BAMs were processed locally. Local Bash was unavailable, so shell syntax checks remain
in the GitHub Actions workflow (`bash -n scripts/*.sh`, applied per file).

No source reads, source assemblies or existing MAG results were modified.
