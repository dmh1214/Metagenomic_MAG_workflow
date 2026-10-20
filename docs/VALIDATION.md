# Validation

Local tests: `python3 -m unittest discover -s tests -v`

Five tests passed:

- Paired-file discovery, joint/grouped selection and missing-mate detection.
- Relative paths and generic groups in samples.tsv.
- Task-array sizes and ordered Slurm dependencies, including single-stage submission.
- Fake sbatch submission of the full 13-job/array plan, without running analysis.
- Namespacing colliding assembly contigs, CheckM2 selection, translation tables 11/4,
  column-aligned DRAM merging, fractional KO aggregation and sample × feature exports.

Python syntax and JSON configuration were parsed locally. The cluster modules, databases,
shell jobs, large-data resource requirements and full biological analyses have not been run
by this validation. Existing successful project scripts supplied the analysis commands.

No source reads, source assemblies or existing MAG results were modified.
