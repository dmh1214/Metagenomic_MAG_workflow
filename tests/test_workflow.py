import csv
import contextlib
import io
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
import helpers
import submit
import table_functions as tables
import collect_dram_annotations as collector
import read_recruitment as recruitment


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT/"tests", prefix="fixture_")
        self.root = Path(self.temp.name)
        assert self.root.resolve().parent == (ROOT/"tests").resolve()

    def tearDown(self):
        self.temp.cleanup()

    def pair(self, folder, name, tag="_R", lane="", extension=".fq.gz"):
        folder.mkdir(parents=True, exist_ok=True)
        for end in (1, 2):
            (folder/f"{name}{tag}{end}{lane}{extension}").write_bytes(b"fixture")

    def test_discovery_and_strategies(self):
        reads = self.root/"reads"
        self.pair(reads/"Low", "L1", "_good_")
        self.pair(reads/"High", "H1", lane="_001", extension=".fastq.gz")
        rows, strategy = submit.prepare_samples(reads)
        self.assertEqual(strategy, "grouped")
        self.assertEqual({row[3] for row in rows}, {"Low", "High"})
        joint, strategy = submit.prepare_samples(reads, "joint")
        self.assertEqual(strategy, "joint")
        self.assertEqual({row[3] for row in joint}, {"all"})
        self.assertEqual({row[0] for row in rows}, {"H1", "L1"})
        (reads/"Low/L1_good_2.fq.gz").unlink()
        with self.assertRaises(ValueError):
            submit.prepare_samples(reads)

    def test_generic_manifest(self):
        self.pair(self.root, "sample", tag="_", extension=".fastq")
        helpers.write_table(self.root/"samples.tsv", ["sample", "r1", "r2", "group"],
                            [("sample", "sample_1.fastq", "sample_2.fastq", "treatment")])
        rows, strategy = submit.prepare_samples(self.root, "grouped")
        self.assertEqual(rows[0][3], "treatment")
        self.assertTrue(Path(rows[0][1]).is_absolute())

    def test_dependency_plan(self):
        plan = submit.create_plan("all", 3, 75, 4, 3)
        jobs = {job["name"]: job for job in plan}
        self.assertEqual(jobs["assembly_mapping"]["array"], "1-225%4")
        self.assertEqual(jobs["binning"]["array"], "1-9%3")
        self.assertEqual(jobs["catalogue_mapping"]["array"], "1-75%4")
        self.assertEqual(set(jobs["tables"]["dependencies"]), {"dram", "catalogue_mapping"})
        self.assertEqual(jobs["recruitment"]["dependencies"], ["catalogue_mapping"])
        names = set()
        for job in plan:
            self.assertTrue(set(job["dependencies"]) <= names)
            self.assertTrue((ROOT/"scripts"/job["script"]).exists())
            names.add(job["name"])
        self.assertEqual([j["name"] for j in submit.create_plan("dram", 1, 20, 4, 3)], ["dram"])
        quantify = submit.create_plan("quantify", 1, 20, 4, 3)
        self.assertNotIn("assembly", {j["name"] for j in quantify})
        self.assertIn("recruitment", {j["name"] for j in quantify})
        self.assertEqual([j["name"] for j in submit.create_plan("recruitment", 1, 20, 4, 3)], ["recruitment"])

    def test_submit_with_fake_slurm(self):
        reads = self.root/"reads"
        self.pair(reads, "sample")
        script_root = self.root/"package/scripts"
        script_root.mkdir(parents=True)
        for source in (ROOT/"scripts").iterdir():
            if source.suffix in (".sh", ".py", ".json"):
                (script_root/source.name).write_bytes(source.read_bytes())
        (script_root.parent/"config.json").write_bytes((ROOT/"config.example.json").read_bytes())
        fake_module = script_root/"submit.py"
        commands = []
        def fake_sbatch(command, text=True):
            commands.append(command)
            return str(7000+len(commands))+"\n"
        arguments = [str(fake_module), "--input", str(reads), "--stage", "all", "--submit"]
        with patch.object(submit, "__file__", str(fake_module)), patch.object(sys, "argv", arguments), \
                patch.object(submit.subprocess, "check_output", side_effect=fake_sbatch), contextlib.redirect_stdout(io.StringIO()):
            submit.main()
        self.assertEqual(len(commands), 14)
        self.assertFalse(any(arg in ("--partition=", "--qos=", "--account=") for command in commands for arg in command))
        dependency = [arg for arg in commands[-2] if arg.startswith("--dependency=")][0]
        self.assertIn("7012", dependency)
        self.assertIn("7009", dependency)
        self.assertIn("--dependency=afterok:7012", commands[-1])

    def test_recruitment_log_denominator(self):
        log = self.root/'mapping.log'
        log.write_text('1,000 reads; of these:\n  1,000 (100.00%) were paired; of these:\n')
        self.assertEqual(recruitment.log_reads(log), 2000)
        log.write_text('100 reads; of these:\n  90 (90.00%) were paired; of these:\n')
        with self.assertRaises(ValueError):
            recruitment.log_reads(log)

    def test_recruitment_with_and_without_metadata(self):
        fasta = self.root/'MAGs'
        fasta.mkdir()
        (fasta/'M.fa').write_text('>selected\nACGT\n')
        bams = self.root/'bam'
        bams.mkdir()
        for end in (1, 2):
            (self.root/f'R{end}.fq').write_text('@a\nAC\n+\nII\n@b\nAC\n+\nII\n')
        samples = self.root/'samples.tsv'
        helpers.write_table(samples, ['sample', 'r1', 'r2', 'group'],
                            [('S', self.root/'R1.fq', self.root/'R2.fq', 'treatment')])
        def fake(*args):
            if '-H' in args:
                return '@SQ\tSN:selected\tLN:4\n@SQ\tSN:unselected\tLN:10\n'
            if '--version' in args:
                return 'samtools mock\n'
            self.assertIn(2308, args)
            bed = Path(args[args.index('-L')+1])
            self.assertEqual(bed.read_text(), 'selected\t0\t4\n')
            return '2' if '-q' in args else '3'
        for mode in ('fastq', 'metadata'):
            if mode == 'metadata':
                (bams/'S.input_reads.json').write_text(json.dumps({'total_clean_reads': 4}))
                # The metadata path must not reread FASTQs.
                (self.root/'R1.fq').unlink()
                (self.root/'R2.fq').unlink()
            out = self.root/mode
            with patch.object(recruitment, 'command', side_effect=fake), contextlib.redirect_stdout(io.StringIO()):
                recruitment.recruit(samples, bams, fasta, out, 2)
            row = helpers.read_table(out/'sample_MAG_read_recruitment.tsv')[0]
            self.assertEqual(float(row['MAG_mapped_percent']), 75)
            self.assertEqual(float(row['MAG_mapped_percent_MAPQ20']), 50)
            self.assertEqual(row['denominator_source'], 'clean_FASTQ_records' if mode == 'fastq' else 'Bowtie2_input_summary')

    def test_catalogue_and_abundance(self):
        run = self.root/"run"
        run.mkdir()
        groups = self.root/"groups.txt"
        groups.write_text("Low\nHigh\n")
        for group in ("Low", "High"):
            folder = run/"groups"/group/"dastool/final_DASTool_bins"
            folder.mkdir(parents=True)
            (folder/"bin.1.fa").write_text(">k127_1\n"+"A"*1000+"\n")
        helpers.merge(run, groups)
        identities = [line.strip() for path in (run/"merged_MAGs").glob("*.fa") for line in path.read_text().splitlines() if line.startswith(">")]
        self.assertEqual(len(set(identities)), 2)
        (run/"checkm2").mkdir()
        helpers.write_table(run/"checkm2/quality_report.tsv", ["Name", "Completeness", "Contamination", "Translation_Table_Used"],
                            [("Low__bin.1", 95, 1, 11), ("High__bin.1", 80, 2, 4)])
        helpers.select(run, 50, 10)
        helpers.finalize(run, run/"passing_MAGs")
        analysis = run/"analysis"
        for table, mag in ((11, "Low__bin.1"), (4, "High__bin.1")):
            folder = analysis/f"DRAM_table{table}"
            folder.mkdir()
            columns = ["", "fasta", "scaffold", "start_position", "end_position", "strandedness", "ko_id", "heme_regulatory_motif_count"]
            row = [mag+"_gene1", mag, mag+"__c1", 1, 1000, 1, "K00001;K00002" if table == 11 else "", 7]
            if table == 11:
                columns.insert(-1, "cazy_ids")
                row.insert(-1, "GH1")
            helpers.write_table(folder/"annotations.tsv", columns, [row])
        collector.collect(analysis, analysis)
        merged = helpers.read_table(analysis/"all_gene_annotations.tsv")
        self.assertEqual(merged[1]["cazy_ids"], "")
        self.assertEqual(merged[1]["heme_regulatory_motif_count"], "7")
        output = run/"tables"
        helpers.tables_prepare(analysis, output)
        genes = helpers.read_table(output/"gene_annotation_compact.tsv")
        header = ["Geneid", "Chr", "Start", "End", "Strand", "Length", "/bam/S1.bam", "/bam/S2.bam"]
        helpers.write_table(output/"featureCounts.raw.tsv", header,
                            [[r["gene_id"], r["contig"], 1, 1000, "+", 1000, 10, 0] for r in genes])
        helpers.write_table(output/"coverm.raw.tsv", ["Genome", "S1.bam Mean", "S1.bam Covered Fraction", "S2.bam Mean", "S2.bam Covered Fraction"],
                            [("Low__bin.1", 2, 1, 0, 0), ("High__bin.1", 1, 0.8, 0, 0)])
        helpers.write_table(self.root/"samples.tsv", ["sample"], [("S1",), ("S2",)])
        helpers.tables_finish(output, self.root/"samples.tsv")
        ko = helpers.read_table(output/"KO_TPM_within_catalogue.tsv")
        self.assertAlmostEqual(sum(float(row["S1"]) for row in ko), 500000)
        sample_gene = helpers.read_table(output/"sample_x_gene_counts.tsv")
        self.assertEqual(float(sample_gene[0][genes[0]["gene_id"]]), 10)
        self.assertTrue((output/"MAG_KO_TPM_long.tsv").exists())


if __name__ == "__main__":
    unittest.main()
