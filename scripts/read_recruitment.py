#!/usr/bin/env python3
"""Count read ends recruited to the final MAG catalogue."""
import argparse
import csv
import gzip
import json
import re
import subprocess
from pathlib import Path

from helpers import fastas, fasta_lines


def log_reads(path):
    text = Path(path).read_text()
    totals = re.findall(r'^\s*([\d,]+) reads; of these:', text, re.M)
    paired = re.findall(r'^\s*([\d,]+) \([\d.]+%\) were paired; of these:', text, re.M)
    if len(totals) != 1 or len(paired) != 1:
        raise ValueError(f'Expected one paired-end Bowtie2 summary: {path}')
    total, pairs = (int(v.replace(',', '')) for v in (totals[0], paired[0]))
    if total != pairs or pairs <= 0:
        raise ValueError(f'Expected exclusively paired-end input: {path}')
    return 2 * pairs


def fastq_reads(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        lines = sum(1 for _ in f)
    if not lines or lines % 4:
        raise ValueError(f'Expected standard four-line FASTQ records: {path}')
    return lines // 4


def catalogue_contigs(folder):
    contigs = {}
    for path in fastas(folder):
        name = None
        for line in fasta_lines(path):
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in contigs:
                    raise ValueError(f'Duplicate catalogue contig: {name}')
                contigs[name] = 0
            elif name is not None:
                contigs[name] += len(line.strip())
    return contigs


def command(*args):
    return subprocess.check_output(list(map(str, args)), text=True)


def recruit(samples, bam_dir, fasta_dir, out, threads):
    out.mkdir(parents=True, exist_ok=False)
    contigs = catalogue_contigs(fasta_dir)
    bed = out / 'catalogue_contigs.bed'
    bed.write_text(''.join(f'{c}\t0\t{n}\n' for c, n in contigs.items()))
    with samples.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    with (out / 'sample_MAG_read_recruitment.tsv').open('x', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'group', 'total_clean_reads', 'MAG_mapped_reads',
                         'MAG_mapped_percent', 'MAG_mapped_reads_MAPQ20',
                         'MAG_mapped_percent_MAPQ20', 'denominator_source'])
        for row in rows:
            sample = row['sample']
            metadata = bam_dir / (sample + '.input_reads.json')
            if metadata.exists():
                info = json.loads(metadata.read_text())
                total, source = info['total_clean_reads'], 'Bowtie2_input_summary'
                if not isinstance(total, int) or total <= 0:
                    raise ValueError(f'Invalid input read count: {metadata}')
            else:
                # Existing BAMs may omit unmapped reads. Count input FASTQs only then.
                r1, r2 = fastq_reads(row['r1']), fastq_reads(row['r2'])
                if r1 != r2:
                    raise ValueError(f'Unequal paired FASTQ counts: {sample}')
                total, source = r1 + r2, 'clean_FASTQ_records'
            bam = bam_dir / (sample + '.bam')
            lengths = {}
            for line in command('samtools', 'view', '-H', bam).splitlines():
                if line.startswith('@SQ\t'):
                    fields = dict(x.split(':', 1) for x in line.split('\t')[1:])
                    lengths[fields['SN']] = int(fields['LN'])
            if any(lengths.get(c) != n for c, n in contigs.items()):
                raise ValueError(f'Catalogue/BAM reference mismatch: {bam}')
            # 2308 excludes unmapped, secondary and supplementary records.
            args = ['samtools', 'view', '-@', threads, '-c', '-F', 2308, '-M', '-L', bed]
            mapped = int(command(*args, bam))
            highq = int(command(*args, '-q', 20, bam))
            if not 0 <= highq <= mapped <= total:
                raise ValueError(f'Invalid recruitment counts: {sample}')
            writer.writerow([sample, row.get('group', ''), total, mapped, 100*mapped/total,
                             highq, 100*highq/total, source])
            f.flush()
            print(f'{sample}: {100*mapped/total:.2f}%; MAPQ20 {100*highq/total:.2f}%', flush=True)
    (out / 'samtools_version.txt').write_text(command('samtools', '--version'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    save = sub.add_parser('save-input')
    save.add_argument('log', type=Path)
    save.add_argument('output', type=Path)
    count = sub.add_parser('count')
    for name in ('samples', 'bam-dir', 'fasta-dir', 'out'):
        count.add_argument('--' + name, type=Path, required=True)
    count.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.mode == 'save-input':
        with args.output.open('x') as f:
            json.dump({'total_clean_reads': log_reads(args.log), 'source_log': str(args.log)}, f)
    else:
        recruit(args.samples, args.bam_dir, args.fasta_dir, args.out, args.threads)


if __name__ == '__main__':
    main()
