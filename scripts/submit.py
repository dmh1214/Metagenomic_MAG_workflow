#!/usr/bin/env python3
"""Submit a full workflow or one section. Default is preview only."""
import argparse
import csv
import datetime
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path


READ1 = re.compile(r"^(.*?)(_R|_good_|_)1(_001)?(\.(?:fq|fastq)(?:\.gz)?)$")


def prepare_samples(input_path, strategy="auto"):
    input_path = Path(input_path).resolve()
    manifest = input_path if input_path.is_file() else input_path/"samples.tsv"
    data = []
    if manifest.is_file():
        with open(manifest, encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                paths = [(manifest.parent/row[col]).resolve() for col in ("r1", "r2")]
                data.append((row["sample"], str(paths[0]), str(paths[1]), row.get("group", "all") or "all"))
    else:
        for r1 in sorted(input_path.rglob("*")):
            match = READ1.fullmatch(r1.name)
            if not match or not r1.is_file():
                continue
            prefix, tag, lane, extension = match.groups()
            r2 = r1.with_name(prefix+tag+"2"+(lane or "")+extension)
            relative = r1.parent.relative_to(input_path)
            group = "__".join(relative.parts) if relative.parts else "all"
            data.append((prefix, str(r1), str(r2), group))
    if strategy == "auto":
        strategy = "grouped" if len({r[3] for r in data}) > 1 else "joint"
    if strategy == "joint":
        data = [(s,r1,r2,"all") for s,r1,r2,g in data]
    if not data or len({r[0] for r in data}) != len(data):
        raise ValueError("Empty or duplicated sample IDs")
    for sample, r1, r2, group in data:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", group):
            raise ValueError("Invalid sample/group ID")
        for path in (r1, r2):
            if not Path(path).is_file() or any(c in path for c in ",\t\r\n"):
                raise ValueError("Missing read or unsupported path: "+path)
    return sorted(data), strategy


def create_plan(stage, groups, samples, map_concurrency, bin_concurrency):
    plan = []
    def add(name, script, resource, dependencies=(), array="", extra=()):
        plan.append(dict(name=name, script=script, resource=resource, dependencies=list(dependencies), array=array, extra=list(extra)))
    if stage in ("all", "assembly", "assembly_continue"):
        add("assembly", "01_assembly.sh", "assembly", array=f"1-{groups}", extra=("continue",) if stage == "assembly_continue" else ())
    if stage in ("all", "binning"):
        add("assembly_index", "02_assembly_index.sh", "index", ("assembly",) if stage == "all" else (), f"1-{groups}")
        add("assembly_mapping", "03_assembly_mapping.sh", "map", ("assembly_index",), f"1-{groups*samples}%{map_concurrency}")
        add("binning", "04_binning.sh", "bin", ("assembly_mapping",), f"1-{groups*3}%{bin_concurrency}")
        add("dastool", "05_dastool.sh", "dastool", ("binning",), f"1-{groups}")
    if stage in ("all", "binning", "quality"):
        add("quality", "06_checkm2.sh", "catalog", ("dastool",) if stage != "quality" else ())
    if stage in ("all", "binning", "catalog"):
        add("catalog", "06b_catalogue.sh", "catalog", ("quality",) if stage != "catalog" else ())
    if stage in ("all", "annotate", "taxonomy"):
        add("taxonomy", "07_taxonomy.sh", "taxonomy", ("catalog",) if stage == "all" else ())
    if stage in ("all", "annotate", "dram"):
        add("dram", "08_dram.sh", "dram", ("catalog",) if stage == "all" else ())
    if stage in ("all", "distill"):
        add("distill", "09_distill.sh", "distill", ("dram",) if stage == "all" else ())
    if stage in ("all", "quantify"):
        add("catalogue_index", "10_catalogue_index.sh", "index", ("catalog",) if stage == "all" else ())
        add("catalogue_mapping", "11_catalogue_mapping.sh", "map", ("catalogue_index",), f"1-{samples}%{map_concurrency}")
    if stage in ("all", "quantify", "tables"):
        deps = ("catalogue_mapping", "dram") if stage == "all" else ("catalogue_mapping",) if stage == "quantify" else ()
        add("tables", "12_tables.sh", "tables", deps)
    if stage in ("all", "quantify", "recruitment"):
        deps = ("catalogue_mapping",) if stage != "recruitment" else ()
        add("recruitment", "13_read_recruitment.sh", "recruitment", deps)
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Clean-read directory, or samples.tsv")
    parser.add_argument("--strategy", choices=("auto", "joint", "grouped"), default="auto")
    parser.add_argument("--stage", default="all", choices=("all", "assembly", "assembly_continue", "binning", "quality", "catalog", "annotate", "taxonomy", "dram", "distill", "quantify", "tables", "recruitment"))
    parser.add_argument("--config", type=Path, default=Path(__file__).resolve().parents[1]/"config.json")
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--assembly-root", type=Path)
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument("--bam-dir", type=Path)
    parser.add_argument("--after", default="", help="Optional upstream Slurm job ID(s), colon separated")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    if not args.config.is_file():
        parser.error("Create config.json first: cp config.example.json config.json")
    config = json.loads(args.config.read_text(encoding="utf-8-sig"))
    samples, strategy = prepare_samples(args.input, args.strategy)
    groups = sorted({r[3] for r in samples})
    source = Path(__file__).resolve().parent
    input_root = args.input.resolve() if args.input.is_dir() else args.input.resolve().parent
    run = (args.run_dir or input_root.parent/"MAG_results"/(input_root.name+"_"+strategy)).resolve()
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    dataset = input_root.parent.name+"_"+input_root.name+"_"+strategy
    log = source.parent/"logs"/dataset/stamp
    log.mkdir(parents=True)
    snapshot = source/".runs"/dataset/stamp
    snapshot.mkdir(parents=True)
    for path in source.iterdir():
        if path.suffix in (".sh", ".py", ".json"):
            shutil.copy2(path, snapshot/path.name)
    (snapshot/"effective_config.json").write_text(json.dumps(config, indent=2))
    with open(snapshot/"samples.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "r1", "r2", "group"])
        writer.writerows(samples)
    (snapshot/"groups.txt").write_text("\n".join(groups)+"\n")
    env = {"RUN": str(run), "ANALYSIS": str((args.analysis_dir or run/"analysis").resolve()),
           "ASSEMBLY_ROOT": str((args.assembly_root or run/"assembly").resolve()),
           "CATALOGUE_BAMS": str((args.bam_dir or run/"catalogue_bam").resolve()),
           "TABLES": str(run/"tables"), "LOG_DIR": str(log), "STRATEGY": strategy,
           "SAMPLE_COUNT": len(samples), "CONDA_SH": os.path.expanduser(os.path.expandvars(config["conda_sh"])), "DASTOOL_ENV": config["dastool_env"],
           "CHECKM2_DB": os.path.expanduser(os.path.expandvars(config["checkm2_db"])), "GTDB_DB": os.path.expanduser(os.path.expandvars(config["gtdb_db"])),
           "DRAM_CONFIG": os.path.expandvars(config["dram_config"]), "ANI": config["ani"],
           "MIN_COMPLETENESS": config["min_completeness"], "MAX_CONTAMINATION": config["max_contamination"]}
    env.update({name.upper()+"_MODULE": value for name,value in config["modules"].items()})
    (snapshot/"run.env.sh").write_text("\n".join(f"{k}={shlex.quote(str(v))}" for k,v in env.items())+"\n")
    if args.submit:
        run.mkdir(parents=True, exist_ok=True)
    plan = create_plan(args.stage, len(groups), len(samples), config["map_concurrency"], config["bin_concurrency"])
    ids = {}
    with open(log/"jobs.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "job_id", "dependencies"])
        for number, job in enumerate(plan, 1):
            resources = config["resources"]
            cpus, memory, wall = resources.get("recruitment", [8, "8G", "24:00:00"]) if job["resource"] == "recruitment" else resources[job["resource"]]
            command = ["sbatch", "--parsable", "--nodes=1", "--ntasks=1",
                       f"--cpus-per-task={cpus}", f"--mem={memory}", f"--time={wall}",
                       f"--job-name=mag_{job['name']}", f"--output={log}/{job['name']}_%A_%a.out",
                       f"--error={log}/{job['name']}_%A_%a.err"]
            for option in ("partition", "qos"):
                if config.get(option):
                    command.append("--"+option+"="+config[option])
            if config.get("account"):
                command.append("--account="+config["account"])
            deps = [ids[name] for name in job["dependencies"]]
            if args.after and not deps:
                deps.append(args.after)
            if deps:
                command.append("--dependency=afterok:"+":".join(deps))
            if job["array"]:
                command.append("--array="+job["array"])
            command.extend([str(snapshot/job["script"]), str(snapshot), *job["extra"]])
            print(shlex.join(command), flush=True)
            if args.submit:
                job_id = subprocess.check_output(command, text=True).strip().split(";")[0]
                if not job_id.isdigit():
                    raise ValueError("Unexpected sbatch output: "+job_id)
            else:
                job_id = "preview_"+str(number)
            ids[job["name"]] = job_id
            writer.writerow([job["name"], job_id, ":".join(deps)])
            handle.flush()
    print(f"{'SUBMITTED' if args.submit else 'PREVIEW ONLY'} | {len(samples)} samples | {len(groups)} assemblies | {strategy}")
    print("Results:", run)
    print("Logs:", log)


if __name__ == "__main__":
    main()
